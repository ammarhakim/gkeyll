#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_gr_maxwell_slope_limiter.h>
#include <gkyl_dg_gr_maxwell_slope_limiter_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_util.h>

gkyl_dg_gr_maxwell_slope_limiter*
gkyl_dg_gr_maxwell_slope_limiter_inew(const struct gkyl_dg_gr_maxwell_slope_limiter_inp *inp)
{
  assert(!inp->use_gpu);  // GPU path not implemented yet
  struct gkyl_dg_gr_maxwell_slope_limiter *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;
  enum gkyl_basis_type b_type = inp->conf_basis->b_type;

  up->cdim = cdim;
  up->conf_grid = *inp->conf_grid;
  up->use_gpu = inp->use_gpu;

  // Default limiter_fac matches the Euler slope-limiter convention:
  // cell_avg(f) = f0 * 1/sqrt(2^cdim), slope basis psi_1 = sqrt(3)/sqrt(2^cdim) * xi,
  // so the cell-slope mode that matches an adjacent-cell average jump is
  // (1/sqrt(3)) * (q[0]_R - q[0]_L). Decrease the factor (e.g. to 0.5/sqrt(3))
  // to increase diffusion; increase it to decrease diffusion.
  up->limiter_fac = (inp->limiter_fac == 0.0) ? 0.5773502691896258 : inp->limiter_fac;

  for (int d = 0; d < cdim; ++d) {
    up->kernels[d] = choose_slope_limiter_kern(d, b_type, cdim, poly_order);
  }
  // Require x-direction kernel to exist (current implementation: 2D ser p=1).
  // Other directions may be NULL on basis/order combinations we haven't
  // implemented yet; we skip those passes.

  for (int i = 0; i < cdim; ++i) {
    up->theta_pole_lo[i] = inp->theta_pole_lo ? inp->theta_pole_lo[i] : 0;
    up->theta_pole_up[i] = inp->theta_pole_up ? inp->theta_pole_up[i] : 0;
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up;

  return up;
}

void
gkyl_dg_gr_maxwell_slope_limiter_advance(struct gkyl_dg_gr_maxwell_slope_limiter *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij,
  const struct gkyl_surf_and_vol_node_arrays *det_h,
  struct gkyl_array *field_con)
{
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  // Apply r-direction limiter, then theta-direction. Per design doc Section
  // 10, this ordering matches the "r-first then theta" convention; for our
  // axisymmetric setup it's the natural choice. q[3] is zeroed conservatively
  // by either pass if it fires.
  for (int dir = 0; dir < cdim; ++dir) {
    slope_limiter_kern_t kfn = up->kernels[dir];
    if (kfn == NULL) continue;  // direction not implemented for this basis/order

    // Which surface nodal array to use for this direction's face metric.
    // Convention: each cell owns its lower-d face nodal arrays at cidx.
    const struct gkyl_array *lapse_surf  = (dir == 0) ? lapse->nodal_arr_surf_x
                                          : (dir == 1) ? lapse->nodal_arr_surf_y
                                          : lapse->nodal_arr_surf_z;
    const struct gkyl_array *shift_surf  = (dir == 0) ? shift->nodal_arr_surf_x
                                          : (dir == 1) ? shift->nodal_arr_surf_y
                                          : shift->nodal_arr_surf_z;
    const struct gkyl_array *h_ij_surf   = (dir == 0) ? h_ij->nodal_arr_surf_x
                                          : (dir == 1) ? h_ij->nodal_arr_surf_y
                                          : h_ij->nodal_arr_surf_z;
    const struct gkyl_array *det_h_surf  = (dir == 0) ? det_h->nodal_arr_surf_x
                                          : (dir == 1) ? det_h->nodal_arr_surf_y
                                          : det_h->nodal_arr_surf_z;

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, conf_range);
    while (gkyl_range_iter_next(&iter)) {
      gkyl_copy_int_arr(cdim, iter.idx, idxc);
      long linc = gkyl_range_idx(conf_range, idxc);

      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);
      idxl[dir] -= 1;
      idxr[dir] += 1;

      long linl = gkyl_range_idx(conf_range, idxl);
      long linr = gkyl_range_idx(conf_range, idxr);

      // Pole-adjacency flags for this cell along this direction. Only meaningful
      // for dir == 1 (theta) in our current setup, but we pass them generically
      // and let the kernel decide.
      int pole_lo = (idxc[dir] == conf_range->lower[dir]) ? up->theta_pole_lo[dir] : 0;
      int pole_up = (idxc[dir] == conf_range->upper[dir]) ? up->theta_pole_up[dir] : 0;

      const double *lapse_d = gkyl_array_cfetch(lapse_surf, linc);
      const double *shift_d = gkyl_array_cfetch(shift_surf, linc);
      const double *h_ij_d  = gkyl_array_cfetch(h_ij_surf,  linc);
      const double *det_h_d = gkyl_array_cfetch(det_h_surf, linc);

      const double *fc_l = gkyl_array_cfetch(field_con, linl);
      const double *fc_c = gkyl_array_cfetch(field_con, linc);
      const double *fc_r = gkyl_array_cfetch(field_con, linr);
      double *fc_out = gkyl_array_fetch(field_con, linc);

      kfn(up->limiter_fac, pole_lo, pole_up,
        lapse_d, shift_d, h_ij_d, det_h_d,
        fc_l, fc_c, fc_r, fc_out);
    }
  }
}

void
gkyl_dg_gr_maxwell_slope_limiter_release(struct gkyl_dg_gr_maxwell_slope_limiter *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
