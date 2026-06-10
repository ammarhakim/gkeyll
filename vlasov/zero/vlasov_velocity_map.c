#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm_io.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_dg_vlasov_divide_Jv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_velocity_map_priv.h>

// Computational coordinates = physical coordinates mapping (default).
static void
vlasov_vmap_identity(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = vc[0];
}

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  return on_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Populate the (host) velocity map arrays from the input 1D mappings by
// constructing a C^1 piecewise cubic representation in each velocity
// dimension and evaluating its derivative at the Gauss-Legendre quadrature
// points needed by volume, surface, and projection operations.
static void
vlasov_velocity_map_populate(const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange,
  int v_poly_order, struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM],
  struct gkyl_array *vmap, struct gkyl_array *jacob_vel, struct gkyl_array *jacob_vel_surf,
  struct gkyl_array *vmap_pgkyl, struct gkyl_array *vmap_avg_pgkyl,
  struct gkyl_array *jacob_vel_gauss)
{
  int vdim = vgrid->ndim;
  struct gkyl_array *v_nodal[3];
  struct gkyl_array *v_cubic[3];
  struct gkyl_dg_basis_op_mem *mem[3];

  // Make 1D cubic basis for constructing C^1 expansion
  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 1, 3);

  // 1D ranges for indexing 1D cubic bases
  struct gkyl_range local[3], local_ext[3];

  // Loop over number of dimensions and construct 1D mappings
  for (int i=0; i<vdim; ++i) {
    double lower[] = { vgrid->lower[i] }, upper[] = { vgrid->upper[i] };
    int cells[] = { vgrid->cells[i] };

    struct gkyl_rect_grid grid_1d;
    gkyl_rect_grid_init(&grid_1d, 1, lower, upper, cells);

    // nodal grid used for constructing physical coordinates
    double nc_lower[] = { lower[0] - 0.5*grid_1d.dx[0] };
    double nc_upper[] = { upper[0] + 0.5*grid_1d.dx[0] };
    int nc_cells[] = { cells[0] + 1 };
    struct gkyl_rect_grid nc_grid;
    gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

    int nghost[GKYL_MAX_CDIM] = { 0, 0 };
    gkyl_create_grid_ranges(&grid_1d, nghost, &local_ext[i], &local[i]);

    struct gkyl_range nc_local, nc_local_ext;
    gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

    v_nodal[i] = gkyl_array_new(GKYL_DOUBLE, 1, cells[0]+1);
    v_cubic[i] = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext[i].volume);
    mem[i] = gkyl_dg_alloc_cubic_1d(cells[0]);
    double xn[1];

    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid_1d, iter.idx, xn);

      double *pn = gkyl_array_fetch(v_nodal[i], nidx);
      inp_vmap[i].eval_vmap(0.0, xn, pn, inp_vmap[i].ctx);
    }

    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem[i], cells[0], grid_1d.dx[0],
      v_nodal[i], v_cubic[i]);
  }

  // initialize the mapping
  vmap_cubic_t vmap_op = choose_vmap_kern(vdim, v_poly_order);
  const double *v_cubic_dir[3]; // 1D cubic in each direction
  int vidx_1D[1]; // 1D index for indexing correct cubic mapping

  struct gkyl_range_iter iter_vmap;
  gkyl_range_iter_init(&iter_vmap, vrange);
  while (gkyl_range_iter_next(&iter_vmap)) {
    long loc_vel = gkyl_range_idx(vrange, iter_vmap.idx);

    for (int i=0; i<vdim; ++i) {
      vidx_1D[0] = iter_vmap.idx[i];
      long loc_vel_1D = gkyl_range_idx(&local[i], vidx_1D);
      v_cubic_dir[i] = gkyl_array_cfetch(v_cubic[i], loc_vel_1D);
    }
    double *vmap_d = gkyl_array_fetch(vmap, loc_vel);
    double *jacob_vel_d = gkyl_array_fetch(jacob_vel, loc_vel);
    double *jacob_vel_surf_d = gkyl_array_fetch(jacob_vel_surf, loc_vel);
    double *vmap_pgkyl_d = gkyl_array_fetch(vmap_pgkyl, loc_vel);
    double *vmap_avg_pgkyl_d = gkyl_array_fetch(vmap_avg_pgkyl, loc_vel);
    double *jacob_vel_gauss_d = gkyl_array_fetch(jacob_vel_gauss, loc_vel);

    vmap_op(vgrid->dx, v_cubic_dir, vmap_d, jacob_vel_d, jacob_vel_surf_d,
      vmap_pgkyl_d, vmap_avg_pgkyl_d, jacob_vel_gauss_d);
  }

  // free temporary memory
  for (int i=0; i<vdim; ++i) {
    gkyl_array_release(v_nodal[i]);
    gkyl_array_release(v_cubic[i]);
    gkyl_dg_basis_op_mem_release(mem[i]);
  }
}

static void
gkyl_vlasov_velocity_map_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_vlasov_velocity_map *vvm = container_of(ref, struct gkyl_vlasov_velocity_map, ref_count);

  gkyl_array_release(vvm->vmap);
  gkyl_array_release(vvm->jacob_vel);
  gkyl_array_release(vvm->jacob_vel_surf);
  gkyl_array_release(vvm->jacob_vel_gauss);

  gkyl_array_release(vvm->vmap_host);
  gkyl_array_release(vvm->jacob_vel_host);
  gkyl_array_release(vvm->jacob_vel_surf_host);
  gkyl_array_release(vvm->jacob_vel_gauss_host);

  gkyl_array_release(vvm->vmap_pgkyl_host);
  gkyl_array_release(vvm->vmap_avg_pgkyl_host);

  gkyl_free(vvm);
}

struct gkyl_vlasov_velocity_map*
gkyl_vlasov_velocity_map_new(const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange,
  const struct gkyl_basis *vel_basis, struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM],
  bool use_gpu)
{
  struct gkyl_vlasov_velocity_map *vvm = gkyl_malloc(sizeof(*vvm));

  int vdim = vgrid->ndim;
  int poly_order = vel_basis->poly_order;

  vvm->grid_vel = *vgrid;
  vvm->local_vel = *vrange;
  vvm->basis_vel = *vel_basis;
  // C^0 linear representation for Serendipity velocity bases is a planned
  // extension; only the C^1 cubic representation is constructed.
  vvm->rep = GKYL_VLASOV_VMAP_C1_CUBIC;

  vvm->is_identity = true;
  for (int v=0; v<vdim; ++v) {
    if (inp_vmap[v].eval_vmap) {
      vvm->is_identity = false;
    }
  }
  vvm->is_mapped = (vel_basis->b_type == GKYL_BASIS_MODAL_TENSOR) || !vvm->is_identity;

  // velocity map is always a C^1 cubic representation in each direction (up to 3V; 3*4=12 components)
  vvm->vmap = mkarr(use_gpu, vdim*4, vrange->volume);
  // velocity-space Jacobian at quadrature points and "surface" quadrature points. Used to compute
  // 1/Jvi nodally in volume and surface operations respectively, with surface operations utilizing
  // more quadrature points to eliminate aliasing errors.
  vvm->jacob_vel = mkarr(use_gpu, vdim*(poly_order+1), vrange->volume);
  vvm->jacob_vel_surf = mkarr(use_gpu, vdim*(poly_order+2), vrange->volume);
  // need special basis sets to get the correct number of coefficients in 2V and 3V for constructing
  // the mapping in post-processing, as well as storing the velocity-space Jacobian at quadrature points.
  struct gkyl_basis vmap_pgkyl_basis, jacob_vel_basis;
  gkyl_cart_modal_serendip(&vmap_pgkyl_basis, vdim, 3);
  gkyl_cart_modal_tensor(&jacob_vel_basis, vdim, poly_order);
  // velocity-space Jacobian at Gaussian quadrature points for projecting distribution functions
  // and dividing out velocity-space Jacobian nodally.
  vvm->jacob_vel_gauss = mkarr(use_gpu, jacob_vel_basis.num_basis, vrange->volume);

  // host-side mirrors for initialization and I/O
  if (use_gpu) {
    vvm->vmap_host = mkarr(false, vvm->vmap->ncomp, vvm->vmap->size);
    vvm->jacob_vel_host = mkarr(false, vvm->jacob_vel->ncomp, vvm->jacob_vel->size);
    vvm->jacob_vel_surf_host = mkarr(false, vvm->jacob_vel_surf->ncomp, vvm->jacob_vel_surf->size);
    vvm->jacob_vel_gauss_host = mkarr(false, vvm->jacob_vel_gauss->ncomp, vvm->jacob_vel_gauss->size);
  }
  else {
    vvm->vmap_host = gkyl_array_acquire(vvm->vmap);
    vvm->jacob_vel_host = gkyl_array_acquire(vvm->jacob_vel);
    vvm->jacob_vel_surf_host = gkyl_array_acquire(vvm->jacob_vel_surf);
    vvm->jacob_vel_gauss_host = gkyl_array_acquire(vvm->jacob_vel_gauss);
  }

  // velocity map for I/O
  vvm->vmap_pgkyl_host = mkarr(false, vdim*vmap_pgkyl_basis.num_basis, vrange->volume);
  vvm->vmap_avg_pgkyl_host = mkarr(false, vdim, vrange->volume);

  // Default unspecified directions to the identity map.
  struct gkyl_vlasov_velocity_map_inp inp[GKYL_MAX_CDIM];
  for (int v=0; v<vdim; ++v) {
    if (inp_vmap[v].eval_vmap) {
      inp[v] = inp_vmap[v];
    }
    else {
      inp[v].eval_vmap = vlasov_vmap_identity;
      inp[v].ctx = 0;
    }
  }

  if (choose_vmap_kern(vdim, poly_order)) {
    vlasov_velocity_map_populate(vgrid, vrange, poly_order, inp,
      vvm->vmap_host, vvm->jacob_vel_host, vvm->jacob_vel_surf_host,
      vvm->vmap_pgkyl_host, vvm->vmap_avg_pgkyl_host, vvm->jacob_vel_gauss_host);

    if (use_gpu) {
      gkyl_array_copy(vvm->vmap, vvm->vmap_host);
      gkyl_array_copy(vvm->jacob_vel, vvm->jacob_vel_host);
      gkyl_array_copy(vvm->jacob_vel_surf, vvm->jacob_vel_surf_host);
      gkyl_array_copy(vvm->jacob_vel_gauss, vvm->jacob_vel_gauss_host);
    }
  }
  else {
    // No C^1 cubic kernel for this (vdim, poly_order) combination; only
    // allowed when the map is the identity (e.g. Serendipity p=1, whose
    // divide/rescale kernels do not reference these arrays).
    assert(vvm->is_identity);
    if (use_gpu) {
      gkyl_array_clear(vvm->vmap, 0.0);
      gkyl_array_clear(vvm->jacob_vel, 0.0);
      gkyl_array_clear(vvm->jacob_vel_surf, 0.0);
      gkyl_array_clear(vvm->jacob_vel_gauss, 0.0);
    }
  }

  vvm->flags = 0;
  if (use_gpu) {
    GKYL_SET_CU_ALLOC(vvm->flags);
  }
  vvm->ref_count = gkyl_ref_count_init(gkyl_vlasov_velocity_map_free);

  return vvm;
}

bool
gkyl_vlasov_velocity_map_is_cu_dev(const struct gkyl_vlasov_velocity_map *vvm)
{
  return GKYL_IS_CU_ALLOC(vvm->flags);
}

void
gkyl_vlasov_velocity_map_write(const struct gkyl_vlasov_velocity_map *vvm,
  struct gkyl_comm *comm, struct gkyl_msgpack_data *mt,
  const char *app_name, const char *species_name, bool write_cell_avg)
{
  int rank;
  gkyl_comm_get_rank(comm, &rank);
  if (rank == 0) {
    const char *fmt_vmap = "%s-%s_vmap.gkyl";
    int sz_vmap = gkyl_calc_strlen(fmt_vmap, app_name, species_name);
    char fileNm_vmap[sz_vmap+1]; // ensures no buffer overflow
    snprintf(fileNm_vmap, sizeof fileNm_vmap, fmt_vmap, app_name, species_name);

    gkyl_grid_sub_array_write(&vvm->grid_vel, &vvm->local_vel,
      mt, vvm->vmap_pgkyl_host, fileNm_vmap);

    if (write_cell_avg) {
      const char *fmt_vmap_avg = "%s-%s_vmap_avg.gkyl";
      int sz_vmap_avg = gkyl_calc_strlen(fmt_vmap_avg, app_name, species_name);
      char fileNm_vmap_avg[sz_vmap_avg+1]; // ensures no buffer overflow
      snprintf(fileNm_vmap_avg, sizeof fileNm_vmap_avg, fmt_vmap_avg, app_name, species_name);

      gkyl_grid_sub_array_write(&vvm->grid_vel, &vvm->local_vel,
        mt, vvm->vmap_avg_pgkyl_host, fileNm_vmap_avg);
    }
  }
}

void
gkyl_vlasov_velocity_map_divide_jacobvel(const struct gkyl_vlasov_velocity_map *vvm,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range,
  const struct gkyl_array *Jf, struct gkyl_array *f_no_J)
{
  gkyl_dg_vlasov_divide_Jv(conf_basis, phase_basis, &vvm->local_vel, phase_range,
    vvm->jacob_vel_gauss, Jf, f_no_J, gkyl_vlasov_velocity_map_is_cu_dev(vvm));
}

void
gkyl_vlasov_velocity_map_rescale_jacobvel(const struct gkyl_vlasov_velocity_map *vvm,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range,
  const struct gkyl_array *f_no_J, struct gkyl_array *Jf)
{
  gkyl_dg_vlasov_rescale_Jv(conf_basis, phase_basis, &vvm->local_vel, phase_range,
    vvm->jacob_vel_gauss, f_no_J, Jf, gkyl_vlasov_velocity_map_is_cu_dev(vvm));
}

struct gkyl_vlasov_velocity_map*
gkyl_vlasov_velocity_map_acquire(const struct gkyl_vlasov_velocity_map *vvm)
{
  gkyl_ref_count_inc(&vvm->ref_count);
  return (struct gkyl_vlasov_velocity_map*) vvm;
}

void
gkyl_vlasov_velocity_map_release(const struct gkyl_vlasov_velocity_map *vvm)
{
  gkyl_ref_count_dec(&vvm->ref_count);
}
