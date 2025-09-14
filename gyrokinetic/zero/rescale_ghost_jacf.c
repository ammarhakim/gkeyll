#include <gkyl_rescale_ghost_jacf.h>
#include <gkyl_rescale_ghost_jacf_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_rescale_ghost_jacf*
gkyl_rescale_ghost_jacf_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_rescale_ghost_jacf *up = gkyl_malloc(sizeof(*up));

  up->dir = dir;
  up->edge = edge;
  up->use_gpu = use_gpu;

  int poly_order = phase_basis->poly_order;
  assert(poly_order == 1); // MF 2024/10/30: Because of the inv_op below.

  // Choose the kernel that does the skin surf from ghost copy
  if (!use_gpu)
    up->kernels = gkyl_malloc(sizeof(struct gkyl_rescale_ghost_jacf_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_rescale_ghost_jacf_kernels));
#endif

  rescale_ghost_jacf_choose_kernel(up->kernels, dir, edge, conf_basis, phase_basis, use_gpu);

  return up;
}

void
gkyl_rescale_ghost_jacf_advance(const struct gkyl_rescale_ghost_jacf *up,
  const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, const struct gkyl_array *jac_sync, struct gkyl_array *jf)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_rescale_ghost_jacf_advance_cu(up,
      conf_skin_r, conf_ghost_r, phase_ghost_r, jac_sync, jf);
    return;
  }
#endif

  int sidx[GKYL_MAX_DIM]; // Skin index.
  int cdim = conf_ghost_r->ndim; 

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; ++d) rem_dir[d] = 1;

  const int num_basis_conf_surf = 8; // MF 2024/11/12: Hardcoded to 2x p2 for now.
  const int num_basis_phase_surf = 24; // MF 2024/11/12: Hardcoded to 2x2v GkHybrid for now.
 
  gkyl_range_iter_init(&conf_iter, conf_ghost_r);
  while (gkyl_range_iter_next(&conf_iter)) {

    // Get skin cell corresponding to this ghost cell.
    gkyl_copy_int_arr(cdim, conf_iter.idx, sidx);
    sidx[up->dir] = up->edge == GKYL_LOWER_EDGE? conf_iter.idx[up->dir]+1 : conf_iter.idx[up->dir]-1; 

    long clinidx_skin  = gkyl_range_idx(conf_skin_r, sidx);
    long clinidx_ghost = gkyl_range_idx(conf_ghost_r, conf_iter.idx);

    const double *jacskin_surf_c = gkyl_array_cfetch(jac_sync, clinidx_skin);
    const double *jacghost_surf_c = gkyl_array_cfetch(jac_sync, clinidx_ghost);

    // Compute the reciprocal of the ghost cell jacobian.
    double jacghost_surf_inv_c[num_basis_conf_surf];
    up->kernels->conf_inv_op(jacghost_surf_c, jacghost_surf_inv_c);

    gkyl_range_deflate(&vel_rng, phase_ghost_r, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      long plinidx_ghost = gkyl_range_idx(&vel_rng, vel_iter.idx);
      double *jf_c = gkyl_array_fetch(jf, plinidx_ghost);

      // Deflate Jf in the ghost cell.
      double jf_surf_c[num_basis_phase_surf];
      up->kernels->deflate_phase_ghost_op(jf_c, jf_surf_c);

      // Multiply the surface Jf by the surface 1/J_ghost and by the surface J_skin.
      up->kernels->conf_phase_mul_op(jacghost_surf_inv_c, jf_surf_c, jf_surf_c);
      up->kernels->conf_phase_mul_op(jacskin_surf_c, jf_surf_c, jf_surf_c);

      // Inflate Jf.
      for (int k=0; k<jf->ncomp; k++)
        jf_c[k] = 0.0;
      up->kernels->inflate_phase_ghost_op(jf_surf_c, jf_c);
    }
  }
}

void
gkyl_rescale_ghost_jacf_release(struct gkyl_rescale_ghost_jacf *up)
{
  // Release memory associated with this updater.
  if (!up->use_gpu)
    gkyl_free(up->kernels);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels);
#endif
  gkyl_free(up);
}
