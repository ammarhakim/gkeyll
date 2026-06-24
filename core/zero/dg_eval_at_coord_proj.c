#include <assert.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eval_at_coord_proj.h>
#include <gkyl_dg_eval_at_coord_proj_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

static int
eval_dirs_to_mask(int num_eval_dirs, const int *eval_dirs)
{
  // Encode eval_dirs as a bitmask: bit d is set iff direction d is evaluated.
  // Used to index the dispatch table as [ndim_do-1][mask-1].
  int mask = 0;
  for (int i=0; i<num_eval_dirs; i++)
    mask |= (1 << eval_dirs[i]);
  return mask;
}

struct gkyl_dg_eval_at_coord_proj*
gkyl_dg_eval_at_coord_proj_new(const struct gkyl_basis *basis_do, const struct gkyl_basis *basis_tar,
  int num_eval_dirs, const int *eval_dirs, bool use_gpu)
{
  int ndim_do  = basis_do->ndim;
  int ndim_tar = basis_tar ? basis_tar->ndim : 0;

  assert(ndim_do >= 1 && ndim_do <= 3);
  assert(num_eval_dirs >= 1 && num_eval_dirs <= ndim_do);
  assert(ndim_tar == ndim_do - num_eval_dirs);
  assert(basis_do->poly_order >= 1 && basis_do->poly_order <= 3);

  struct gkyl_dg_eval_at_coord_proj *up = gkyl_malloc(sizeof(*up));

  up->use_gpu       = use_gpu;
  up->ndim_do       = ndim_do;
  up->ndim_tar      = ndim_tar;
  up->num_basis_do  = basis_do->num_basis;
  up->num_basis_tar = (ndim_tar > 0 && basis_tar)? basis_tar->num_basis : 1;
  up->num_eval_dirs = num_eval_dirs;

  for (int i=0; i<num_eval_dirs; i++)
    up->eval_dirs[i] = eval_dirs[i];

  for (int d=0; d<GKYL_MAX_DIM; d++)
    up->is_eval[d] = false;

  for (int i=0; i<num_eval_dirs; i++)
    up->is_eval[eval_dirs[i]] = true;

  up->dir_mask = eval_dirs_to_mask(num_eval_dirs, eval_dirs);

  int poly_order = basis_do->poly_order;

  switch (basis_do->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      up->kernel = choose_ser_eval_at_coord_kern(ndim_do, up->dir_mask, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      up->kernel = choose_ten_eval_at_coord_kern(ndim_do, up->dir_mask, poly_order);
      break;
    default:
      assert(false);
      break;
  }
  assert(up->kernel);

  return up;
}

void
gkyl_dg_eval_at_coord_proj_advance(struct gkyl_dg_eval_at_coord_proj *up, const double *eval_coords,
  const struct gkyl_rect_grid *grid, const bool *pick_lower, const int *known_index,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_dg_eval_at_coord_proj_advance_cu(up, eval_coords, grid, pick_lower, known_index,
      rng_do, rng_tar, fdo, ftar);
    return;
  }
#endif

  // Build full ndim_do-dimensional point for gkyl_rect_grid_find_cell.
  // Evaluated directions use eval_coords[i]; non-evaluated dirs use the
  // center of the lowest cell in rng_do for that direction.
  double point[GKYL_MAX_DIM];
  int eval_ctr = 0;
  for (int d=0; d<up->ndim_do; d++) {
    if (up->is_eval[d])
      point[d] = eval_coords[eval_ctr++];
    else
      point[d] = grid->lower[d] + (rng_do->lower[d] - 0.5) * grid->dx[d];
  }

  int cell_idx[GKYL_MAX_DIM];
  gkyl_rect_grid_find_cell(grid, point, pick_lower, known_index, cell_idx);

  // Convert comp eval_coords to logical coords.
  double eval_coords_log[GKYL_MAX_DIM];
  for (int i=0; i<up->num_eval_dirs; i++) {
    int d = up->eval_dirs[i];
    double xc_d = grid->lower[d] + (cell_idx[d] - 0.5) * grid->dx[d];
    eval_coords_log[i] = 2.0 * (eval_coords[i] - xc_d) / grid->dx[d];
    // Force to be between -1 and 1, to account for floating point precision
    // causing results outside of [-1,1].
    eval_coords_log[i] = GKYL_MAX2(-1.0, GKYL_MIN2(eval_coords_log[i],1.0) );
  }

  int idx_do[GKYL_MAX_DIM] = {0};
  int ncomp = fdo->ncomp / up->num_basis_do;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, rng_tar);
  while (gkyl_range_iter_next(&iter)) {
    eval_at_coord_get_idx_do(up->is_eval, up->ndim_do, iter.idx, cell_idx, idx_do);

    long linidx_do  = gkyl_range_idx(rng_do,  idx_do);
    long linidx_tar = gkyl_range_idx(rng_tar, iter.idx);

    const double *fdo_c  = gkyl_array_cfetch(fdo,  linidx_do);
    double       *ftar_c = gkyl_array_fetch(ftar, linidx_tar);

    for (int n=0; n<ncomp; n++)
      up->kernel(eval_coords_log, fdo_c+n*up->num_basis_do, ftar_c+n*up->num_basis_tar);
  }
}

void
gkyl_dg_eval_at_coord_proj_release(struct gkyl_dg_eval_at_coord_proj *up)
{
  gkyl_free(up);
}
