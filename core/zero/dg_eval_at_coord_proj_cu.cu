/* -*- c++ -*- */

extern "C" {
#include <gkyl_array_ops_priv.h>
#include <dg_eval_at_coord_proj.h>
#include <dg_eval_at_coord_proj_priv.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
}

__global__ void
dg_eval_at_coord_proj_range_cu_kernel(struct dg_eval_at_coord_proj up_val,
  dg_evproj_struct_double_t eval_coords_log, dg_evproj_struct_int_t cell_idx,
  struct gkyl_array *ftar, const struct gkyl_array *fdo,
  int ncomp, struct gkyl_range rng_tar, struct gkyl_range rng_do)
{
  int idx_tar[GKYL_MAX_DIM];
  int idx_do[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < rng_tar.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&rng_tar, linc1, idx_tar);

    eval_at_coord_get_idx_do(up_val.is_eval, up_val.ndim_do,
      idx_tar, cell_idx.c, idx_do);

    long start_do  = gkyl_range_idx(&rng_do,  idx_do);
    long start_tar = gkyl_range_idx(&rng_tar, idx_tar);

    const double *fdo_c  = (const double *) gkyl_array_cfetch(fdo,  start_do);
    double       *ftar_c = (double *)       gkyl_array_fetch(ftar, start_tar);

    for (int n = 0; n < ncomp; n++)
      up_val.kernel(eval_coords_log.c, fdo_c+n*up_val.num_basis_do, ftar_c+n*up_val.num_basis_tar);
  }
}

void
dg_eval_at_coord_proj_advance_cu(struct dg_eval_at_coord_proj *up, const double *eval_coords,
  const struct gkyl_rect_grid *grid, const bool *pick_lower, const int *known_index,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar)
{
  int ncomp = fdo->ncomp / up->num_basis_do;

  // Build full ndim_do-dimensional point for gkyl_rect_grid_find_cell.
  double point[GKYL_MAX_DIM];
  int eval_ctr = 0;
  for (int d=0; d<up->ndim_do; d++) {
    if (up->is_eval[d])
      point[d] = eval_coords[eval_ctr++];
    else
      point[d] = grid->lower[d] + (rng_do->lower[d] - 0.5) * grid->dx[d];
  }

  dg_evproj_struct_int_t cell_idx = {0};
  gkyl_rect_grid_find_cell(grid, point, pick_lower, known_index, cell_idx.c);

  // Convert comp eval_coords to logical coords.
  dg_evproj_struct_double_t eval_coords_log = {0};
  for (int i=0; i<up->num_eval_dirs; i++) {
    int d = up->eval_dirs[i];
    double xc_d = grid->lower[d] + (cell_idx.c[d] - 0.5) * grid->dx[d];
    eval_coords_log.c[i] = 2.0 * (eval_coords[i] - xc_d) / grid->dx[d];
  }

  int nblocks = rng_tar->nblocks;
  int nthreads = rng_tar->nthreads;

  dg_eval_at_coord_proj_range_cu_kernel<<<nblocks, nthreads>>>(*up, eval_coords_log, cell_idx,
    ftar->on_dev, fdo->on_dev, ncomp, *rng_tar, *rng_do);
}
