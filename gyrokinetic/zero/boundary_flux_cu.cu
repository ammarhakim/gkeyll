/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_boundary_flux.h>
#include <gkyl_boundary_flux_priv.h>
}

struct gkyl_boundary_flux*
gkyl_boundary_flux_cu_dev_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_dg_eqn *equation, double skip_cell_threshold)
{
  struct gkyl_boundary_flux *up = (struct gkyl_boundary_flux*) gkyl_malloc(sizeof(struct gkyl_boundary_flux));

  up->dir = dir;
  up->edge = edge;
  up->grid = *grid;
  up->skin_r = *skin_r;
  up->ghost_r = *ghost_r;
  up->use_gpu = true;

  if (skip_cell_threshold > 0.0)
    up->skip_cell_threshold = skip_cell_threshold * pow(sqrt(2.0), grid->ndim);
  else
    up->skip_cell_threshold = -1.0;

  struct gkyl_dg_eqn *eqn = gkyl_dg_eqn_acquire(equation);
  up->equation = eqn->on_dev;
  
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_boundary_flux *up_cu = (struct gkyl_boundary_flux*) gkyl_cu_malloc(sizeof(struct gkyl_boundary_flux));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_boundary_flux), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;
  
  up->equation = eqn;
  
  return up;
}

__global__ static void
gkyl_boundary_flux_advance_cu_ker(const struct gkyl_boundary_flux *up,
  const struct gkyl_array* fIn, struct gkyl_array* fluxOut)
{
  int idx_g[GKYL_MAX_DIM], idx_s[GKYL_MAX_DIM];
  double xc_g[GKYL_MAX_DIM], xc_s[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < up->ghost_r.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&up->ghost_r, tid, idx_g);

    gkyl_copy_int_arr(up->ghost_r.ndim, idx_g, idx_s);
    idx_s[up->dir] = up->edge == GKYL_LOWER_EDGE? idx_g[up->dir]+1 : idx_g[up->dir]-1;

    gkyl_rect_grid_cell_center(&up->grid, idx_g, xc_g);
    gkyl_rect_grid_cell_center(&up->grid, idx_s, xc_s);

    long linidx_g = gkyl_range_idx(&up->ghost_r, idx_g);
    long linidx_s = gkyl_range_idx(&up->skin_r, idx_s);

    const double* fg_c = (const double*) gkyl_array_cfetch(fIn, linidx_g);
    const double* fs_c = (const double*) gkyl_array_cfetch(fIn, linidx_s);
    double *fluxOut_g = (double*) gkyl_array_fetch(fluxOut, linidx_g);

    if (fabs(fs_c[0]) < up->skip_cell_threshold && fabs(fg_c[0]) < up->skip_cell_threshold)
    {
      for (int d=0; d<fluxOut->ncomp; ++d) {
        fluxOut_g[d] = 0.0;
      }
    }
    else {
      up->equation->boundary_surf_term(up->equation, up->dir, xc_s, xc_g,
        up->grid.dx, up->grid.dx, idx_s, idx_g, up->edge == GKYL_LOWER_EDGE? -1 : 1,
        fs_c, fg_c, fluxOut_g);
    }
  }
}

void
gkyl_boundary_flux_advance_cu(struct gkyl_boundary_flux *up,
  const struct gkyl_array *fIn, struct gkyl_array *fluxOut)
{
  int nblocks = up->ghost_r.nblocks, nthreads = up->ghost_r.nthreads;
  
  gkyl_boundary_flux_advance_cu_ker<<<nblocks, nthreads>>>(up->on_dev,
    fIn->on_dev, fluxOut->on_dev);
}
