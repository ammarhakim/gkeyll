/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_return_flux_gyrokinetic.h>
#include <gkyl_bc_return_flux_gyrokinetic_priv.h>
}

__global__ static void
gkyl_bc_return_flux_gyrokinetic_advance_cu_ker(int dir, enum gkyl_edge_loc edge, int disp_dir, int disp_cellsD2,
  const struct gkyl_range bflux_ghost_r, const struct gkyl_range skin_r, 
  const struct gkyl_array *bflux, struct gkyl_array *rhs)
{
  int idx[GKYL_MAX_DIM];
  int gidx[GKYL_MAX_DIM]; // Ghost index shifted along displacement direction.

  for (unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < skin_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&skin_r, linc, idx);

    // Fill ghost index.
    gkyl_copy_int_arr(skin_r.ndim, idx, gidx);
    gidx[dir] = edge == GKYL_LOWER_EDGE? gidx[dir]-1 : gidx[dir]+1;
    // Shift the index in the displacement direction.
    int disp_idx = idx[disp_dir];
    gidx[disp_dir] = disp_idx < disp_cellsD2 ? disp_idx+disp_cellsD2 : disp_idx-disp_cellsD2;

    long linidx_ghost = gkyl_range_idx(&bflux_ghost_r, gidx);
    long linidx_skin = gkyl_range_idx(&skin_r, idx);

    const double *bflux_d = (const double*) gkyl_array_cfetch(bflux, linidx_ghost);
    double *rhs_d = (double*) gkyl_array_fetch(rhs, linidx_skin);

    for (int i=0; i<rhs->ncomp; i++)
      rhs_d[i] += bflux_d[i];
  }
}

void
gkyl_bc_return_flux_gyrokinetic_advance_cu(const struct gkyl_bc_return_flux_gyrokinetic *up,
  const struct gkyl_array *bflux, struct gkyl_array *rhs)
{
  int nblocks = up->skin_r->nblocks, nthreads = up->skin_r->nthreads;

  gkyl_bc_return_flux_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>(up->dir, up->edge,
    up->disp_dir, up->disp_dir_num_cellsD2, *up->bflux_ghost_r, *up->skin_r, bflux->on_dev, rhs->on_dev);
}
