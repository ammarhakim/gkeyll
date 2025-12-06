#include <gkyl_bc_return_flux_gyrokinetic.h>
#include <gkyl_bc_return_flux_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_bc_return_flux_gyrokinetic*
gkyl_bc_return_flux_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, int displacement_dir,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *bflux_ghost_r,
  const struct gkyl_range *skin_r, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_return_flux_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->dir = dir;
  up->edge = edge;
  up->disp_dir = displacement_dir;
  up->disp_dir_num_cellsD2 = grid->cells[displacement_dir]/2.0;
  up->bflux_ghost_r = bflux_ghost_r;
  up->skin_r = skin_r;
  up->use_gpu = use_gpu;

  return up;
}

void
gkyl_bc_return_flux_gyrokinetic_advance(const struct gkyl_bc_return_flux_gyrokinetic *up,
  const struct gkyl_array *bflux, struct gkyl_array *rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_return_flux_gyrokinetic_advance_cu(up, bflux, rhs);
    return;
  }
#endif

  int gidx[GKYL_MAX_DIM]; // Ghost index.
  int disp_cellsD2 = up->disp_dir_num_cellsD2; // Number of cells in displacement direction / 2.

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, up->skin_r);
  while (gkyl_range_iter_next(&iter)) {

    // Fill ghost index.
    gkyl_copy_int_arr(up->skin_r->ndim, iter.idx, gidx);
    gidx[up->dir] = up->edge == GKYL_LOWER_EDGE? gidx[up->dir]-1 : gidx[up->dir]+1;
    // Shift the index in the displacement direction
    int disp_idx = iter.idx[up->disp_dir];
    gidx[up->disp_dir] = disp_idx < disp_cellsD2 ? disp_idx+disp_cellsD2 : disp_idx-disp_cellsD2;

    long linidx_ghost = gkyl_range_idx(up->bflux_ghost_r, gidx);
    long linidx_skin = gkyl_range_idx(up->skin_r, iter.idx);

    const double *bflux_d = (const double*) gkyl_array_cfetch(bflux, linidx_ghost);
    double *rhs_d = (double*) gkyl_array_fetch(rhs, linidx_skin);

    for (int i=0; i<rhs->ncomp; i++)
      rhs_d[i] += bflux_d[i];
  }
}

void gkyl_bc_return_flux_gyrokinetic_release(struct gkyl_bc_return_flux_gyrokinetic *up)
{
  // Release memory associated with this updater.
  gkyl_free(up);
}
