#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>

#include <assert.h>

static void
bc_twistshift_filter_enabled(struct gkyl_dg_lowpass_filter *filt_up,
  struct gkyl_array *GKYL_RESTRICT finout, struct gkyl_array *GKYL_RESTRICT fbuff)
{
  gkyl_dg_lowpass_filter_advance(filt_up, finout, fbuff);
  gkyl_array_copy(finout, fbuff);
}

static void
bc_twistshift_filter_disabled(struct gkyl_dg_lowpass_filter *filt_up,
  struct gkyl_array *GKYL_RESTRICT finout, struct gkyl_array *GKYL_RESTRICT fbuff)
{
  // Do nothing.
}

static void
bc_twistshift_refine_enabled(struct gkyl_dg_interpolate *refine,
  struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  gkyl_dg_interpolate_advance(refine, fdo, ftar);
}

static void
bc_twistshift_refine_disabled(struct gkyl_dg_interpolate *refine,
  struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  gkyl_array_copy(ftar, fdo);
}

static void
bc_twistshift_coarsen_enabled(struct gkyl_dg_interpolate *coarsen,
  struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  gkyl_dg_interpolate_advance(coarsen, fdo, ftar);
}

static void
bc_twistshift_coarsen_disabled(struct gkyl_dg_interpolate *coarsen,
  struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  gkyl_array_copy(ftar, fdo);
}

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_new(const struct gkyl_bc_twistshift_inp *inp)
{
  struct gkyl_bc_twistshift *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = inp->use_gpu;
  up->filter_half_width = inp->filter_half_width;
  up->filter_cutoff_wavelength = inp->filter_cutoff_wavelength;
  up->upsample_factor = inp->upsample_factor > 1 ? inp->upsample_factor : 1;

  // Ensure that we use filter if we upsample.
  if (up->upsample_factor > 1) {
    assert(up->filter_half_width > 1);
    assert(up->filter_cutoff_wavelength > 0.0);
  }

  // Grid and update range for the twist-shift.
  const int ndim = inp->bcdir_ext_update_r.ndim;
  int fine_cells[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++) fine_cells[d] = inp->grid.cells[d];
  fine_cells[inp->shear_dir] *= up->upsample_factor;
  gkyl_rect_grid_init(&up->ts_grid, ndim, inp->grid.lower, inp->grid.upper, fine_cells);

  struct gkyl_range fine_ext, fine_local;
  gkyl_create_grid_ranges(&up->ts_grid, inp->num_ghost, &fine_ext, &fine_local);
  int flo[GKYL_MAX_DIM], fup[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++) 
  { 
    flo[d] = fine_local.lower[d]; 
    fup[d] = fine_local.upper[d]; 
  }
  flo[inp->bc_dir] = fine_ext.lower[inp->bc_dir];
  fup[inp->bc_dir] = fine_ext.upper[inp->bc_dir];
  gkyl_sub_range_init(&up->ts_update_r, &fine_ext, flo, fup);

  up->ffine = gkyl_array_new(GKYL_DOUBLE, inp->basis.num_basis, fine_ext.volume);

  // The pure twist-shift updater.
  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir = inp->bc_dir,
    .shift_dir = inp->shift_dir,
    .shear_dir = inp->shear_dir,
    .edge = inp->edge,
    .cdim = inp->cdim,
    .bcdir_ext_update_r = up->ts_update_r,
    .num_ghost = inp->num_ghost,
    .basis = inp->basis,
    .grid = up->ts_grid,
    .shift_func = inp->shift_func,
    .shift_func_ctx = inp->shift_func_ctx,
    .shift_dg = inp->shift_dg,
    .use_gpu = inp->use_gpu,
    .shift_poly_order = inp->shift_poly_order,
  };
  up->ts = gkyl_twistshift_dg_new(&tsinp);

  // Ghost plane the twist-shift fills, on the ts grid.
  if (inp->edge == GKYL_LOWER_EDGE)
    gkyl_range_shorten_from_above(&up->ghost_r, &up->ts_update_r, inp->bc_dir, inp->num_ghost[inp->bc_dir]);
  else
    gkyl_range_shorten_from_below(&up->ghost_r, &up->ts_update_r, inp->bc_dir, inp->num_ghost[inp->bc_dir]);

  // Optional low-pass filter along shear_dir.
  up->filter = NULL;
  up->filt_buff = NULL;
  up->filter_func = bc_twistshift_filter_disabled;
  up->refine = NULL;
  up->coarsen = NULL;
  up->refine_func = bc_twistshift_refine_disabled;
  up->coarsen_func = bc_twistshift_coarsen_disabled;
  if (up->filter_half_width > 0) {
    up->filt_buff = gkyl_array_new(GKYL_DOUBLE, inp->basis.num_basis, up->ffine->size);
    up->filter = gkyl_dg_lowpass_filter_new(inp->shear_dir, up->filter_half_width,
      up->filter_cutoff_wavelength, &inp->basis, &up->ts_grid, &up->ghost_r, inp->use_gpu);
      up->filter_func = bc_twistshift_filter_enabled;
    // Optional refine/coarsen interpolators along shear_dir.
    if (up->upsample_factor > 1) {
      if (inp->edge == GKYL_LOWER_EDGE)
        gkyl_range_shorten_from_above(&up->coarse_ghost_r, &inp->bcdir_ext_update_r, inp->bc_dir, inp->num_ghost[inp->bc_dir]);
      else
        gkyl_range_shorten_from_below(&up->coarse_ghost_r, &inp->bcdir_ext_update_r, inp->bc_dir, inp->num_ghost[inp->bc_dir]);
      up->refine = gkyl_dg_interpolate_new(inp->cdim, &inp->basis, &inp->grid, &up->ts_grid,
        &up->coarse_ghost_r, &up->ghost_r, inp->num_ghost, inp->use_gpu);
      up->coarsen = gkyl_dg_interpolate_new(inp->cdim, &inp->basis, &up->ts_grid, &inp->grid,
        &up->ghost_r, &up->coarse_ghost_r, inp->num_ghost, inp->use_gpu);
      up->refine_func = bc_twistshift_refine_enabled;
      up->coarsen_func = bc_twistshift_coarsen_enabled;
    }
  }

  return up;
}

void
gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar)
{
  // Refine the data to the fine grid.
  up->refine_func(up->refine, fdo, up->ffine);
  // Apply twist-shift.
  gkyl_twistshift_dg_advance(up->ts, up->ffine, up->ffine);
  // Apply the low-pass filter.
  up->filter_func(up->filter, up->ffine, up->filt_buff);
  // Move back the data to the original grid.
  up->coarsen_func(up->coarsen, up->ffine, ftar);
}

void
gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up)
{
  gkyl_twistshift_dg_release(up->ts);
  gkyl_array_release(up->ffine);
  if (up->filter_half_width > 0) {
    gkyl_array_release(up->filt_buff);
    gkyl_dg_lowpass_filter_release(up->filter);
    if (up->upsample_factor > 1) {
      gkyl_dg_interpolate_release(up->refine);
      gkyl_dg_interpolate_release(up->coarsen);
    }
  }
  gkyl_free(up);
}
