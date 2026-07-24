#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_dg_lowpass_filter_priv.h>
#include <gkyl_alloc.h>

struct gkyl_dg_lowpass_filter*
gkyl_dg_lowpass_filter_new(int dir, int half_width, double cutoff_wavelength,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_dg_lowpass_filter *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  assert(!up->use_gpu); // GPU implementation pending.
  up->ndim = basis->ndim;
  up->dir = dir;
  up->half_width = half_width;
  up->num_basis = basis->num_basis;
  up->grid = *grid;
  up->range = *range;

  // Perform some basic checks:
  assert(grid->ndim == range->ndim);
  assert(0 <= dir && dir < grid->ndim);
  assert(half_width > 0);

  // Normalized cutoff frequency in cycles per cell.
  double fc = grid->dx[dir]/cutoff_wavelength;
  assert(0.0 < fc && fc <= 0.5);

  up->weights = gkyl_malloc((2*half_width+1)*sizeof(double));
  dg_lpf_calc_weights(half_width, fc, up->weights);

  return up;
}

void
gkyl_dg_lowpass_filter_advance(gkyl_dg_lowpass_filter *up,
  struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  assert(fdo != ftar); // Stencil operation, can't be done in-place.
  assert(fdo->ncomp == ftar->ncomp);
  assert(fdo->size == ftar->size);
  assert(fdo->ncomp == up->num_basis);

  int dir = up->dir;
  int M = up->half_width;
  int num_basis = up->num_basis;
  int idx_do[GKYL_MAX_DIM];

  // Loop over the target range.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->range);
  while (gkyl_range_iter_next(&iter)) {

    long linidx_tar = gkyl_range_idx(&up->range, iter.idx);
    double *ftar_c = gkyl_array_fetch(ftar, linidx_tar);

    // Truncate the stencil at the boundaries of the range.
    int koff_lo = GKYL_MAX2(-M, up->range.lower[dir]-iter.idx[dir]);
    int koff_up = GKYL_MIN2( M, up->range.upper[dir]-iter.idx[dir]);

    gkyl_copy_int_arr(up->range.ndim, iter.idx, idx_do);

    for (int c=0; c<num_basis; c++)
      ftar_c[c] = 0.0;

    // Loop over the donor cells contributing to this target cell.
    double wsum = 0.0;
    for (int k=koff_lo; k<koff_up+1; k++) {
      idx_do[dir] = iter.idx[dir]+k;

      long linidx_do = gkyl_range_idx(&up->range, idx_do);
      const double *fdo_c = gkyl_array_cfetch(fdo, linidx_do);

      double w = up->weights[k+M];
      for (int c=0; c<num_basis; c++)
        ftar_c[c] += w*fdo_c[c];
      wsum += w;
    }

    // Renormalize truncated stencils to preserve constants at boundaries.
    for (int c=0; c<num_basis; c++)
      ftar_c[c] /= wsum;
  }
}

void
gkyl_dg_lowpass_filter_release(gkyl_dg_lowpass_filter *up)
{
  // Release memory associated with this updater.
  gkyl_free(up->weights);
  gkyl_free(up);
}
