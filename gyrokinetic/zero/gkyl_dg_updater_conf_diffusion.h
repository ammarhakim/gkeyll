#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

typedef struct gkyl_dg_updater_conf_diffusion gkyl_dg_updater_conf_diffusion;

struct gkyl_dg_updater_conf_diffusion_tm {
  double diffusion_tm;
};

/**
 * Create a conservative configuration-space diffusion updater. The
 * coefficient must be represented in the same DG basis as the input field.
 * The implementation supports 1x-3x serendipity p1 updates on the host.
 * Set zero_flux_flags[d] and zero_flux_flags[d+cdim] for homogeneous
 * diffusive-flux conditions at the lower and upper boundaries in direction d;
 * leave them zero on periodic directions.
 */
struct gkyl_dg_updater_conf_diffusion* gkyl_dg_updater_conf_diffusion_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const struct gkyl_range *conf_range, const struct gkyl_array *diffusion_tensor,
  const struct gkyl_array *jacobgeo_inv,
  const int zero_flux_flags[2*GKYL_MAX_CDIM], bool use_gpu);

/** Increment rhs with div(D grad f) and cflrate with its stability rate. */
void gkyl_dg_updater_conf_diffusion_advance(
  struct gkyl_dg_updater_conf_diffusion *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn,
  struct gkyl_array *GKYL_RESTRICT cflrate,
  struct gkyl_array *GKYL_RESTRICT rhs);

struct gkyl_dg_updater_conf_diffusion_tm
gkyl_dg_updater_conf_diffusion_get_tm(
  const struct gkyl_dg_updater_conf_diffusion *up);

void gkyl_dg_updater_conf_diffusion_release(
  struct gkyl_dg_updater_conf_diffusion *up);
