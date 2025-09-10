#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_deflate_geo gkyl_deflate_geo;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */


gkyl_deflate_geo* gkyl_deflate_geo_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, bool use_gpu);

/**
 * Advance deflate_geo (compute the derived_geo coefficients).
 *
 * @param up deflate_geo updater object.
 * @param crange Config-space range.
 * @param gFld field containing DG rep of the metric coefficients
 * @param jFld output field where jacobian will be placed
 */


void gkyl_deflate_geo_advance(const gkyl_deflate_geo *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_geo_release(gkyl_deflate_geo* up);
