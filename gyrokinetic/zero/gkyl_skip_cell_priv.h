#pragma once

// Private header for skip_cell object, not for direct use in user code.

#include <gkyl_skip_cell.h>
#include <gkyl_ref_count.h>

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
void
gkyl_skip_cell_free(const struct gkyl_ref_count *ref);

#ifdef GKYL_HAVE_CUDA

/**
 * Create a new skip cell object on CUDA device.
 *
 * @param skip_cell_ho Host-side skip cell object.
 * @return New skip cell object on device.
 */
struct gkyl_skip_cell*
gkyl_skip_cell_cu_dev_new(struct gkyl_skip_cell *skip_cell_ho);

/**
 * CUDA device function to update skip cell mask on GPU.
 *
 * @param skip_cell Skip cell object.
 * @param distf Distribution function array.
 */
void
gkyl_skip_cell_advance_cu(struct gkyl_skip_cell *skip_cell, const struct gkyl_array *distf);

#endif
