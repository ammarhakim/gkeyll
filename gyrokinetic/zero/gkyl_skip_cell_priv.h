#pragma once

// Private header for skip_cell object, not for direct use in user code.

#include <gkyl_skip_cell.h>
#include <gkyl_ref_count.h>


#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to update skip cell mask on GPU.
 *
 * @param skip_cell Skip cell object.
 * @param distf Distribution function array.
 */
void
gkyl_skip_cell_advance_cu(struct gkyl_skip_cell *skip_cell, const struct gkyl_array *distf);

#endif
