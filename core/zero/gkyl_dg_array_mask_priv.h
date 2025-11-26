#pragma once

// Private header for dg_array_mask object, not for direct use in user code.

#include <gkyl_dg_array_mask.h>
#include <gkyl_ref_count.h>

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
void
gkyl_dg_array_mask_free(const struct gkyl_ref_count *ref);

#ifdef GKYL_HAVE_CUDA

/**
 * Create a new dg_array_mask object on CUDA device.
 *
 * @param mask_ho Host-side dg_array_mask object.
 * @return New dg_array_mask object on device.
 */
struct gkyl_dg_array_mask*
gkyl_dg_array_mask_cu_dev_new(struct gkyl_dg_array_mask *mask_ho);

/**
 * CUDA device function to update dg_array_mask on GPU.
 *
 * @param mask dg_array_mask object.
 * @param arr_to_mask Array to mask.
 */
void
gkyl_dg_array_mask_advance_cu(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask);

#endif
