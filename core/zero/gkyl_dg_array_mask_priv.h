#pragma once

// Private header for dg_array_mask object, not for direct use in user code.

#include <gkyl_dg_array_mask.h>
#include <gkyl_ref_count.h>

/**
 * Skip cell object definition.
 */
struct gkyl_dg_array_mask {
  enum gkyl_dg_array_mask_types type; // Type of mask operation.
  struct gkyl_array *mask; // Mask array (1.0 is true, -1.0 is false).
  double val_threshold; // Threshold for marking cells as masked.
  struct gkyl_range phase_rng; // Phase-space range.
  bool use_gpu; // Flag indicating GPU usage.
  
  uint32_t flags;
  struct gkyl_dg_array_mask *on_dev; // Pointer to device object.
  
  struct gkyl_ref_count ref_count; // Reference counter.
};

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
void
gkyl_dg_array_mask_free(const struct gkyl_ref_count *ref);

GKYL_CU_DH
static inline bool
gkyl_dg_array_mask_eval_ker(struct gkyl_dg_array_mask *mask, long lidx)
{
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    return false;
  }
  const double *mask_c = (const double *) gkyl_array_cfetch(mask->mask, lidx);
  return *mask_c > 0; // Returns true if the mask is true.
}

GKYL_CU_DH
static inline bool
gkyl_dg_array_mask_eval_idx_ker(struct gkyl_dg_array_mask *mask, const int* idx)
{
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    return false;
  }
  long linidx = gkyl_range_idx(&mask->phase_rng, idx);
  const double *mask_c = (const double *) gkyl_array_cfetch(mask->mask, linidx);
  return *mask_c > 0; // Returns true if the mask is true.
}

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