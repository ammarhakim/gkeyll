#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object definition
typedef struct gkyl_dg_array_mask gkyl_dg_array_mask;

enum gkyl_dg_array_mask_types {
  GKYL_DG_ARRAY_MASK_NONE = 0, // No mask applied.
  GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD, // Mask applied based on the 0th component of the array.
  GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD, // Mask applied based on the 0th component of the array.
  GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD, // Applies the mask based on the 0th component of the input array, which is spatially dependent.
  GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD, // Applies the mask based on the 0th component of the input array, which is spatially dependent.
  GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL, // Applies the mask based on the 0th component of the input array, which is spatially dependent.
  GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL, // Applies the mask based on the 0th component of the input array, which is spatially dependent.
};

// Input structure for creating a mask object.
struct gkyl_dg_array_mask_inp {
  enum gkyl_dg_array_mask_types type; // Type of mask.
  double val_threshold; // Threshold for marking cells as masked.
  double frac_threshold;
  struct gkyl_range phase_rng; // Phase-space range.
  struct gkyl_range config_rng; // Configuration-space range.
  struct gkyl_range vel_rng; // Velocity-space range.
  bool use_gpu; // Flag indicating GPU usage.
};

/**
 * Create a new mask object.
 *
 * @param mask_inp Input parameters for mask.
 * @return New mask object.
 */
struct gkyl_dg_array_mask*
gkyl_dg_array_mask_new(struct gkyl_dg_array_mask_inp mask_inp);

/**
 * Update the mask array.
 *
 * @param mask Mask updater.
 * @param arr_to_mask Array which the mask is applied on top of.
 */
void
gkyl_dg_array_mask_advance(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask);

/**
 * Evaluate if the conditional mask is true at a given cell.
 *
 * @param mask Mask object.
 * @param idx Linear index of the cell to evaluate.
 * @return True if the mask is true at the cell, false otherwise.
 */
bool gkyl_dg_array_mask_eval(struct gkyl_dg_array_mask *mask, long lidx);


/** 
 * Evaluate if the conditional mask is true at a given multi-dimensional index.
 *
 * @param mask Mask object.
 * @param idx Multi-dimensional index array.
 * @return True if the mask is true at the index, false otherwise.
 */
bool gkyl_dg_array_mask_eval_idx(struct gkyl_dg_array_mask *mask, const int* idx);

/**
 * mask = mask * arr_to_multiply
 * Scale the dg_array_mask by an array, modifying the mask.
 *`
 * @param mask Mask object.
 * @param arr_to_multiply Array to be multiplied by the mask.
 */
void
gkyl_dg_array_mask_scale_by_cell(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_multiply);


/**
 * Acquire a reference to the mask object.
 * 
 * @param mask Mask object to acquire.
 */

struct gkyl_dg_array_mask*
gkyl_dg_array_mask_acquire(struct gkyl_dg_array_mask *mask);


/**
 * Get the underlying mask array.
 *
 * @param mask Mask object.
 * @return Pointer to the mask array.
 */
const struct gkyl_array*
gkyl_dg_array_mask_get_mask(const struct gkyl_dg_array_mask *mask);

/**
 * Get device pointer to mask object.
 *
 * @param mask Mask object.
 */
struct gkyl_dg_array_mask*
gkyl_dg_array_mask_get_dev_ptr(struct gkyl_dg_array_mask *mask);


/**
 * Release memory associated with mask object.
 *
 * @param mask Mask object to release.
 */
void
gkyl_dg_array_mask_release(struct gkyl_dg_array_mask *mask);
