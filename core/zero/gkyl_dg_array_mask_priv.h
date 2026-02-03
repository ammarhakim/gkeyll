#pragma once

// Private header for dg_array_mask object, not for direct use in user code.

#include <gkyl_dg_array_mask.h>
#include <gkyl_ref_count.h>

/**
  * Skip cell object definition.
  */
struct gkyl_dg_array_mask {
  enum gkyl_dg_array_mask_types type;     // Type of mask operation.
  bool default_value;                     // Default value for mask (true/false) if no masking is applied. Defaults to false (-1.0).
  struct gkyl_array *mask_arr;            // Mask array (1.0 is true, -1.0 is false).
  double threshold;                       // Threshold for marking cells as masked. Scaled absolute value for *_THRESHOLD types, fraction for *_FRAC_THRESHOLD types.
  double frac_threshold;                  // Fractional threshold, which is 0-1 for global fraction mask.
  const struct gkyl_range *mask_rng;      // Pointer to range over which mask is applied (phase_rng for kinetic, conf_rng for fluid).
  const struct gkyl_range *mask_rng_ext;  // Pointer to extended range for mask allocation.
  const struct gkyl_range *conf_rng;      // Configuration-space range.
  const struct gkyl_range *conf_rng_ext;  // Extended configuration-space range.
  const struct gkyl_range *vel_rng;       // Velocity-space range.
  bool use_gpu;                           // Flag indicating GPU usage.

  struct gkyl_array *local_max_arr; // Pre-allocated config-space array for spatial fractional masks.
  double *global_max;               // Pre-allocated array for global reduction results (fractional threshold masks).

  // Function pointer for advance method (CPU), set at init time based on mask type.
  void (*advance_func)(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask);
  void (*advance_threshold_func)(struct gkyl_dg_array_mask *mask, const double global_max);

  // Function pointer for advance method (GPU), set at init time based on mask type.
  void (*advance_func_cu)(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask);
  
  // Function that evaluates the mask.
  bool (*eval_idx_func)(struct gkyl_dg_array_mask *mask, const int *idx);

  // Function pointer for scale_by_cell method, set at init time based on mask type.
  void (*scale_by_cell_func)(struct gkyl_dg_array_mask *mask,
    const struct gkyl_array *arr_to_multiply);

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
void gkyl_dg_array_mask_free(const struct gkyl_ref_count *ref);

GKYL_CU_DH
static bool
eval_idx_ker_disabled(struct gkyl_dg_array_mask *mask, const int *idx)
{
  return mask->default_value;
}

GKYL_CU_DH
static bool
eval_idx_ker_enabled(struct gkyl_dg_array_mask *mask, const int *idx)
{
  long linidx = gkyl_range_idx(mask->mask_rng, idx);
  const double *mask_c = (const double *)gkyl_array_cfetch(mask->mask_arr, linidx);
  return *mask_c > 0; // Returns true if the mask is true.
}

#ifdef GKYL_HAVE_CUDA

/**
  * Create a new dg_array_mask object on CUDA device.
  *
  * @param mask_ho Host-side dg_array_mask object.
  * @return New dg_array_mask object on device.
  */
struct gkyl_dg_array_mask *
gkyl_dg_array_mask_cu_dev_new(struct gkyl_dg_array_mask *mask_ho);

/**
  * CUDA device function to update dg_array_mask on GPU.
  *
  * @param mask dg_array_mask object.
  * @param arr_to_mask Array to mask.
  */
void gkyl_dg_array_mask_advance_cu(struct gkyl_dg_array_mask *mask,
  const struct gkyl_array *arr_to_mask);

#endif
