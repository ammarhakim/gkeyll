#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>

void
gkyl_dg_array_mask_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_array_mask *mask = container_of(ref, struct gkyl_dg_array_mask, ref_count);
  
  if (mask->mask)
    gkyl_array_release(mask->mask);
  
  if (mask->local_max_arr)
    gkyl_array_release(mask->local_max_arr);
  
  if (mask->global_max)
    gkyl_free(mask->global_max);
  
  if (GKYL_IS_CU_ALLOC(mask->flags))
    gkyl_cu_free(mask->on_dev);
  
  gkyl_free(mask);
}

struct gkyl_dg_array_mask*
gkyl_dg_array_mask_new(struct gkyl_dg_array_mask_inp mask_inp)
{
  struct gkyl_dg_array_mask *mask = gkyl_malloc(sizeof(*mask));

  mask->type = mask_inp.type;
  mask->default_value = mask_inp.default_value;
  mask->use_gpu = mask_inp.use_gpu;
  mask->val_threshold = 0.0;
  mask->mask = NULL;
  mask->local_max_arr = NULL;
  mask->global_max = NULL;
  mask->flags = 0;
  GKYL_CLEAR_CU_ALLOC(mask->flags);
  mask->ref_count = gkyl_ref_count_init(gkyl_dg_array_mask_free);
  mask->on_dev = mask; // CPU mask points to itself

  if (mask->type != GKYL_DG_ARRAY_MASK_NONE) {
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD ||
      mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD) {
        mask->val_threshold = mask_inp.val_threshold * pow(sqrt(2.0), mask_inp.phase_rng.ndim);
    } else if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD ||
               mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD) {
      // Threshold will be set during advance based on global max value.
      mask->frac_threshold = mask_inp.frac_threshold;
      // Pre-allocate array for global reduction
      mask->global_max = (double*) gkyl_malloc(sizeof(double) * mask_inp.phase_rng.ndim);
    } else if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL ||
               mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL) {
      // Threshold will be set during advance based on global max value.
      mask->frac_threshold = mask_inp.frac_threshold;
    }
    
    // Initialize the mask array on host.
    mask->phase_rng = mask_inp.phase_rng;
    mask->conf_rng = mask_inp.config_rng;
    mask->vel_rng = mask_inp.vel_rng;
    mask->mask = gkyl_array_new(GKYL_DOUBLE, 1, mask_inp.phase_rng.volume);
    gkyl_array_clear(mask->mask, -1.0); // Initialize all cells to false for safety.
  }

  struct gkyl_dg_array_mask *mask_out = mask;
#ifdef GKYL_HAVE_CUDA
  if (mask->use_gpu) {
    mask_out = gkyl_dg_array_mask_cu_dev_new(mask);
    gkyl_dg_array_mask_release(mask);
  }
#endif

  return mask_out;
}

void
gkyl_dg_array_mask_advance(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask)
{
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    return;
  }

  assert(arr_to_mask->type == GKYL_DOUBLE);
  assert(mask->mask->size == arr_to_mask->size);

#ifdef GKYL_HAVE_CUDA
  if (mask->use_gpu) {
    gkyl_dg_array_mask_advance_cu(mask, arr_to_mask);
    return;
  }
#endif

  // Iterate over phase space and update mask.
  // Could be improved in the future by using nodal values or quadrature points
  
  if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD ||
      mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD) {

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &mask->phase_rng);

    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&mask->phase_rng, iter.idx);
      const double *arr_to_mask_c = gkyl_array_cfetch(arr_to_mask, linidx);
      double *mask_c = gkyl_array_fetch(mask->mask, linidx);
      
      double abs_val = fabs(arr_to_mask_c[0]);
      if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD) {
        *mask_c = (abs_val < mask->val_threshold) ? 1.0 : -1.0;
      }
      else { // GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD
        *mask_c = (abs_val > mask->val_threshold) ? 1.0 : -1.0;
      }
    }
  } else if (mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD ||
             mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD) {

    // First find the global max value of the 0th component over configuration space
    gkyl_array_reduce(mask->global_max, arr_to_mask, GKYL_MAX);
    double global_max_c0 = mask->global_max[0];
    mask->val_threshold = mask->frac_threshold * global_max_c0;

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &mask->phase_rng);

    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&mask->phase_rng, iter.idx);
      const double *arr_to_mask_c = gkyl_array_cfetch(arr_to_mask, linidx);
      double *mask_c = gkyl_array_fetch(mask->mask, linidx);
      
      double abs_val = fabs(arr_to_mask_c[0]);
      if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD) {
        *mask_c = (abs_val < mask->val_threshold) ? 1.0 : -1.0;
      }
      else { // GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD
        *mask_c = (abs_val > mask->val_threshold) ? 1.0 : -1.0;
      }
    }
  } else if (mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL ||
             mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL) {

    struct gkyl_range_iter iter_vel, iter_conf;
    gkyl_range_iter_init(&iter_conf, &mask->conf_rng);

    // For each configuration space cell, find the max over velocity space
    while (gkyl_range_iter_next(&iter_conf)) {
      // Find max in velocity space for this configuration space cell
      double local_max = -DBL_MAX;
      
      gkyl_range_iter_init(&iter_vel, &mask->vel_rng); // Reset before FIRST loop
      while (gkyl_range_iter_next(&iter_vel)) {
        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < mask->conf_rng.ndim; d++) {
          pidx[d] = iter_conf.idx[d];
        }
        for (int d = 0; d < mask->vel_rng.ndim; d++) {
          pidx[mask->conf_rng.ndim + d] = iter_vel.idx[d];
        }
        long linidx_phase = gkyl_range_idx(&mask->phase_rng, pidx);
        const double *arr_to_mask_c = gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double abs_val = fabs(arr_to_mask_c[0]);
        if (abs_val > local_max) {
          local_max = abs_val;
        }
      }
      
      // Now compute the mask in velocity space for this configuration space cell
      mask->val_threshold = mask->frac_threshold * local_max;
      
      gkyl_range_iter_init(&iter_vel, &mask->vel_rng); // Reset before SECOND loop
      while (gkyl_range_iter_next(&iter_vel)) {
        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < mask->conf_rng.ndim; d++) {
          pidx[d] = iter_conf.idx[d];
        }
        for (int d = 0; d < mask->vel_rng.ndim; d++) {
          pidx[mask->conf_rng.ndim + d] = iter_vel.idx[d];
        }
        long linidx_phase = gkyl_range_idx(&mask->phase_rng, pidx);
        const double *arr_to_mask_c = gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double abs_val = fabs(arr_to_mask_c[0]);
        
        double *mask_c = gkyl_array_fetch(mask->mask, linidx_phase);
        if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL) {
          *mask_c = (abs_val < mask->val_threshold) ? 1.0 : -1.0;
        }
        else if (mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL) {
          *mask_c = (abs_val > mask->val_threshold) ? 1.0 : -1.0;
        }
      }
    }
  }
}

void
gkyl_dg_array_mask_scale_by_cell(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_multiply)
{
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    return;
  }
  gkyl_array_scale_by_cell(mask->mask, arr_to_multiply);
}

struct gkyl_dg_array_mask*
gkyl_dg_array_mask_acquire(struct gkyl_dg_array_mask *mask)
{
  gkyl_ref_count_inc(&mask->ref_count);
  return (struct gkyl_dg_array_mask*) mask;
}

bool
gkyl_dg_array_mask_eval(struct gkyl_dg_array_mask *mask, long lidx)
{
  return gkyl_dg_array_mask_eval_ker(mask, lidx);
}

bool
gkyl_dg_array_mask_eval_idx(struct gkyl_dg_array_mask *mask, const int* idx)
{
  return gkyl_dg_array_mask_eval_idx_ker(mask, idx);
}

struct gkyl_dg_array_mask*
gkyl_dg_array_mask_get_dev_ptr(struct gkyl_dg_array_mask *mask)
{
  return mask->on_dev;
} 

const struct gkyl_array*
gkyl_dg_array_mask_get_mask(const struct gkyl_dg_array_mask *mask)
{
  return mask->mask;
}

void
gkyl_dg_array_mask_release(struct gkyl_dg_array_mask *mask)
{
  if (mask)
    gkyl_ref_count_dec(&mask->ref_count);
}
