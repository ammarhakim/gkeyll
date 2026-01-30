#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>

// Static helper: Apply less-than threshold mask over phase range
static void
apply_mask_less_than(struct gkyl_array *mask_arr, const struct gkyl_array *arr_to_mask,
  const struct gkyl_range *phase_rng, double threshold)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_rng);

  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(phase_rng, iter.idx);
    const double *arr_c = gkyl_array_cfetch(arr_to_mask, linidx);
    double *mask_c = gkyl_array_fetch(mask_arr, linidx);
    
    double abs_val = fabs(arr_c[0]);
    *mask_c = (abs_val < threshold) ? 1.0 : -1.0;
  }
}

// Static helper: Apply greater-than threshold mask over phase range
static void
apply_mask_greater_than(struct gkyl_array *mask_arr, const struct gkyl_array *arr_to_mask,
  const struct gkyl_range *phase_rng, double threshold)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_rng);

  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(phase_rng, iter.idx);
    const double *arr_c = gkyl_array_cfetch(arr_to_mask, linidx);
    double *mask_c = gkyl_array_fetch(mask_arr, linidx);
    
    double abs_val = fabs(arr_c[0]);
    *mask_c = (abs_val > threshold) ? 1.0 : -1.0;
  }
}

// Static helper: Find max value in velocity space for a given config cell
static double
find_local_max_in_vel_space(const struct gkyl_array *arr_to_mask,
  const struct gkyl_range *conf_rng, const struct gkyl_range *vel_rng,
  const struct gkyl_range *phase_rng, const int *conf_idx)
{
  double local_max = -DBL_MAX;
  
  struct gkyl_range_iter iter_vel;
  gkyl_range_iter_init(&iter_vel, vel_rng);
  
  while (gkyl_range_iter_next(&iter_vel)) {
    int pidx[GKYL_MAX_DIM];
    for (int d = 0; d < conf_rng->ndim; d++) {
      pidx[d] = conf_idx[d];
    }
    for (int d = 0; d < vel_rng->ndim; d++) {
      pidx[conf_rng->ndim + d] = iter_vel.idx[d];
    }
    
    long linidx_phase = gkyl_range_idx(phase_rng, pidx);
    const double *arr_c = gkyl_array_cfetch(arr_to_mask, linidx_phase);
    double abs_val = fabs(arr_c[0]);
    
    if (abs_val > local_max) {
      local_max = abs_val;
    }
  }
  
  return local_max;
}

// Static helper: Apply spatial fractional mask (less-than) for a given config cell
static void
apply_spatial_mask_less_than(struct gkyl_array *mask_arr, const struct gkyl_array *arr_to_mask,
  const struct gkyl_range *conf_rng, const struct gkyl_range *vel_rng,
  const struct gkyl_range *phase_rng, const int *conf_idx, double threshold)
{
  struct gkyl_range_iter iter_vel;
  gkyl_range_iter_init(&iter_vel, vel_rng);
  
  while (gkyl_range_iter_next(&iter_vel)) {
    int pidx[GKYL_MAX_DIM];
    for (int d = 0; d < conf_rng->ndim; d++) {
      pidx[d] = conf_idx[d];
    }
    for (int d = 0; d < vel_rng->ndim; d++) {
      pidx[conf_rng->ndim + d] = iter_vel.idx[d];
    }
    
    long linidx_phase = gkyl_range_idx(phase_rng, pidx);
    const double *arr_c = gkyl_array_cfetch(arr_to_mask, linidx_phase);
    double *mask_c = gkyl_array_fetch(mask_arr, linidx_phase);
    
    double abs_val = fabs(arr_c[0]);
    *mask_c = (abs_val < threshold) ? 1.0 : -1.0;
  }
}

// Static helper: Apply spatial fractional mask (greater-than) for a given config cell
static void
apply_spatial_mask_greater_than(struct gkyl_array *mask_arr, const struct gkyl_array *arr_to_mask,
  const struct gkyl_range *conf_rng, const struct gkyl_range *vel_rng,
  const struct gkyl_range *phase_rng, const int *conf_idx, double threshold)
{
  struct gkyl_range_iter iter_vel;
  gkyl_range_iter_init(&iter_vel, vel_rng);
  
  while (gkyl_range_iter_next(&iter_vel)) {
    int pidx[GKYL_MAX_DIM];
    for (int d = 0; d < conf_rng->ndim; d++) {
      pidx[d] = conf_idx[d];
    }
    for (int d = 0; d < vel_rng->ndim; d++) {
      pidx[conf_rng->ndim + d] = iter_vel.idx[d];
    }
    
    long linidx_phase = gkyl_range_idx(phase_rng, pidx);
    const double *arr_c = gkyl_array_cfetch(arr_to_mask, linidx_phase);
    double *mask_c = gkyl_array_fetch(mask_arr, linidx_phase);
    
    double abs_val = fabs(arr_c[0]);
    *mask_c = (abs_val > threshold) ? 1.0 : -1.0;
  }
}

void
gkyl_dg_array_mask_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_array_mask *mask = container_of(ref, struct gkyl_dg_array_mask, ref_count);
  
  if (mask->mask)
    gkyl_array_release(mask->mask);
  
  if (mask->local_max_arr)
    gkyl_array_release(mask->local_max_arr);
  
  if (mask->global_max) {
    if (GKYL_IS_CU_ALLOC(mask->flags))
      gkyl_cu_free(mask->global_max);
    else
      gkyl_free(mask->global_max);
  }
  
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
  mask->threshold = 0.0;
  mask->mask = NULL;
  mask->local_max_arr = NULL;
  mask->global_max = NULL;
  mask->flags = 0;
  GKYL_CLEAR_CU_ALLOC(mask->flags);
  mask->ref_count = gkyl_ref_count_init(gkyl_dg_array_mask_free);
  mask->on_dev = mask; // CPU mask points to itself

  if (mask->type != GKYL_DG_ARRAY_MASK_NONE) {
    // Store all ranges from input as pointers.
    mask->phase_rng = mask_inp.phase_rng;
    mask->phase_rng_ext = mask_inp.phase_rng_ext;
    mask->conf_rng = mask_inp.conf_rng;
    mask->conf_rng_ext = mask_inp.conf_rng_ext;
    mask->vel_rng = mask_inp.vel_rng;
    
    // Determine mask range: use phase_rng for kinetic species, conf_rng for fluid species.
    // A NULL or zero-volume phase_rng indicates a fluid species.
    bool is_kinetic = (mask_inp.phase_rng && mask_inp.phase_rng->volume > 0);
    if (is_kinetic) {
      mask->mask_rng = mask_inp.phase_rng;
      mask->mask_rng_ext = mask_inp.phase_rng_ext;
    } else {
      mask->mask_rng = mask_inp.conf_rng;
      mask->mask_rng_ext = mask_inp.conf_rng_ext;
    }
    
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS ||
      mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER) {
        mask->threshold = mask_inp.threshold * pow(sqrt(2.0), mask->mask_rng->ndim);
    } else if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC ||
               mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC) {
      // Threshold will be set during advance based on global max value.
      mask->threshold = mask_inp.threshold;
      // Pre-allocate array for global reduction
      mask->global_max = (double*) gkyl_malloc(sizeof(double) * mask->mask_rng_ext->ndim);
    } else if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF ||
               mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF) {
      // Threshold will be set during advance based on global max value.
      mask->threshold = mask_inp.threshold;
    }
    
    // Initialize the mask array on host.
    mask->mask = gkyl_array_new(GKYL_DOUBLE, 1, mask->mask_rng_ext->volume);
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

  // Apply mask based on type using static helper functions
  if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS) {
    apply_mask_less_than(mask->mask, arr_to_mask, mask->mask_rng, mask->threshold);
  }
  else if (mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER) {
    apply_mask_greater_than(mask->mask, arr_to_mask, mask->mask_rng, mask->threshold);
  }
  else if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC ||
           mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC) {
    // Find global max and compute threshold
    gkyl_array_reduce(mask->global_max, arr_to_mask, GKYL_MAX);
    double frac_threshold = mask->threshold * mask->global_max[0];
    
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC) {
      apply_mask_less_than(mask->mask, arr_to_mask, mask->mask_rng, frac_threshold);
    } else {
      apply_mask_greater_than(mask->mask, arr_to_mask, mask->mask_rng, frac_threshold);
    }
  }
  else if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF ||
           mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF) {
    // For each config cell, find local max in velocity space and apply mask.
    // Note: _FRAC_CONF types are only valid for kinetic species (phase_rng must have volume).
    struct gkyl_range_iter iter_conf;
    gkyl_range_iter_init(&iter_conf, mask->conf_rng);

    while (gkyl_range_iter_next(&iter_conf)) {
      // Find max in velocity space for this configuration cell
      double local_max = find_local_max_in_vel_space(arr_to_mask,
        mask->conf_rng, mask->vel_rng, mask->mask_rng, iter_conf.idx);
      
      // Compute threshold and apply mask for this config cell
      double frac_threshold = mask->threshold * local_max;
      
      if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF) {
        apply_spatial_mask_less_than(mask->mask, arr_to_mask,
          mask->conf_rng, mask->vel_rng, mask->mask_rng, iter_conf.idx, frac_threshold);
      } else {
        apply_spatial_mask_greater_than(mask->mask, arr_to_mask,
          mask->conf_rng, mask->vel_rng, mask->mask_rng, iter_conf.idx, frac_threshold);
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
