/* -*- c++ -*- */
extern "C" {
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_reduce.h>
#include <assert.h>
#include <float.h>
}

#ifdef GKYL_HAVE_CUDA

struct gkyl_dg_array_mask* 
gkyl_dg_array_mask_cu_dev_new(struct gkyl_dg_array_mask *mask_ho)
{
  struct gkyl_dg_array_mask *mask = (struct gkyl_dg_array_mask *) gkyl_malloc(sizeof(*mask));

  mask->type = mask_ho->type;
  mask->val_threshold = mask_ho->val_threshold;
  mask->frac_threshold = mask_ho->frac_threshold;
  mask->use_gpu = mask_ho->use_gpu;
  mask->phase_rng = mask_ho->phase_rng;
  mask->phase_rng_ext = mask_ho->phase_rng_ext;
  mask->conf_rng = mask_ho->conf_rng;
  mask->conf_rng_ext = mask_ho->conf_rng_ext;
  mask->vel_rng = mask_ho->vel_rng;
  mask->mask = NULL;
  mask->local_max_arr = NULL;
  mask->global_max = NULL;

  mask->flags = 0;
  GKYL_SET_CU_ALLOC(mask->flags);
  mask->ref_count = gkyl_ref_count_init(gkyl_dg_array_mask_free);

  // For NONE type, don't allocate mask array
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    // Initialize the device object.
    struct gkyl_dg_array_mask *mask_cu = (struct gkyl_dg_array_mask*) gkyl_cu_malloc(sizeof(*mask_cu));
    gkyl_cu_memcpy(mask_cu, mask, sizeof(struct gkyl_dg_array_mask), GKYL_CU_MEMCPY_H2D);
    mask->on_dev = mask_cu;
  }
  else {
    struct gkyl_array *mask_array = gkyl_array_cu_dev_new(GKYL_DOUBLE, mask_ho->mask->ncomp, mask_ho->mask->size);
    gkyl_array_copy(mask_array, mask_ho->mask);
    mask->mask = mask_array->on_dev;
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL ||
             mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL) {
      mask->local_max_arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, mask->conf_rng_ext.volume);
    }
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD ||
        mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD) {
      mask->global_max = (double*) gkyl_cu_malloc(sizeof(double));
    }

    // Initialize the device object.
    struct gkyl_dg_array_mask *mask_cu = (struct gkyl_dg_array_mask*) gkyl_cu_malloc(sizeof(*mask_cu));
    gkyl_cu_memcpy(mask_cu, mask, sizeof(struct gkyl_dg_array_mask), GKYL_CU_MEMCPY_H2D);
    mask->on_dev = mask_cu;
    
    // The returned object should store host pointer to gkyl_array.
    mask->mask = mask_array;
  }
  
  return mask;
}

// CUDA kernel for LESS_THAN_THRESHOLD masks
__global__ void
gkyl_dg_array_mask_less_than_kernel(struct gkyl_range phase_rng, 
  const struct gkyl_array *arr_to_mask, struct gkyl_array *mask, double f_threshold)
{
  int idx[GKYL_MAX_DIM];
  
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
       tid < phase_rng.volume; tid += blockDim.x*gridDim.x) {
    
    gkyl_sub_range_inv_idx(&phase_rng, tid, idx);
    
    if (gkyl_range_contains_idx(&phase_rng, idx)) {
      long linidx = gkyl_range_idx(&phase_rng, idx);
      const double *distf_c = (const double*) gkyl_array_cfetch(arr_to_mask, linidx);
      double *mask_c = (double*)gkyl_array_fetch(mask, linidx);
      
      double abs_val = fabs(distf_c[0]);
      *mask_c = (abs_val < f_threshold) ? 1.0 : -1.0;
    }
  }
}

// CUDA kernel for GREATER_THAN_THRESHOLD masks
__global__ void
gkyl_dg_array_mask_greater_than_kernel(struct gkyl_range phase_rng,
  const struct gkyl_array *arr_to_mask, struct gkyl_array *mask, double f_threshold)
{
  int idx[GKYL_MAX_DIM];
  
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
       tid < phase_rng.volume; tid += blockDim.x*gridDim.x) {
    
    gkyl_sub_range_inv_idx(&phase_rng, tid, idx);
    
    if (gkyl_range_contains_idx(&phase_rng, idx)) {
      long linidx = gkyl_range_idx(&phase_rng, idx);
      const double *distf_c = (const double*) gkyl_array_cfetch(arr_to_mask, linidx);
      double *mask_c = (double*)gkyl_array_fetch(mask, linidx);
      
      double abs_val = fabs(distf_c[0]);
      *mask_c = (abs_val > f_threshold) ? 1.0 : -1.0;
    }
  }
}

// CUDA kernel to find max value in velocity space for each configuration cell
__global__ void
gkyl_dg_array_mask_find_local_max_kernel(struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  struct gkyl_range phase_rng, const struct gkyl_array *arr_to_mask, struct gkyl_array *local_max_arr)
{
  int conf_idx[GKYL_MAX_DIM];
  
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
       tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    
    gkyl_sub_range_inv_idx(&conf_rng, tid, conf_idx);
    
    if (gkyl_range_contains_idx(&conf_rng, conf_idx)) {
      double local_max = -DBL_MAX;
      
      // Iterate over velocity space for this config cell
      for (unsigned long vel_tid = 0; vel_tid < vel_rng.volume; vel_tid++) {
        int vel_idx[GKYL_MAX_DIM];
        gkyl_sub_range_inv_idx(&vel_rng, vel_tid, vel_idx);
        
        // Build phase space index
        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < conf_rng.ndim; d++) {
          pidx[d] = conf_idx[d];
        }
        for (int d = 0; d < vel_rng.ndim; d++) {
          pidx[conf_rng.ndim + d] = vel_idx[d];
        }
        
        long linidx_phase = gkyl_range_idx(&phase_rng, pidx);
        const double *arr_c = (const double*) gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double abs_val = fabs(arr_c[0]);
        
        if (abs_val > local_max) {
          local_max = abs_val;
        }
      }
      
      long conf_linidx = gkyl_range_idx(&conf_rng, conf_idx);
      double *local_max_c = (double*) gkyl_array_fetch(local_max_arr, conf_linidx);
      local_max_c[0] = local_max;
    }
  }
}

// CUDA kernel to apply spatial fractional mask (LESS_THAN)
__global__ void
gkyl_dg_array_mask_spatial_frac_less_than_kernel(struct gkyl_range conf_rng, 
  struct gkyl_range vel_rng, struct gkyl_range phase_rng, const struct gkyl_array *arr_to_mask,
  struct gkyl_array *mask, const struct gkyl_array *local_max_arr, double frac_threshold)
{
  int conf_idx[GKYL_MAX_DIM];
  
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
       tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    
    gkyl_sub_range_inv_idx(&conf_rng, tid, conf_idx);
    
    if (gkyl_range_contains_idx(&conf_rng, conf_idx)) {
      long conf_linidx = gkyl_range_idx(&conf_rng, conf_idx);
      const double *local_max_c = (const double*) gkyl_array_cfetch(local_max_arr, conf_linidx);
      double local_threshold = frac_threshold * local_max_c[0];
      
      // Apply mask to all velocity cells in this config cell
      for (unsigned long vel_tid = 0; vel_tid < vel_rng.volume; vel_tid++) {
        int vel_idx[GKYL_MAX_DIM];
        gkyl_sub_range_inv_idx(&vel_rng, vel_tid, vel_idx);
        
        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < conf_rng.ndim; d++) {
          pidx[d] = conf_idx[d];
        }
        for (int d = 0; d < vel_rng.ndim; d++) {
          pidx[conf_rng.ndim + d] = vel_idx[d];
        }
        
        long linidx_phase = gkyl_range_idx(&phase_rng, pidx);
        const double *arr_c = (const double*) gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double *mask_c = (double*)gkyl_array_fetch(mask, linidx_phase);
        
        double abs_val = fabs(arr_c[0]);
        *mask_c = (abs_val < local_threshold) ? 1.0 : -1.0;
      }
    }
  }
}

// CUDA kernel to apply spatial fractional mask (GREATER_THAN)
__global__ void
gkyl_dg_array_mask_spatial_frac_greater_than_kernel(struct gkyl_range conf_rng,
  struct gkyl_range vel_rng, struct gkyl_range phase_rng, const struct gkyl_array *arr_to_mask,
  struct gkyl_array *mask, const struct gkyl_array *local_max_arr, double frac_threshold)
{
  int conf_idx[GKYL_MAX_DIM];
  
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
       tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    
    gkyl_sub_range_inv_idx(&conf_rng, tid, conf_idx);
    
    if (gkyl_range_contains_idx(&conf_rng, conf_idx)) {
      long conf_linidx = gkyl_range_idx(&conf_rng, conf_idx);
      const double *local_max_c = (const double*) gkyl_array_cfetch(local_max_arr, conf_linidx);
      double local_threshold = frac_threshold * local_max_c[0];
      
      // Apply mask to all velocity cells in this config cell
      for (unsigned long vel_tid = 0; vel_tid < vel_rng.volume; vel_tid++) {
        int vel_idx[GKYL_MAX_DIM];
        gkyl_sub_range_inv_idx(&vel_rng, vel_tid, vel_idx);
        
        int pidx[GKYL_MAX_DIM];
        for (int d = 0; d < conf_rng.ndim; d++) {
          pidx[d] = conf_idx[d];
        }
        for (int d = 0; d < vel_rng.ndim; d++) {
          pidx[conf_rng.ndim + d] = vel_idx[d];
        }
        
        long linidx_phase = gkyl_range_idx(&phase_rng, pidx);
        const double *arr_c = (const double*) gkyl_array_cfetch(arr_to_mask, linidx_phase);
        double *mask_c = (double*)gkyl_array_fetch(mask, linidx_phase);
        
        double abs_val = fabs(arr_c[0]);
        *mask_c = (abs_val > local_threshold) ? 1.0 : -1.0;
      }
    }
  }
}

// Host function to launch CUDA kernel.
void
gkyl_dg_array_mask_advance_cu(struct gkyl_dg_array_mask *mask, const struct gkyl_array *arr_to_mask)
{
  // Do nothing for NONE type
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    return;
  }
  
  int nblocks = mask->mask->nblocks;
  int nthreads = mask->mask->nthreads;
  
  // Simple threshold masks - launch specialized kernel directly
  if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD) {
    gkyl_dg_array_mask_less_than_kernel<<<nblocks, nthreads>>>(
      mask->phase_rng, arr_to_mask->on_dev, mask->mask->on_dev, mask->val_threshold);
    return;
  }
  
  if (mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD) {
    gkyl_dg_array_mask_greater_than_kernel<<<nblocks, nthreads>>>(
      mask->phase_rng, arr_to_mask->on_dev, mask->mask->on_dev, mask->val_threshold);
    return;
  }
  
  // Global fractional threshold masks - compute max, then launch kernel
  if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD ||
      mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD) {
    gkyl_array_reduce(mask->global_max, arr_to_mask, GKYL_MAX);
    // Copy result from device to host
    double global_max_c0;
    gkyl_cu_memcpy(&global_max_c0, mask->global_max, sizeof(double), GKYL_CU_MEMCPY_D2H);
    double threshold = mask->frac_threshold * global_max_c0;
    
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD) {
      gkyl_dg_array_mask_less_than_kernel<<<nblocks, nthreads>>>(
        mask->phase_rng, arr_to_mask->on_dev, mask->mask->on_dev, threshold);
    } else {
      gkyl_dg_array_mask_greater_than_kernel<<<nblocks, nthreads>>>(
        mask->phase_rng, arr_to_mask->on_dev, mask->mask->on_dev, threshold);
    }
    return;
  }
  
  // Spatial fractional threshold masks - two-phase GPU approach
  if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL ||
      mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL) {
    
    // Phase 1: Find max in velocity space for each configuration cell
    int conf_nblocks = mask->conf_rng.nblocks;
    gkyl_dg_array_mask_find_local_max_kernel<<<conf_nblocks, nthreads>>>(
      mask->conf_rng, mask->vel_rng, mask->phase_rng, arr_to_mask->on_dev, mask->local_max_arr->on_dev);
    
    // Phase 2: Apply mask based on local thresholds
    if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL) {
      gkyl_dg_array_mask_spatial_frac_less_than_kernel<<<conf_nblocks, nthreads>>>(
        mask->conf_rng, mask->vel_rng, mask->phase_rng, arr_to_mask->on_dev,
        mask->mask->on_dev, mask->local_max_arr->on_dev, mask->frac_threshold);
    }
    else {
      gkyl_dg_array_mask_spatial_frac_greater_than_kernel<<<conf_nblocks, nthreads>>>(
        mask->conf_rng, mask->vel_rng, mask->phase_rng, arr_to_mask->on_dev,
        mask->mask->on_dev, mask->local_max_arr->on_dev, mask->frac_threshold);
    }
    
    return;
  }
}

#endif
