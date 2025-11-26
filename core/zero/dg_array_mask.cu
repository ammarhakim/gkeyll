/* -*- c++ -*- */
extern "C" {
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_alloc_flags_priv.h>
#include <assert.h>
}

#ifdef GKYL_HAVE_CUDA

struct gkyl_dg_array_mask* 
gkyl_dg_array_mask_cu_dev_new(struct gkyl_dg_array_mask *mask_ho)
{
  struct gkyl_dg_array_mask *mask = (struct gkyl_dg_array_mask *) gkyl_malloc(sizeof(*mask));

  mask->type = mask_ho->type;
  mask->val_threshold = mask_ho->val_threshold;
  mask->use_gpu = mask_ho->use_gpu;
  mask->phase_rng = mask_ho->phase_rng;

  mask->flags = 0;
  GKYL_SET_CU_ALLOC(mask->flags);
  mask->ref_count = gkyl_ref_count_init(gkyl_dg_array_mask_free);

  // For NONE type, don't allocate mask array
  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    mask->mask = NULL;
    
    // Initialize the device object.
    struct gkyl_dg_array_mask *mask_cu = (struct gkyl_dg_array_mask*) gkyl_cu_malloc(sizeof(*mask_cu));
    gkyl_cu_memcpy(mask_cu, mask, sizeof(struct gkyl_dg_array_mask), GKYL_CU_MEMCPY_H2D);
    mask->on_dev = mask_cu;
  }
  else {
    struct gkyl_array *mask_array = gkyl_array_cu_dev_new(GKYL_DOUBLE, mask_ho->mask->ncomp, mask_ho->mask->size);
    gkyl_array_copy(mask_array, mask_ho->mask);
    mask->mask = mask_array->on_dev;
    
    // Initialize the device object.
    struct gkyl_dg_array_mask *mask_cu = (struct gkyl_dg_array_mask*) gkyl_cu_malloc(sizeof(*mask_cu));
    gkyl_cu_memcpy(mask_cu, mask, sizeof(struct gkyl_dg_array_mask), GKYL_CU_MEMCPY_H2D);
    mask->on_dev = mask_cu;
    
    // The returned object should store host pointer to gkyl_array.
    mask->mask = mask_array;
  }
  
  return mask;
}

// CUDA kernel to update skip cell mask.
__global__ void
gkyl_dg_array_mask_advance_kernel(int mask_type,
  struct gkyl_range phase_rng, const struct gkyl_array *arr_to_mask, 
  struct gkyl_array *mask, double f_threshold)
{
  int idx[GKYL_MAX_DIM];
  
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
       tid < phase_rng.volume; tid += blockDim.x*gridDim.x) {
    
    // Convert linear index to multi-dimensional index.
    gkyl_sub_range_inv_idx(&phase_rng, tid, idx);
    
    // Check if we're inside the range.
    if (gkyl_range_contains_idx(&phase_rng, idx)) {
      long linidx = gkyl_range_idx(&phase_rng, idx);
      
      // Fetch the distribution function at this cell.
      const double *distf_c = (const double*) gkyl_array_cfetch(arr_to_mask, linidx);
      
      // Fetch the mask value for this cell.
      double *mask_c = (double*)gkyl_array_fetch(mask, linidx);
      
      // Evaluate mask condition based on type.
      double abs_val = fabs(distf_c[0]);
      
      if (mask_type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD) {
        *mask_c = (abs_val < f_threshold) ? 1.0 : -1.0;
      }
      else if (mask_type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD) {
        *mask_c = (abs_val > f_threshold) ? 1.0 : -1.0;
      }
      else {
        *mask_c = -1.0;
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
  
  gkyl_dg_array_mask_advance_kernel<<<nblocks, nthreads>>>(
    mask->type, mask->phase_rng, arr_to_mask->on_dev, mask->mask->on_dev,
    mask->val_threshold);
}

#endif
