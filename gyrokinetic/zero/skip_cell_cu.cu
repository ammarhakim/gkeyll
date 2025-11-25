/* -*- c++ -*- */

extern "C" {
#include <gkyl_skip_cell.h>
#include <gkyl_skip_cell_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <assert.h>
}

#ifdef GKYL_HAVE_CUDA

// CUDA kernel to update skip cell mask.
__global__ void
gkyl_skip_cell_advance_kernel(struct gkyl_range phase_rng,
  const struct gkyl_array *distf, struct gkyl_array *mask,
  double f_threshold)
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
      const double *distf_c = (const double*) gkyl_array_cfetch(distf, linidx);
      
      // Fetch the boolean mask value for this cell.
      bool *skip_c = (bool*)gkyl_array_fetch(mask, linidx);
      
      // Mark cell as skippable if distribution function is below threshold.
      if (fabs(distf_c[0]) < f_threshold) {
        *skip_c = true;
      } else {
        *skip_c = false;
      }
    }
  }
}

// Host function to launch CUDA kernel.
void
gkyl_skip_cell_advance_cu(struct gkyl_skip_cell *skip_cell, const struct gkyl_array *distf)
{
  int nblocks = skip_cell->mask->nblocks;
  int nthreads = skip_cell->mask->nthreads;
  
  gkyl_skip_cell_advance_kernel<<<nblocks, nthreads>>>(
    skip_cell->phase_rng, distf->on_dev, skip_cell->mask->on_dev,
    skip_cell->f_threshold);
}

#endif
