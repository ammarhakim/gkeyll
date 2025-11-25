/* -*- c++ -*- */
extern "C" {
#include <gkyl_velocity_map.h>
#include <gkyl_velocity_map_priv.h>
#include <gkyl_skip_cell.h>
#include <gkyl_skip_cell_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_alloc_flags_priv.h>
#include <assert.h>
}

#ifdef GKYL_HAVE_CUDA

struct gkyl_skip_cell* 
gkyl_skip_cell_new_cu_dev(struct gkyl_skip_cell *skip_cell_ho)
{
  struct gkyl_skip_cell *skip_cell = (struct gkyl_skip_cell *) gkyl_malloc(sizeof(*skip_cell));

  skip_cell->type = skip_cell_ho->type;
  skip_cell->f_threshold = skip_cell_ho->f_threshold;
  skip_cell->use_gpu = skip_cell_ho->use_gpu;
  skip_cell->phase_rng = skip_cell_ho->phase_rng;

  // Copy the host-side initialized array to the device.
  struct gkyl_array *mask = gkyl_array_cu_dev_new(GKYL_DOUBLE, skip_cell_ho->mask->ncomp, skip_cell_ho->mask->size);
  gkyl_array_copy(mask, skip_cell_ho->mask);

  skip_cell->mask = mask->on_dev;

  skip_cell->flags = 0;
  GKYL_SET_CU_ALLOC(skip_cell->flags);
  skip_cell->ref_count = gkyl_ref_count_init(gkyl_skip_cell_free);

  // Initialize the device object.
  struct gkyl_skip_cell *skip_cell_cu = (struct gkyl_skip_cell*) gkyl_cu_malloc(sizeof(*skip_cell_cu));
  gkyl_cu_memcpy(skip_cell_cu, skip_cell, sizeof(struct gkyl_skip_cell), GKYL_CU_MEMCPY_H2D);
  skip_cell->on_dev = skip_cell_cu;

  // The returned object should store host pointer to gkyl_array.
  skip_cell->mask = mask;

  return skip_cell;
}

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
        *skip_c = 1.0;
      } else {
        *skip_c = 0.0;
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
