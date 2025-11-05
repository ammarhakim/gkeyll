#include <gkyl_skip_cell.h>
#include <gkyl_skip_cell_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>

#include <assert.h>
#include <float.h>

// Free function for reference counting.
static void
skip_cell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_skip_cell *skip_cell = container_of(ref, struct gkyl_skip_cell, ref_count);
  
  // Release the boolean array.
  if (skip_cell->booleans)
    gkyl_array_release(skip_cell->booleans);
  
  gkyl_free(skip_cell);
}

struct gkyl_skip_cell*
gkyl_skip_cell_new(struct gkyl_skip_cell_inp skip_cell_inp, struct gkyl_range phase_rng, bool use_gpu)
{
  // Allocate space for new skip cell object.
  struct gkyl_skip_cell *skip_cell = gkyl_malloc(sizeof(*skip_cell));

  if (skip_cell_inp.threshold > 0.0)
    skip_cell->skip_cell_threshold = skip_cell_inp.threshold * pow(sqrt(2.0), phase_rng.ndim);
  else
    skip_cell->skip_cell_threshold = -DBL_MAX;

  skip_cell->use_gpu = use_gpu;
  skip_cell->phase_rng = phase_rng;
  
  // Initialize the boolean mask array.
  // Size is the total number of cells in phase space.
  if (!use_gpu) {
    skip_cell->booleans = gkyl_array_new(GKYL_BOOL, 1, phase_rng.volume);
  }
#ifdef GKYL_HAVE_CUDA
  else {
    skip_cell->booleans = gkyl_array_cu_dev_new(GKYL_BOOL, 1, phase_rng.volume);
  }
#endif
  
  // Initialize all cells to false (not skipped).
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip_cell->phase_rng);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&skip_cell->phase_rng, iter.idx);
    bool *skip_c = gkyl_array_fetch(skip_cell->booleans, linidx);
    *skip_c = false;
  }

  // Initialize reference counter.
  skip_cell->ref_count = gkyl_ref_count_init(skip_cell_free);
  
  return skip_cell;
}

void
gkyl_skip_cell_advance(struct gkyl_skip_cell *skip_cell, const struct gkyl_array *distf)
{
#ifdef GKYL_HAVE_CUDA
  if (skip_cell->use_gpu) {
    gkyl_skip_cell_advance_cu(skip_cell, distf);
    return;
  }
#endif

  // Host implementation: iterate over phase space and update boolean mask.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip_cell->phase_rng);
  
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&skip_cell->phase_rng, iter.idx);

    // Fetch the distribution function at this cell.
    const double *distf_c = gkyl_array_cfetch(distf, linidx);
    
    // Fetch the boolean mask value for this cell.
    bool *skip_c = gkyl_array_fetch(skip_cell->booleans, linidx);
    
    // Mark cell as skippable if distribution function is below threshold.
    // Using fabs to check absolute value of the cell-averaged distribution function.
    if (fabs(distf_c[0]) < skip_cell->skip_cell_threshold) {
      *skip_c = true;
    } else {
      *skip_c = false;
    }
  }
}

struct gkyl_skip_cell*
gkyl_skip_cell_acquire(struct gkyl_skip_cell *skip_cell)
{
  gkyl_ref_count_inc(&skip_cell->ref_count);
  return (struct gkyl_skip_cell*) skip_cell;
}

void
gkyl_skip_cell_release(struct gkyl_skip_cell *skip_cell)
{
  if (skip_cell)
    gkyl_ref_count_dec(&skip_cell->ref_count);
}
