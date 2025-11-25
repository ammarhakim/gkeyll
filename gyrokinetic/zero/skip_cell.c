#include <gkyl_skip_cell.h>
#include <gkyl_skip_cell_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>

#include <assert.h>
#include <float.h>

// Free function for reference counting.
void
gkyl_skip_cell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_skip_cell *skip_cell = container_of(ref, struct gkyl_skip_cell, ref_count);
  
  // Release the mask array.
  if (skip_cell->mask)
    gkyl_array_release(skip_cell->mask);
  
  if (GKYL_IS_CU_ALLOC(skip_cell->flags))
    gkyl_cu_free(skip_cell->on_dev);
  
  gkyl_free(skip_cell);
}

struct gkyl_skip_cell*
gkyl_skip_cell_new(struct gkyl_skip_cell_inp skip_cell_inp, struct gkyl_range phase_rng, bool use_gpu)
{
  // Allocate space for new skip cell object.
  struct gkyl_skip_cell *skip_cell = gkyl_malloc(sizeof(*skip_cell));

  skip_cell->type = skip_cell_inp.type;

  if (skip_cell->type == GKYL_GK_SKIP_CELL_F_THRESHOLD) {
    skip_cell->f_threshold = skip_cell_inp.threshold * pow(sqrt(2.0), phase_rng.ndim);
  }

  skip_cell->use_gpu = use_gpu;
  skip_cell->phase_rng = phase_rng;
  
  // Initialize the mask array on host.
  skip_cell->mask = gkyl_array_new(GKYL_DOUBLE, 1, phase_rng.volume);
  gkyl_array_clear(skip_cell->mask, false);

  skip_cell->flags = 0;
  GKYL_CLEAR_CU_ALLOC(skip_cell->flags);
  skip_cell->ref_count = gkyl_ref_count_init(gkyl_skip_cell_free);
  skip_cell->on_dev = skip_cell; // CPU skip_cell points to itself

  struct gkyl_skip_cell *skip_cell_out = skip_cell;
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    skip_cell_out = gkyl_skip_cell_new_cu_dev(skip_cell);
    gkyl_skip_cell_release(skip_cell);
  }
#endif
  
  return skip_cell_out;
}

void
gkyl_skip_cell_advance(struct gkyl_skip_cell *skip_cell, const struct gkyl_array *distf)
{
  if (skip_cell->type == GKYL_GK_SKIP_CELL_NONE) {
    return;
  }

#ifdef GKYL_HAVE_CUDA
  if (skip_cell->use_gpu) {
    gkyl_skip_cell_advance_cu(skip_cell, distf);
    return;
  }
#endif

  // Host implementation: iterate over phase space and update mask.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip_cell->phase_rng);
  
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&skip_cell->phase_rng, iter.idx);

    // Fetch the distribution function at this cell.
    const double *distf_c = gkyl_array_cfetch(distf, linidx);
    
    // Fetch the mask value for this cell.
    double *skip_c = gkyl_array_fetch(skip_cell->mask, linidx);
    
    // Mark cell as skippable if distribution function is below threshold.
    // Using fabs to check absolute value of the cell-averaged distribution function.
    if (fabs(distf_c[0]) < skip_cell->f_threshold) {
      *skip_c = 1.0;
    } else {
      *skip_c = 0.0;
    }
  }
}

void
gkyl_skip_cell_invert_mask(struct gkyl_skip_cell *skip_cell)
{
  gkyl_array_shiftc(skip_cell->mask, -1.0, 0);
  gkyl_array_scale(skip_cell->mask, -1.0);
}

bool
gkyl_skip_cell_is_cu_dev(const struct gkyl_skip_cell* skip_cell)
{
  return GKYL_IS_CU_ALLOC(skip_cell->flags);
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
