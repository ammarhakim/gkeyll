#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>

void
gkyl_dg_array_mask_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_array_mask *mask = container_of(ref, struct gkyl_dg_array_mask, ref_count);
  
  if (mask->mask)
    gkyl_array_release(mask->mask);
  
  if (GKYL_IS_CU_ALLOC(mask->flags))
    gkyl_cu_free(mask->on_dev);
  
  gkyl_free(mask);
}

struct gkyl_dg_array_mask*
gkyl_dg_array_mask_new(struct gkyl_dg_array_mask_inp mask_inp)
{
  struct gkyl_dg_array_mask *mask = gkyl_malloc(sizeof(*mask));

  mask->type = mask_inp.type;
  mask->use_gpu = mask_inp.use_gpu;
  mask->val_threshold = 0.0;
  mask->mask = NULL;
  mask->flags = 0;
  GKYL_CLEAR_CU_ALLOC(mask->flags);
  mask->ref_count = gkyl_ref_count_init(gkyl_dg_array_mask_free);
  mask->on_dev = mask; // CPU mask points to itself

  if (mask->type == GKYL_DG_ARRAY_MASK_NONE) {
    return mask;
  }
  else {
  if (mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD ||
    mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD) {
      mask->val_threshold = mask_inp.val_threshold * pow(sqrt(2.0), mask_inp.phase_rng.ndim);
    }
    
    // Initialize the mask array on host.
    mask->phase_rng = mask_inp.phase_rng;    
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
    else if (mask->type == GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD) {
      *mask_c = (abs_val > mask->val_threshold) ? 1.0 : -1.0;
    }
    else {
      *mask_c = -1.0;
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

void
gkyl_dg_array_mask_release(struct gkyl_dg_array_mask *mask)
{
  if (mask)
    gkyl_ref_count_dec(&mask->ref_count);
}
