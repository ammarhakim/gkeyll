#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_array_mask.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <math.h>

// Helper function to create test arrays
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                 : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Test basic mask creation and initialization
void test_mask_new()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 1e-10,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  TEST_CHECK( mask != NULL );
  TEST_CHECK( mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD );
  TEST_CHECK( mask->use_gpu == false );
  TEST_CHECK( mask->mask != NULL );
  TEST_CHECK( mask->mask->size == range.volume );
  TEST_CHECK( mask->mask->ncomp == 1 );
  
  // Check that mask is initialized to -1.0 (false)
  double *mask_d = mask->mask->data;
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
  }
  
  gkyl_dg_array_mask_release(mask);
}

// Test mask with NONE type
void test_mask_none_type()
{
  int shape[] = {5, 10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_NONE,
    .val_threshold = 0.0,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  TEST_CHECK( mask != NULL );
  TEST_CHECK( mask->type == GKYL_DG_ARRAY_MASK_NONE );
  
  // Create a test array
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr, 10.0);
  
  // Advance should do nothing for NONE type
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Mask should still be at initial -1.0 values
  double *mask_d = mask->mask->data;
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with threshold
void test_mask_advance_threshold()
{
  int shape[] = {10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with varying values
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  double *arr_d = arr->data;
  for (unsigned i = 0; i < arr->size; ++i) {
    // Values alternate above and below threshold
    arr_d[i] = (i % 2 == 0) ? 0.1 : 1.0;
  }
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Check mask values
  double *mask_d = mask->mask->data;
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    if (fabs(arr_d[i]) < expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) ); // True (masked)
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) ); // False (not masked)
    }
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values below threshold
void test_mask_advance_all_below()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values below threshold
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr, 1e-12);
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // All cells should be masked (1.0)
  double *mask_d = mask->mask->data;
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values above threshold
void test_mask_advance_all_above()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values above threshold
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr, 100.0);
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // No cells should be masked (-1.0)
  double *mask_d = mask->mask->data;
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with negative values
void test_mask_advance_negative_values()
{
  int shape[] = {8};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with negative and positive values
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  double *arr_d = arr->data;
  double vals[] = {0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0};
  for (unsigned i = 0; i < arr->size; ++i) {
    arr_d[i] = vals[i];
  }
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Check mask values (uses fabs, so negative values treated as positive)
  double *mask_d = mask->mask->data;
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    if (fabs(arr_d[i]) < expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
    }
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with GREATER_THAN threshold
void test_mask_advance_greater_than_threshold()
{
  int shape[] = {10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with varying values
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  double *arr_d = arr->data;
  for (unsigned i = 0; i < arr->size; ++i) {
    // Values alternate above and below threshold
    arr_d[i] = (i % 2 == 0) ? 0.1 : 1.0;
  }
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Check mask values (opposite of LESS_THAN)
  double *mask_d = mask->mask->data;
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    if (fabs(arr_d[i]) > expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) ); // True (masked)
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) ); // False (not masked)
    }
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values above threshold
void test_mask_advance_greater_than_all_above()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values above threshold
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr, 100.0);
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // All cells should be masked (1.0) since values > threshold
  double *mask_d = mask->mask->data;
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values below threshold
void test_mask_advance_greater_than_all_below()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values below threshold
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr, 1e-12);
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // No cells should be masked (-1.0) since values < threshold
  double *mask_d = mask->mask->data;
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with negative values
void test_mask_advance_greater_than_negative_values()
{
  int shape[] = {8};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with negative and positive values
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  double *arr_d = arr->data;
  double vals[] = {0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0};
  for (unsigned i = 0; i < arr->size; ++i) {
    arr_d[i] = vals[i];
  }
  
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Check mask values (uses fabs, so negative values treated as positive)
  double *mask_d = mask->mask->data;
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    if (fabs(arr_d[i]) > expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
    }
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask eval inline function
void test_mask_eval()
{
  int shape[] = {6};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Manually set mask values
  double *mask_d = mask->mask->data;
  mask_d[0] = 1.0;  // True
  mask_d[1] = -1.0; // False
  mask_d[2] = 1.0;  // True
  mask_d[3] = -1.0; // False
  mask_d[4] = 0.5;  // Should be true (> 0)
  mask_d[5] = -0.5; // Should be false (<= 0)
  
  // Test eval function
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 0) == true );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 1) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 2) == true );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 3) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 4) == true );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 5) == false );
  
  gkyl_dg_array_mask_release(mask);
}

// Test mask eval with NONE type
void test_mask_eval_none_type()
{
  int shape[] = {4};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_NONE,
    .val_threshold = 0.0,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Manually set mask values
  double *mask_d = mask->mask->data;
  mask_d[0] = 1.0;
  mask_d[1] = -1.0;
  mask_d[2] = 1.0;
  mask_d[3] = -1.0;
  
  // For NONE type, eval should always return false
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 0) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 1) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 2) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 3) == false );
  
  gkyl_dg_array_mask_release(mask);
}

// Test mask scale_by_cell
void test_mask_scale_by_cell()
{
  int shape[] = {5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 0.5,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Set mask values
  double *mask_d = mask->mask->data;
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    mask_d[i] = (i % 2 == 0) ? 1.0 : -1.0;
  }
  
  // Save original mask values
  double mask_orig[mask->mask->size];
  for (unsigned i = 0; i < mask->mask->size; ++i) {
    mask_orig[i] = mask_d[i];
  }
  
  // Create array to multiply (1 component per cell)
  struct gkyl_array *arr = mkarr(false, 1, range.volume);
  double *arr_d = arr->data;
  for (unsigned i = 0; i < arr->size; ++i) {
    arr_d[i] = 2.0 + i;
  }
  
  // Scale mask by array
  gkyl_dg_array_mask_scale_by_cell(mask, arr);
  
  // Check results: mask should be multiplied by arr values
  for (unsigned cell = 0; cell < range.volume; ++cell) {
    double arr_val = arr_d[cell];
    double expected = mask_orig[cell] * arr_val;
    TEST_CHECK( gkyl_compare(mask_d[cell], expected, 1e-14) );
  }
  
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask acquire and release
void test_mask_acquire_release()
{
  int shape[] = {3, 3};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 1.0,
    .phase_rng = range,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Acquire reference
  struct gkyl_dg_array_mask *mask2 = gkyl_dg_array_mask_acquire(mask);
  TEST_CHECK( mask2 == mask );
  
  // Release both references
  gkyl_dg_array_mask_release(mask2);
  gkyl_dg_array_mask_release(mask);
}

// Test threshold scaling with different dimensions
void test_mask_threshold_scaling()
{
  // Test 1D
  int shape1d[] = {10};
  struct gkyl_range range1d;
  gkyl_range_init_from_shape(&range1d, 1, shape1d);
  
  double threshold = 1.0;
  struct gkyl_dg_array_mask_inp mask_inp1d = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range1d,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask1d = gkyl_dg_array_mask_new(mask_inp1d);
  double expected1d = threshold * sqrt(2.0);
  TEST_CHECK( gkyl_compare(mask1d->val_threshold, expected1d, 1e-14) );
  gkyl_dg_array_mask_release(mask1d);
  
  // Test 2D
  int shape2d[] = {5, 5};
  struct gkyl_range range2d;
  gkyl_range_init_from_shape(&range2d, 2, shape2d);
  
  struct gkyl_dg_array_mask_inp mask_inp2d = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range2d,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask2d = gkyl_dg_array_mask_new(mask_inp2d);
  double expected2d = threshold * pow(sqrt(2.0), 2);
  TEST_CHECK( gkyl_compare(mask2d->val_threshold, expected2d, 1e-14) );
  gkyl_dg_array_mask_release(mask2d);
  
  // Test 3D
  int shape3d[] = {3, 3, 3};
  struct gkyl_range range3d;
  gkyl_range_init_from_shape(&range3d, 3, shape3d);
  
  struct gkyl_dg_array_mask_inp mask_inp3d = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range3d,
    .use_gpu = false
  };
  
  struct gkyl_dg_array_mask *mask3d = gkyl_dg_array_mask_new(mask_inp3d);
  double expected3d = threshold * pow(sqrt(2.0), 3);
  TEST_CHECK( gkyl_compare(mask3d->val_threshold, expected3d, 1e-14) );
  gkyl_dg_array_mask_release(mask3d);
}

#ifdef GKYL_HAVE_CUDA

// Test basic mask creation and initialization on GPU
void test_cu_mask_new()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 1e-10,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  TEST_CHECK( mask != NULL );
  TEST_CHECK( mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD );
  TEST_CHECK( mask->use_gpu == true );
  TEST_CHECK( mask->mask != NULL );
  TEST_CHECK( mask->mask->size == range.volume );
  TEST_CHECK( mask->mask->ncomp == 1 );
  
  // Copy to host and check initialization to -1.0
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  double *mask_d = mask_ho->data;
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
  }
  
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with threshold on GPU
void test_cu_mask_advance_threshold()
{
  int shape[] = {10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array on host
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  double *arr_d = arr_ho->data;
  for (unsigned i = 0; i < arr_ho->size; ++i) {
    // Values alternate above and below threshold
    arr_d[i] = (i % 2 == 0) ? 0.1 : 1.0;
  }
  
  // Copy to device
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  // Run advance on GPU
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // Check mask values
  double *mask_d = mask_ho->data;
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    if (fabs(arr_d[i]) < expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) ); // True (masked)
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) ); // False (not masked)
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values below threshold on GPU
void test_cu_mask_advance_all_below()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values below threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 1e-12);
  
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // All cells should be masked (1.0)
  double *mask_d = mask_ho->data;
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values above threshold on GPU
void test_cu_mask_advance_all_above()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values above threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 100.0);
  
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // No cells should be masked (-1.0)
  double *mask_d = mask_ho->data;
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with negative values on GPU
void test_cu_mask_advance_negative_values()
{
  int shape[] = {8};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with negative and positive values
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  double *arr_d = arr_ho->data;
  double vals[] = {0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0};
  for (unsigned i = 0; i < arr_ho->size; ++i) {
    arr_d[i] = vals[i];
  }
  
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // Check mask values (uses fabs, so negative values treated as positive)
  double *mask_d = mask_ho->data;
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    if (fabs(arr_d[i]) < expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with GREATER_THAN threshold on GPU
void test_cu_mask_advance_greater_than_threshold()
{
  int shape[] = {10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array on host
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  double *arr_d = arr_ho->data;
  for (unsigned i = 0; i < arr_ho->size; ++i) {
    // Values alternate above and below threshold
    arr_d[i] = (i % 2 == 0) ? 0.1 : 1.0;
  }
  
  // Copy to device
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  // Run advance on GPU
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // Check mask values (opposite of LESS_THAN)
  double *mask_d = mask_ho->data;
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    if (fabs(arr_d[i]) > expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) ); // True (masked)
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) ); // False (not masked)
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values above threshold on GPU
void test_cu_mask_advance_greater_than_all_above()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values above threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 100.0);
  
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // All cells should be masked (1.0) since values > threshold
  double *mask_d = mask_ho->data;
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values below threshold on GPU
void test_cu_mask_advance_greater_than_all_below()
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values below threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 1e-12);
  
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // No cells should be masked (-1.0) since values < threshold
  double *mask_d = mask_ho->data;
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with negative values on GPU
void test_cu_mask_advance_greater_than_negative_values()
{
  int shape[] = {8};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with negative and positive values
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  double *arr_d = arr_ho->data;
  double vals[] = {0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0};
  for (unsigned i = 0; i < arr_ho->size; ++i) {
    arr_d[i] = vals[i];
  }
  
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  gkyl_dg_array_mask_advance(mask, arr_dev);
  
  // Copy mask back to host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask->mask);
  
  // Check mask values (uses fabs, so negative values treated as positive)
  double *mask_d = mask_ho->data;
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    if (fabs(arr_d[i]) > expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[i], 1.0, 1e-14) );
    } else {
      TEST_CHECK( gkyl_compare(mask_d[i], -1.0, 1e-14) );
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask scale_by_cell on GPU
void test_cu_mask_scale_by_cell()
{
  int shape[] = {5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 0.5,
    .phase_rng = range,
    .use_gpu = true
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Set mask values on host
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  double *mask_ho_d = mask_ho->data;
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    mask_ho_d[i] = (i % 2 == 0) ? 1.0 : -1.0;
  }
  
  // Copy mask to device
  gkyl_array_copy(mask->mask, mask_ho);
  
  // Save original mask values
  double mask_orig[mask_ho->size];
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    mask_orig[i] = mask_ho_d[i];
  }
  
  // Create array to multiply (1 component per cell)
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  double *arr_ho_d = arr_ho->data;
  for (unsigned i = 0; i < arr_ho->size; ++i) {
    arr_ho_d[i] = 2.0 + i;
  }
  
  // Copy to device
  struct gkyl_array *arr_dev = mkarr(true, 1, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  // Scale mask by array on GPU
  gkyl_dg_array_mask_scale_by_cell(mask, arr_dev);
  
  // Copy mask back to host
  gkyl_array_copy(mask_ho, mask->mask);
  
  // Check results: mask should be multiplied by arr values
  for (unsigned cell = 0; cell < range.volume; ++cell) {
    double arr_val = arr_ho_d[cell];
    double expected = mask_orig[cell] * arr_val;
    TEST_CHECK( gkyl_compare(mask_ho_d[cell], expected, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

#endif

TEST_LIST = {
  { "mask_new", test_mask_new },
  { "mask_none_type", test_mask_none_type },
  { "mask_advance_threshold", test_mask_advance_threshold },
  { "mask_advance_all_below", test_mask_advance_all_below },
  { "mask_advance_all_above", test_mask_advance_all_above },
  { "mask_advance_negative_values", test_mask_advance_negative_values },
  { "mask_advance_greater_than_threshold", test_mask_advance_greater_than_threshold },
  { "mask_advance_greater_than_all_above", test_mask_advance_greater_than_all_above },
  { "mask_advance_greater_than_all_below", test_mask_advance_greater_than_all_below },
  { "mask_advance_greater_than_neg_vals", test_mask_advance_greater_than_negative_values },
  { "mask_eval", test_mask_eval },
  { "mask_eval_none_type", test_mask_eval_none_type },
  { "mask_scale_by_cell", test_mask_scale_by_cell },
  { "mask_acquire_release", test_mask_acquire_release },
  { "mask_threshold_scaling", test_mask_threshold_scaling },
#ifdef GKYL_HAVE_CUDA
  { "cu_mask_new", test_cu_mask_new },
  { "cu_mask_advance_threshold", test_cu_mask_advance_threshold },
  { "cu_mask_advance_all_below", test_cu_mask_advance_all_below },
  { "cu_mask_advance_all_above", test_cu_mask_advance_all_above },
  { "cu_mask_advance_negative_values", test_cu_mask_advance_negative_values },
  { "cu_mask_advance_greater_than_threshold", test_cu_mask_advance_greater_than_threshold },
  { "cu_mask_advance_greater_than_all_above", test_cu_mask_advance_greater_than_all_above },
  { "cu_mask_advance_greater_than_all_below", test_cu_mask_advance_greater_than_all_below },
  { "cu_mask_advance_greater_than_neg_vals", test_cu_mask_advance_greater_than_negative_values },
  { "cu_mask_scale_by_cell", test_cu_mask_scale_by_cell },
#endif
  { NULL, NULL },
};
