// Test for dg_array_mask CUDA implementation.
#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_array_mask.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <math.h>

void test_cu_mask_new();
void test_cu_mask_advance_threshold();
void test_cu_mask_advance_all_below();
void test_cu_mask_advance_all_above();
void test_cu_mask_advance_negative_values();
void test_cu_mask_advance_greater_than_threshold();
void test_cu_mask_advance_greater_than_all_above();
void test_cu_mask_advance_greater_than_all_below();
void test_cu_mask_advance_greater_than_negative_values();
void test_cu_mask_scale_by_cell();

// Helper function to create test arrays
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                 : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
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
  double *mask_d = mask_ho->data;
  for (unsigned i = 0; i < mask_ho->size; ++i) {
    mask_d[i] = (i % 2 == 0) ? 1.0 : -1.0;
  }
  
  // Copy mask to device
  gkyl_array_copy(mask->mask, mask_ho);
  
  // Create array to multiply (3 components per cell)
  struct gkyl_array *arr_ho = mkarr(false, 3, range.volume);
  double *arr_d = arr_ho->data;
  for (unsigned i = 0; i < arr_ho->size; ++i) {
    arr_d[i] = 2.0 + i;
  }
  
  // Save original values
  double arr_orig[arr_ho->size];
  for (unsigned i = 0; i < arr_ho->size; ++i) {
    arr_orig[i] = arr_d[i];
  }
  
  // Copy to device
  struct gkyl_array *arr_dev = mkarr(true, 3, range.volume);
  gkyl_array_copy(arr_dev, arr_ho);
  
  // Scale by mask on GPU
  gkyl_dg_array_mask_scale_by_cell(mask, arr_dev);
  
  // Copy result back to host
  gkyl_array_copy(arr_ho, arr_dev);
  
  // Check results: each cell's components should be multiplied by mask value
  for (unsigned cell = 0; cell < range.volume; ++cell) {
    double mask_val = mask_d[cell];
    for (int comp = 0; comp < 3; ++comp) {
      int idx = cell * 3 + comp;
      double expected = arr_orig[idx] * mask_val;
      TEST_CHECK( gkyl_compare(arr_d[idx], expected, 1e-14) );
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr_dev);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

#endif

TEST_LIST = {
#ifdef GKYL_HAVE_CUDA
  { "cu_mask_new", test_cu_mask_new },
  { "cu_mask_advance_threshold", test_cu_mask_advance_threshold },
  { "cu_mask_advance_all_below", test_cu_mask_advance_all_below },
  { "cu_mask_advance_all_above", test_cu_mask_advance_all_above },
  { "cu_mask_advance_negative_values", test_cu_mask_advance_negative_values },
  { "cu_mask_advance_greater_than_threshold", test_cu_mask_advance_greater_than_threshold },
  { "cu_mask_advance_greater_than_all_above", test_cu_mask_advance_greater_than_all_above },
  { "cu_mask_advance_greater_than_all_below", test_cu_mask_advance_greater_than_all_below },
  { "cu_mask_advance_greater_than_negative_values", test_cu_mask_advance_greater_than_negative_values },
  { "cu_mask_scale_by_cell", test_cu_mask_scale_by_cell },
#endif
  { NULL, NULL },
};
