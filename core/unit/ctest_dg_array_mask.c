#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_array_mask.h>
#include <gkyl_dg_array_mask_priv.h>
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
void test_mask_new(bool use_gpu)
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 1e-10,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  TEST_CHECK( mask != NULL );
  TEST_CHECK( mask->type == GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD );
  TEST_CHECK( mask->use_gpu == use_gpu );
  TEST_CHECK( mask->mask != NULL );
  TEST_CHECK( mask->mask->size == range.volume );
  TEST_CHECK( mask->mask->ncomp == 1 );
  
  // Check that mask is initialized to -1.0 (false)
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);

  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK( gkyl_compare(mask_d[0], -1.0, 1e-14) );
  }
  
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask with NONE type
void test_mask_none_type(bool use_gpu)
{
  struct gkyl_dg_array_mask_inp mask_inp = {
    .use_gpu = use_gpu
  };
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);

  TEST_CHECK( mask->type == GKYL_DG_ARRAY_MASK_NONE );
  TEST_CHECK( mask->use_gpu == use_gpu );
  TEST_CHECK( mask->val_threshold == 0.0 );

  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with threshold
void test_mask_advance_threshold(bool use_gpu)
{
  int shape[] = {10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array on host with varying values
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  for (unsigned i = 0; i < range.volume; ++i) {
    // Values alternate above and below threshold
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = (i % 2 == 0) ? 0.0 : 10.0;
  }
  
  // Copy to device if needed and run advance
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // Check mask values
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) < expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[0], 1.0, 1e-14) ); // True (masked)
    } else {
      TEST_CHECK( gkyl_compare(mask_d[0], -1.0, 1e-14) ); // False (not masked)
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values below threshold
void test_mask_advance_all_below(bool use_gpu)
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values below threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 1e-12);
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);

  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // All cells should be masked (1.0)
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK( gkyl_compare(mask_d[0], 1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with all values above threshold
void test_mask_advance_all_above(bool use_gpu)
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values above threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 100.0);
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // No cells should be masked (-1.0)
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK( gkyl_compare(mask_d[0], -1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with negative values
void test_mask_advance_negative_values(bool use_gpu)
{
  int shape[] = {8};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with negative and positive values
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  double vals[] = {0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0};
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = vals[i];
  }
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // Check mask values (uses fabs, so negative values treated as positive)
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) < expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[0], 1.0, 1e-14) );
    } else {
      TEST_CHECK( gkyl_compare(mask_d[0], -1.0, 1e-14) );
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask advance with GREATER_THAN threshold
void test_mask_advance_greater_than_threshold(bool use_gpu)
{
  int shape[] = {10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with varying values
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  for (unsigned i = 0; i < range.volume; ++i) {
    // Values alternate above and below threshold
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = (i % 2 == 0) ? 0.1 : 1.0;
  }
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // Check mask values (opposite of LESS_THAN)
  double expected_threshold = threshold * sqrt(2.0); // 1D scaling
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    if (fabs(arr_d[0]) > expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[0], 1.0, 1e-14) ); // True (masked)
    } else {
      TEST_CHECK( gkyl_compare(mask_d[0], -1.0, 1e-14) ); // False (not masked)
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values above threshold
void test_mask_advance_greater_than_all_above(bool use_gpu)
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 1e-15;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values above threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 100.0);
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // All cells should be masked (1.0) since values > threshold
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK( gkyl_compare(mask_d[0], 1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with all values below threshold
void test_mask_advance_greater_than_all_below(bool use_gpu)
{
  int shape[] = {5, 5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  double threshold = 10.0;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with all values below threshold
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  gkyl_array_clear(arr_ho, 1e-12);
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // No cells should be masked (-1.0) since values < threshold
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    TEST_CHECK( gkyl_compare(mask_d[0], -1.0, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test GREATER_THAN mask with negative values
void test_mask_advance_greater_than_negative_values(bool use_gpu)
{
  int shape[] = {8};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Create test array with negative and positive values
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  double vals[] = {0.1, -0.1, 1.0, -1.0, 0.3, -0.3, 2.0, -2.0};
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = vals[i];
  }
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Copy mask back to host for verification
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // Check mask values (uses fabs, so negative values treated as positive)
  double expected_threshold = threshold * sqrt(2.0);
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *arr_d = gkyl_array_cfetch(arr_ho, i);
    const double *mask_d = gkyl_array_cfetch(mask_ho, i);
    if (fabs(arr_d[0]) > expected_threshold) {
      TEST_CHECK( gkyl_compare(mask_d[0], 1.0, 1e-14) );
    } else {
      TEST_CHECK( gkyl_compare(mask_d[0], -1.0, 1e-14) );
    }
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask eval inline function
void test_mask_eval(bool use_gpu)
{
  int shape[] = {6};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  double threshold = 0.5;
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = threshold,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Use advance to set mask values via a test array
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  // Set up values so that mask will be: true, false, true, false, true, false
  double *arr_d;
  arr_d = gkyl_array_fetch(arr_ho, 0);
  arr_d[0] = 0.1;   // below threshold -> true
  arr_d = gkyl_array_fetch(arr_ho, 1);
  arr_d[0] = 1.0;   // above threshold -> false
  arr_d = gkyl_array_fetch(arr_ho, 2);
  arr_d[0] = 0.2;   // below threshold -> true
  arr_d = gkyl_array_fetch(arr_ho, 3);
  arr_d[0] = 2.0;   // above threshold -> false
  arr_d = gkyl_array_fetch(arr_ho, 4);
  arr_d[0] = 0.3;   // below threshold -> true
  arr_d = gkyl_array_fetch(arr_ho, 5);
  arr_d[0] = 5.0;   // above threshold -> false

  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Test eval function
  // This part does not work on the GPU because the boolean is evaluated and
  // returned on the GPU, but needs a host evaluation for the TEST_CHECK
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 0) == true );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 1) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 2) == true );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 3) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 4) == true );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 5) == false );
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_dg_array_mask_release(mask);
}

// Test mask eval with NONE type
void test_mask_eval_none_type(bool use_gpu)
{
  struct gkyl_dg_array_mask_inp mask_inp = {
    .use_gpu = use_gpu
  };
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
    
  // For NONE type, eval should always return false regardless of index
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 0) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 1) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 2) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 3) == false );
  TEST_CHECK( gkyl_dg_array_mask_eval(mask, 100) == false ); // Even invalid indices should be safe
  
  gkyl_dg_array_mask_release(mask);
}

// Test mask scale_by_cell
void test_mask_scale_by_cell(bool use_gpu)
{
  int shape[] = {5};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 0.5,
    .phase_rng = range,
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask = gkyl_dg_array_mask_new(mask_inp);
  
  // Set mask values
  // Create an array with alternating 1.0 and -1.0, then use mask_advance to set mask values
  struct gkyl_array *arr_ho = mkarr(false, 1, range.volume);
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);
    arr_d[0] = (i % 2 == 0) ? 0.0 : 1.0; // 0.0 will be masked (mask=1.0), 1.0 will not (mask=-1.0)
  }
  
  struct gkyl_array *arr = mkarr(use_gpu, 1, range.volume);
  gkyl_array_copy(arr, arr_ho);
  gkyl_dg_array_mask_advance(mask, arr);
  
  // Save original mask values
  const struct gkyl_array *mask_arr = gkyl_dg_array_mask_get_mask(mask);
  struct gkyl_array *mask_orig_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_orig_ho, mask_arr);
  double mask_orig[range.volume];
  for (unsigned i = 0; i < range.volume; ++i) {
    const double *mask_d = gkyl_array_cfetch(mask_orig_ho, i);
    mask_orig[i] = mask_d[0];
  }
  
  // Create array to multiply (1 component per cell)
  gkyl_array_clear(arr_ho, 0.0);
  for (unsigned i = 0; i < range.volume; ++i) {
    double *arr_d = gkyl_array_fetch(arr_ho, i);    
    arr_d[0] = 2.0 + i;
  }
  
  gkyl_array_copy(arr, arr_ho);
  
  // Scale mask by array
  gkyl_dg_array_mask_scale_by_cell(mask, arr);
  
  // Copy mask back to host for verification
  struct gkyl_array *mask_ho = mkarr(false, 1, range.volume);
  gkyl_array_copy(mask_ho, mask_arr);
  
  // Check results: mask should be multiplied by arr values
  for (unsigned cell = 0; cell < range.volume; ++cell) {
    const double *arr_val = gkyl_array_cfetch(arr_ho, cell);
    double expected = mask_orig[cell] * arr_val[0];
    const double *mask_d = gkyl_array_cfetch(mask_ho, cell);
    TEST_CHECK( gkyl_compare(mask_d[0], expected, 1e-14) );
  }
  
  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(mask_orig_ho);
  gkyl_array_release(mask_ho);
  gkyl_dg_array_mask_release(mask);
}

// Test mask acquire and release
void test_mask_acquire_release(bool use_gpu)
{
  int shape[] = {3, 3};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_dg_array_mask_inp mask_inp = {
    .type = GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD,
    .val_threshold = 1.0,
    .phase_rng = range,
    .use_gpu = use_gpu
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
void test_mask_threshold_scaling(bool use_gpu)
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
    .use_gpu = use_gpu
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
    .use_gpu = use_gpu
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
    .use_gpu = use_gpu
  };
  
  struct gkyl_dg_array_mask *mask3d = gkyl_dg_array_mask_new(mask_inp3d);
  double expected3d = threshold * pow(sqrt(2.0), 3);
  TEST_CHECK( gkyl_compare(mask3d->val_threshold, expected3d, 1e-14) );
  gkyl_dg_array_mask_release(mask3d);
}

// CPU test wrappers
void test_mask_new_ho() {
  test_mask_new(false);
}

void test_mask_none_type_ho() {
  test_mask_none_type(false);
}

void test_mask_advance_threshold_ho() {
  test_mask_advance_threshold(false);
}

void test_mask_advance_all_below_ho() {
  test_mask_advance_all_below(false);
}

void test_mask_advance_all_above_ho() {
  test_mask_advance_all_above(false);
}

void test_mask_advance_negative_values_ho() {
  test_mask_advance_negative_values(false);
}

void test_mask_advance_greater_than_threshold_ho() {
  test_mask_advance_greater_than_threshold(false);
}

void test_mask_advance_greater_than_all_above_ho() {
  test_mask_advance_greater_than_all_above(false);
}

void test_mask_advance_greater_than_all_below_ho() {
  test_mask_advance_greater_than_all_below(false);
}

void test_mask_advance_greater_than_negative_values_ho() {
  test_mask_advance_greater_than_negative_values(false);
}

void test_mask_eval_ho() {
  test_mask_eval(false);
}

void test_mask_eval_none_type_ho() {
  test_mask_eval_none_type(false);
}

void test_mask_scale_by_cell_ho() {
  test_mask_scale_by_cell(false);
}

void test_mask_acquire_release_ho() {
  test_mask_acquire_release(false);
}

void test_mask_threshold_scaling_ho() {
  test_mask_threshold_scaling(false);
}

#ifdef GKYL_HAVE_CUDA

// GPU test wrappers
void test_mask_new_dev() {
  test_mask_new(true);
}

void test_mask_none_type_dev() {
  test_mask_none_type(true);
}

void test_mask_advance_threshold_dev() {
  test_mask_advance_threshold(true);
}

void test_mask_advance_all_below_dev() {
  test_mask_advance_all_below(true);
}

void test_mask_advance_all_above_dev() {
  test_mask_advance_all_above(true);
}

void test_mask_advance_negative_values_dev() {
  test_mask_advance_negative_values(true);
}

void test_mask_advance_greater_than_threshold_dev() {
  test_mask_advance_greater_than_threshold(true);
}

void test_mask_advance_greater_than_all_above_dev() {
  test_mask_advance_greater_than_all_above(true);
}

void test_mask_advance_greater_than_all_below_dev() {
  test_mask_advance_greater_than_all_below(true);
}

void test_mask_advance_greater_than_negative_values_dev() {
  test_mask_advance_greater_than_negative_values(true);
}

void test_mask_eval_none_type_dev() {
  test_mask_eval_none_type(true);
}

void test_mask_scale_by_cell_dev() {
  test_mask_scale_by_cell(true);
}

void test_mask_acquire_release_dev() {
  test_mask_acquire_release(true);
}

void test_mask_threshold_scaling_dev() {
  test_mask_threshold_scaling(true);
}

#endif

TEST_LIST = {
  { "mask_new", test_mask_new_ho },
  { "mask_none_type", test_mask_none_type_ho },
  { "mask_advance_threshold", test_mask_advance_threshold_ho },
  { "mask_advance_all_below", test_mask_advance_all_below_ho },
  { "mask_advance_all_above", test_mask_advance_all_above_ho },
  { "mask_advance_negative_values", test_mask_advance_negative_values_ho },
  { "mask_advance_greater_than_threshold", test_mask_advance_greater_than_threshold_ho },
  { "mask_advance_greater_than_all_above", test_mask_advance_greater_than_all_above_ho },
  { "mask_advance_greater_than_all_below", test_mask_advance_greater_than_all_below_ho },
  { "mask_advance_greater_than_neg_vals", test_mask_advance_greater_than_negative_values_ho },
  { "mask_eval", test_mask_eval_ho },
  { "mask_eval_none_type", test_mask_eval_none_type_ho },
  { "mask_scale_by_cell", test_mask_scale_by_cell_ho },
  { "mask_acquire_release", test_mask_acquire_release_ho },
  { "mask_threshold_scaling", test_mask_threshold_scaling_ho },
#ifdef GKYL_HAVE_CUDA
  { "cu_mask_new", test_mask_new_dev },
  { "cu_mask_none_type", test_mask_none_type_dev },
  { "cu_mask_advance_threshold", test_mask_advance_threshold_dev },
  { "cu_mask_advance_all_below", test_mask_advance_all_below_dev },
  { "cu_mask_advance_all_above", test_mask_advance_all_above_dev },
  { "cu_mask_advance_negative_values", test_mask_advance_negative_values_dev },
  { "cu_mask_advance_greater_than_threshold", test_mask_advance_greater_than_threshold_dev },
  { "cu_mask_advance_greater_than_all_above", test_mask_advance_greater_than_all_above_dev },
  { "cu_mask_advance_greater_than_all_below", test_mask_advance_greater_than_all_below_dev },
  { "cu_mask_advance_greater_than_neg_vals", test_mask_advance_greater_than_negative_values_dev },
  { "cu_mask_eval_none_type", test_mask_eval_none_type_dev },
  { "cu_mask_scale_by_cell", test_mask_scale_by_cell_dev },
  { "cu_mask_acquire_release", test_mask_acquire_release_dev },
  { "cu_mask_threshold_scaling", test_mask_threshold_scaling_dev },
#endif
  { NULL, NULL },
};


