#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>
#include <time.h>

void test_reduce_min_max()
{
  int ncomp = 3, ncells = 200;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, ncells);

  for (size_t i=0; i<arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (size_t c=0; c<ncomp; ++c)
      d[c] = 0.5*i + 0.1*c;
  }
  
  double amin[ncomp], amax[ncomp];
  gkyl_array_reduce(&amin[0], arr, GKYL_MIN);
  gkyl_array_reduce(&amax[0], arr, GKYL_MAX);

  for (size_t c=0; c<ncomp; ++c) {
    TEST_CHECK( amin[c] == 0.1*c );
    TEST_CHECK( amax[c] == 0.5*(ncells-1) + 0.1*c );
  }

  gkyl_array_release(arr);
}

void test_reduce_min_max_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  int count = -1000;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = count+0.5;
    d[1] = count+1.5;
    d[2] = count+2.5;

    count += 1;
  }

  double amin[3], amax[3];
  gkyl_array_reduce_range(amin, arr, GKYL_MIN, &range);

  TEST_CHECK( amin[0] == -999.5 );
  TEST_CHECK( amin[1] == -998.5 );
  TEST_CHECK( amin[2] == -997.5 );

  gkyl_array_reduce_range(amax, arr, GKYL_MAX, &range);
  
  TEST_CHECK( amax[0] == -800.5 );
  TEST_CHECK( amax[1] == -799.5 );
  TEST_CHECK( amax[2] == -798.5 );

  gkyl_array_release(arr);
}

void test_reduce_sum_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  int ncomp = 3;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = 0.5;
    d[1] = 1.5;
    d[2] = 2.5;
  }

  double asum[ncomp];
  gkyl_array_reduce_range(asum, arr, GKYL_SUM, &range);

  TEST_CHECK( asum[0] == 0.5*range.volume );
  TEST_CHECK( asum[1] == 1.5*range.volume );
  TEST_CHECK( asum[2] == 2.5*range.volume );

  gkyl_array_release(arr);
}

void test_reduce_abs_max()
{
  int ncomp = 3, ncells = 200;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, ncells);

  for (int i=0; i<arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (int c=0; c<ncomp; ++c)
      d[c] = -0.5*i + 0.1*c;
  }
  
  double absmax[ncomp];
  gkyl_array_reduce(absmax, arr, GKYL_ABS_MAX);

  for (int c=0; c<ncomp; ++c) {
    TEST_CHECK( absmax[c] == fabs(-0.5*(ncells-1)+0.1*c) );
    TEST_MSG( "Got: %.9e | Expected: %.9e\n",absmax[c], fabs(-0.5*(ncells-1)+0.1*c) );
  }

  gkyl_array_release(arr);
}

void test_reduce_sq_sum_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = 0.5;
    d[1] = 1.5;
    d[2] = 2.5;
  }

  double asum[3];
  gkyl_array_reduce_range(asum, arr, GKYL_SQ_SUM, &range);

  TEST_CHECK( asum[0] == 0.5*0.5*range.volume );
  TEST_CHECK( asum[1] == 1.5*1.5*range.volume );
  TEST_CHECK( asum[2] == 2.5*2.5*range.volume );

  gkyl_array_release(arr);
}

void test_reduce_weighted_min_max()
{
  int ncomp = 3, ncells = 200;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, ncells);
  struct gkyl_array *brr = gkyl_array_new(GKYL_DOUBLE, ncomp, ncells);

  for (size_t i=0; i<arr->size; ++i) {
    double *arr_c = gkyl_array_fetch(arr, i);
    double *brr_c = gkyl_array_fetch(brr, i);
    for (size_t c=0; c<ncomp; ++c) {
      arr_c[c] = 0.5*i + 0.1*c;
      brr_c[c] = 0.25*i + 0.4*c;
    }
  }
  
  double amin[ncomp], amax[ncomp];
  gkyl_array_reduce_weighted(amin, arr, brr, GKYL_MIN);
  gkyl_array_reduce_weighted(amax, arr, brr, GKYL_MAX);

  for (size_t c=0; c<ncomp; ++c) {
    TEST_CHECK( amin[c] == (0.4*c)*(0.1*c) );
    TEST_CHECK( amax[c] == (0.25*(ncells-1) + 0.4*c)*(0.5*(ncells-1) + 0.1*c) );
  }

  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

void test_reduce_weighted_sq_sum_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  int ncomp = 3;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *brr = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *arr_c = gkyl_array_fetch(arr, loc);
    double *brr_c = gkyl_array_fetch(brr, loc);
    arr_c[0] = 0.5;
    arr_c[1] = 1.5;
    arr_c[2] = 2.5;

    brr_c[0] = 0.05;
    brr_c[1] = 0.15;
    brr_c[2] = 0.25;
  }

  double asum[ncomp];
  gkyl_array_reduce_weighted_range(asum, arr, brr, GKYL_SQ_SUM, &range);

  TEST_CHECK( gkyl_compare(asum[0], 0.05*0.5*0.5*range.volume, 1e-12) );
  TEST_CHECK( gkyl_compare(asum[1], 0.15*1.5*1.5*range.volume, 1e-12) );
  TEST_CHECK( gkyl_compare(asum[2], 0.25*2.5*2.5*range.volume, 1e-12) );

  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

void test_cu_array_reduce_max()
{
  unsigned long numComp = 3, numCells = 10;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_ho = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_ho_d = a1_ho->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_ho_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1, a1_ho);
  // create a device and host arrays to store reduction
  double* a1max = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  double* a1max_ho = gkyl_malloc(numComp*sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max, a1, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max_ho, a1max, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max_ho[k], (double)(numCells-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_ho);
  gkyl_array_release(a1);
  gkyl_array_release(a1_ho);

}

void test_cu_array_reduce_abs_max()
{
  unsigned long numComp = 3, numCells = 10;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_d = a1->data;
  for (int i=0; i<numCells; ++i) {
    for (int k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = -i+k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1absmax = gkyl_malloc(numComp*sizeof(double));
  double* a1absmax_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1absmax_cu, a1_cu, GKYL_ABS_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1absmax, a1absmax_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1absmax[k], fabs(-(double)(numCells-1)+(double)k*0.10), 1e-14) );
  }

  gkyl_free(a1absmax);
  gkyl_cu_free(a1absmax_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

}

void test_cu_array_reduce_max_big()
{

  unsigned long numComp = 12, numCells = 1000;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max_cu, a1_cu, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    // Make print statements to manually check the comparison
    printf("a1max[%d] = %g\n", k, a1max[k]);
    printf("numCells-1+(double)k*0.10 = %g\n", (double)(numCells-1)+(double)k*0.10);
    TEST_CHECK( gkyl_compare(a1max[k], (double)(numCells-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

}

void test_cu_array_reduce_max_range_1d()
{
  unsigned long numComp = 1, numCells = 10;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  int lower[] = {1}, upper[] = {numCells};
  int sublower[] = {3}, subupper[] = {7};
  struct gkyl_range range, subrange;
  gkyl_range_init(&range, 1, lower, upper);
  gkyl_sub_range_init(&subrange, &range, sublower, subupper);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);
  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(upper[0]-1)+(double)k*0.10, 1e-14) );
  }

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &subrange);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(subupper[0]-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_reduce_max_range_2d()
{
  int cells[] = {8, 10};
  int ghost[] = {1, 0};
  double lower[] = {0., -1.};
  double upper[] = {1., 1.};

  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i=0; i<numCells; ++i) {
    b1_d[i] = (double) i+1000;
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)numCells-1-i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_correct = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  double* b1max = gkyl_malloc(1*sizeof(double));
  double* b1max_cu = (double*) gkyl_cu_malloc(1*sizeof(double));
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double)*1);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(b1max_cu, b1_cu, GKYL_MAX, &range_ext);
  // Copy to host and check values.
  gkyl_cu_memcpy(b1max, b1max_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( gkyl_compare(b1max[0], (double)(1000+numCells-1), 1e-14) );

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);

  gkyl_array_reduce_range(a1max_correct, a1, GKYL_MAX, &range);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], a1max_correct[k], 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  gkyl_cu_free(b1max);
  gkyl_cu_free(b1max_cu);   
  gkyl_array_release(b1);
  gkyl_array_release(b1_cu);
}

void
test_cu_array_reduce_max_range_timer(int NX, int NY, int VX, int VY)
{
  int cells[] = {NX, NY, VX, VY};
  int ghost[] = {1, 1, 0, 0};
  double lower[] = {0.0, 0.0, 0.0, 0.0};
  double upper[] = {1.0, 1.0, 1.0, 1.0};

  int ndim = 4;
  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i=0; i<numCells; ++i) {
    b1_d[i] = (double) i+1000;
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)numCells-1-i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_correct = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));
    
  double* b1max = gkyl_malloc(1*sizeof(double));
  double* b1max_cu = (double*) gkyl_cu_malloc(1*sizeof(double));
  
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double)*1);

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);

  struct timespec tm = gkyl_wall_clock();

  for (int i=0; i<100; ++i)
    gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);
  double red_tm = gkyl_time_diff_now_sec(tm);

  printf("100 reductions on (%d,%d,%d,%d) took %g sec\n", NX, NY, VX, VY,
    red_tm);

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  gkyl_cu_free(b1max);
  gkyl_cu_free(b1max_cu);   
  gkyl_array_release(b1);
  gkyl_array_release(b1_cu);
}

void
test_cu_array_reduce_max_range_timer_32x32x40x40()
{
  test_cu_array_reduce_max_range_timer(32, 32, 40, 40);
}

void
test_cu_array_reduce_max_range_timer_32x32x32x32()
{
  test_cu_array_reduce_max_range_timer(32, 32, 32, 32);
}

void test_reduce_sum_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  int ncomp = 3;
  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *arr_d = gkyl_array_fetch(arr, loc);
    arr_d[0] = 0.5;
    arr_d[1] = 1.5;
    arr_d[2] = 2.5;
  }
  gkyl_array_copy(arr, arr_ho);

  double* asum = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  gkyl_array_reduce_range(asum, arr, GKYL_SUM, &range);

  double asum_ho[ncomp];
  gkyl_cu_memcpy(asum_ho, asum, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);

  TEST_CHECK( asum_ho[0] == 0.5*range.volume );
  TEST_CHECK( asum_ho[1] == 1.5*range.volume );
  TEST_CHECK( asum_ho[2] == 2.5*range.volume );

  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
  gkyl_cu_free(asum);
}

void test_cu_array_reduce_weighted_max()
{
  unsigned long numComp = 3, numCells = 10;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_ho = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1_ho = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_ho_d = a1_ho->data;
  double *b1_ho_d = b1_ho->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_ho_d[i*numComp+k] = i+k*0.10;
      b1_ho_d[i*numComp+k] = 0.1*i+k*0.01;
  }}
  gkyl_array_copy(a1, a1_ho);
  gkyl_array_copy(b1, b1_ho);
  // create a device and host arrays to store reduction
  double* a1max = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  double* a1max_ho = gkyl_malloc(numComp*sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max, a1, b1, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max_ho, a1max, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max_ho[k], (0.1*(numCells-1)+k*0.01)*((numCells-1)+k*0.10), 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_ho);
  gkyl_array_release(a1);
  gkyl_array_release(a1_ho);

}

void test_reduce_weighted_sum_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  int ncomp = 3;
  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *brr = gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, range.volume);
  struct gkyl_array *brr_ho = gkyl_array_new(GKYL_DOUBLE, ncomp, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *arr_d = gkyl_array_fetch(arr, loc);
    arr_d[0] = 0.5;
    arr_d[1] = 1.5;
    arr_d[2] = 2.5;

    double *brr_d = gkyl_brray_fetch(brr, loc);
    brr_d[0] = 0.05;
    brr_d[1] = 0.15;
    brr_d[2] = 0.25;
  }
  gkyl_array_copy(arr, arr_ho);
  gkyl_array_copy(brr, brr_ho);

  double* asum = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  gkyl_array_reduce_weighted_range(asum, arr, brr, GKYL_SUM, &range);

  double asum_ho[ncomp];
  gkyl_cu_memcpy(asum_ho, asum, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);

  TEST_CHECK( asum_ho[0] == 0.05*0.5*0.5*range.volume );
  TEST_CHECK( asum_ho[1] == 0.15*1.5*1.5*range.volume );
  TEST_CHECK( asum_ho[2] == 0.25*2.5*2.5*range.volume );

  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
  gkyl_array_release(brr);
  gkyl_array_release(brr_ho);
  gkyl_cu_free(asum);
}

#endif

TEST_LIST = {
  { "array_reduce_min_max", test_reduce_min_max },
  { "array_reduce_min_max_range", test_reduce_min_max_range },
  { "array_reduce_sum_range", test_reduce_sum_range },
  { "array_reduce_abs_max", test_reduce_abs_max },
  { "array_reduce_sq_sum_range", test_reduce_sq_sum_range },
  { "array_reduce_weighted_min_max", test_reduce_weighted_min_max },
  { "array_reduce_weighted_sq_sum_range", test_reduce_weighted_sq_sum_range },
#ifdef GKYL_HAVE_CUDA
  { "cu_array_reduce_max", test_cu_array_reduce_max },
  { "cu_array_reduce_abs_max", test_cu_array_reduce_abs_max },
  { "cu_array_reduce_max_big", test_cu_array_reduce_max_big },
  { "cu_array_reduce_max_range_1d", test_cu_array_reduce_max_range_1d },
  { "cu_array_reduce_max_range_2d", test_cu_array_reduce_max_range_2d },
  { "cu_array_reduce_max_range_timer_32x32x40x40", test_cu_array_reduce_max_range_timer_32x32x40x40  },
  { "cu_array_reduce_max_range_timer_32x32x32x32", test_cu_array_reduce_max_range_timer_32x32x32x32  },  
  { "cu_array_reduce_weighted_max", test_cu_array_reduce_weighted_max },
#endif
  { NULL, NULL },
};
