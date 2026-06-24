#include <acutest.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_differentiate.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct fin_ctx {
  double fac_m0; // For cell average.
  double fac_m1[GKYL_MAX_DIM]; // For linear monomials.
  double fac_m2[GKYL_MAX_DIM]; // For quadratic monomials in 1 variable.
  double fac_m3[3]; // For monomials in 2 variables.
  double fac_m4; // For monomials in 3 variables.
};

void fin_1x_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct fin_ctx *params = ctx;
  double fac0 = params->fac_m0;
  double fac1 = params->fac_m1[0];

  double x = xn[0];
  fout[0] = fac0 + fac1*x;
}

void
test_dg_differentiate_1x(int poly_order, bool use_gpu)
{
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {6};

  // Parameters for projected function.
  struct fin_ctx inp_params = {0};
  inp_params.fac_m0 = 2.0;
  inp_params.fac_m1[0] = 0.3;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project the input function.
  struct gkyl_proj_on_basis *proj_fin = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, fin_1x_func, &inp_params);

  // Input field array.
  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size)
                                     : gkyl_array_acquire(fin);

  gkyl_proj_on_basis_advance(proj_fin, 0.0, &local, fin_ho);
  gkyl_array_copy(fin, fin_ho);

  struct gkyl_array *derf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *derf_ho = use_gpu? mkarr(false, derf->ncomp, derf->size)
                                     : gkyl_array_acquire(derf);

  // Differentiate input field.
  int diff_dir = 0;
  int diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    m = 0;
    ref_val = inp_params.fac_m1[0]*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, derf_c[m], ref_val);

    m = 1;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, derf_c[m], ref_val);
  }

  gkyl_proj_on_basis_release(proj_fin);
  gkyl_array_release(fin);
  gkyl_array_release(fin_ho);
}

void fin_2x_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct fin_ctx *params = ctx;
  double fac0 = params->fac_m0;
  double fac1x = params->fac_m1[0];
  double fac1y = params->fac_m1[1];
  double fac3xy = params->fac_m3[0];

  double x = xn[0], y = xn[1];
  fout[0] = fac0 + fac1x*x + fac1y*y + fac3xy*x*y;
}

void
test_dg_differentiate_2x(int poly_order, bool use_gpu)
{
  double lower[] = {-M_PI, -2.0}, upper[] = {M_PI, 2.0};
  int cells[] = {6, 4};

  // Parameters for projected function.
  struct fin_ctx inp_params = {0};
  inp_params.fac_m0 = 2.0;
  inp_params.fac_m1[0] = 0.3;
  inp_params.fac_m1[1] = 1.7;
  inp_params.fac_m3[0] = 0.;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project the input function.
  struct gkyl_proj_on_basis *proj_fin = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, fin_2x_func, &inp_params);

  // Input field array.
  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size)
                                     : gkyl_array_acquire(fin);

  gkyl_proj_on_basis_advance(proj_fin, 0.0, &local, fin_ho);
  gkyl_array_copy(fin, fin_ho);

  struct gkyl_array *derf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *derf_ho = use_gpu? mkarr(false, derf->ncomp, derf->size)
                                     : gkyl_array_acquire(derf);
  
  // Differentiate input field along x.
  int diff_dir = 0;
  int diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    m = 0;
    ref_val = inp_params.fac_m1[0]*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 1;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 2;
    ref_val = inp_params.fac_m3[0]*2.0/sqrt(3.0);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 3;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);
  }

  // Differentiate input field along y.
  diff_dir = 1;
  diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    m = 0;
    ref_val = inp_params.fac_m1[1]*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 1;
    ref_val = inp_params.fac_m3[0]*2.0/sqrt(3.0);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 2;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 3;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);
  }

  gkyl_proj_on_basis_release(proj_fin);
  gkyl_array_release(fin);
  gkyl_array_release(fin_ho);
}

void fin_3x_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct fin_ctx *params = ctx;
  double fac0 = params->fac_m0;
  double fac1x = params->fac_m1[0];
  double fac1y = params->fac_m1[1];
  double fac1z = params->fac_m1[2];
  double fac3xy = params->fac_m3[0];
  double fac3xz = params->fac_m3[1];
  double fac3yz = params->fac_m3[2];
  double fac4 = params->fac_m4;

  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = fac0 + fac1x*x + fac1y*y + fac1z*z + fac3xy*x*y + fac3xz*x*z + fac3yz*y*z + fac4*x*y*z;
}

void
test_dg_differentiate_3x(int poly_order, bool use_gpu)
{
  double lower[] = {-M_PI, -2.0, 1.0}, upper[] = {M_PI, 2.0, 3.5};
  int cells[] = {6, 4, 4};

  // Parameters for projected function.
  struct fin_ctx inp_params = {0};
  inp_params.fac_m0 = 2.0;
  inp_params.fac_m1[0] = 0.3;
  inp_params.fac_m1[1] = 1.7;
  inp_params.fac_m1[2] = 0.9;
  inp_params.fac_m3[0] = 0.;
  inp_params.fac_m3[1] = 0.;
  inp_params.fac_m3[2] = 0.;
  inp_params.fac_m4 = 0.;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project the input function.
  struct gkyl_proj_on_basis *proj_fin = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, fin_3x_func, &inp_params);

  // Input field array.
  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size)
                                     : gkyl_array_acquire(fin);

  gkyl_proj_on_basis_advance(proj_fin, 0.0, &local, fin_ho);
  gkyl_array_copy(fin, fin_ho);

  struct gkyl_array *derf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *derf_ho = use_gpu? mkarr(false, derf->ncomp, derf->size)
                                     : gkyl_array_acquire(derf);
  
  // Differentiate input field along x.
  int diff_dir = 0;
  int diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    m = 0;
    ref_val = inp_params.fac_m1[0]*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 1;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 2;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 3;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);
  }

  // Differentiate input field along y.
  diff_dir = 1;
  diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    m = 0;
    ref_val = inp_params.fac_m1[1]*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 1;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 2;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 3;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);
  }

  // Differentiate input field along z.
  diff_dir = 2;
  diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    m = 0;
    ref_val = inp_params.fac_m1[2]*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 1;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 2;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 3;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);
  }

  gkyl_proj_on_basis_release(proj_fin);
  gkyl_array_release(fin);
  gkyl_array_release(fin_ho);
}

void test_dg_differentiate_1x_p1_ho() {
  test_dg_differentiate_1x(1, false);
}

void test_dg_differentiate_2x_p1_ho() {
  test_dg_differentiate_2x(1, false);
}

void test_dg_differentiate_3x_p1_ho() {
  test_dg_differentiate_3x(1, false);
}

#ifdef GKYL_HAVE_CUDA
void test_dg_differentiate_1x_p1_cu() {
  test_dg_differentiate_1x(1, true);
}

void test_dg_differentiate_2x_p1_cu() {
  test_dg_differentiate_2x(1, true);
}

void test_dg_differentiate_3x_p1_cu() {
  test_dg_differentiate_3x(1, true);
}
#endif

TEST_LIST = {
  { "test_dg_differentiate_1x_p1_ho", test_dg_differentiate_1x_p1_ho },
  { "test_dg_differentiate_2x_p1_ho", test_dg_differentiate_2x_p1_ho },
  { "test_dg_differentiate_3x_p1_ho", test_dg_differentiate_3x_p1_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_dg_differentiate_1x_p1_cu", test_dg_differentiate_1x_p1_cu },
  { "test_dg_differentiate_2x_p1_cu", test_dg_differentiate_2x_p1_cu },
  { "test_dg_differentiate_3x_p1_cu", test_dg_differentiate_3x_p1_cu },
#endif
  { NULL, NULL },
};
