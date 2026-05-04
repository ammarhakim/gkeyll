#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_positivity_shift_zhang_shu.h>
#include <gkyl_positivity_shift_zhang_shu_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <acutest.h>

#include <float.h>
#include <math.h>
#include <string.h>

// Helper function to compare two arrays within a tolerance
static int
gkyl_array_compare(const struct gkyl_array *arr1, const struct gkyl_array *arr2, double tol)
{
  assert(arr1->size == arr2->size);
  assert(arr1->ncomp == arr2->ncomp);
  
  for (long i = 0; i < arr1->size; ++i) {
    const double *d1 = gkyl_array_cfetch(arr1, i);
    const double *d2 = gkyl_array_cfetch(arr2, i);
    for (int c = 0; c < arr1->ncomp; ++c) {
      if (!gkyl_compare(d1[c], d2[c], tol)) {
        return 0;
      }
    }
  }
  return 1;
}

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
set_modal_from_quad(const struct gkyl_basis *basis, const double *quad_vals, double *modal_vals)
{
  memset(modal_vals, 0, sizeof(double)*basis->num_basis);
  for (int k=0; k<basis->num_basis; ++k)
    basis->quad_nodal_to_modal(quad_vals, modal_vals, k);
}

static double
cell_min_from_modal(const struct gkyl_basis *basis, const double *modal_vals)
{
  const int nbasis = basis->num_basis;
  double quad_vals[nbasis];
  double fmin = DBL_MAX;

  for (int k=0; k<nbasis; ++k) {
    quad_vals[k] = 0.0;
    basis->modal_to_quad_nodal(modal_vals, quad_vals, k);
    fmin = fmin < quad_vals[k] ? fmin : quad_vals[k];
  }

  return fmin;
}

static void
fill_case(double *quad_vals, int nbasis, double base, double neg_val)
{
  for (int k=0; k<nbasis; ++k)
    quad_vals[k] = base;
  quad_vals[0] = neg_val;
  if (nbasis > 1)
    quad_vals[1] = 0.5*base;
  if (nbasis > 2)
    quad_vals[2] = 1.4*base;
}

static void
run_case(int cdim)
{
  const int vdim = 2;
  const int ndim = cdim + vdim;
  int cells[GKYL_MAX_DIM] = { 2, 2, 2, 2, 2, 2 };
  double lower[GKYL_MAX_DIM] = {0.0}, upper[GKYL_MAX_DIM] = {1.0};

  struct gkyl_basis basis;
  gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[GKYL_MAX_DIM] = {0};
  struct gkyl_range range, range_ext;
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  struct gkyl_array *f = mkarr(basis.num_basis, range_ext.volume);
  struct gkyl_array *f_ref = mkarr(basis.num_basis, range_ext.volume);

  const int nbasis = basis.num_basis;
  double quad_vals[nbasis];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    double *fc = gkyl_array_fetch(f, lidx);
    double *fr = gkyl_array_fetch(f_ref, lidx);

    fill_case(quad_vals, nbasis, 1.0 + 0.1*lidx, -0.2 - 0.05*lidx);
    set_modal_from_quad(&basis, quad_vals, fc);
    memcpy(fr, fc, sizeof(double)*nbasis);
  }

  struct gkyl_positivity_shift_zhang_shu *up = gkyl_positivity_shift_zhang_shu_new(basis);
  bool success = gkyl_positivity_shift_zhang_shu_advance(up, &range, f);
  TEST_CHECK(success);

  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    const double *fc = gkyl_array_cfetch(f, lidx);
    const double *fr = gkyl_array_cfetch(f_ref, lidx);

    TEST_CHECK(gkyl_compare(fc[0], fr[0], 1e-14));

    double fbar = up->cellav_fac * fr[0];
    double fmin = cell_min_from_modal(&basis, fr);
    TEST_CHECK(fmin < 0.0);

    double expected_theta = fbar / (fbar - fmin);
    for (int k=1; k<nbasis; ++k)
      TEST_CHECK(gkyl_compare(fc[k], expected_theta*fr[k], 1e-12));

    TEST_CHECK(cell_min_from_modal(&basis, fc) >= -1e-13);
  }

  struct gkyl_array *f_twice = mkarr(basis.num_basis, range_ext.volume);
  gkyl_array_copy(f_twice, f);
  success = gkyl_positivity_shift_zhang_shu_advance(up, &range, f_twice);
  TEST_CHECK(success);
  TEST_CHECK(gkyl_array_compare(f_twice, f, 1e-14));

  gkyl_array_release(f_twice);
  gkyl_positivity_shift_zhang_shu_release(up);
  gkyl_array_release(f_ref);
  gkyl_array_release(f);
}

static void
test_1x2v()
{
  run_case(1);
}

static void
test_2x2v()
{
  run_case(2);
}

static void
test_nonpositive_mean()
{
  const int cdim = 1, vdim = 2, ndim = cdim + vdim;
  int cells[GKYL_MAX_DIM] = { 1, 1, 1, 1, 1, 1 };
  double lower[GKYL_MAX_DIM] = {0.0}, upper[GKYL_MAX_DIM] = {1.0};

  struct gkyl_basis basis;
  gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[GKYL_MAX_DIM] = {0};
  struct gkyl_range range, range_ext;
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  struct gkyl_array *f = mkarr(basis.num_basis, range_ext.volume);
  double quad_vals[basis.num_basis];
  for (int k=0; k<basis.num_basis; ++k)
    quad_vals[k] = -1.0;

  double *fc = gkyl_array_fetch(f, 0);
  set_modal_from_quad(&basis, quad_vals, fc);
  double before[basis.num_basis];
  memcpy(before, fc, sizeof(before));

  struct gkyl_positivity_shift_zhang_shu *up = gkyl_positivity_shift_zhang_shu_new(basis);
  bool success = gkyl_positivity_shift_zhang_shu_advance(up, &range, f);

  TEST_CHECK(!success);
  for (int k=0; k<basis.num_basis; ++k)
    TEST_CHECK(gkyl_compare(fc[k], before[k], 1e-14));

  gkyl_positivity_shift_zhang_shu_release(up);
  gkyl_array_release(f);
}

TEST_LIST = {
  { "test_1x2v", test_1x2v },
  { "test_2x2v", test_2x2v },
  { "test_nonpositive_mean", test_nonpositive_mean },
  { NULL, NULL },
};