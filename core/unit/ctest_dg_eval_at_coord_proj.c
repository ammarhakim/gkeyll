#include <math.h>
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_dg_eval_at_coord_proj.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // Allocate array (filled with zeros).
  return use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
    : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
eval_f_2d(double t, const double *xn, double *restrict fout, void *ctx)
{
  // f(x,y) = 1 + 2x + 3y + 4xy.
  fout[0] = 1.0 + 2.0*xn[0] + 3.0*xn[1] + 4.0*xn[0]*xn[1];
}

static void
eval_f_1d(double t, const double *xn, double *restrict fout, void *ctx)
{
  // f(x0, y) = (1 + 2*x0) + (3 + 4*x0)*y, with x0 passed through ctx.
  double x0 = *(const double *)ctx;
  fout[0] = (1.0 + 2.0*x0) + (3.0 + 4.0*x0)*xn[0];
}

static void
test_2x_to_1x(int poly_order, bool use_gpu)
{
  // Project a 2D DG field onto 1D by evaluating in x at a fixed
  // computational coordinate x0. Verified by comparing against a direct 1D
  // projection of f(x0, y).

  // 2D donor grid.
  double lower_do[] = {0.0, 0.0}, upper_do[] = {1.0, 1.0};
  int cells_do[] = {2, 2};

  double x0 = 0.75; // lies in the second x-cell [0.5, 1.0]

  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);
  struct gkyl_rect_grid grid_do;
  gkyl_rect_grid_init(&grid_do, ndim_do, lower_do, upper_do, cells_do);

  // 1D target grid spans the y-direction.
  double lower_tar[] = {lower_do[1]}, upper_tar[] = {upper_do[1]};
  int cells_tar[] = {cells_do[1]};
  int ndim_tar = sizeof(cells_tar)/sizeof(cells_tar[0]);
  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim_tar, lower_tar, upper_tar, cells_tar);

  // Basis.
  struct gkyl_basis basis_do, basis_tar;
  gkyl_cart_modal_serendip(&basis_do, ndim_do, poly_order);
  gkyl_cart_modal_serendip(&basis_tar, ndim_tar, poly_order);

  // Ranges.
  int ghost_do[] = {1, 1, 1, 1, 1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &local_ext_do, &local_do);

  int ghost_tar[] = {ghost_do[1]};
  struct gkyl_range local_tar, local_ext_tar;
  gkyl_create_grid_ranges(&grid_tar, ghost_tar, &local_ext_tar, &local_tar);

  // Project f(x,y).
  struct gkyl_array *fdo = mkarr(use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_ho = use_gpu? mkarr(false, fdo->ncomp, fdo->size)
                                     : gkyl_array_acquire(fdo);
  
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid_do, &basis_do,
    poly_order+1, 1, eval_f_2d, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo_ho);
  gkyl_array_copy(fdo, fdo_ho);

  // Project f(x0,y) directly onto the 1D basis.
  struct gkyl_array *fref_ho = mkarr(false, basis_tar.num_basis, local_ext_tar.volume);
  gkyl_proj_on_basis *proj_tar = gkyl_proj_on_basis_new(&grid_tar, &basis_tar,
    poly_order+1, 1, eval_f_1d, &x0);
  gkyl_proj_on_basis_advance(proj_tar, 0.0, &local_tar, fref_ho);

  // Create the updater: evaluate in x (dir 0), keep y (1D target).
  int eval_dirs[] = {0};
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    &basis_do, &basis_tar, num_eval_dirs, eval_dirs, use_gpu);

  // Apply updater at x = x0.
  struct gkyl_array *ftar = mkarr(use_gpu, basis_tar.num_basis, local_ext_tar.volume);
  struct gkyl_array *ftar_ho = use_gpu? mkarr(false, ftar->ncomp, ftar->size)
                                      : gkyl_array_acquire(ftar);
  double eval_coords[] = {x0};
  bool pick_lower[] = {false, false};
  int known_index[] = {-1, -1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid_do, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Check answer.
  gkyl_array_copy(ftar_ho, ftar);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_tar);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local_tar, iter.idx);
    const double *ftar_c = gkyl_array_cfetch(ftar_ho, linidx);
    const double *fref_c = gkyl_array_cfetch(fref_ho, linidx);
    for (int k = 0; k < basis_tar.num_basis; k++) {
      TEST_CHECK(gkyl_compare(fref_c[k], ftar_c[k], 1e-14));
      TEST_MSG(" k:%d | Expected %.6e | Got %6.e\n", k, fref_c[k], ftar_c[k]);
    }
  }

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_proj_on_basis_release(proj_tar);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_ho);
  gkyl_array_release(ftar);
  gkyl_array_release(ftar_ho);
  gkyl_array_release(fref_ho);
}

static void
test_2x_to_scalar(int poly_order, bool use_gpu)
{
  // Evaluate a 2D DG field at a fixed (x0, y0), reducing to a scalar.
  // Verified by comparing against the analytic value f(x0, y0).

  // 2D donor grid.
  double lower_do[] = {0.0, 0.0}, upper_do[] = {1.0, 1.0};
  int cells_do[] = {2, 2};

  double x0 = 0.75, y0 = 0.25; // x0 in [0.5, 1.0], y0 in [0.0, 0.5]

  int ndim_do = sizeof(cells_do)/sizeof(cells_do[0]);
  struct gkyl_rect_grid grid_do;
  gkyl_rect_grid_init(&grid_do, ndim_do, lower_do, upper_do, cells_do);

  // Donor basis and range.
  struct gkyl_basis basis_do;
  gkyl_cart_modal_serendip(&basis_do, ndim_do, poly_order);

  int ghost_do[] = {1, 1, 1, 1, 1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &local_ext_do, &local_do);

  // Project f(x,y) onto the donor basis.
  struct gkyl_array *fdo = mkarr(use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_ho = use_gpu? mkarr(false, fdo->ncomp, fdo->size)
                                     : gkyl_array_acquire(fdo);
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid_do, &basis_do,
    poly_order+1, 1, eval_f_2d, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo_ho);
  gkyl_array_copy(fdo, fdo_ho);

  // Target range with 1 cell.
  int lower_tar[] = {1}, upper_tar[] = {1};
  struct gkyl_range local_tar;
  gkyl_range_init(&local_tar, 1, lower_tar, upper_tar);

  struct gkyl_array *ftar = mkarr(use_gpu, 1, local_tar.volume);
  struct gkyl_array *ftar_ho = use_gpu? mkarr(false, ftar->ncomp, ftar->size)
                                      : gkyl_array_acquire(ftar);

  // Evaluate and project field.
  int eval_dirs[] = {0, 1};
  int num_eval_dirs = sizeof(eval_dirs)/sizeof(eval_dirs[0]);
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    &basis_do, NULL, num_eval_dirs, eval_dirs, use_gpu);

  // Apply updater at (x0, y0).
  double eval_coords[] = {x0, y0};
  bool pick_lower[] = {false, false};
  int known_index[] = {-1, -1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid_do, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Check the results.
  double fref;
  eval_f_2d(0.0, (const double[2]){x0,y0}, &fref, 0);
  // Apply normalization for p0.
  fref *= 1.0/sqrt(2.0);

  gkyl_array_copy(ftar_ho, ftar);
  const double *ftar_c = gkyl_array_cfetch(ftar_ho, 0);
  TEST_CHECK(gkyl_compare(fref, ftar_c[0], 1e-14));
  TEST_MSG(" Expected %.6e | Got %.6e\n", fref, ftar_c[0]);

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_ho);
  gkyl_array_release(ftar);
  gkyl_array_release(ftar_ho);
}

void test_dg_evproj_2x_to_1x_p1_ho(void) {
  test_2x_to_1x(1, false);
}

void test_dg_evproj_2x_to_1x_p2_ho(void) {
  test_2x_to_1x(2, false);
}

void test_dg_evproj_2x_to_scalar_p1_ho(void) {
  test_2x_to_scalar(1, false);
}

void test_dg_evproj_2x_to_scalar_p2_ho(void) {
  test_2x_to_scalar(2, false);
}

#ifdef GKYL_HAVE_CUDA
void test_dg_evproj_2x_to_1x_p1_dev(void) {
  test_2x_to_1x(1, true);
}

void test_dg_evproj_2x_to_1x_p2_dev(void) {
  test_2x_to_1x(2, true);
}

void test_dg_evproj_2x_to_scalar_p1_dev(void) {
  test_2x_to_scalar(1, true);
}

void test_dg_evproj_2x_to_scalar_p2_dev(void) {
  test_2x_to_scalar(2, true);
}
#endif

TEST_LIST = {
  { "test_dg_evproj_2x_to_1x_p1_ho", test_dg_evproj_2x_to_1x_p1_ho },
  { "test_dg_evproj_2x_to_1x_p2_ho", test_dg_evproj_2x_to_1x_p2_ho },
  { "test_dg_evproj_2x_to_scalar_p1_ho", test_dg_evproj_2x_to_scalar_p1_ho },
  { "test_dg_evproj_2x_to_scalar_p2_ho", test_dg_evproj_2x_to_scalar_p2_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_dg_evproj_2x_to_1x_p1_dev", test_dg_evproj_2x_to_1x_p1_dev },
  { "test_dg_evproj_2x_to_1x_p2_dev", test_dg_evproj_2x_to_1x_p2_dev },
  { "test_dg_evproj_2x_to_scalar_p1_dev", test_dg_evproj_2x_to_scalar_p1_dev },
  { "test_dg_evproj_2x_to_scalar_p2_dev", test_dg_evproj_2x_to_scalar_p2_dev },
#endif
  { NULL, NULL },
};
