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
test_2x_to_1x(int poly_order)
{
  // Project a 2D DG field onto 1D by evaluating in x at a fixed
  // computational coordinate x0. Verified by comparing against a direct 1D
  // projection of f(x0, y).

  // 2D donor grid.
  double lower_do[] = {0.0, 0.0}, upper_do[] = {1.0, 1.0};
  int cells_do[] = {2, 2};

  double x0 = 0.75; // lies in the second x-cell [0.5, 1.0]

  int ndim = sizeof(cells_do)/sizeof(cells_do[0]);
  struct gkyl_rect_grid grid2;
  gkyl_rect_grid_init(&grid2, ndim, lower_do, upper_do, cells_do);

  // 1D target grid spans the y-direction.
  double lower_tar[] = {lower_do[1]}, upper_tar[] = {upper_do[1]};
  int cells_tar[] = {cells_do[1]};
  struct gkyl_rect_grid grid1;
  gkyl_rect_grid_init(&grid1, 1, lower_tar, upper_tar, cells_tar);

  struct gkyl_basis basis_do, basis_tar;
  gkyl_cart_modal_serendip(&basis_do, ndim, poly_order);
  gkyl_cart_modal_serendip(&basis_tar, 1, poly_order);

  // Ranges.
  int ghost_do[] = {1, 1};
  struct gkyl_range local_do, local_ext_do;
  gkyl_create_grid_ranges(&grid2, ghost_do, &local_ext_do, &local_do);

  int ghost_tar[] = {ghost_do[1]};
  struct gkyl_range local_tar, local_ext_tar;
  gkyl_create_grid_ranges(&grid1, ghost_tar, &local_ext_tar, &local_tar);

  // Project f(x,y).
  struct gkyl_array *fdo = gkyl_array_new(GKYL_DOUBLE, basis_do.num_basis, local_ext_do.volume);
  gkyl_proj_on_basis *proj_do = gkyl_proj_on_basis_new(&grid2, &basis_do,
    poly_order+1, 1, eval_f_2d, NULL);
  gkyl_proj_on_basis_advance(proj_do, 0.0, &local_do, fdo);

  // Project f(x0,y) directly onto the 1D basis.
  struct gkyl_array *fref = gkyl_array_new(GKYL_DOUBLE, basis_tar.num_basis, local_ext_tar.volume);
  gkyl_proj_on_basis *proj_tar = gkyl_proj_on_basis_new(&grid1, &basis_tar,
    poly_order+1, 1, eval_f_1d, &x0);
  gkyl_proj_on_basis_advance(proj_tar, 0.0, &local_tar, fref);

  // Create the updater: evaluate in x (dir 0), keep y (1D target).
  int eval_dirs[] = {0};
  struct gkyl_dg_eval_at_coord_proj *up = gkyl_dg_eval_at_coord_proj_new(
    &basis_do, &basis_tar, 1, eval_dirs, false);

  // Apply updater at x = x0.
  struct gkyl_array *ftar = gkyl_array_new(GKYL_DOUBLE, basis_tar.num_basis, local_ext_tar.volume);
  double eval_coords[] = {x0};
  bool pick_lower[] = {false, false};
  int known_index[] = {-1, -1};
  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid2, pick_lower, known_index,
    &local_do, &local_tar, fdo, ftar);

  // Check answer.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_tar);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local_tar, iter.idx);
    const double *ftar_c = gkyl_array_cfetch(ftar, linidx);
    const double *fref_c = gkyl_array_cfetch(fref, linidx);
    for (int k = 0; k < basis_tar.num_basis; k++)
      TEST_CHECK(gkyl_compare(ftar_c[k], fref_c[k], 1e-14));
  }

  gkyl_dg_eval_at_coord_proj_release(up);
  gkyl_proj_on_basis_release(proj_do);
  gkyl_proj_on_basis_release(proj_tar);
  gkyl_array_release(fdo);
  gkyl_array_release(ftar);
  gkyl_array_release(fref);
}

void test_2x_to_1x_p1(void) { test_2x_to_1x(1); }
void test_2x_to_1x_p2(void) { test_2x_to_1x(2); }

TEST_LIST = {
  { "test_2x_to_1x_p1", test_2x_to_1x_p1 },
  { "test_2x_to_1x_p2", test_2x_to_1x_p2 },
  { NULL, NULL },
};
