#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_dg_reduce_dir.h>
#include <gkyl_array_dg_reduce_dir_priv.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <math.h>

// 1D test function: f(z) = -z^2 + 1
// Maximum at z=0 with value 1, minimum at endpoints.
static void
test_func_1d_parabola(double t, const double *xn, double *fout, void *ctx)
{
  double z = xn[0];
  fout[0] = -z*z + 1.0;
}

// 1D test function: f(z) = sin(z) on [0, pi]
// Maximum at z=pi/2 with value 1.
static void
test_func_1d_sin(double t, const double *xn, double *fout, void *ctx)
{
  double z = xn[0];
  fout[0] = sin(z);
}

// 2D test function: f(psi, z) = psi * (-z^2 + 1)
// Maximum along z is at z=0 for each psi, with value psi.
static void
test_func_2d_parabola(double t, const double *xn, double *fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  fout[0] = psi * (-z*z + 1.0);
}

// 2D mirror-like function: f(psi, z) = B0(psi) * (1 + (R-1)*sin^2(pi*z/L))
// Minimum at z=0 for all psi.
static void
test_func_2d_mirror(double t, const double *xn, double *fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  double L = 2.0;
  double B0 = 1.0 + 0.1*psi;
  double R = 4.0;
  double sinval = sin(M_PI * z / L);
  fout[0] = B0 * (1.0 + (R - 1.0) * sinval * sinval);
}

// Test 1D reduction with MAX operation.
void
test_1d_reduce_max(int poly_order)
{
  // Grid: z in [-1, 1].
  double lower[] = {-1.0};
  double upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_1d_parabola, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create reducer.
  struct gkyl_array_dg_reduce_dir_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .reduce_dir = 0,
    .op = GKYL_REDUCE_OP_MAX,
    .use_gpu = false,
  };
  struct gkyl_array_dg_reduce_dir *reducer = gkyl_array_dg_reduce_dir_new(&inp);

  // Compute reduction.
  gkyl_array_dg_reduce_dir_advance(reducer, f);

  // Check results: maximum of -z^2+1 on [-1,1] is 1 at z=0.
  const struct gkyl_array *vals = gkyl_array_dg_reduce_dir_get_vals(reducer);
  const struct gkyl_array *coords = gkyl_array_dg_reduce_dir_get_coords(reducer);
  
  const double *val = gkyl_array_cfetch(vals, 0);
  const double *coord = gkyl_array_cfetch(coords, 0);
  
  TEST_CHECK(gkyl_compare_double(val[0], 1.0, 1e-14));
  TEST_CHECK(gkyl_compare_double(fabs(coord[0]), 0.0, 1e-14));

  gkyl_array_release(f);
  gkyl_array_dg_reduce_dir_release(reducer);
}

// Test 1D reduction with MIN operation.
void
test_1d_reduce_min(int poly_order)
{
  // Grid: z in [-1, 1].
  double lower[] = {-1.0};
  double upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_1d_parabola, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create reducer.
  struct gkyl_array_dg_reduce_dir_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .reduce_dir = 0,
    .op = GKYL_REDUCE_OP_MIN,
    .use_gpu = false,
  };
  struct gkyl_array_dg_reduce_dir *reducer = gkyl_array_dg_reduce_dir_new(&inp);

  // Compute reduction.
  gkyl_array_dg_reduce_dir_advance(reducer, f);

  // Check results: minimum of -z^2+1 on [-1,1] is 0 at z=±1.
  const struct gkyl_array *vals = gkyl_array_dg_reduce_dir_get_vals(reducer);
  const struct gkyl_array *coords = gkyl_array_dg_reduce_dir_get_coords(reducer);
  
  const double *val = gkyl_array_cfetch(vals, 0);
  const double *coord = gkyl_array_cfetch(coords, 0);
  
  TEST_CHECK(gkyl_compare_double(val[0], 0.0, 1e-14));
  TEST_CHECK(gkyl_compare_double(fabs(coord[0]), 1.0, 1e-14));  // Either -1 or 1.

  gkyl_array_release(f);
  gkyl_array_dg_reduce_dir_release(reducer);
}

// Test 2D reduction along z direction with MAX operation.
void
test_2d_reduce_max(int poly_order)
{
  // Grid: psi in [0.5, 2.0], z in [-1, 1].
  double lower[] = {0.5, -1.0};
  double upper[] = {2.0, 1.0};
  int cells[] = {8, 16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_2d_parabola, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create reducer (reduce along z, which is direction 1).
  struct gkyl_array_dg_reduce_dir_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .reduce_dir = 1,
    .op = GKYL_REDUCE_OP_MAX,
    .use_gpu = false,
  };
  struct gkyl_array_dg_reduce_dir *reducer = gkyl_array_dg_reduce_dir_new(&inp);

  // Compute reduction.
  gkyl_array_dg_reduce_dir_advance(reducer, f);

  // Check results.
  const struct gkyl_basis *out_basis = gkyl_array_dg_reduce_dir_get_basis(reducer);
  const struct gkyl_range *out_range = gkyl_array_dg_reduce_dir_get_range(reducer);
  const struct gkyl_rect_grid *out_grid = gkyl_array_dg_reduce_dir_get_grid(reducer);
  
  // Access nodal arrays directly for testing.
  const struct gkyl_array *vals_nodal = reducer->out_vals_nodal;
  const struct gkyl_array *coords_nodal = reducer->out_coords_nodal;
  const struct gkyl_range *out_nrange = &reducer->out_nrange;

  // Check each nodal point.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, out_nrange);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(out_nrange, iter.idx);
    
    const double *val_nodal = gkyl_array_cfetch(vals_nodal, linidx);
    const double *coord_nodal = gkyl_array_cfetch(coords_nodal, linidx);
    
    // Compute physical psi coordinate at this nodal point.
    // For p=1: node 0 at lower bound, node ncells at upper bound.
    int num_cells_psi = cells[0];
    double dpsi = (upper[0] - lower[0]) / num_cells_psi;
    double psi_phys = lower[0] + iter.idx[0] * dpsi;
    
    // Expected: max of psi*(-z^2+1) over z is psi at z=0.
    double expected_val = psi_phys;
    double expected_coord = 0.0;

    TEST_CHECK(gkyl_compare_double(val_nodal[0], expected_val, 1e-14));
    TEST_CHECK(gkyl_compare_double(coord_nodal[0], expected_coord, 1e-14));
  }

  gkyl_array_release(f);
  gkyl_array_dg_reduce_dir_release(reducer);
}

// Test 2D reduction with mirror-like function (finding minimum bmag).
void
test_2d_reduce_min_mirror(int poly_order)
{
  // Grid: psi in [0.0, 1.0], z in [-1, 1].
  double lower[] = {0.0, -1.0};
  double upper[] = {1.0, 1.0};
  int cells[] = {4, 16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project mirror function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_2d_mirror, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create reducer (reduce along z, find minimum).
  struct gkyl_array_dg_reduce_dir_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .reduce_dir = 1,
    .op = GKYL_REDUCE_OP_MIN,
    .use_gpu = false,
  };
  struct gkyl_array_dg_reduce_dir *reducer = gkyl_array_dg_reduce_dir_new(&inp);

  // Compute reduction.
  gkyl_array_dg_reduce_dir_advance(reducer, f);

  // Access nodal arrays directly for testing.
  const struct gkyl_array *vals_nodal = reducer->out_vals_nodal;
  const struct gkyl_array *coords_nodal = reducer->out_coords_nodal;
  const struct gkyl_range *out_nrange = &reducer->out_nrange;

  // Check each nodal point.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, out_nrange);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(out_nrange, iter.idx);
    
    const double *val_nodal = gkyl_array_cfetch(vals_nodal, linidx);
    const double *coord_nodal = gkyl_array_cfetch(coords_nodal, linidx);
    
    // Compute physical psi coordinate.
    int num_cells_psi = cells[0];
    double dpsi = (upper[0] - lower[0]) / num_cells_psi;
    double psi_phys = lower[0] + iter.idx[0] * dpsi;
    
    // Expected: minimum of B0*(1+(R-1)*sin^2(pi*z/L)) is B0 at z=0.
    double B0 = 1.0 + 0.1*psi_phys;
    double expected_val = B0;
    double expected_coord = 0.0;

    TEST_CHECK(gkyl_compare_double(val_nodal[0], expected_val, 1e-14));
    TEST_CHECK(gkyl_compare_double(coord_nodal[0], expected_coord, 1e-14));
  }

  gkyl_array_release(f);
  gkyl_array_dg_reduce_dir_release(reducer);
}

// Test eval_at_extremum functionality.
void
test_2d_eval_at_extremum(int poly_order)
{
  // Grid: psi in [0.5, 2.0], z in [-1, 1].
  double lower[] = {0.5, -1.0};
  double upper[] = {2.0, 1.0};
  int cells[] = {8, 16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_2d_parabola, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create reducer.
  struct gkyl_array_dg_reduce_dir_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .reduce_dir = 1,
    .op = GKYL_REDUCE_OP_MAX,
    .use_gpu = false,
  };
  struct gkyl_array_dg_reduce_dir *reducer = gkyl_array_dg_reduce_dir_new(&inp);

  // Compute reduction.
  gkyl_array_dg_reduce_dir_advance(reducer, f);

  // Evaluate f at the extremum coordinates.
  const struct gkyl_range *out_range_ext = gkyl_array_dg_reduce_dir_get_range_ext(reducer);
  const struct gkyl_basis *out_basis = gkyl_array_dg_reduce_dir_get_basis(reducer);
  
  struct gkyl_array *f_at_max = gkyl_array_new(GKYL_DOUBLE, out_basis->num_basis, out_range_ext->volume);
  gkyl_array_dg_reduce_dir_eval_at_extremum(reducer, f, f_at_max);

  // The value of f at its maximum should equal the maximum value.
  const struct gkyl_array *vals = gkyl_array_dg_reduce_dir_get_vals(reducer);
  
  // Compare at cell centers.
  double xc_log[1] = {0.0};
  const struct gkyl_range *out_range = gkyl_array_dg_reduce_dir_get_range(reducer);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, out_range);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(out_range, iter.idx);
    
    const double *max_val = gkyl_array_cfetch(vals, linidx);
    const double *eval_val = gkyl_array_cfetch(f_at_max, linidx);
    
    double max_at_center = out_basis->eval_expand(xc_log, max_val);
    double eval_at_center = out_basis->eval_expand(xc_log, eval_val);

    TEST_CHECK(gkyl_compare_double(eval_at_center, max_at_center, 1e-14));
  }

  gkyl_array_release(f);
  gkyl_array_release(f_at_max);
  gkyl_array_dg_reduce_dir_release(reducer);
}

void test_1d_max_p1() { test_1d_reduce_max(1); }
void test_1d_min_p1() { test_1d_reduce_min(1); }
void test_2d_max_p1() { test_2d_reduce_max(1); }
void test_2d_min_mirror_p1() { test_2d_reduce_min_mirror(1); }
void test_2d_eval_at_extremum_p1() { test_2d_eval_at_extremum(1); }

TEST_LIST = {
  {"test_1d_max_p1", test_1d_max_p1},
  {"test_1d_min_p1", test_1d_min_p1},
  {"test_2d_max_p1", test_2d_max_p1},
  {"test_2d_min_mirror_p1", test_2d_min_mirror_p1},
  {"test_2d_eval_at_extremum_p1", test_2d_eval_at_extremum_p1},
  {NULL, NULL},
};
