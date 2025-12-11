#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_dg_find_peaks.h>
#include <gkyl_array_dg_find_peaks_priv.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <math.h>
#include <stdio.h>

// 1D test function with multiple peaks: f(z) = cos(2*pi*z/L) 
// Has maxima at z=0, z=L and minimum at z=L/2.
static void
test_func_1d_cos(double t, const double *xn, double *fout, void *ctx)
{
  double z = xn[0];
  double L = 2.0;  // Period.
  fout[0] = cos(2.0 * M_PI * z / L);
}

// 1D test function that looks like a mirror bmag profile:
// f(z) = B0 * (1 + (R-1)*sin^2(pi*z/L))
// Has minimum at z=0, maxima at z=-L/2 and z=L/2 (mirror throats).
static void
test_func_1d_mirror(double t, const double *xn, double *fout, void *ctx)
{
  double z = xn[0];
  double L = 2.0;       // Half-length.
  double B0 = 1.0;      // Minimum B.
  double R = 4.0;       // Mirror ratio.
  double sinval = sin(M_PI * z / L);
  fout[0] = B0 * (1.0 + (R - 1.0) * sinval * sinval);
}

// 2D test function: f(psi, z) = (1 + 0.1*psi) * cos(2*pi*z/L)
// The peaks vary slightly with psi.
static void
test_func_2d_cos(double t, const double *xn, double *fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  double L = 2.0;
  fout[0] = (1.0 + 0.1*psi) * cos(2.0 * M_PI * z / L);
}

// 2D mirror-like function: peaks at z = +/- z_m(psi).
static void
test_func_2d_mirror(double t, const double *xn, double *fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  double L = 2.0;
  double B0 = 1.0 + 0.1*psi;  // Varies with psi.
  double R = 4.0;
  double sinval = sin(M_PI * z / L);
  fout[0] = B0 * (1.0 + (R - 1.0) * sinval * sinval);
}

// Complex test function: f(z) = (sin^2(z) + 0.1) * exp(-z^2/100)
// Has 9 peaks on [-5, 5]: edges at z=±5, local extrema at z≈±3π/2, ±π, ±π/2, 0.
static void
test_func_1d_complex(double t, const double *xn, double *fout, void *ctx)
{
  double z = xn[0];
  double sinz = sin(z);
  fout[0] = (sinz * sinz + 0.1) * exp(-z * z / 100.0);
}

// 2D complex test function: f(psi, z) = (sin^2(z) + 0.1) * exp(-z^2/100) * psi
// Peaks scale linearly with psi.
static void
test_func_2d_complex(double t, const double *xn, double *fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  double sinz = sin(z);
  fout[0] = (sinz * sinz + 0.1) * exp(-z * z / 100.0) * psi;
}

// Test function to project onto peaks: g(psi, z) = z^2 * psi^2
static void
test_func_quadratic_2d(double t, const double *xn, double *fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  fout[0] = z * z * psi * psi;
}

// 1D version: g(z) = z^2
static void
test_func_quadratic_1d(double t, const double *xn, double *fout, void *ctx)
{
  double z = xn[0];
  fout[0] = z * z;
}

// Test 1D peak finding with cos function.
void
test_1d_find_peaks_cos(int poly_order)
{
  // Grid: z in [-1, 1] (one period of cos(2*pi*z/2)).
  double lower[] = {-1.0};
  double upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  // Ranges.
  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_1d_cos, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder.
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 0,
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);

  // Compute peaks.
  gkyl_array_dg_find_peaks_advance(peaks, f);

  // Check results: cos(pi*z) on [-1,1] has EDGE_LO at z=-1, LOCAL_MAX at z=0, EDGE_HI at z=1.
  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 3);

  struct {
    enum gkyl_peak_type type;
    double z_expected;
  } expected_peaks[] = {
    {GKYL_PEAK_EDGE_LO,    -1.0,},
    {GKYL_PEAK_LOCAL_MAX,   0.0,},
    {GKYL_PEAK_EDGE_HI,     1.0,},
  };

  for (int p = 0; p < 3 && p < num_peaks; p++) {
    enum gkyl_peak_type ptype = gkyl_array_dg_find_peaks_get_type(peaks, p);
    const struct gkyl_array *vals = gkyl_array_dg_find_peaks_get_vals(peaks, p);
    const struct gkyl_array *coords = gkyl_array_dg_find_peaks_get_coords(peaks, p);
    
    const double *val = gkyl_array_cfetch(vals, 0);
    const double *coord = gkyl_array_cfetch(coords, 0);
    
    double z = coord[0];
    double expected_val[1];
    test_func_1d_cos(0.0, &z, expected_val, NULL);

    TEST_CHECK(ptype == expected_peaks[p].type);
    TEST_CHECK(fabs(coord[0] - expected_peaks[p].z_expected) < 0.1);
    TEST_CHECK(gkyl_compare_double(val[0], expected_val[0], 0.2));
  }

  gkyl_array_release(f);
  gkyl_array_dg_find_peaks_release(peaks);
}

// Test 1D peak finding with mirror-like function.
void
test_1d_find_peaks_mirror(int poly_order)
{
  // Grid: z in [-1, 1].
  double lower[] = {-1.0};
  double upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  // Ranges.
  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_1d_mirror, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder.
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 0,
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);

  // Compute peaks.
  gkyl_array_dg_find_peaks_advance(peaks, f);

  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 3);

  for (int p = 0; p < num_peaks; p++) {
    enum gkyl_peak_type ptype = gkyl_array_dg_find_peaks_get_type(peaks, p);
    const struct gkyl_array *vals = gkyl_array_dg_find_peaks_get_vals(peaks, p);
    const struct gkyl_array *coords = gkyl_array_dg_find_peaks_get_coords(peaks, p);
    
    const double *val = gkyl_array_cfetch(vals, 0);
    const double *coord = gkyl_array_cfetch(coords, 0);
    
    // Check specific peaks.
    if (ptype == GKYL_PEAK_EDGE_LO) {
      TEST_CHECK(gkyl_compare_double(val[0], 4.0, 1e-15));
      TEST_CHECK(fabs(coord[0] - (-1.0)) < 1e-15);
    }
    else if (ptype == GKYL_PEAK_LOCAL_MIN) {
      TEST_CHECK(gkyl_compare_double(val[0], 1.0, 1e-15));
      TEST_CHECK(fabs(coord[0]) < 1e-15);
    }
    else if (ptype == GKYL_PEAK_EDGE_HI) {
      TEST_CHECK(gkyl_compare_double(val[0], 4.0, 1e-15));
      TEST_CHECK(fabs(coord[0] - 1.0) < 1e-15);
    }
  }

  gkyl_array_release(f);
  gkyl_array_dg_find_peaks_release(peaks);
}

// Test 2D peak finding.
void
test_2d_find_peaks(int poly_order)
{
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

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_2d_mirror, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder (search along z, which is direction 1).
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 1,  // Search along z.
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);

  // Compute peaks.
  gkyl_array_dg_find_peaks_advance(peaks, f);

  // Check results: Mirror function should have 3 peaks along z.
  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 3);

  const struct gkyl_basis *out_basis = gkyl_array_dg_find_peaks_get_basis(peaks);
  const struct gkyl_range *out_range = gkyl_array_dg_find_peaks_get_range(peaks);
  const struct gkyl_rect_grid *out_grid = gkyl_array_dg_find_peaks_get_grid(peaks);

  // Check that values and coordinates are reasonable for each peak.
  for (int p = 0; p < num_peaks; p++) {
    enum gkyl_peak_type ptype = gkyl_array_dg_find_peaks_get_type(peaks, p);
    const struct gkyl_array *vals = gkyl_array_dg_find_peaks_get_vals(peaks, p);
    const struct gkyl_array *coords = gkyl_array_dg_find_peaks_get_coords(peaks, p);
    
    double xc_log[1] = {0.0};
    
    // Check first and last psi cells.
    for (int cell_idx = out_range->lower[0]; cell_idx <= out_range->upper[0]; 
         cell_idx += (out_range->upper[0] - out_range->lower[0])) {
      long linidx = gkyl_range_idx(out_range, (int[]){cell_idx});
      const double *val_d = gkyl_array_cfetch(vals, linidx);
      const double *coord_d = gkyl_array_cfetch(coords, linidx);

      double val_at_center = out_basis->eval_expand(xc_log, val_d);
      double coord_at_center = out_basis->eval_expand(xc_log, coord_d);
      double psi_phys = out_grid->lower[0] + (cell_idx - 0.5) * out_grid->dx[0];
      
      // Compute expected value at detected coordinate.
      double xn[2] = {psi_phys, coord_at_center};
      double expected_val[1];
      test_func_2d_mirror(0.0, xn, expected_val, NULL);
      
      // Check value matches analytical function.
      TEST_CHECK(gkyl_compare_double(val_at_center, expected_val[0], 1e-15));
      
      // Check that coordinate matches expected peak location.
      if (ptype == GKYL_PEAK_EDGE_LO) {
        TEST_CHECK(fabs(coord_at_center - (-1.0)) < 1e-15);
      } else if (ptype == GKYL_PEAK_LOCAL_MIN) {
        TEST_CHECK(fabs(coord_at_center) < 1e-15);
      } else if (ptype == GKYL_PEAK_EDGE_HI) {
        TEST_CHECK(fabs(coord_at_center - 1.0) < 1e-15);
      }
    }
  }

  gkyl_array_release(f);
  gkyl_array_dg_find_peaks_release(peaks);
}

// Test 1D peak finding with complex oscillatory function.
void
test_1d_find_peaks_complex(int poly_order)
{
  double lower[] = {-2.0*M_PI};
  double upper[] = {2.0*M_PI};
  int cells[] = {64};  // Need fine resolution to capture oscillations.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_1d_complex, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder.
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 0,
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);

  // Compute peaks.
  gkyl_array_dg_find_peaks_advance(peaks, f);

  // Check results.
  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  
  TEST_CHECK(num_peaks == 9);

  // Define expected peak locations and types.
  struct {
    enum gkyl_peak_type type;
    double z_expected;
  } expected_peaks[] = {
    {GKYL_PEAK_EDGE_LO,    -2.0*M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -3.0*M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,  -M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   0.0,},
    {GKYL_PEAK_LOCAL_MAX,   M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   M_PI,},
    {GKYL_PEAK_LOCAL_MAX,   3.0*M_PI/2.0,},
    {GKYL_PEAK_EDGE_HI,     2.0*M_PI,},
  };

  for (int p = 0; p < num_peaks; p++) {
    enum gkyl_peak_type ptype = gkyl_array_dg_find_peaks_get_type(peaks, p);
    const struct gkyl_array *vals = gkyl_array_dg_find_peaks_get_vals(peaks, p);
    const struct gkyl_array *coords = gkyl_array_dg_find_peaks_get_coords(peaks, p);
    
    const double *val = gkyl_array_cfetch(vals, 0);
    const double *coord = gkyl_array_cfetch(coords, 0);
    
    double z = coord[0];
    double expected_val[1];
    test_func_1d_complex(0.0, &z, expected_val, NULL);

    TEST_CHECK(ptype == expected_peaks[p].type);
    TEST_CHECK(fabs(coord[0] - expected_peaks[p].z_expected) < 1e-15);
    double rel_error = fabs(val[0] - expected_val[0]) / fabs(expected_val[0]);
    TEST_CHECK(rel_error < 1e-15);
  }

  gkyl_array_release(f);
  gkyl_array_dg_find_peaks_release(peaks);
}

// Test 2D peak finding with complex oscillatory function.
void
test_2d_find_peaks_complex(int poly_order)
{
  // Grid: psi in [0.5, 2.0], z in [-5, 5].
  double lower[] = {0.5, -2.0*M_PI};
  double upper[] = {2.0, 2.0*M_PI};
  int cells[] = {16, 64};
  int ndim = 2;
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Ranges.
  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function onto basis.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_2d_complex, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder (search along z, which is direction 1).
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 1,  // Search along z.
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);

  // Compute peaks.
  gkyl_array_dg_find_peaks_advance(peaks, f);

  // Check results.
  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 9);

  const struct gkyl_basis *out_basis = gkyl_array_dg_find_peaks_get_basis(peaks);
  const struct gkyl_range *out_range = gkyl_array_dg_find_peaks_get_range(peaks);
  const struct gkyl_rect_grid *out_grid = gkyl_array_dg_find_peaks_get_grid(peaks);

  // Define expected peak locations and types (same as 1D).
  struct {
    enum gkyl_peak_type type;
    double z_expected;
  } expected_peaks[] = {
    {GKYL_PEAK_EDGE_LO,    -2.0*M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -3.0*M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,  -M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   0.0,},
    {GKYL_PEAK_LOCAL_MAX,   M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   M_PI,},
    {GKYL_PEAK_LOCAL_MAX,   3.0*M_PI/2.0,},
    {GKYL_PEAK_EDGE_HI,     2.0*M_PI,},
  };

  // Get node locations for output basis.
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, out_basis->ndim, out_basis->num_basis);
  out_basis->node_list(gkyl_array_fetch(nodes, 0));

  // Check peak types and validate values at each psi cell.
  for (int p = 0; p < num_peaks; p++) {
    enum gkyl_peak_type ptype = gkyl_array_dg_find_peaks_get_type(peaks, p);
    TEST_CHECK(ptype == expected_peaks[p].type);
    
    const struct gkyl_array *vals = gkyl_array_dg_find_peaks_get_vals(peaks, p);
    const struct gkyl_array *coords = gkyl_array_dg_find_peaks_get_coords(peaks, p);
    
    // Check each psi cell.
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, out_range);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(out_range, iter.idx);
      
      const double *val_d = gkyl_array_cfetch(vals, linidx);
      const double *coord_d = gkyl_array_cfetch(coords, linidx);
      
      // Get cell center for physical psi coordinate.
      double xc_out[1];
      gkyl_rect_grid_cell_center(out_grid, (int[]){iter.idx[0]}, xc_out);
      double psi_phys = xc_out[0];
      
      // Evaluate at each nodal point in this cell.
      for (int n = 0; n < out_basis->num_basis; n++) {
        const double *nod_log = gkyl_array_cfetch(nodes, n);
        double val_at_node = out_basis->eval_expand(nod_log, val_d);
        double z_at_node = out_basis->eval_expand(nod_log, coord_d);
        
        // Compute physical psi coordinate at this node.
        // dx/2 away from the center is the nodal location.
        double nod_phys[1];
        nod_phys[0] = xc_out[0] + nod_log[0] * out_grid->dx[0]/2.0;
        
        // Compute expected value at detected coordinates.
        double xn[2] = {nod_phys[0], z_at_node};
        double expected_val[1];
        test_func_2d_complex(0.0, xn, expected_val, NULL);

        TEST_CHECK(fabs(z_at_node - expected_peaks[p].z_expected) < 1e-15);
        double rel_error = fabs(val_at_node - expected_val[0]) / fabs(expected_val[0]);
        TEST_CHECK(rel_error < 1e-15);
      }
    }
  }

  gkyl_array_release(nodes);
  gkyl_array_release(f);
  gkyl_array_dg_find_peaks_release(peaks);
}

// Test 1D project_on_peaks with complex function.
void
test_1d_project_on_peaks(int poly_order)
{
  double lower[] = {-2.0*M_PI};
  double upper[] = {2.0*M_PI};
  int cells[] = {64};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function for peak finding.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_1d_complex, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Project quadratic function to evaluate at peaks.
  struct gkyl_array *g = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_quadratic_1d, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, g);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder.
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 0,
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);
  gkyl_array_dg_find_peaks_advance(peaks, f);

  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 9);

  // Allocate output arrays for projected values.
  const struct gkyl_range *out_range_ext = gkyl_array_dg_find_peaks_get_range_ext(peaks);
  const struct gkyl_basis *out_basis = gkyl_array_dg_find_peaks_get_basis(peaks);
  struct gkyl_array *g_at_peaks[GKYL_DG_FIND_PEAKS_MAX];
  for (int p = 0; p < num_peaks; p++) {
    g_at_peaks[p] = gkyl_array_new(GKYL_DOUBLE, out_basis->num_basis, out_range_ext->volume);
  }
  gkyl_array_dg_find_peaks_project_on_peaks(peaks, g, g_at_peaks);

  // Verify that g evaluated at each peak matches analytical values.
  // For 1D->0D, output is p=0, so the value is already the cell average.
  // The cell average of a p=0 expansion is value / sqrt(volume), and
  // for a 1D cell with volume=1, it's just the value / sqrt(1.0) = value.
  struct {
    enum gkyl_peak_type type;
    double z_expected;
  } expected_peaks[] = {
    {GKYL_PEAK_EDGE_LO,    -2.0*M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -3.0*M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,  -M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   0.0,},
    {GKYL_PEAK_LOCAL_MAX,   M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   M_PI,},
    {GKYL_PEAK_LOCAL_MAX,   3.0*M_PI/2.0,},
    {GKYL_PEAK_EDGE_HI,     2.0*M_PI,},
  };
  for (int p = 0; p < num_peaks; p++) {
    const double *g_val = gkyl_array_cfetch(g_at_peaks[p], 0);
    double z = expected_peaks[p].z_expected;
    double expected = z * z;
    TEST_CHECK(gkyl_compare_double(g_val[0], expected, 1e-12));
  }

  for (int p = 0; p < num_peaks; p++) {
    gkyl_array_release(g_at_peaks[p]);
  }
  gkyl_array_release(f);
  gkyl_array_release(g);
  gkyl_array_dg_find_peaks_release(peaks);
}

// Test 2D project_on_peaks with complex function.
void
test_2d_project_on_peaks(int poly_order)
{
  double lower[] = {0.5, -2.0*M_PI};
  double upper[] = {2.0, 2.0*M_PI};
  int cells[] = {16, 64};
  int ndim = 2;
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function for peak finding.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_2d_complex, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Project quadratic function to evaluate at peaks: g(psi, z) = z^2 * psi^2
  struct gkyl_array *g = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_quadratic_2d, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, g);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder (search along z, which is direction 1).
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 1,  // Search along z.
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);
  gkyl_array_dg_find_peaks_advance(peaks, f);

  // Check results.
  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 9);

  const struct gkyl_basis *out_basis = gkyl_array_dg_find_peaks_get_basis(peaks);
  const struct gkyl_range *out_range = gkyl_array_dg_find_peaks_get_range(peaks);
  const struct gkyl_range *out_range_ext = gkyl_array_dg_find_peaks_get_range_ext(peaks);
  const struct gkyl_rect_grid *out_grid = gkyl_array_dg_find_peaks_get_grid(peaks);

  // Allocate output arrays for projected values.
  struct gkyl_array *g_at_peaks[GKYL_DG_FIND_PEAKS_MAX];
  for (int p = 0; p < num_peaks; p++) {
    g_at_peaks[p] = gkyl_array_new(GKYL_DOUBLE, out_basis->num_basis, out_range_ext->volume);
  }
  gkyl_array_dg_find_peaks_project_on_peaks(peaks, g, g_at_peaks);

  // Define expected peak locations (same as before).
  double expected_z_peaks[] = {
    -2.0*M_PI, -3.0*M_PI/2.0, -M_PI, -M_PI/2.0, 0.0,
    M_PI/2.0, M_PI, 3.0*M_PI/2.0, 2.0*M_PI
  };

  // Get node locations for output basis.
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, out_basis->ndim, out_basis->num_basis);
  out_basis->node_list(gkyl_array_fetch(nodes, 0));

  // Verify that g evaluated at each peak matches analytical values.
  for (int p = 0; p < num_peaks; p++) {
    const struct gkyl_array *coords = gkyl_array_dg_find_peaks_get_coords(peaks, p);
    
    // Check each psi cell.
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, out_range);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(out_range, iter.idx);
      
      const double *g_val_d = gkyl_array_cfetch(g_at_peaks[p], linidx);
      const double *coord_d = gkyl_array_cfetch(coords, linidx);
      
      // Get cell center for physical psi coordinate.
      double xc_out[1];
      gkyl_rect_grid_cell_center(out_grid, (int[]){iter.idx[0]}, xc_out);
      
      // Evaluate at each nodal point in this cell.
      for (int n = 0; n < out_basis->num_basis; n++) {
        const double *nod_log = gkyl_array_cfetch(nodes, n);
        double g_at_node = out_basis->eval_expand(nod_log, g_val_d);
        double z_at_node = out_basis->eval_expand(nod_log, coord_d);
        
        // Compute physical psi coordinate at this node.
        double nod_phys[1];
        nod_phys[0] = xc_out[0] + nod_log[0] * out_grid->dx[0]/2.0;
        double psi = nod_phys[0];
        
        // Analytical value: g(psi, z) = z^2 * psi^2
        double expected = z_at_node * z_at_node * psi * psi;
        TEST_CHECK(fabs(z_at_node - expected_z_peaks[p]) < 1e-15);
        TEST_CHECK(gkyl_compare_double(g_at_node, expected, 1e-15));
      }
    }
  }

  // Clean up.
  gkyl_array_release(nodes);
  for (int p = 0; p < num_peaks; p++) {
    gkyl_array_release(g_at_peaks[p]);
  }
  gkyl_array_release(f);
  gkyl_array_release(g);
  gkyl_array_dg_find_peaks_release(peaks);
}


// Test 1D project_on_peak_idx with complex function.
void
test_1d_project_on_peak_idx(int poly_order)
{
  double lower[] = {-2.0*M_PI};
  double upper[] = {2.0*M_PI};
  int cells[] = {64};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function for peak finding.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_1d_complex, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Project quadratic function to evaluate at peaks.
  struct gkyl_array *g = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_quadratic_1d, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, g);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder.
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 0,
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);
  gkyl_array_dg_find_peaks_advance(peaks, f);

  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 9);

  // Allocate output arrays for projected values.
  const struct gkyl_range *out_range_ext = gkyl_array_dg_find_peaks_get_range_ext(peaks);
  const struct gkyl_basis *out_basis = gkyl_array_dg_find_peaks_get_basis(peaks);
  struct gkyl_array *g_at_peaks = gkyl_array_new(GKYL_DOUBLE, out_basis->num_basis, out_range_ext->volume);

  int chosen_idx = 1;
  gkyl_array_dg_find_peaks_project_on_peak_idx(peaks, g, chosen_idx, g_at_peaks);

  // Verify that g evaluated at each peak matches analytical values.
  // For 1D->0D, output is p=0, so the value is already the cell average.
  // The cell average of a p=0 expansion is value / sqrt(volume), and
  // for a 1D cell with volume=1, it's just the value / sqrt(1.0) = value.
  struct {
    enum gkyl_peak_type type;
    double z_expected;
  } expected_peaks[] = {
    {GKYL_PEAK_EDGE_LO,    -2.0*M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -3.0*M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,  -M_PI,},
    {GKYL_PEAK_LOCAL_MAX,  -M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   0.0,},
    {GKYL_PEAK_LOCAL_MAX,   M_PI/2.0,},
    {GKYL_PEAK_LOCAL_MIN,   M_PI,},
    {GKYL_PEAK_LOCAL_MAX,   3.0*M_PI/2.0,},
    {GKYL_PEAK_EDGE_HI,     2.0*M_PI,},
  };
  
  const double *g_val = gkyl_array_cfetch(g_at_peaks, 0);
  double z = expected_peaks[chosen_idx].z_expected;
  double expected = z * z;
  TEST_CHECK(gkyl_compare_double(g_val[0], expected, 1e-12));

  gkyl_array_release(g_at_peaks);
  gkyl_array_release(f);
  gkyl_array_release(g);
  gkyl_array_dg_find_peaks_release(peaks);
}

// Test 2D project_on_peak_idx with complex function.
void
test_2d_project_on_peak_idx(int poly_order)
{
  double lower[] = {0.5, -2.0*M_PI};
  double upper[] = {2.0, 2.0*M_PI};
  int cells[] = {16, 64};
  int ndim = 2;
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project test function for peak finding.
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_2d_complex, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, f);
  gkyl_eval_on_nodes_release(ev);

  // Project quadratic function to evaluate at peaks: g(psi, z) = z^2 * psi^2
  struct gkyl_array *g = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  ev = gkyl_eval_on_nodes_new(&grid, &basis, 1, test_func_quadratic_2d, NULL);
  gkyl_eval_on_nodes_advance(ev, 0.0, &local, g);
  gkyl_eval_on_nodes_release(ev);

  // Create peak finder (search along z, which is direction 1).
  struct gkyl_array_dg_find_peaks_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &local_ext,
    .search_dir = 1,  // Search along z.
    .use_gpu = false,
  };
  struct gkyl_array_dg_find_peaks *peaks = gkyl_array_dg_find_peaks_new(&inp, f);
  gkyl_array_dg_find_peaks_advance(peaks, f);

  // Check results.
  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peaks);
  TEST_CHECK(num_peaks == 9);

  const struct gkyl_basis *out_basis = gkyl_array_dg_find_peaks_get_basis(peaks);
  const struct gkyl_range *out_range = gkyl_array_dg_find_peaks_get_range(peaks);
  const struct gkyl_range *out_range_ext = gkyl_array_dg_find_peaks_get_range_ext(peaks);
  const struct gkyl_rect_grid *out_grid = gkyl_array_dg_find_peaks_get_grid(peaks);

  // Allocate output arrays for projected values.
  struct gkyl_array *g_at_peaks[GKYL_DG_FIND_PEAKS_MAX];
  for (int p = 0; p < num_peaks; p++) {
    g_at_peaks[p] = gkyl_array_new(GKYL_DOUBLE, out_basis->num_basis, out_range_ext->volume);
  }
  gkyl_array_dg_find_peaks_project_on_peaks(peaks, g, g_at_peaks);

  // Define expected peak locations (same as before).
  double expected_z_peaks[] = {
    -2.0*M_PI, -3.0*M_PI/2.0, -M_PI, -M_PI/2.0, 0.0,
    M_PI/2.0, M_PI, 3.0*M_PI/2.0, 2.0*M_PI
  };

  // Get node locations for output basis.
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, out_basis->ndim, out_basis->num_basis);
  out_basis->node_list(gkyl_array_fetch(nodes, 0));

  // Verify that g evaluated at each peak matches analytical values.
  for (int p = 0; p < num_peaks; p++) {
    const struct gkyl_array *coords = gkyl_array_dg_find_peaks_get_coords(peaks, p);
    
    // Check each psi cell.
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, out_range);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(out_range, iter.idx);
      
      const double *g_val_d = gkyl_array_cfetch(g_at_peaks[p], linidx);
      const double *coord_d = gkyl_array_cfetch(coords, linidx);
      
      // Get cell center for physical psi coordinate.
      double xc_out[1];
      gkyl_rect_grid_cell_center(out_grid, (int[]){iter.idx[0]}, xc_out);
      
      // Evaluate at each nodal point in this cell.
      for (int n = 0; n < out_basis->num_basis; n++) {
        const double *nod_log = gkyl_array_cfetch(nodes, n);
        double g_at_node = out_basis->eval_expand(nod_log, g_val_d);
        double z_at_node = out_basis->eval_expand(nod_log, coord_d);
        
        // Compute physical psi coordinate at this node.
        double nod_phys[1];
        nod_phys[0] = xc_out[0] + nod_log[0] * out_grid->dx[0]/2.0;
        double psi = nod_phys[0];
        
        // Analytical value: g(psi, z) = z^2 * psi^2
        double expected = z_at_node * z_at_node * psi * psi;
        TEST_CHECK(fabs(z_at_node - expected_z_peaks[p]) < 1e-15);
        TEST_CHECK(gkyl_compare_double(g_at_node, expected, 1e-15));
      }
    }
  }

  // Clean up.
  gkyl_array_release(nodes);
  for (int p = 0; p < num_peaks; p++) {
    gkyl_array_release(g_at_peaks[p]);
  }
  gkyl_array_release(f);
  gkyl_array_release(g);
  gkyl_array_dg_find_peaks_release(peaks);
}












void test_1d_cos_p1() { test_1d_find_peaks_cos(1); }
void test_1d_mirror_p1() { test_1d_find_peaks_mirror(1); }
void test_1d_complex_p1() { test_1d_find_peaks_complex(1); }
void test_2d_p1() { test_2d_find_peaks(1); }
void test_2d_complex_p1() { test_2d_find_peaks_complex(1); }
void test_1d_project_p1() { test_1d_project_on_peaks(1); }
void test_2d_project_p1() { test_2d_project_on_peaks(1); }
void test_1d_project_idx_p1() { test_1d_project_on_peak_idx(1); }
void test_2d_project_idx_p1() { test_2d_project_on_peak_idx(1); }

TEST_LIST = {
  {"test_1d_cos_p1", test_1d_cos_p1},
  {"test_1d_mirror_p1", test_1d_mirror_p1},
  {"test_1d_complex_p1", test_1d_complex_p1},
  {"test_2d_p1", test_2d_p1},
  {"test_2d_complex_p1", test_2d_complex_p1},
  {"test_1d_project_p1", test_1d_project_p1},
  {"test_2d_project_p1", test_2d_project_p1},
  {"test_1d_project_idx_p1", test_1d_project_idx_p1},
  // {"test_2d_project_idx_p1", test_2d_project_idx_p1},
  {NULL, NULL},
};
