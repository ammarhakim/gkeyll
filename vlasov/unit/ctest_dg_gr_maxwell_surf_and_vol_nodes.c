/*
 * Unit tests for dg_gr_maxwell_surf_and_vol_nodes.
 *
 * What is being tested:
 * 1. The private node builders in gkyl_dg_gr_maxwell_surf_and_vol_nodes_priv.h.
 *    These routines generate the Gauss-Legendre volume nodes and the lower-face
 *    surface nodes used by the updater. The test checks their coordinates
 *    directly against the expected tensor-product construction in 1x, 2x, and 3x.
 *
 * 2. The public updater construction path with the identity coordinate mapping.
 *    The tests build the updater in 1x, 2x, and 3x, then verify that
 *    gkyl_dg_gr_maxwell_surf_and_vol_nodes_fetch_node returns the same basis
 *    nodes supplied by basis.node_list().
 *
 * 3. The nodal array allocation layout produced by
 *    gkyl_surf_and_vol_node_arrays_new. For each dimension we check the expected
 *    number of stored components for the volume and surface arrays.
 *
 * 4. The main advance routine,
 *    gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance. The test uses a known
 *    analytic function with three return values, evaluates it on every volume
 *    node and every lower-face surface node, and compares the stored values to
 *    the exact analytic result in each cell.
 *
 * Analytic strategy:
 * - The updater evaluates the callback at physical coordinates obtained from the
 *   cell center plus the local node location scaled by dx/2.
 * - Because this test is intentionally identity-only, the physical coordinate is
 *   just that computational coordinate.
 * - For each cell and each node, the test reconstructs the exact evaluation
 *   point, computes the expected analytic values, and checks the flattened nodal
 *   storage entry-by-entry in component-major order.
 *
 * Coverage by test case:
 * - private_node_builders:
 *   Verifies the raw node-generation helpers in 1x, 2x, and 3x.
 * - surf_and_vol_nodes_1x_identity:
 *   Checks volume nodes and x-surface nodes in 1x.
 * - surf_and_vol_nodes_2x_identity:
 *   Checks volume nodes and x/y lower-face surface nodes in 2x.
 * - surf_and_vol_nodes_3x_identity:
 *   Checks volume nodes and x/y/z lower-face surface nodes in 3x.
 *
 * This test does not exercise mapped coordinates on purpose. Its role is to
 * validate the current identity-coordinate implementation and the exact nodal
 * values stored by the updater.
 */

#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes_priv.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <stdbool.h>


struct eval_ctx {
  int ndim;
};

static void
calc_expected_eval(int ndim, double tm, const double *x, double *fout)
{
  double x0 = ndim > 0 ? x[0] : 0.0;
  double x1 = ndim > 1 ? x[1] : 0.0;
  double x2 = ndim > 2 ? x[2] : 0.0;

  fout[0] = 0.5 + 0.25*tm + 1.1*x0 - 0.7*x1 + 0.9*x2;
  fout[1] = -0.2 + x0*x0 + 0.5*x0*x1 + 0.25*x1*x1 + 0.75*x2;
  fout[2] = 1.0 - 0.5*tm + (x0 - 0.25)*(ndim > 1 ? (x1 + 0.5) : 1.0)
    + (ndim > 2 ? 0.2*x2*x2 : 0.0);
}

static void
analytic_eval(double tm, const double *x, double *fout, void *ctx)
{
  struct eval_ctx *ev = ctx;
  calc_expected_eval(ev->ndim, tm, x, fout);
}

static void
expected_vol_node(int ndim, int polyorder, int node, double *eta)
{
  const double *gl_nodes = gkyl_gauss_ordinates[polyorder + 1];
  int nn = polyorder + 1;

  switch (ndim) {
    case 1:
      eta[0] = gl_nodes[node];
      break;
    case 2:
      eta[0] = gl_nodes[node/nn];
      eta[1] = gl_nodes[node%nn];
      break;
    case 3:
      eta[0] = gl_nodes[node/(nn*nn)];
      eta[1] = gl_nodes[(node/nn)%nn];
      eta[2] = gl_nodes[node%nn];
      break;
  }
}

static void
expected_surf_node(int ndim, int dir, int polyorder, int node, double *eta)
{
  const double *gl_nodes = gkyl_gauss_ordinates[polyorder + 1];
  int nn = polyorder + 1;

  switch (ndim) {
    case 1:
      eta[0] = -1.0;
      break;
    case 2:
      if (dir == 0) {
        eta[0] = -1.0;
        eta[1] = gl_nodes[node];
      }
      else {
        eta[0] = gl_nodes[node];
        eta[1] = -1.0;
      }
      break;
    case 3:
      if (dir == 0) {
        eta[0] = -1.0;
        eta[1] = gl_nodes[node/nn];
        eta[2] = gl_nodes[node%nn];
      }
      else if (dir == 1) {
        eta[0] = gl_nodes[node/nn];
        eta[1] = -1.0;
        eta[2] = gl_nodes[node%nn];
      }
      else {
        eta[0] = gl_nodes[node/nn];
        eta[1] = gl_nodes[node%nn];
        eta[2] = -1.0;
      }
      break;
  }
}

static void
assert_close(double expected, double actual, double tol, const char *kind,
  const int *idx, int node, int component)
{
  bool ok = gkyl_compare(expected, actual, tol);
  TEST_CHECK(ok);
  if (!ok) {
    TEST_MSG("%s mismatch at cell (%d,%d,%d), node %d, component %d",
      kind, idx[0], idx[1], idx[2], node, component);
    TEST_MSG("Expected %.15e", expected);
    TEST_MSG("Produced %.15e", actual);
  }
}

static void
assert_node_array_matches(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *nodes, const struct gkyl_array *actual, double tm, const char *kind)
{
  double xc[GKYL_MAX_DIM], xphys[GKYL_MAX_DIM], expected[3];
  const long num_nodes = nodes->size;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(range, iter.idx);
    const double *f = gkyl_array_cfetch(actual, lidx);
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    for (long node=0; node<nodes->size; ++node) {
      const double *eta = gkyl_array_cfetch(nodes, node);
      for (int d=0; d<grid->ndim; ++d) {
        xphys[d] = 0.5*grid->dx[d]*eta[d] + xc[d];
      }

      calc_expected_eval(grid->ndim, tm, xphys, expected);
      for (int c=0; c<3; ++c)
        assert_close(expected[c], f[node + num_nodes*c], 1e-14, kind, iter.idx, node, c);
    }
  }
}

static void
test_private_node_builders(void)
{
  int cases[][2] = { { 1, 1 }, { 2, 2 }, { 3, 2 } };

  for (int n=0; n<3; ++n) {
    int ndim = cases[n][0];
    int polyorder = cases[n][1];
    int nn = polyorder + 1;
    long num_nodes_vol = ndim == 1 ? nn : (ndim == 2 ? nn*nn : nn*nn*nn);
    long num_nodes_surf = ndim == 1 ? 1 : (ndim == 2 ? nn : nn*nn);

    struct gkyl_array *vol_nodes = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes_vol);
    struct gkyl_array *surf_x = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes_surf);
    struct gkyl_array *surf_y = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes_surf);
    struct gkyl_array *surf_z = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes_surf);

    vol_node_list(vol_nodes, ndim, polyorder);
    surf_node_list(surf_x, ndim, 0, polyorder);
    if (ndim > 1) surf_node_list(surf_y, ndim, 1, polyorder);
    if (ndim > 2) surf_node_list(surf_z, ndim, 2, polyorder);

    for (long i=0; i<num_nodes_vol; ++i) {
      double expected[GKYL_MAX_DIM] = { 0.0 };
      const double *actual = gkyl_array_cfetch(vol_nodes, i);

      expected_vol_node(ndim, polyorder, i, expected);
      for (int d=0; d<ndim; ++d)
        assert_close(expected[d], actual[d], 1e-14, "vol-node", (int[3]) { 0, 0, 0 }, i, d);
    }

    for (long i=0; i<num_nodes_surf; ++i) {
      double expected[GKYL_MAX_DIM] = { 0.0 };
      const double *actual_x = gkyl_array_cfetch(surf_x, i);

      expected_surf_node(ndim, 0, polyorder, i, expected);
      for (int d=0; d<ndim; ++d)
        assert_close(expected[d], actual_x[d], 1e-14, "surf-x-node", (int[3]) { 0, 0, 0 }, i, d);

      if (ndim > 1) {
        const double *actual_y = gkyl_array_cfetch(surf_y, i);
        expected_surf_node(ndim, 1, polyorder, i, expected);
        for (int d=0; d<ndim; ++d)
          assert_close(expected[d], actual_y[d], 1e-14, "surf-y-node", (int[3]) { 0, 0, 0 }, i, d);
      }

      if (ndim > 2) {
        const double *actual_z = gkyl_array_cfetch(surf_z, i);
        expected_surf_node(ndim, 2, polyorder, i, expected);
        for (int d=0; d<ndim; ++d)
          assert_close(expected[d], actual_z[d], 1e-14, "surf-z-node", (int[3]) { 0, 0, 0 }, i, d);
      }
    }

    gkyl_array_release(vol_nodes);
    gkyl_array_release(surf_x);
    gkyl_array_release(surf_y);
    gkyl_array_release(surf_z);
  }
}

static void
test_surf_and_vol_nodes_1x_identity(void)
{
  double lower[] = { -1.5 }, upper[] = { 0.5 };
  int cells[] = { 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 1, 2);

  struct eval_ctx ev = { .ndim = 1 };
  gkyl_dg_gr_maxwell_surf_and_vol_nodes *up =
    gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(&grid, &basis, 3, 2, analytic_eval, &ev);

  double basis_nodes[3];
  basis.node_list(basis_nodes);
  for (int i=0; i<3; ++i) {
    double *node = gkyl_dg_gr_maxwell_surf_and_vol_nodes_fetch_node(up, i);
    assert_close(basis_nodes[i], node[0], 1e-14, "basis-node", (int[3]) { 0, 0, 0 }, i, 0);
  }

  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_surf_and_vol_node_arrays *nodal =
    gkyl_surf_and_vol_node_arrays_new(up, range.volume, false);

  TEST_CHECK(nodal->use_gpu == false);
  TEST_CHECK(nodal->ndim == 1);
  TEST_CHECK(nodal->nodal_arr_vol == nodal->nodal_arr_vol_host);
  TEST_CHECK(nodal->nodal_arr_surf_x == nodal->nodal_arr_surf_x_host);
  TEST_CHECK(nodal->nodal_arr_vol->ncomp == 9);
  TEST_CHECK(nodal->nodal_arr_surf_x->ncomp == 3);

  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(up, 0.125, &range, nodal);

  struct gkyl_array *vol_nodes = gkyl_array_new(GKYL_DOUBLE, 1, 3);
  struct gkyl_array *surf_x = gkyl_array_new(GKYL_DOUBLE, 1, 1);
  vol_node_list(vol_nodes, 1, 2);
  surf_node_list(surf_x, 1, 0, 2);

  assert_node_array_matches(&grid, &range, vol_nodes, nodal->nodal_arr_vol_host, 0.125, "vol");
  assert_node_array_matches(&grid, &range, surf_x, nodal->nodal_arr_surf_x_host, 0.125, "surf-x");

  gkyl_array_release(vol_nodes);
  gkyl_array_release(surf_x);
  gkyl_surf_and_vol_node_arrays_release(nodal);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(up);
}

static void
test_surf_and_vol_nodes_2x_identity(void)
{
  double lower[] = { -1.0, 0.25 }, upper[] = { 1.0, 2.25 };
  int cells[] = { 2, 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 2);

  struct eval_ctx ev = { .ndim = 2 };
  gkyl_dg_gr_maxwell_surf_and_vol_nodes *up =
    gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(&grid, &basis, 3, 2, analytic_eval, &ev);

  double basis_nodes[18];
  basis.node_list(basis_nodes);
  for (int i=0; i<basis.num_basis; ++i) {
    double *node = gkyl_dg_gr_maxwell_surf_and_vol_nodes_fetch_node(up, i);
    for (int d=0; d<2; ++d)
      assert_close(basis_nodes[2*i + d], node[d], 1e-14, "basis-node", (int[3]) { 0, 0, 0 }, i, d);
  }

  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_surf_and_vol_node_arrays *nodal =
    gkyl_surf_and_vol_node_arrays_new(up, range.volume, false);

  TEST_CHECK(nodal->nodal_arr_vol->ncomp == 27);
  TEST_CHECK(nodal->nodal_arr_surf_x->ncomp == 9);
  TEST_CHECK(nodal->nodal_arr_surf_y->ncomp == 9);

  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(up, 0.25, &range, nodal);

  struct gkyl_array *vol_nodes = gkyl_array_new(GKYL_DOUBLE, 2, 9);
  struct gkyl_array *surf_x = gkyl_array_new(GKYL_DOUBLE, 2, 3);
  struct gkyl_array *surf_y = gkyl_array_new(GKYL_DOUBLE, 2, 3);
  vol_node_list(vol_nodes, 2, 2);
  surf_node_list(surf_x, 2, 0, 2);
  surf_node_list(surf_y, 2, 1, 2);

  assert_node_array_matches(&grid, &range, vol_nodes, nodal->nodal_arr_vol_host, 0.25, "vol");
  assert_node_array_matches(&grid, &range, surf_x, nodal->nodal_arr_surf_x_host, 0.25, "surf-x");
  assert_node_array_matches(&grid, &range, surf_y, nodal->nodal_arr_surf_y_host, 0.25, "surf-y");

  gkyl_array_release(vol_nodes);
  gkyl_array_release(surf_x);
  gkyl_array_release(surf_y);
  gkyl_surf_and_vol_node_arrays_release(nodal);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(up);
}

static void
test_surf_and_vol_nodes_3x_identity(void)
{
  double lower[] = { -2.0, 0.5, 1.0 }, upper[] = { 0.0, 2.5, 2.5 };
  int cells[] = { 2, 1, 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 3, 2);

  struct eval_ctx ev = { .ndim = 3 };
  gkyl_dg_gr_maxwell_surf_and_vol_nodes *up =
    gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(&grid, &basis, 3, 2, analytic_eval, &ev);

  double basis_nodes[81];
  basis.node_list(basis_nodes);
  for (int i=0; i<basis.num_basis; ++i) {
    double *node = gkyl_dg_gr_maxwell_surf_and_vol_nodes_fetch_node(up, i);
    for (int d=0; d<3; ++d)
      assert_close(basis_nodes[3*i + d], node[d], 1e-14, "basis-node", (int[3]) { 0, 0, 0 }, i, d);
  }

  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_surf_and_vol_node_arrays *nodal =
    gkyl_surf_and_vol_node_arrays_new(up, range.volume, false);

  TEST_CHECK(nodal->nodal_arr_vol->ncomp == 81);
  TEST_CHECK(nodal->nodal_arr_surf_x->ncomp == 27);
  TEST_CHECK(nodal->nodal_arr_surf_y->ncomp == 27);
  TEST_CHECK(nodal->nodal_arr_surf_z->ncomp == 27);

  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(up, 0.375, &range, nodal);

  struct gkyl_array *vol_nodes = gkyl_array_new(GKYL_DOUBLE, 3, 27);
  struct gkyl_array *surf_x = gkyl_array_new(GKYL_DOUBLE, 3, 9);
  struct gkyl_array *surf_y = gkyl_array_new(GKYL_DOUBLE, 3, 9);
  struct gkyl_array *surf_z = gkyl_array_new(GKYL_DOUBLE, 3, 9);
  vol_node_list(vol_nodes, 3, 2);
  surf_node_list(surf_x, 3, 0, 2);
  surf_node_list(surf_y, 3, 1, 2);
  surf_node_list(surf_z, 3, 2, 2);

  assert_node_array_matches(&grid, &range, vol_nodes, nodal->nodal_arr_vol_host, 0.375, "vol");
  assert_node_array_matches(&grid, &range, surf_x, nodal->nodal_arr_surf_x_host, 0.375, "surf-x");
  assert_node_array_matches(&grid, &range, surf_y, nodal->nodal_arr_surf_y_host, 0.375, "surf-y");
  assert_node_array_matches(&grid, &range, surf_z, nodal->nodal_arr_surf_z_host, 0.375, "surf-z");

  gkyl_array_release(vol_nodes);
  gkyl_array_release(surf_x);
  gkyl_array_release(surf_y);
  gkyl_array_release(surf_z);
  gkyl_surf_and_vol_node_arrays_release(nodal);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(up);
}

TEST_LIST = {
  { "private_node_builders", test_private_node_builders },
  { "surf_and_vol_nodes_1x_identity", test_surf_and_vol_nodes_1x_identity },
  { "surf_and_vol_nodes_2x_identity", test_surf_and_vol_nodes_2x_identity },
  { "surf_and_vol_nodes_3x_identity", test_surf_and_vol_nodes_3x_identity },
  { NULL, NULL },
};
