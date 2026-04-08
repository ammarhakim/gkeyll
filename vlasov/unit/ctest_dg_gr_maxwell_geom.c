#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static inline void
log_to_comp(int ndim, const double *eta, const double *dx, const double *xc, double *xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d] + xc[d];
}


void
test_ks_r_theta_2x_geom_p1()
{

  bool use_gpu = false;

  // GR paramters
  // populate the geometry context struct with mass and spin
  struct gkyl_dg_gr_maxwell_geom_ctx ctx = {
    .mass_bh = 1.0,
    .spin_bh = 0.0
  };
  double Mass = ctx.mass_bh;
  double Spin = ctx.spin_bh;

  // poylnomial order for this test is fixed to 1
  int poly_order = 1;

  double lower[] = {1.0, 0.0}, upper[] = {1.5, M_PI};
  int cells[] = {1, 4};
  int cdim = 2;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis confBasis;
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  // Create the fields (JD^i, JB^i) by setting the bare minimum required fields
  // in gkyl_vlasov_app, gkyl_vm input. 
  
  // Create the geometry alpha, \beta^i, h_ij, and J_c nodally 

  // Record which configured field boundaries use the theta-pole BC.
  enum gkyl_triad_preset_geom_type triad_preset_geom_type = GKYL_TRIAD_GR_KERR_SCHILD_RTHETA;
  struct gkyl_surf_and_vol_node_arrays *lapse, *shift, *h_ij, *det_h;

  // Evaluation of geometry at surface and volume nodal points.
  // Lapse - \alpha in the ADM split
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* lapse_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 1, poly_order, gkyl_dg_gr_maxwell_preset_lapse(triad_preset_geom_type), &ctx);
  lapse = gkyl_surf_and_vol_node_arrays_new(lapse_proj, confLocal_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(lapse_proj, 0.0, &confLocal_ext, lapse);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(lapse_proj);

  // shift - \beta^i components in the ADM split (contravariant)
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* shift_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 3, poly_order, gkyl_dg_gr_maxwell_preset_shift(triad_preset_geom_type), &ctx);
  shift = gkyl_surf_and_vol_node_arrays_new(shift_proj, confLocal_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(shift_proj, 0.0, &confLocal_ext, shift);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(shift_proj);

  // h_ij - Covariant components of the spatial metric (assumed to allways be a upper 
  // triangular matrix of 6 unique elements)
  // Allocate arrays for specified metric inverse
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* h_ij_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 6, poly_order, gkyl_dg_gr_maxwell_preset_h_ij(triad_preset_geom_type), &ctx);
  h_ij = gkyl_surf_and_vol_node_arrays_new(h_ij_proj, confLocal_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(h_ij_proj, 0.0, &confLocal_ext, h_ij);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(h_ij_proj);

  // Allocate arrays for the metric determinant (computed from J = sqrt(det(h_ij)))
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* det_h_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 1, poly_order, gkyl_dg_gr_maxwell_preset_det_h(triad_preset_geom_type), &ctx);
  det_h = gkyl_surf_and_vol_node_arrays_new(det_h_proj, confLocal_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(det_h_proj, 0.0, &confLocal_ext, det_h);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(det_h_proj);
  

  // Number of nodes (surf and volume)
  int num_surf_nodes = (poly_order + 1);
  int num_vol_nodes = (poly_order + 1) * (poly_order + 1);
  const double *gl_nodes = gkyl_gauss_ordinates[poly_order + 1];

  // Check the geometry has been correctly built against analytic expressions
  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    double xc[GKYL_MAX_DIM], xphys[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&confGrid, iter.idx, xc);

    // 1. Compare the surface-x alpha 
    const double *lapse_surf_x_d = gkyl_array_cfetch(lapse->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {
      double expected = 0;

      double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
      log_to_comp(cdim, eta, confGrid.dx, xc, xphys);

      // For each index, grab the expected value, for surf_x alpha, they happen to all be the same as 
      // alpha(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
      double alpha_lower_x_surf = 1.0 / sqrt( 1.0 + 2.0 * Mass / xphys[0]);
      // if (iter.idx[1] == 1) expected = alpha_lower_x_surf;
      // if (iter.idx[1] == 2) expected = alpha_lower_x_surf;
      // if (iter.idx[1] == 3) expected = alpha_lower_x_surf;
      // if (iter.idx[1] == 4) expected = alpha_lower_x_surf;
      expected = alpha_lower_x_surf;
      //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,lapse_surf_x_d[m], expected);
      TEST_CHECK( gkyl_compare_double(lapse_surf_x_d[m], expected, 1e-12) );
    }

    // 2. Compare the surface-y alpha
    const double *lapse_surf_y_d = gkyl_array_cfetch(lapse->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {
      double expected = 0;

      // Compute the radius we are evaluating at
      double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
      log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
      double r_val = xphys[0];

      // For each index, grab the expected value, for surf_y alpha, they happen to all be the same as 
      // alpha(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
      double alpha_lower_y_surf = 1.0 / sqrt( 1.0 + 2.0 * Mass / r_val);
      expected = alpha_lower_y_surf;
      //printf("(r=%1.1f), Value (Surf-y alpha) [%d]: %1.16e, Reference: %1.16e\n",r_val, m,lapse_surf_y_d[m], expected);
      TEST_CHECK( gkyl_compare_double(lapse_surf_y_d[m], expected, 1e-12) );
    }

    // 3. Compare the volume alpha
    const double *lapse_vol_d = gkyl_array_cfetch(lapse->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {
      double expected = 0;

      // Compute the radius we are evaluating at
      int theta_indx = m % num_surf_nodes;
      int r_indx = m / num_surf_nodes;
      double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
      log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
      double r_val = xphys[0];

      // For each index, grab the expected value, for vol alpha, they happen to all be the same as 
      // alpha(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
      double alpha_vol = 1.0 / sqrt( 1.0 + 2.0 * Mass / r_val);
      expected = alpha_vol;
      //printf("(r=%1.1f), Value (vol alpha) [%d]: %1.16e, Reference: %1.16e\n",r_val, m,lapse_vol_d[m], expected);
      TEST_CHECK( gkyl_compare_double(lapse_vol_d[m], expected, 1e-12) );
    }




    // 4. Compare the surface-x beta^i 
    const double *shift_surf_x_d = gkyl_array_cfetch(shift->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // Shift has 3 components
      double NC = 3; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];

        // For each index, grab the expected value, for surf_x beta, they happen to all be the same as 
        // beta(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
        double beta_lower_x_surf;
        if (n == 0) beta_lower_x_surf = ( 2.0 * Mass * r_val ) / ( r_val * r_val + 2.0 * Mass * r_val );
        if (n == 1) beta_lower_x_surf = 0.0;
        if (n == 2) beta_lower_x_surf = 0.0;
        expected = beta_lower_x_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,shift_surf_x_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(shift_surf_x_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }


    // 5. Compare the surface-y beta^i 
    const double *shift_surf_y_d = gkyl_array_cfetch(shift->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // Shift has 3 components
      double NC = 3; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];

        // For each index, grab the expected value, for surf_y beta, they happen to all be the same as 
        // beta(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
        double beta_lower_y_surf;
        if (n == 0) beta_lower_y_surf = ( 2.0 * Mass * r_val ) / ( r_val * r_val + 2.0 * Mass * r_val );
        if (n == 1) beta_lower_y_surf = 0.0;
        if (n == 2) beta_lower_y_surf = 0.0;
        expected = beta_lower_y_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,shift_surf_y_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(shift_surf_y_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 6. Compare the vol beta^i 
    const double *shift_vol_d = gkyl_array_cfetch(shift->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {

      // Shift has 3 components
      double NC = 3; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        int theta_indx = m % num_surf_nodes;
        int r_indx = m / num_surf_nodes;
        double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];

        // For each index, grab the expected value, for vol beta, they happen to all be the same as 
        // beta(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
        double beta_vol;
        if (n == 0) beta_vol = ( 2.0 * Mass * r_val ) / ( r_val * r_val + 2.0 * Mass * r_val );
        if (n == 1) beta_vol = 0.0;
        if (n == 2) beta_vol = 0.0;
        expected = beta_vol;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,shift_vol_d[m + n*num_vol_nodes], expected);
        TEST_CHECK( gkyl_compare_double(shift_vol_d[m + n*num_vol_nodes], expected, 1e-12) );
      }
    }



    // 7. Compare the surface-x h_ij
    const double *h_ij_surf_x_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // h_ij has 3 components
      double NC = 6; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_x h_ij
        double h_ij_lower_x_surf;
        if (n == 0) h_ij_lower_x_surf = (1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_lower_x_surf = 0.0;
        if (n == 2) h_ij_lower_x_surf = 0.0;
        if (n == 3) h_ij_lower_x_surf = r_val * r_val;
        if (n == 4) h_ij_lower_x_surf = 0.0;
        if (n == 5) h_ij_lower_x_surf = r_val * r_val * sin(theta_val) * sin(theta_val);
        expected = h_ij_lower_x_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_surf_x_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_surf_x_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }


    // 8. Compare the surface-y h_ij 
    const double *h_ij_surf_y_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // h_ij has 3 components
      double NC = 6; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_y h_ij
        double h_ij_lower_y_surf;
        if (n == 0) h_ij_lower_y_surf = (1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_lower_y_surf = 0.0;
        if (n == 2) h_ij_lower_y_surf = 0.0;
        if (n == 3) h_ij_lower_y_surf = r_val * r_val;
        if (n == 4) h_ij_lower_y_surf = 0.0;
        if (n == 5) h_ij_lower_y_surf = r_val * r_val * sin(theta_val) * sin(theta_val);
        expected = h_ij_lower_y_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_surf_y_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_surf_y_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 9. Compare the vol h_ij
    const double *h_ij_vol_d = gkyl_array_cfetch(h_ij->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {

      // h_ij has 3 components
      double NC = 6; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        int theta_indx = m % num_surf_nodes;
        int r_indx = m / num_surf_nodes;
        double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for vol h_ij
        double h_ij_vol;
        if (n == 0) h_ij_vol = (1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_vol = 0.0;
        if (n == 2) h_ij_vol = 0.0;
        if (n == 3) h_ij_vol = r_val * r_val;
        if (n == 4) h_ij_vol = 0.0;
        if (n == 5) h_ij_vol = r_val * r_val * sin(theta_val) * sin(theta_val);
        expected = h_ij_vol;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_vol_d[m + n*num_vol_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_vol_d[m + n*num_vol_nodes], expected, 1e-12) );
      }
    }


    // 10. Compare the surface-x det_h
    const double *det_h_surf_x_d = gkyl_array_cfetch(det_h->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // det_h has 3 components
      double NC = 1; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_x det_h
        double det_h_lower_x_surf = r_val * r_val * sin(theta_val) * sqrt( 1.0 + 2.0 * Mass / r_val);
        expected = det_h_lower_x_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,det_h_surf_x_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(det_h_surf_x_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }


    // 11. Compare the surface-y det_h 
    const double *det_h_surf_y_d = gkyl_array_cfetch(det_h->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // det_h has 3 components
      double NC = 1; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_y det_h
        double det_h_lower_y_surf = r_val * r_val * sin(theta_val) * sqrt( 1.0 + 2.0 * Mass / r_val);
        expected = det_h_lower_y_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,det_h_surf_y_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(det_h_surf_y_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 12. Compare the vol det_h
    const double *det_h_vol_d = gkyl_array_cfetch(det_h->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {

      // det_h has 3 components
      double NC = 1; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        int theta_indx = m % num_surf_nodes;
        int r_indx = m / num_surf_nodes;
        double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for vol det_h
        double det_h_vol = r_val * r_val * sin(theta_val) * sqrt( 1.0 + 2.0 * Mass / r_val);
        expected = det_h_vol;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,det_h_vol_d[m + n*num_vol_nodes], expected);
        TEST_CHECK( gkyl_compare_double(det_h_vol_d[m + n*num_vol_nodes], expected, 1e-12) );
      }
    }

  }

  // Release vmg data
  gkyl_surf_and_vol_node_arrays_release(lapse);
  gkyl_surf_and_vol_node_arrays_release(shift);
  gkyl_surf_and_vol_node_arrays_release(h_ij);
  gkyl_surf_and_vol_node_arrays_release(det_h);

}

void test_ks_r_theta_2x_p1() { test_ks_r_theta_2x_geom_p1(); }

TEST_LIST = {
  { "test_ks_r_theta_2x_p1", test_ks_r_theta_2x_p1},
  {NULL, NULL}
};
