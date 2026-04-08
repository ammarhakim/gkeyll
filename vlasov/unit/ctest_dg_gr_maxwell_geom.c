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
#include <gkyl_dg_gr_maxwell_divide_Jc.h>
#include <gkyl_dg_gr_maxwell_divide_Jc_priv.h>

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

void field_init_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct gkyl_dg_gr_maxwell_geom_ctx *app = ctx;
  double r = xn[0]; 
  double theta = xn[1];
  double massBH = app->mass_bh;
  double spinBH = app->spin_bh;
  double B_0 = 1.0; // Magnitude of the fields

  // Fields
  double Dr = 0.0; // Total electric field (r-direction).
  double Dtheta = 0.0; // Total electric field (theta-direction).
  double Dphi = - 2.0 * massBH * B_0 / ( r * r * sqrt( 1.0 + 2.0 * massBH / r ) );  // Total electric field (phi-direction).

  double Br = - B_0 * cos(theta) / ( sqrt( 1.0 + 2.0 * massBH / r ) );   // Total magnetic field (r-direction).
  double Btheta = B_0 * sin(theta) / ( r * sqrt( 1.0 + 2.0 * massBH / r ) ); // Total magnetic field (theta-direction).
  double Bphi = 0.0; // Total magnetic field (phi-direction).

  // Must return conserved variables
  double rho = sqrt(r * r + spinBH * spinBH * cos(theta) * cos(theta) );
  double metric_det = rho * sqrt(2*massBH*r + rho * rho) * sin(theta);
  
  // Compute Jc * D^i and Jc * B^i
  double JDr = metric_det * Dr;
  double JDtheta = metric_det * Dtheta;
  double JDphi = metric_det * Dphi;
  double JBr = metric_det * Br;
  double JBtheta = metric_det * Btheta;
  double JBphi = metric_det * Bphi;

  // Hand off the conserved varaibles (J * Q^\xi)
  fout[0] = JDr;
  fout[1] = JDtheta;
  fout[2] = JDphi;
  fout[3] = JBr;
  fout[4] = JBtheta;
  fout[5] = JBphi;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

static const double Dphi_ref[4][4] = {
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00},
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00},
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00},
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00}
};

static const double Br_ref[4][4] = {
  {-1.1141258010329544e+00, -3.9838967223318700e-02, +1.0646084038672876e-01, +3.8068321609655564e-03},
  {-4.6148601697763836e-01, -1.6501840534835825e-02, +2.5701920472327799e-01, +9.1905058326813213e-03},
  {+4.6148601697763814e-01, +1.6501840534835860e-02, +2.5701920472327799e-01, +9.1905058326813283e-03},
  {+1.1141258010329544e+00, +3.9838967223318700e-02, +1.0646084038672876e-01, +3.8068321609655564e-03}
};

static const double Btheta_ref[4][4] = {
  {+3.7263287721296423e-01, -2.9826465965693272e-02, +2.0753349447565991e-01, -1.6611499114132122e-02},
  {+8.9961534595364645e-01, -7.2007458652036282e-02, +8.5963188058500117e-02, -6.8807082244221884e-03},
  {+8.9961534595364645e-01, -7.2007458652036255e-02, -8.5963188058500034e-02, +6.8807082244221884e-03},
  {+3.7263287721296434e-01, -2.9826465965693265e-02, -2.0753349447565989e-01, +1.6611499114132122e-02}
};

static const double JDphi_ref[4][4] = {
  {-1.4915589107194009e+00, +0.0000000000000000e+00, -8.3070617727859664e-01, +0.0000000000000000e+00},
  {-3.6009417513372184e+00, -3.3306690738754696e-16, -3.4408976497590316e-01, +0.0000000000000000e+00},
  {-3.6009417513372184e+00, -2.2204460492503131e-16, +3.4408976497590282e-01, +0.0000000000000000e+00},
  {-1.4915589107194016e+00, +0.0000000000000000e+00, +8.3070617727859641e-01, -2.7755575615628914e-17}
};

static const double JBr_ref[4][4] = {
  {-1.0064417291126031e+00, -2.2936950127375039e-01, -4.9045618369336014e-01, -1.1177566171621522e-01},
  {-1.0064417291126033e+00, -2.2936950127375047e-01, +4.9045618369336008e-01, +1.1177566171621517e-01},
  {+1.0064417291126031e+00, +2.2936950127375044e-01, +4.9045618369336008e-01, +1.1177566171621517e-01},
  {+1.0064417291126035e+00, +2.2936950127375053e-01, -4.9045618369335997e-01, -1.1177566171621522e-01}
};

static const double JBtheta_ref[4][4] = {
  {+4.5544074017426023e-01, +5.2589766787906389e-02, +3.8720225028423161e-01, +4.4710264686486062e-02},
  {+2.0445592598257387e+00, +2.3608536780690687e-01, +3.8720225028423155e-01, +4.4710264686485979e-02},
  {+2.0445592598257387e+00, +2.3608536780690687e-01, -3.8720225028423133e-01, -4.4710264686486034e-02},
  {+4.5544074017426045e-01, +5.2589766787906451e-02, -3.8720225028423166e-01, -4.4710264686486104e-02}
};

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


  // (TEST I) Geometry objects projected nodally 
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

      // det_h has 1 component
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

      // det_h has 1 component
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

      // det_h has 1 component
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

  // (TEST II) division by J_c
  // Allocate the field arrays in memory
  struct gkyl_array *field_no_J_con, *field_with_J_con;
  field_no_J_con = mkarr(8*confBasis.num_basis, confLocal_ext.volume);
  field_with_J_con = mkarr(8*confBasis.num_basis, confLocal_ext.volume);

  // Project the inital conditions (Kerr-Schild Coordinates, Schwarzschild Wald solution)
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1,
     8, field_init_func, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &confLocal_ext, field_with_J_con);
  gkyl_proj_on_basis_release(proj);

  // Divide out the spatial Jacobian
  gkyl_dg_gr_maxwell_divide_Jc(&confBasis, &confLocal, det_h->nodal_arr_vol, 
    field_with_J_con, field_no_J_con, false);

  // Test the fields before and after division by J_c (mode by mode).
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&confLocal, iter.idx);

    const double *em      = gkyl_array_cfetch(field_with_J_con, lidx);
    const double *em_no_J = gkyl_array_cfetch(field_no_J_con, lidx);

    // with J
    const double *JDr     = &em[0*confBasis.num_basis];
    const double *JDtheta = &em[1*confBasis.num_basis];
    const double *JDphi   = &em[2*confBasis.num_basis];
    const double *JBr     = &em[3*confBasis.num_basis];
    const double *JBtheta = &em[4*confBasis.num_basis];
    const double *JBphi   = &em[5*confBasis.num_basis];

    // without J
    const double *Dr      = &em_no_J[0*confBasis.num_basis];
    const double *Dtheta  = &em_no_J[1*confBasis.num_basis];
    const double *Dphi    = &em_no_J[2*confBasis.num_basis];
    const double *Br      = &em_no_J[3*confBasis.num_basis];
    const double *Btheta  = &em_no_J[4*confBasis.num_basis];
    const double *Bphi    = &em_no_J[5*confBasis.num_basis];

    int jloc = iter.idx[1] - confLocal.lower[1]; // 0..3 for Nz1 = 4

    const double *Dphi_mode_ref    = Dphi_ref[jloc];
    const double *Br_mode_ref      = Br_ref[jloc];
    const double *Btheta_mode_ref  = Btheta_ref[jloc];

    const double *JDphi_mode_ref   = JDphi_ref[jloc];
    const double *JBr_mode_ref     = JBr_ref[jloc];
    const double *JBtheta_mode_ref = JBtheta_ref[jloc];

    // Loop over all bases
    for (int k = 0; k < confBasis.num_basis; ++k) {

      // printf("\n---------\n");    
      // printf("(Dphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx, k,Dphi[k],   Dphi_mode_ref[k]);
      // printf("(Br)   Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx, k,Br[k],     Br_mode_ref[k]);
      // printf("(Bphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx, k,Btheta[k], Btheta_mode_ref[k]);  
      // printf("(JDphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx,k,JDphi[k],   JDphi_mode_ref[k]);
      // printf("(JBr)   Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx,k,JBr[k],     JBr_mode_ref[k]);
      // printf("(JBphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx,k,JBtheta[k], JBtheta_mode_ref[k]);

      // no J
      TEST_CHECK(gkyl_compare_double(Dr[k],     0.0,                1e-12));
      TEST_CHECK(gkyl_compare_double(Dtheta[k], 0.0,                1e-12));
      TEST_CHECK(gkyl_compare_double(Dphi[k],   Dphi_mode_ref[k],   1e-12));
      TEST_CHECK(gkyl_compare_double(Br[k],     Br_mode_ref[k],     1e-12));
      TEST_CHECK(gkyl_compare_double(Btheta[k], Btheta_mode_ref[k], 1e-12));
      TEST_CHECK(gkyl_compare_double(Bphi[k],   0.0,                1e-12));

      // with J
      TEST_CHECK(gkyl_compare_double(JDr[k],     0.0,                 1e-12));
      TEST_CHECK(gkyl_compare_double(JDtheta[k], 0.0,                 1e-12));
      TEST_CHECK(gkyl_compare_double(JDphi[k],   JDphi_mode_ref[k],   1e-12));
      TEST_CHECK(gkyl_compare_double(JBr[k],     JBr_mode_ref[k],     1e-12));
      TEST_CHECK(gkyl_compare_double(JBtheta[k], JBtheta_mode_ref[k], 1e-12));
      TEST_CHECK(gkyl_compare_double(JBphi[k],   0.0,                 1e-12));
    }
  }
    


  // (TEST III) Configuration space flux and maximum eigenvalue calculation
  // Compute the x-surface
  int theta_pole_lo[3] = {0,1,0};
  struct gkyl_array *conf_flux_surf = mkarr(cdim*8*num_surf_nodes, confLocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, confLocal_ext.volume);
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM]; 
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx);
    long cidx = gkyl_range_idx(&confLocal, idx);

    // Grab the cell center location for NC bracket calculation 
    double xcC[GKYL_MAX_DIM], xcR[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&confGrid, idx, xcC);

    const double *field_no_J_con_c = gkyl_array_cfetch(field_no_J_con, cidx); 
    double *cflrate_d = gkyl_array_fetch(cflrate, cidx);
    double *flux = gkyl_array_fetch(conf_flux_surf, cidx); 

    // Each cell owns *lower* fluxes in each configuration-space direction. 
    // So we need the distribution function in our current cell, and the cell
    // one index lower in each direction. If we are at the lower configuration-space
    // edge, we call ghost cells
    for (int dir = 0; dir<cdim; ++dir) {

      // Select the geometry elements based on direction
      const double *lapse_d = 0;
      const double *shift_d = 0;
      const double *h_ij_d = 0;
      const double *det_h_d = 0;

      if (dir == 0) {
        lapse_d = gkyl_array_cfetch(lapse->nodal_arr_surf_x, cidx);
        shift_d = gkyl_array_cfetch(shift->nodal_arr_surf_x, cidx);
        h_ij_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_x, cidx);
        det_h_d = gkyl_array_cfetch(det_h->nodal_arr_surf_x, cidx);
      } else if (dir == 1) {
        lapse_d = gkyl_array_cfetch(lapse->nodal_arr_surf_y, cidx);
        shift_d = gkyl_array_cfetch(shift->nodal_arr_surf_y, cidx);
        h_ij_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_y, cidx);
        det_h_d = gkyl_array_cfetch(det_h->nodal_arr_surf_y, cidx);
      } else {
        lapse_d = gkyl_array_cfetch(lapse->nodal_arr_surf_z, cidx);
        shift_d = gkyl_array_cfetch(shift->nodal_arr_surf_z, cidx);
        h_ij_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_z, cidx);
        det_h_d = gkyl_array_cfetch(det_h->nodal_arr_surf_z, cidx);
      }

      // Create an index for the left cell (which may be a ghost cell) 
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir]-1;
      long cidx_l = gkyl_range_idx(&confLocal, idx_l); 
      const double *field_no_J_con_l = gkyl_array_cfetch(field_no_J_con, cidx_l);

      // For Points not along the domain-edge in theta, compute the left hand surface 
      // conf-flux.
      int theta_pole = 0;

      // If at the left edge in the theta dir
      if(idx[dir] == confLocal.lower[dir] && theta_pole_lo[dir]) {
        theta_pole = 1;
      }
      

      // 2 quadrature nodes on a 2x surface, 6 field components.
      double alpha_quad_x[2] = {0.0}; 
      double flux_l_x[12] = {0.0};
      double flux_r_x[12] = {0.0};
      double alpha_quad_y[2] = {0.0}; 
      double flux_l_y[12] = {0.0};
      double flux_r_y[12] = {0.0};

      // Compute the fluxes in the first two directions
      if (dir == 0) { 
        dg_gr_maxwell_alpha_quad_x_2x_ser_p1(xcC, confGrid.dx, theta_pole,
          lapse_d, shift_d, h_ij_d, det_h_d, field_no_J_con_l, field_no_J_con_c, flux_l_x, flux_r_x, alpha_quad_x);
          
        double cflrate = lax_flux_nodal_to_modal_x_2x_ser_p1(confGrid.dx, det_h_d, flux_l_x, flux_r_x, alpha_quad_x,
          field_no_J_con_l, field_no_J_con_c, flux); 
      }
      else if (dir == 1) { 
        dg_gr_maxwell_alpha_quad_y_2x_ser_p1(xcC, confGrid.dx, theta_pole,
          lapse_d, shift_d, h_ij_d, det_h_d, field_no_J_con_l, field_no_J_con_c, flux_l_y, flux_r_y, alpha_quad_y);
          
        double cflrate = lax_flux_nodal_to_modal_y_2x_ser_p1(confGrid.dx, det_h_d, flux_l_y, flux_r_y, alpha_quad_y,
          field_no_J_con_l, field_no_J_con_c, flux); 
      }
    }
  }

  // Release vmg data
  gkyl_surf_and_vol_node_arrays_release(lapse);
  gkyl_surf_and_vol_node_arrays_release(shift);
  gkyl_surf_and_vol_node_arrays_release(h_ij);
  gkyl_surf_and_vol_node_arrays_release(det_h);
  
  
  // Release field arrays
  gkyl_array_release(field_with_J_con);
  gkyl_array_release(field_no_J_con);
  gkyl_array_release(conf_flux_surf);
  gkyl_array_release(cflrate);

}

void test_ks_r_theta_2x_p1() { test_ks_r_theta_2x_geom_p1(); }

TEST_LIST = {
  { "test_ks_r_theta_2x_p1", test_ks_r_theta_2x_p1},
  {NULL, NULL}
};
