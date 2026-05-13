// Unit tests for the modular tetrad-basis GR Euler equation
// (wv_gr_euler_tetrad_mod). Strategy mirrors ctest_wv_gr_euler_mod.c:
// every test sets up the packed wv_gr_euler_tetrad equation (71-comp state
// with spacetime baked in) and the mod-tetrad equation (5-comp hydro state
// plus an external spacetime-products array) from the same spacetime
// callbacks, then verifies they produce equivalent results.
//
// The tetrad-specific arithmetic lives in the flat/correction flux split,
// so the per-cell flux equivalence checks compare BOTH the flat SR flux
// AND the corrected GR flux directly.
//
// Any divergence here means the flat-flux / correction port drifted away
// from the packed tetrad path and the regression-test acceptance gate will
// not hold.

#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_euler_tetrad_mod.h>
#include <gkyl_wv_gr_euler_tetrad_mod_priv.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>

// Helper: pack the packed-layout spacetime block (q[5..66]) and the mod
// products row (prods[0..NCOMP_BASE)) from the spacetime callbacks at the
// physical point (x,y,z). Identical to ctest_wv_gr_euler_mod.c — the tetrad
// equation reads the same spacetime layout.
static void
fill_spacetime(struct gkyl_gr_spacetime *spacetime, double x, double y, double z,
  double q[71], double *prods)
{
  double spatial_det, lapse;
  bool in_excision_region;
  double *shift = gkyl_malloc(sizeof(double[3]));
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }
  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));
  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++)
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }
  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z,
    pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z,
    pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z,
    pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z,
    pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

  q[5] = lapse;
  q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];
  q[9]  = spatial_metric[0][0]; q[10] = spatial_metric[0][1]; q[11] = spatial_metric[0][2];
  q[12] = spatial_metric[1][0]; q[13] = spatial_metric[1][1]; q[14] = spatial_metric[1][2];
  q[15] = spatial_metric[2][0]; q[16] = spatial_metric[2][1]; q[17] = spatial_metric[2][2];
  q[18] = extrinsic_curvature[0][0]; q[19] = extrinsic_curvature[0][1]; q[20] = extrinsic_curvature[0][2];
  q[21] = extrinsic_curvature[1][0]; q[22] = extrinsic_curvature[1][1]; q[23] = extrinsic_curvature[1][2];
  q[24] = extrinsic_curvature[2][0]; q[25] = extrinsic_curvature[2][1]; q[26] = extrinsic_curvature[2][2];
  q[27] = in_excision_region ? -1.0 : 1.0;
  q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
  q[31] = shift_der[0][0]; q[32] = shift_der[0][1]; q[33] = shift_der[0][2];
  q[34] = shift_der[1][0]; q[35] = shift_der[1][1]; q[36] = shift_der[1][2];
  q[37] = shift_der[2][0]; q[38] = shift_der[2][1]; q[39] = shift_der[2][2];
  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        q[40 + 9*k + 3*i + j] = spatial_metric_der[k][i][j];
      }
    }
  }
  q[67] = 0.0;
  q[68] = x; q[69] = y; q[70] = z;

  for (int k = 0; k < 62; k++) prods[k] = q[5 + k];

  gkyl_gr_euler_tetrad_inv_spatial_metric(q, &inv_spatial_metric);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_INV_GIJ + 3*i + j] = inv_spatial_metric[i][j];
  prods[GKYL_GR_SP_SPATIAL_DET] = spatial_det;
  prods[GKYL_GR_SP_CELLCTR + 0] = x;
  prods[GKYL_GR_SP_CELLCTR + 1] = y;
  prods[GKYL_GR_SP_CELLCTR + 2] = z;

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]);
    gkyl_free(shift_der[i]);
    gkyl_free(inv_spatial_metric[i]);
    for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]);
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(extrinsic_curvature);
  gkyl_free(shift_der);
  gkyl_free(spatial_metric_der);
  gkyl_free(inv_spatial_metric);
  gkyl_free(shift);
  gkyl_free(lapse_der);
}

static void
check_rotated_spacetime(const double q_local_packed[71],
  const double prod_local[GKYL_GR_SP_NCOMP_BASE])
{
  for (int k = 0; k < 62; k++)
    TEST_CHECK( gkyl_compare(prod_local[k], q_local_packed[5 + k], 1e-12) );
  TEST_CHECK( gkyl_compare(prod_local[GKYL_GR_SP_CELLCTR + 0], q_local_packed[68], 1e-12) );
  TEST_CHECK( gkyl_compare(prod_local[GKYL_GR_SP_CELLCTR + 1], q_local_packed[69], 1e-12) );
  TEST_CHECK( gkyl_compare(prod_local[GKYL_GR_SP_CELLCTR + 2], q_local_packed[70], 1e-12) );
}

// Rotation / per-cell math equivalence test runner. Compares mod-tetrad vs
// packed-tetrad for prim_vars, max_abs_speed, check_inv, rotation round-trip,
// the flat-space SR flux, the GR flux correction, flux_jump, and locally-
// rotated products.
static void
run_rotation_equivalence(struct gkyl_gr_spacetime *spacetime, int half_extent,
  double step)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_wv_eqn *gr_packed = gkyl_wv_gr_euler_tetrad_new(
    gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  int lower[1] = { 0 };
  int upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  struct gkyl_wv_eqn *gr_mod = gkyl_wv_gr_euler_tetrad_mod_new(
    gas_gamma, &conf_range, false);
  gkyl_gr_euler_tetrad_mod_set_auxfields(gr_mod,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  TEST_CHECK( gr_mod->num_equations == 5 );
  TEST_CHECK( gr_packed->num_equations == 71 );
  TEST_CHECK( gr_mod->set_interface_idx_func != NULL );
  TEST_CHECK( gr_mod->set_cell_idx_func != NULL );

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };
  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };
  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  struct wv_gr_euler_tetrad_mod *grm =
    container_of(gr_mod, struct wv_gr_euler_tetrad_mod, eqn);

  for (int x_ind = -half_extent; x_ind < half_extent + 1; x_ind++) {
    for (int y_ind = -half_extent; y_ind < half_extent + 1; y_ind++) {
      double x = step * x_ind;
      double y = step * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double q[71];
      double *prods_row = gkyl_array_fetch(prods, 0);
      fill_spacetime(spacetime, x, y, 0.0, q, prods_row);

      if (q[27] < 0.0) continue;

      double spatial_det = prods_row[GKYL_GR_SP_SPATIAL_DET];

      double vel[3] = { u, v, w };
      double v_sq = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          v_sq += prods_row[GKYL_GR_SP_GIJ + 3*i + j] * vel[i] * vel[j];

      double W = 1.0 / sqrt(1.0 - v_sq);
      double h = 1.0 + (p / rho) * (gas_gamma / (gas_gamma - 1.0));

      q[0] = sqrt(spatial_det) * rho * W;
      q[1] = sqrt(spatial_det) * rho * h * (W*W) * u;
      q[2] = sqrt(spatial_det) * rho * h * (W*W) * v;
      q[3] = sqrt(spatial_det) * rho * h * (W*W) * w;
      q[4] = sqrt(spatial_det) * ((rho * h * (W*W)) - p - (rho * W));

      double q_mod[5] = { q[0], q[1], q[2], q[3], q[4] };

      int idx_l[1] = { 0 };
      int idx_r[1] = { 0 };

      // prim_vars equivalence (packed tetrad vs mod tetrad).
      double v_packed[71], v_mod[5];
      gkyl_gr_euler_tetrad_prim_vars(gas_gamma, q, v_packed);
      gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, q_mod, prods_row, v_mod);
      for (int i = 0; i < 5; i++)
        TEST_CHECK( gkyl_compare(v_mod[i], v_packed[i], 1e-12) );

      gr_mod->set_cell_idx_func(gr_mod, idx_l);
      double maxs_packed = gr_packed->max_speed_func(gr_packed, q);
      double maxs_mod    = gr_mod->max_speed_func(gr_mod, q_mod);
      TEST_CHECK( gkyl_compare(maxs_mod, maxs_packed, 1e-12) );

      gr_mod->set_cell_idx_func(gr_mod, idx_l);
      bool inv_packed = gr_packed->check_inv_func(gr_packed, q);
      bool inv_mod    = gr_mod->check_inv_func(gr_mod, q_mod);
      TEST_CHECK( inv_packed == inv_mod );

      gr_mod->set_interface_idx_func(gr_mod, idx_l, idx_r);

      for (int d = 0; d < 3; d++) {
        double q_local_packed[71], q_global_back[71];
        gr_packed->rotate_to_local_func(gr_packed,
          tau1[d], tau2[d], norm[d], q, q_local_packed);
        gr_packed->rotate_to_global_func(gr_packed,
          tau1[d], tau2[d], norm[d], q_local_packed, q_global_back);
        for (int i = 0; i < 71; i++)
          TEST_CHECK( gkyl_compare(q[i], q_global_back[i], 1e-12) );

        gr_mod->set_interface_idx_func(gr_mod, idx_l, idx_r);

        double q_local_mod_l[5], q_local_mod_r[5], q_global_mod[5];
        gr_mod->rotate_to_local_func(gr_mod,
          tau1[d], tau2[d], norm[d], q_mod, q_local_mod_l);
        gr_mod->rotate_to_local_func(gr_mod,
          tau1[d], tau2[d], norm[d], q_mod, q_local_mod_r);
        gr_mod->rotate_to_global_func(gr_mod,
          tau1[d], tau2[d], norm[d], q_local_mod_l, q_global_mod);

        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(q_mod[i], q_global_mod[i], 1e-12) );

        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(q_local_mod_l[i], q_local_packed[i], 1e-12) );

        // Flat-space SR flux equivalence — the meat of the tetrad split.
        double flux_sr_packed[71], flux_sr_mod[5];
        gkyl_gr_euler_tetrad_flux(gas_gamma, q_local_packed, flux_sr_packed);
        gkyl_gr_euler_tetrad_mod_flux(gas_gamma, q_local_mod_l, grm->prodl_local, flux_sr_mod);
        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(flux_sr_mod[i], flux_sr_packed[i], 1e-12) );

        // GR correction step. Both consume the flat flux + their own
        // locally-rotated spacetime block and produce the curved-space flux.
        double flux_gr_packed[71], flux_gr_mod[5];
        gkyl_gr_euler_tetrad_flux_correction(gas_gamma,
          q_local_packed, flux_sr_packed, flux_gr_packed);
        gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma,
          q_local_mod_l, grm->prodl_local, flux_sr_mod, flux_gr_mod);
        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(flux_gr_mod[i], flux_gr_packed[i], 1e-12) );

        // flux_jump_func via the function-pointer interface. Single-cell
        // test: jump is identically 0 on both sides.
        double fj_packed[71], fj_mod[5];
        double amax_packed = gr_packed->flux_jump(gr_packed,
          q_local_packed, q_local_packed, fj_packed);
        double amax_mod    = gr_mod->flux_jump(gr_mod,
          q_local_mod_l, q_local_mod_l, fj_mod);
        TEST_CHECK( gkyl_compare(amax_mod, amax_packed, 1e-12) );
        for (int i = 0; i < 5; i++) {
          TEST_CHECK( gkyl_compare(fj_mod[i],     0.0, 1e-12) );
          TEST_CHECK( gkyl_compare(fj_packed[i],  0.0, 1e-12) );
        }

        check_rotated_spacetime(q_local_packed, grm->prodl_local);
        check_rotated_spacetime(q_local_packed, grm->prodr_local);

        TEST_CHECK( gkyl_compare(grm->prodl_local[GKYL_GR_SP_SPATIAL_DET],
                                 prods_row[GKYL_GR_SP_SPATIAL_DET], 1e-12) );
        double check[3][3] = {{0}};
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
              check[i][j] +=
                grm->prodl_local[GKYL_GR_SP_GIJ + 3*i + k] *
                grm->prodl_local[GKYL_GR_SP_INV_GIJ + 3*k + j];
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            TEST_CHECK( gkyl_compare(check[i][j], (i == j) ? 1.0 : 0.0, 1e-10) );

        double w_packed[71], q_packed_cycle[71];
        gr_packed->cons_to_riem(gr_packed, q_local_packed, q_local_packed, w_packed);
        gr_packed->riem_to_cons(gr_packed, q_local_packed, w_packed, q_packed_cycle);
        for (int i = 0; i < 71; i++)
          TEST_CHECK( gkyl_compare(q_local_packed[i], q_packed_cycle[i], 1e-12) );

        double w_mod[5], q_mod_cycle[5];
        gr_mod->cons_to_riem(gr_mod, q_local_mod_l, q_local_mod_l, w_mod);
        gr_mod->riem_to_cons(gr_mod, q_local_mod_l, w_mod, q_mod_cycle);
        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(q_local_mod_l[i], q_mod_cycle[i], 1e-12) );
      }
    }
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(gr_mod);
  gkyl_wv_eqn_release(gr_packed);
}

static void
test_gr_euler_tetrad_mod_construction()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  int lower[1] = { 0 };
  int upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_euler_tetrad_mod_new(5.0/3.0, &conf_range, false);
  TEST_CHECK( eqn->type == GKYL_EQN_GR_EULER_TETRAD_MOD );
  TEST_CHECK( eqn->num_equations == 5 );
  TEST_CHECK( eqn->num_diag == 5 );
  TEST_CHECK( eqn->num_waves == 2 );      // HLL default
  TEST_CHECK( eqn->set_interface_idx_func != NULL );
  TEST_CHECK( eqn->set_cell_idx_func != NULL );
  TEST_CHECK( eqn->rotate_to_local_func != NULL );
  TEST_CHECK( eqn->rotate_to_global_func != NULL );
  TEST_CHECK( eqn->cons_to_diag != NULL );

  TEST_CHECK( gkyl_compare(gkyl_wv_gr_euler_tetrad_mod_gas_gamma(eqn), 5.0/3.0, 1e-15) );

  struct gkyl_wv_eqn *eqn_lax = gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = 1.4, .conf_range = conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_LAX, .use_gpu = false,
    });
  TEST_CHECK( eqn_lax->num_waves == 2 );

  struct gkyl_wv_eqn *eqn_roe = gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = 1.4, .conf_range = conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_ROE, .use_gpu = false,
    });
  TEST_CHECK( eqn_roe->num_waves == 3 );

  gkyl_wv_eqn_release(eqn);
  gkyl_wv_eqn_release(eqn_lax);
  gkyl_wv_eqn_release(eqn_roe);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_rotation_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_rotation_equivalence(spacetime, 10, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_rotation_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_rotation_equivalence(spacetime, 10, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_rotation_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_rotation_equivalence(spacetime, 10, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

// Riemann-solver equivalence: identical structure to the regular mod
// test, but constructs packed-tetrad and mod-tetrad equations and compares
// hydro waves/qfluct between them.
static void
run_riemann_equivalence(struct gkyl_gr_spacetime *spacetime,
  enum gkyl_wv_gr_euler_tetrad_rp rp_type, int num_waves)
{
  double gas_gamma = 5.0 / 3.0;

  struct gkyl_wv_eqn *gr_packed = gkyl_wv_gr_euler_tetrad_inew(
    &(struct gkyl_wv_gr_euler_tetrad_inp){
      .gas_gamma = gas_gamma,
      .spacetime_gauge = GKYL_STATIC_GAUGE,
      .reinit_freq = 0,
      .spacetime = spacetime,
      .rp_type = rp_type,
      .use_gpu = false,
    });

  int lower[1] = { 0 }, upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);

  struct gkyl_wv_eqn *gr_mod = gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = gas_gamma,
      .conf_range = conf_range,
      .rp_type = rp_type,
      .use_gpu = false,
    });
  gkyl_gr_euler_tetrad_mod_set_auxfields(gr_mod,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  TEST_CHECK( gr_mod->num_waves == num_waves );
  TEST_CHECK( gr_packed->num_waves == num_waves );

  struct wv_gr_euler_tetrad_mod *grm = container_of(gr_mod,
    struct wv_gr_euler_tetrad_mod, eqn);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  for (int x_ind = -5; x_ind < 6; x_ind++) {
    for (int y_ind = -5; y_ind < 6; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double ql_packed[71], qr_packed[71];
      double *prods_row = gkyl_array_fetch(prods, 0);
      fill_spacetime(spacetime, x, y, 0.0, ql_packed, prods_row);
      for (int i = 5; i < 71; i++) qr_packed[i] = ql_packed[i];

      if (ql_packed[27] < 0.0) continue;

      double rho_l = 1.0, ul = 0.10, vl = 0.20, wl = 0.30, pl = 1.5;
      double rho_r = 0.5, ur = 0.05, vr = 0.10, wr = 0.15, pr = 0.7;

      double spatial_det = prods_row[GKYL_GR_SP_SPATIAL_DET];
      double vell[3] = { ul, vl, wl }, velr[3] = { ur, vr, wr };
      double vsq_l = 0.0, vsq_r = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          vsq_l += prods_row[GKYL_GR_SP_GIJ + 3*i + j] * vell[i] * vell[j];
          vsq_r += prods_row[GKYL_GR_SP_GIJ + 3*i + j] * velr[i] * velr[j];
        }
      double Wl = 1.0 / sqrt(1.0 - vsq_l);
      double Wr = 1.0 / sqrt(1.0 - vsq_r);
      double hl = 1.0 + (pl / rho_l) * (gas_gamma / (gas_gamma - 1.0));
      double hr = 1.0 + (pr / rho_r) * (gas_gamma / (gas_gamma - 1.0));

      ql_packed[0] = sqrt(spatial_det) * rho_l * Wl;
      ql_packed[1] = sqrt(spatial_det) * rho_l * hl * (Wl*Wl) * ul;
      ql_packed[2] = sqrt(spatial_det) * rho_l * hl * (Wl*Wl) * vl;
      ql_packed[3] = sqrt(spatial_det) * rho_l * hl * (Wl*Wl) * wl;
      ql_packed[4] = sqrt(spatial_det) * ((rho_l * hl * (Wl*Wl)) - pl - (rho_l * Wl));

      qr_packed[0] = sqrt(spatial_det) * rho_r * Wr;
      qr_packed[1] = sqrt(spatial_det) * rho_r * hr * (Wr*Wr) * ur;
      qr_packed[2] = sqrt(spatial_det) * rho_r * hr * (Wr*Wr) * vr;
      qr_packed[3] = sqrt(spatial_det) * rho_r * hr * (Wr*Wr) * wr;
      qr_packed[4] = sqrt(spatial_det) * ((rho_r * hr * (Wr*Wr)) - pr - (rho_r * Wr));

      double ql_mod[5] = { ql_packed[0], ql_packed[1], ql_packed[2], ql_packed[3], ql_packed[4] };
      double qr_mod[5] = { qr_packed[0], qr_packed[1], qr_packed[2], qr_packed[3], qr_packed[4] };

      double ql_packed_local[71], qr_packed_local[71];
      gr_packed->rotate_to_local_func(gr_packed, tau1, tau2, norm, ql_packed, ql_packed_local);
      gr_packed->rotate_to_local_func(gr_packed, tau1, tau2, norm, qr_packed, qr_packed_local);

      int idx_l[1] = { 0 }, idx_r[1] = { 0 };
      gr_mod->set_interface_idx_func(gr_mod, idx_l, idx_r);
      double ql_mod_local[5], qr_mod_local[5];
      gr_mod->rotate_to_local_func(gr_mod, tau1, tau2, norm, ql_mod, ql_mod_local);
      gr_mod->rotate_to_local_func(gr_mod, tau1, tau2, norm, qr_mod, qr_mod_local);

      double delta_packed[71], delta_mod[5];
      for (int i = 0; i < 71; i++) delta_packed[i] = qr_packed_local[i] - ql_packed_local[i];
      for (int i = 0; i <  5; i++) delta_mod[i]    = qr_mod_local[i]    - ql_mod_local[i];

      double waves_packed[3 * 71];
      double speeds_packed[3];
      double waves_mod[3 * 5];
      double speeds_mod[3];

      double maxs_packed = gr_packed->waves_func(gr_packed,
        GKYL_WV_HIGH_ORDER_FLUX, delta_packed,
        ql_packed_local, qr_packed_local, 1.0, 1.0,
        waves_packed, speeds_packed);
      double maxs_mod    = gr_mod->waves_func(gr_mod,
        GKYL_WV_HIGH_ORDER_FLUX, delta_mod,
        ql_mod_local, qr_mod_local, 1.0, 1.0,
        waves_mod, speeds_mod);

      TEST_CHECK( gkyl_compare(maxs_mod, maxs_packed, 1e-12) );

      for (int m = 0; m < num_waves; m++) {
        TEST_CHECK( gkyl_compare(speeds_mod[m], speeds_packed[m], 1e-12) );
        for (int i = 0; i < 5; i++) {
          double wp = waves_packed[m * 71 + i];
          double wm = waves_mod[m * 5 + i];
          TEST_CHECK( gkyl_compare(wm, wp, 1e-12) );
        }
      }

      double amdq_packed[71], apdq_packed[71];
      double amdq_mod[5],     apdq_mod[5];
      gr_packed->qfluct_func(gr_packed, GKYL_WV_HIGH_ORDER_FLUX,
        ql_packed_local, qr_packed_local, 1.0, 1.0,
        waves_packed, speeds_packed, amdq_packed, apdq_packed);
      gr_mod->qfluct_func(gr_mod, GKYL_WV_HIGH_ORDER_FLUX,
        ql_mod_local, qr_mod_local, 1.0, 1.0,
        waves_mod, speeds_mod, amdq_mod, apdq_mod);

      for (int i = 0; i < 5; i++) {
        TEST_CHECK( gkyl_compare(amdq_mod[i], amdq_packed[i], 1e-12) );
        TEST_CHECK( gkyl_compare(apdq_mod[i], apdq_packed[i], 1e-12) );
      }

      if (rp_type != WV_GR_EULER_TETRAD_RP_ROE) {
        for (int i = 0; i < 5; i++) {
          double sum = 0.0;
          for (int m = 0; m < num_waves; m++) sum += waves_mod[m * 5 + i];
          TEST_CHECK( gkyl_compare(sum, qr_mod_local[i] - ql_mod_local[i], 1e-11) );
        }
      }
      (void)grm;
    }
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(gr_mod);
  gkyl_wv_eqn_release(gr_packed);
}

// Roe-specific consistency checks for the tetrad mod equation. We test the
// four classical Roe properties:
//   (1) Wave-sum consistency:   ∑_k w_k        = qR_local - qL_local
//   (2) Flux-jump consistency:  ∑_k s_k · w_k  = f_SR(qR_local) - f_SR(qL_local)
//   (3) Fluctuation balance:    amdq + apdq    = ∑_k s_k · w_k
//   (4) Eigenvalue ordering:    s[0] ≤ s[1] ≤ s[2]
//
// Empirically only (3) and (4) hold to floating-point precision. (1) and (2)
// fail by O(1) magnitudes on every spacetime — including pure Minkowski,
// where any genuine SR Roe construction would satisfy both exactly. This
// pins down that the "Roe" entry in wv_gr_euler{,_tetrad}.c is in fact a
// non-conservative approximate Riemann solver derived from the SR Roe
// algebra, NOT a true Roe solver. We assert (3) and (4) as hard checks and
// emit the (1)/(2) residuals as informational TEST_MSG so that any future
// improvement that turns this into a real Roe solver will visibly drive the
// residuals down. See the worst-case residual values printed in the test
// report for the current state.
static void
run_roe_properties(struct gkyl_gr_spacetime *spacetime)
{
  double gas_gamma = 5.0 / 3.0;

  int lower[1] = { 0 }, upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);

  struct gkyl_wv_eqn *gr_mod = gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = gas_gamma,
      .conf_range = conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_ROE,
      .use_gpu = false,
    });
  gkyl_gr_euler_tetrad_mod_set_auxfields(gr_mod,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  TEST_CHECK( gr_mod->num_waves == 3 );

  struct wv_gr_euler_tetrad_mod *grm = container_of(gr_mod,
    struct wv_gr_euler_tetrad_mod, eqn);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  // Tracked across the (x,y) sweep; emitted as informational TEST_MSG
  // after the loop completes. These properties do NOT hold for the current
  // approximate-Roe implementation (see header comment).
  double max_wsum_res = 0.0;
  double max_fj_res   = 0.0;

  for (int x_ind = -5; x_ind < 6; x_ind++) {
    for (int y_ind = -5; y_ind < 6; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double q_seed[71];
      double *prods_row = gkyl_array_fetch(prods, 0);
      fill_spacetime(spacetime, x, y, 0.0, q_seed, prods_row);

      if (q_seed[27] < 0.0) continue;

      // Two distinct hydro states at the same spacetime cell.
      double rho_l = 1.0, ul = 0.10, vl = 0.20, wl = 0.30, pl = 1.5;
      double rho_r = 0.5, ur = 0.05, vr = 0.10, wr = 0.15, pr = 0.7;

      double spatial_det = prods_row[GKYL_GR_SP_SPATIAL_DET];
      double vell[3] = { ul, vl, wl }, velr[3] = { ur, vr, wr };
      double vsq_l = 0.0, vsq_r = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          vsq_l += prods_row[GKYL_GR_SP_GIJ + 3*i + j] * vell[i] * vell[j];
          vsq_r += prods_row[GKYL_GR_SP_GIJ + 3*i + j] * velr[i] * velr[j];
        }
      double Wl = 1.0 / sqrt(1.0 - vsq_l);
      double Wr = 1.0 / sqrt(1.0 - vsq_r);
      double hl = 1.0 + (pl / rho_l) * (gas_gamma / (gas_gamma - 1.0));
      double hr = 1.0 + (pr / rho_r) * (gas_gamma / (gas_gamma - 1.0));

      double ql_mod[5], qr_mod[5];
      ql_mod[0] = sqrt(spatial_det) * rho_l * Wl;
      ql_mod[1] = sqrt(spatial_det) * rho_l * hl * (Wl*Wl) * ul;
      ql_mod[2] = sqrt(spatial_det) * rho_l * hl * (Wl*Wl) * vl;
      ql_mod[3] = sqrt(spatial_det) * rho_l * hl * (Wl*Wl) * wl;
      ql_mod[4] = sqrt(spatial_det) * ((rho_l * hl * (Wl*Wl)) - pl - (rho_l * Wl));

      qr_mod[0] = sqrt(spatial_det) * rho_r * Wr;
      qr_mod[1] = sqrt(spatial_det) * rho_r * hr * (Wr*Wr) * ur;
      qr_mod[2] = sqrt(spatial_det) * rho_r * hr * (Wr*Wr) * vr;
      qr_mod[3] = sqrt(spatial_det) * rho_r * hr * (Wr*Wr) * wr;
      qr_mod[4] = sqrt(spatial_det) * ((rho_r * hr * (Wr*Wr)) - pr - (rho_r * Wr));

      int idx_l[1] = { 0 }, idx_r[1] = { 0 };
      gr_mod->set_interface_idx_func(gr_mod, idx_l, idx_r);

      double ql_local[5], qr_local[5];
      gr_mod->rotate_to_local_func(gr_mod, tau1, tau2, norm, ql_mod, ql_local);
      gr_mod->rotate_to_local_func(gr_mod, tau1, tau2, norm, qr_mod, qr_local);

      double delta[5];
      for (int i = 0; i < 5; i++) delta[i] = qr_local[i] - ql_local[i];

      double waves[3 * 5], speeds[3];
      double maxs = gr_mod->waves_func(gr_mod, GKYL_WV_HIGH_ORDER_FLUX,
        delta, ql_local, qr_local, 1.0, 1.0, waves, speeds);

      // ---- (1) Wave-sum consistency: ∑ w_k = ∆q ----
      // True for any complete eigenvector basis with correctly-computed
      // wave amplitudes; the fix to a1, a2, a3 in wv_gr_euler_tetrad_mod.c
      // restored this identity (it failed by O(1) with the pre-fix amps).
      for (int i = 0; i < 5; i++) {
        double sum = 0.0;
        for (int k = 0; k < 3; k++) sum += waves[k * 5 + i];
        double res = fabs(sum - delta[i]);
        if (res > max_wsum_res) max_wsum_res = res;
      }

      // ---- (2) Flux-jump consistency: ∑ s_k · w_k = √γ · ∆f_SR ----
      // SR Roe satisfies A_Roe · ∆q_SR = ∆f_SR by construction. In the
      // tetrad framework the conserved state vector carries an explicit
      // √γ prefactor (D = √γ·ρW, τ = √γ·(ρhW²−p−ρW)), so the Roe
      // decomposition operates on √γ-scaled deltas:
      //   ∆q_GR = √γ · ∆q_SR  ⇒  ∑ s_k · w_k_GR = √γ · ∆f_SR.
      // gkyl_gr_euler_tetrad_mod_flux returns the bare flat-space SR flux
      // (no √γ), so we scale it by √γ here before the comparison. In
      // Minkowski √γ=1 and the scaling is a no-op.
      double fl_sr[5], fr_sr[5];
      gkyl_gr_euler_tetrad_mod_flux(gas_gamma, ql_local, grm->prodl_local, fl_sr);
      gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qr_local, grm->prodr_local, fr_sr);
      double sqrt_det_test = sqrt(grm->prodl_local[GKYL_GR_SP_SPATIAL_DET]);
      double df_sr[5];
      for (int i = 0; i < 5; i++) df_sr[i] = sqrt_det_test * (fr_sr[i] - fl_sr[i]);

      // Informational only — the implementation does NOT satisfy this.
      double sw[5] = {0};
      for (int k = 0; k < 3; k++)
        for (int i = 0; i < 5; i++)
          sw[i] += speeds[k] * waves[k * 5 + i];
      for (int i = 0; i < 5; i++) {
        double res = fabs(sw[i] - df_sr[i]);
        if (res > max_fj_res) max_fj_res = res;
      }

      // ---- (3) Fluctuation balance: amdq + apdq = ∑ s_k · w_k ----
      // qfluct_roe assembles amdq from negative-speed parts and apdq from
      // positive-speed parts of each wave; together they must reconstruct
      // ∑ s_k · w_k component-by-component.
      double amdq[5], apdq[5];
      gr_mod->qfluct_func(gr_mod, GKYL_WV_HIGH_ORDER_FLUX,
        ql_local, qr_local, 1.0, 1.0, waves, speeds, amdq, apdq);
      for (int i = 0; i < 5; i++) {
        TEST_CHECK( gkyl_compare(amdq[i] + apdq[i], sw[i], 1e-13) );
      }

      // ---- (4) Eigenvalue ordering: s[0] ≤ s[1] ≤ s[2] ----
      // The Roe wavespeeds for GR Euler are left-acoustic, contact, and
      // right-acoustic. Strict-inequality breakdowns at sonic points are not
      // expected for our test states (subsonic, distinct primitives), so a
      // monotone-non-decreasing check is appropriate.
      TEST_CHECK( speeds[0] <= speeds[1] );
      TEST_CHECK( speeds[1] <= speeds[2] );

      // max_abs_speed returned by waves_func must dominate every wavespeed.
      double max_abs = 0.0;
      for (int k = 0; k < 3; k++) if (fabs(speeds[k]) > max_abs) max_abs = fabs(speeds[k]);
      TEST_CHECK( maxs + 1e-12 >= max_abs );

      (void)grm;
    }
  }

  // With the corrected wave amplitudes AND eigenvectors (per the
  // Eulderink-Mellema derivation in wv_gr_euler_roe_derivation.md) both
  // Roe identities hold to floating-point precision. The flux-jump
  // residual involves a long pipeline of sqrt/divide ops on Roe-averaged
  // primitives, so 1e-9 is the achievable tolerance; wave-sum (which is
  // pure linear algebra) tightens to 1e-10.
  TEST_CHECK_( max_wsum_res < 1e-10,
    "Roe wave-sum residual:   max |∑ w_k − ∆q|        = %.3e", max_wsum_res );
  TEST_CHECK_( max_fj_res < 1e-9,
    "Roe flux-jump residual:  max |∑ s_k·w_k − ∆f_SR| = %.3e", max_fj_res );

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(gr_mod);
}

// Full-tetrad-chain validation test. The "tetrad mod" code does NOT
// explicitly perform a Lorentz transformation to a tetrad frame; it
// operates on the curved-coord-frame conserved variables directly,
// using an SR-form Roe eigenstructure plus a flux-factorization
// (flat-formula × geometric correction). This test asks the question:
// is the code's approach mathematically equivalent to an explicit
// tetrad transformation?
//
// Path A: feed q_GR (curved-frame, with √γ) to gr_mod->waves_func and
//          read out the waves.
//
// Path B: explicit tetrad chain ---
//   1. Build orthonormal triad ε^i_a from γ_ij via Cholesky:
//        γ = L L^T  ⇒  ε = (L^{-1})^T  ⇒  γ_ij ε^i_a ε^j_b = δ_ab
//   2. Transform q_GR → q_tet (strip √γ from D and τ; apply E = L^T to
//      momentum to get tetrad-frame momentum). Now q_tet looks like a
//      flat-space SR-conserved state.
//   3. Run pure SR Roe (same equation object, but pointed at a
//      Minkowski-flat prods array so prim_vars and waves_func use the
//      flat metric internally) → w_tet.
//   4. Transform w_tet → w_back (apply ε = (L^{-1})^T to momentum
//      components, multiply hydro slots by √γ).
//
// Assertion: |w_back − w_A| < 1e-9 component-by-component.
//
// In Minkowski (γ = I, √γ = 1) both paths are pure identity and the
// test is trivial. In curved spacetime, if they disagree, the code's
// "tetrad mod" is NOT actually doing a tetrad transformation — it's a
// different (possibly approximate) Riemann solver in curved-coord
// variables. The disagreement magnitude tells us how much.

// Cholesky decomposition of a 3×3 symmetric positive-definite matrix.
// Returns lower-triangular L such that g = L L^T.
static void
cholesky_3x3(const double g[3][3], double L[3][3])
{
  L[0][0] = sqrt(g[0][0]);
  L[1][0] = g[1][0] / L[0][0];
  L[1][1] = sqrt(g[1][1] - L[1][0]*L[1][0]);
  L[2][0] = g[2][0] / L[0][0];
  L[2][1] = (g[2][1] - L[2][0]*L[1][0]) / L[1][1];
  L[2][2] = sqrt(g[2][2] - L[2][0]*L[2][0] - L[2][1]*L[2][1]);
  L[0][1] = L[0][2] = L[1][2] = 0.0;
}

// Invert a 3×3 lower-triangular matrix.
static void
invert_lower_3x3(const double L[3][3], double Linv[3][3])
{
  Linv[0][0] = 1.0 / L[0][0];
  Linv[1][1] = 1.0 / L[1][1];
  Linv[2][2] = 1.0 / L[2][2];
  Linv[1][0] = -L[1][0] * Linv[0][0] / L[1][1];
  Linv[2][1] = -L[2][1] * Linv[1][1] / L[2][2];
  Linv[2][0] = -(L[2][0] * Linv[0][0] + L[2][1] * Linv[1][0]) / L[2][2];
  Linv[0][1] = Linv[0][2] = Linv[1][2] = 0.0;
}

// Apply matrix M (3×3) to a 3-vector v: out = M v.
static void
mat_vec_3x3(const double M[3][3], const double v[3], double out[3])
{
  for (int i = 0; i < 3; i++) {
    out[i] = M[i][0]*v[0] + M[i][1]*v[1] + M[i][2]*v[2];
  }
}

// Build a Minkowski prods row (lapse=1, shift=0, γ_ij=δ_ij, etc.).
static void
fill_minkowski_prods(double *prods)
{
  for (int k = 0; k < GKYL_GR_SP_NCOMP_BASE; k++) prods[k] = 0.0;
  prods[GKYL_GR_SP_LAPSE] = 1.0;
  for (int i = 0; i < 3; i++) {
    prods[GKYL_GR_SP_GIJ + 3*i + i]     = 1.0;
    prods[GKYL_GR_SP_INV_GIJ + 3*i + i] = 1.0;
  }
  prods[GKYL_GR_SP_SPATIAL_DET] = 1.0;
  prods[GKYL_GR_SP_EXCISION]    = 1.0;  // not excised
}

static void
run_full_tetrad_chain_validation(struct gkyl_gr_spacetime *spacetime)
{
  double gas_gamma = 5.0 / 3.0;

  int lower[1] = { 0 }, upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_array *prods_curved = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  struct gkyl_array *prods_mink = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);

  // Curved-frame equation: receives q_GR.
  struct gkyl_wv_eqn *gr_curved = gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = gas_gamma, .conf_range = conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_ROE, .use_gpu = false,
    });
  gkyl_gr_euler_tetrad_mod_set_auxfields(gr_curved,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods_curved });

  // Tetrad-frame equation: receives q_tet, will see Minkowski via prods_mink.
  struct gkyl_wv_eqn *gr_mink = gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = gas_gamma, .conf_range = conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_ROE, .use_gpu = false,
    });
  gkyl_gr_euler_tetrad_mod_set_auxfields(gr_mink,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods_mink });

  // Pre-fill Minkowski prods (lapse=1, γ=I, √γ=1).
  double *prods_mink_row = gkyl_array_fetch(prods_mink, 0);
  fill_minkowski_prods(prods_mink_row);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double max_chain_diff = 0.0;

  for (int x_ind = -5; x_ind < 6; x_ind++) {
    for (int y_ind = -5; y_ind < 6; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      // Fill the curved prods at this point (skip excised cells).
      double q_seed[71];
      double *prods_curved_row = gkyl_array_fetch(prods_curved, 0);
      fill_spacetime(spacetime, x, y, 0.0, q_seed, prods_curved_row);
      if (q_seed[27] < 0.0) continue;

      // ----- Build a Riemann problem on q_GR (curved frame) -----
      double sqrt_det = sqrt(prods_curved_row[GKYL_GR_SP_SPATIAL_DET]);
      double g_ij[3][3];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          g_ij[i][j] = prods_curved_row[GKYL_GR_SP_GIJ + 3*i + j];

      double rho_l = 1.0, ul = 0.10, vl = 0.20, wl = 0.30, pl = 1.5;
      double rho_r = 0.5, ur = 0.05, vr = 0.10, wr = 0.15, pr = 0.7;

      double vell[3] = { ul, vl, wl }, velr[3] = { ur, vr, wr };
      double vsq_l = 0.0, vsq_r = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          vsq_l += g_ij[i][j] * vell[i] * vell[j];
          vsq_r += g_ij[i][j] * velr[i] * velr[j];
        }
      double Wl = 1.0 / sqrt(1.0 - vsq_l);
      double Wr = 1.0 / sqrt(1.0 - vsq_r);
      double hl = 1.0 + (pl / rho_l) * (gas_gamma / (gas_gamma - 1.0));
      double hr = 1.0 + (pr / rho_r) * (gas_gamma / (gas_gamma - 1.0));

      // Curved-frame conserved (S_i is lowered with γ_ij).
      double Sl_lower[3], Sr_lower[3];
      for (int i = 0; i < 3; i++) {
        Sl_lower[i] = 0.0;
        Sr_lower[i] = 0.0;
        for (int j = 0; j < 3; j++) {
          Sl_lower[i] += g_ij[i][j] * vell[j];
          Sr_lower[i] += g_ij[i][j] * velr[j];
        }
        Sl_lower[i] *= sqrt_det * rho_l * hl * (Wl*Wl);
        Sr_lower[i] *= sqrt_det * rho_r * hr * (Wr*Wr);
      }
      double ql_GR[5] = {
        sqrt_det * rho_l * Wl,
        Sl_lower[0], Sl_lower[1], Sl_lower[2],
        sqrt_det * ((rho_l * hl * (Wl*Wl)) - pl - (rho_l * Wl))
      };
      double qr_GR[5] = {
        sqrt_det * rho_r * Wr,
        Sr_lower[0], Sr_lower[1], Sr_lower[2],
        sqrt_det * ((rho_r * hr * (Wr*Wr)) - pr - (rho_r * Wr))
      };

      // ===== Path A: code's curved-frame Roe =====
      int idx_l[1] = { 0 }, idx_r[1] = { 0 };
      gr_curved->set_interface_idx_func(gr_curved, idx_l, idx_r);
      double ql_loc_A[5], qr_loc_A[5];
      gr_curved->rotate_to_local_func(gr_curved, tau1, tau2, norm, ql_GR, ql_loc_A);
      gr_curved->rotate_to_local_func(gr_curved, tau1, tau2, norm, qr_GR, qr_loc_A);
      double delta_A[5], waves_A[3 * 5], speeds_A[3];
      for (int i = 0; i < 5; i++) delta_A[i] = qr_loc_A[i] - ql_loc_A[i];
      gr_curved->waves_func(gr_curved, GKYL_WV_HIGH_ORDER_FLUX,
        delta_A, ql_loc_A, qr_loc_A, 1.0, 1.0, waves_A, speeds_A);

      // ===== Path B: explicit tetrad chain =====
      // Step 1: build triad ε^i_a from γ_ij via Cholesky.
      //   γ = L L^T, then ε = (L^{-1})^T  ⇒  γ_ij ε^i_a ε^j_b = δ_ab.
      double L[3][3], Linv[3][3];
      cholesky_3x3(g_ij, L);
      invert_lower_3x3(L, Linv);
      // ε^i_a = (L^{-1})^T at row i, column a: ε[i][a] = Linv[a][i].
      // E^a_i = ε^{-1} = L^T at row a, column i: E[a][i] = L[i][a].
      double eps[3][3], E[3][3];
      for (int i = 0; i < 3; i++) for (int a = 0; a < 3; a++) {
        eps[i][a] = Linv[a][i];
        E[a][i]   = L[i][a];
      }

      // Step 2: transform q_GR → q_tet at the cell center.
      // D_tet = D_GR / √γ;   τ_tet = τ_GR / √γ
      // S^a_tet = E^a_i (γ^ij S_j_GR / √γ)
      //   In curved-coord with S_i_GR = √γ · ρhW² · v_i  (covariant),
      //   the SR S^i = γ^ij · S_j_GR/√γ = ρhW² v^i (contravariant).
      //   The tetrad S^a_tet = E^a_i · ρhW² v^i = ρhW² v^a.
      double inv_g[3][3];
      // Compute γ^ij from γ_ij directly (assumes 3×3 SPD; consistent
      // with what cholesky_3x3 wants).
      {
        double det_g = g_ij[0][0]*(g_ij[1][1]*g_ij[2][2] - g_ij[1][2]*g_ij[2][1])
                     - g_ij[0][1]*(g_ij[1][0]*g_ij[2][2] - g_ij[1][2]*g_ij[2][0])
                     + g_ij[0][2]*(g_ij[1][0]*g_ij[2][1] - g_ij[1][1]*g_ij[2][0]);
        inv_g[0][0] =  (g_ij[1][1]*g_ij[2][2] - g_ij[1][2]*g_ij[2][1]) / det_g;
        inv_g[0][1] = -(g_ij[0][1]*g_ij[2][2] - g_ij[0][2]*g_ij[2][1]) / det_g;
        inv_g[0][2] =  (g_ij[0][1]*g_ij[1][2] - g_ij[0][2]*g_ij[1][1]) / det_g;
        inv_g[1][0] = inv_g[0][1];
        inv_g[1][1] =  (g_ij[0][0]*g_ij[2][2] - g_ij[0][2]*g_ij[2][0]) / det_g;
        inv_g[1][2] = -(g_ij[0][0]*g_ij[1][2] - g_ij[0][2]*g_ij[1][0]) / det_g;
        inv_g[2][0] = inv_g[0][2];
        inv_g[2][1] = inv_g[1][2];
        inv_g[2][2] =  (g_ij[0][0]*g_ij[1][1] - g_ij[0][1]*g_ij[1][0]) / det_g;
      }

      // Raise momentum, then apply E.
      double Sl_upper[3], Sr_upper[3], Sl_tet[3], Sr_tet[3];
      double Sl_lower_SR[3], Sr_lower_SR[3];
      for (int i = 0; i < 3; i++) {
        Sl_lower_SR[i] = ql_GR[i+1] / sqrt_det;
        Sr_lower_SR[i] = qr_GR[i+1] / sqrt_det;
      }
      for (int i = 0; i < 3; i++) {
        Sl_upper[i] = inv_g[i][0]*Sl_lower_SR[0] + inv_g[i][1]*Sl_lower_SR[1] + inv_g[i][2]*Sl_lower_SR[2];
        Sr_upper[i] = inv_g[i][0]*Sr_lower_SR[0] + inv_g[i][1]*Sr_lower_SR[1] + inv_g[i][2]*Sr_lower_SR[2];
      }
      mat_vec_3x3(E, Sl_upper, Sl_tet);
      mat_vec_3x3(E, Sr_upper, Sr_tet);

      double ql_tet[5] = {
        ql_GR[0] / sqrt_det,
        Sl_tet[0], Sl_tet[1], Sl_tet[2],
        ql_GR[4] / sqrt_det
      };
      double qr_tet[5] = {
        qr_GR[0] / sqrt_det,
        Sr_tet[0], Sr_tet[1], Sr_tet[2],
        qr_GR[4] / sqrt_det
      };

      // Step 3: run pure SR Roe (via gr_mink with Minkowski prods).
      gr_mink->set_interface_idx_func(gr_mink, idx_l, idx_r);
      double ql_loc_B[5], qr_loc_B[5];
      gr_mink->rotate_to_local_func(gr_mink, tau1, tau2, norm, ql_tet, ql_loc_B);
      gr_mink->rotate_to_local_func(gr_mink, tau1, tau2, norm, qr_tet, qr_loc_B);
      double delta_B[5], waves_tet[3 * 5], speeds_tet[3];
      for (int i = 0; i < 5; i++) delta_B[i] = qr_loc_B[i] - ql_loc_B[i];
      gr_mink->waves_func(gr_mink, GKYL_WV_HIGH_ORDER_FLUX,
        delta_B, ql_loc_B, qr_loc_B, 1.0, 1.0, waves_tet, speeds_tet);

      // Step 4: transform w_tet → w_back.
      //   w_back[0] (D)   = √γ · w_tet[0]
      //   w_back[i+1] (S_i, lowered, with √γ): the tetrad wave momentum
      //     components live in tetrad-frame slots Sx_tet, Sy_tet, Sz_tet.
      //     Map back: S_i_GR_wave = √γ · γ_ij · ε^j_a · w_tet[a+1]
      //                          = √γ · (γ ε)^a_i · w_tet[a+1]
      //                          = √γ · L_i^a · w_tet[a+1]   (since γ ε = L)
      //   w_back[4] (τ)   = √γ · w_tet[4]
      double waves_back[3 * 5];
      for (int k = 0; k < 3; k++) {
        const double *w_tet = &waves_tet[k * 5];
        double *w_back      = &waves_back[k * 5];
        w_back[0] = sqrt_det * w_tet[0];
        for (int i = 0; i < 3; i++) {
          // S_i_back = √γ · L_i^a w_tet[a+1]
          w_back[i+1] = sqrt_det * (L[i][0]*w_tet[1] + L[i][1]*w_tet[2] + L[i][2]*w_tet[3]);
        }
        w_back[4] = sqrt_det * w_tet[4];
      }

      // ----- Compare A vs B -----
      // Note: rotate_to_local on q_GR rotates the lowered momentum
      // (Sx, Sy, Sz) in Cartesian fashion; same operation on q_tet
      // rotates the tetrad-frame momentum. The triad construction at
      // the cell center is Cartesian-rotation-invariant (we built it
      // from γ_ij which is the same on both sides), so we can compare
      // waves at the rotated frame slot-by-slot.
      //
      // For an honest comparison we need both paths to share the same
      // pre-rotation interpretation of the slots; they do, modulo the
      // tetrad transform on slots 1..3.
      for (int i = 0; i < 5 * 3; i++) {
        double d = fabs(waves_back[i] - waves_A[i]);
        if (d > max_chain_diff) max_chain_diff = d;
      }
    }
  }

  // Minkowski: trivially equal. Curved: pass iff the code's "tetrad mod"
  // is mathematically equivalent to an explicit tetrad transformation.
  TEST_CHECK_( max_chain_diff < 1e-9,
    "Tetrad-chain residual: max |w_chain − w_code| = %.3e", max_chain_diff );

  gkyl_array_release(prods_curved);
  gkyl_array_release(prods_mink);
  gkyl_wv_eqn_release(gr_curved);
  gkyl_wv_eqn_release(gr_mink);
}

void
test_gr_euler_tetrad_mod_full_chain_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_full_tetrad_chain_validation(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_full_chain_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_full_tetrad_chain_validation(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_full_chain_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_full_tetrad_chain_validation(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_roe_properties_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_roe_properties(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_roe_properties_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_roe_properties(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_roe_properties_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_roe_properties(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_lax_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_lax_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_lax_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_hll_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_hll_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_hll_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_roe_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_ROE, 3);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_roe_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_ROE, 3);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_tetrad_mod_roe_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_TETRAD_RP_ROE, 3);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_euler_tetrad_mod_construction",           test_gr_euler_tetrad_mod_construction },
  { "gr_euler_tetrad_mod_rotation_minkowski",     test_gr_euler_tetrad_mod_rotation_minkowski },
  { "gr_euler_tetrad_mod_rotation_schwarzschild", test_gr_euler_tetrad_mod_rotation_schwarzschild },
  { "gr_euler_tetrad_mod_rotation_kerr",          test_gr_euler_tetrad_mod_rotation_kerr },
  { "gr_euler_tetrad_mod_lax_minkowski",          test_gr_euler_tetrad_mod_lax_minkowski },
  { "gr_euler_tetrad_mod_lax_schwarzschild",      test_gr_euler_tetrad_mod_lax_schwarzschild },
  { "gr_euler_tetrad_mod_lax_kerr",               test_gr_euler_tetrad_mod_lax_kerr },
  { "gr_euler_tetrad_mod_hll_minkowski",          test_gr_euler_tetrad_mod_hll_minkowski },
  { "gr_euler_tetrad_mod_hll_schwarzschild",      test_gr_euler_tetrad_mod_hll_schwarzschild },
  { "gr_euler_tetrad_mod_hll_kerr",               test_gr_euler_tetrad_mod_hll_kerr },
  { "gr_euler_tetrad_mod_roe_minkowski",          test_gr_euler_tetrad_mod_roe_minkowski },
  { "gr_euler_tetrad_mod_roe_schwarzschild",      test_gr_euler_tetrad_mod_roe_schwarzschild },
  { "gr_euler_tetrad_mod_roe_kerr",               test_gr_euler_tetrad_mod_roe_kerr },
  { "gr_euler_tetrad_mod_roe_properties_minkowski",     test_gr_euler_tetrad_mod_roe_properties_minkowski },
  { "gr_euler_tetrad_mod_roe_properties_schwarzschild", test_gr_euler_tetrad_mod_roe_properties_schwarzschild },
  { "gr_euler_tetrad_mod_roe_properties_kerr",          test_gr_euler_tetrad_mod_roe_properties_kerr },
  { "gr_euler_tetrad_mod_full_chain_minkowski",         test_gr_euler_tetrad_mod_full_chain_minkowski },
  { "gr_euler_tetrad_mod_full_chain_schwarzschild",     test_gr_euler_tetrad_mod_full_chain_schwarzschild },
  { "gr_euler_tetrad_mod_full_chain_kerr",              test_gr_euler_tetrad_mod_full_chain_kerr },
  { NULL, NULL },
};
