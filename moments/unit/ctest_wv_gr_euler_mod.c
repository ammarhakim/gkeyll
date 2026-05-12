// Unit tests for the modular GR Euler equation object (wv_gr_euler_mod).
// Strategy: every test sets up both the packed wv_gr_euler equation (71-comp
// state with spacetime baked in) and the mod equation (5-comp hydro state
// plus an external spacetime-products array) from the same spacetime
// callbacks at the same physical point, then verifies that the mod equation
// produces results equivalent to the packed equation on the components they
// share.
//
// This is the primary correctness mechanism for the refactor: any divergence
// here means the rotation/Riemann math drifted between the two paths and
// the <1e-12 acceptance gate in Phase A will not hold.

#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_mod.h>
#include <gkyl_wv_gr_euler_mod_priv.h>
#include <gkyl_wv_gr_euler_priv.h>

// Helper: pack the packed-layout spacetime block (q[5..66]) and the mod
// products row (prods[0..NCOMP_BASE)) from the spacetime callbacks at the
// physical point (x,y,z). The two layouts are deliberately aligned: index k
// in prods (0..61) corresponds to packed q[5+k]; coords sit at q[68..70] and
// prods[GKYL_GR_SP_CELLCTR..]; the derived inv-metric and determinant blocks
// in prods have no packed analog and are computed here.
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

  // Pack packed-q spacetime block (q[5..66]) and coords (q[67..70]).
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

  // Pack the mod products row. Indices 0..61 of prods coincide with q[5..66].
  for (int k = 0; k < 62; k++) prods[k] = q[5 + k];

  // Derived: inverse spatial metric and determinant. The packed equation
  // recomputes inv γ_ij inline at every flux/source call; the mod pipeline
  // precomputes it once.
  gkyl_gr_euler_inv_spatial_metric(q, &inv_spatial_metric);
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

// Helper: check that packed q_local[5..66] and the mod-equation's per-side
// scratch buffer prodX_local[0..61] hold the same locally-rotated spacetime
// data, plus the coordinate-rotation match. Tolerance is tight (1e-12) since
// the two paths execute the same arithmetic on the same input.
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

// Common runner: for a given spacetime, sweep over an (x,y) grid and check
// rotation equivalence between packed and mod equations at every point and
// every face direction. The hydro state is the same fixed values used by the
// existing packed tests so test points stay inside the physically valid
// regime for all three spacetimes.
static void
run_rotation_equivalence(struct gkyl_gr_spacetime *spacetime, int half_extent,
  double step)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_wv_eqn *gr_euler_packed =
    gkyl_wv_gr_euler_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  // Single-cell range used for cidx lookups inside the mod equation.
  int lower[1] = { 0 };
  int upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  struct gkyl_wv_eqn *gr_euler_mod =
    gkyl_wv_gr_euler_mod_new(gas_gamma, &conf_range, false);
  gkyl_gr_euler_mod_set_auxfields(gr_euler_mod,
    (struct gkyl_wv_gr_euler_mod_auxfields){ .prods = prods });

  TEST_CHECK( gr_euler_mod->num_equations == 5 );
  TEST_CHECK( gr_euler_packed->num_equations == 71 );
  TEST_CHECK( gr_euler_mod->set_interface_idx_func != NULL );
  TEST_CHECK( gr_euler_mod->set_cell_idx_func != NULL );

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

  struct wv_gr_euler_mod *grm =
    container_of(gr_euler_mod, struct wv_gr_euler_mod, eqn);

  for (int x_ind = -half_extent; x_ind < half_extent + 1; x_ind++) {
    for (int y_ind = -half_extent; y_ind < half_extent + 1; y_ind++) {
      double x = step * x_ind;
      double y = step * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double q[71];
      double *prods_row = gkyl_array_fetch(prods, 0);
      fill_spacetime(spacetime, x, y, 0.0, q, prods_row);

      // Skip excision-region cells: the packed flow zeros the state inside
      // them and rotations there are not meaningful.
      if (q[27] < 0.0) continue;

      double spatial_det = prods_row[GKYL_GR_SP_SPATIAL_DET];
      double lapse = prods_row[GKYL_GR_SP_LAPSE];
      (void)lapse;  // not directly used below

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

      // Mod hydro state is the same 5-element prefix of q.
      double q_mod[5] = { q[0], q[1], q[2], q[3], q[4] };

      int idx_l[1] = { 0 };
      int idx_r[1] = { 0 };

      // ---- Per-cell equivalence checks (global frame, no rotation) ----
      // prim_vars recovery.
      double v_packed[71], v_mod[5];
      gkyl_gr_euler_prim_vars(gas_gamma, q, v_packed);
      gkyl_gr_euler_mod_prim_vars(gas_gamma, q_mod, prods_row, v_mod);
      for (int i = 0; i < 5; i++)
        TEST_CHECK( gkyl_compare(v_mod[i], v_packed[i], 1e-12) );

      // max_speed equivalence via the function-pointer interface (which is
      // what wave_prop actually calls). set_cell_idx_func is required for
      // the mod path.
      gr_euler_mod->set_cell_idx_func(gr_euler_mod, idx_l);
      double maxs_packed = gr_euler_packed->max_speed_func(gr_euler_packed, q);
      double maxs_mod    = gr_euler_mod->max_speed_func(gr_euler_mod, q_mod);
      TEST_CHECK( gkyl_compare(maxs_mod, maxs_packed, 1e-12) );

      // check_inv equivalence (both should return true for these test inputs).
      gr_euler_mod->set_cell_idx_func(gr_euler_mod, idx_l);
      bool inv_packed = gr_euler_packed->check_inv_func(gr_euler_packed, q);
      bool inv_mod    = gr_euler_mod->check_inv_func(gr_euler_mod, q_mod);
      TEST_CHECK( inv_packed == inv_mod );

      // ---- Interface-level equivalence checks ----
      // Mod uses both LEFT and RIGHT slots; in this single-cell test both
      // point at cidx 0. We sequence set_interface_idx_func -> rotate(left)
      // -> rotate(right) the same way wave_prop does.
      gr_euler_mod->set_interface_idx_func(gr_euler_mod, idx_l, idx_r);

      for (int d = 0; d < 3; d++) {
        // Packed: rotate to local, then back to global; should be identity.
        double q_local_packed[71], q_global_back[71];
        gr_euler_packed->rotate_to_local_func(gr_euler_packed,
          tau1[d], tau2[d], norm[d], q, q_local_packed);
        gr_euler_packed->rotate_to_global_func(gr_euler_packed,
          tau1[d], tau2[d], norm[d], q_local_packed, q_global_back);
        for (int i = 0; i < 71; i++)
          TEST_CHECK( gkyl_compare(q[i], q_global_back[i], 1e-12) );

        // Mod: rotate-to-local twice (once for each side of the interface);
        // first call populates prodl_local, second populates prodr_local.
        // The parity is reset on each set_interface_idx_func call above.
        gr_euler_mod->set_interface_idx_func(gr_euler_mod, idx_l, idx_r);

        double q_local_mod_l[5], q_local_mod_r[5], q_global_mod[5];
        gr_euler_mod->rotate_to_local_func(gr_euler_mod,
          tau1[d], tau2[d], norm[d], q_mod, q_local_mod_l);
        gr_euler_mod->rotate_to_local_func(gr_euler_mod,
          tau1[d], tau2[d], norm[d], q_mod, q_local_mod_r);
        gr_euler_mod->rotate_to_global_func(gr_euler_mod,
          tau1[d], tau2[d], norm[d], q_local_mod_l, q_global_mod);

        // Hydro round-trip identity.
        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(q_mod[i], q_global_mod[i], 1e-12) );

        // Hydro rotation matches packed.
        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(q_local_mod_l[i], q_local_packed[i], 1e-12) );

        // Flux equivalence in the local frame: rotate q → flux → compare.
        // Both paths consume the locally-rotated spacetime (packed via
        // q_local_packed[5..66], mod via grm->prodl_local).
        double flux_packed[71], flux_mod[5];
        gkyl_gr_euler_flux(gas_gamma, q_local_packed, flux_packed);
        gkyl_gr_euler_mod_flux(gas_gamma, q_local_mod_l, grm->prodl_local, flux_mod);
        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(flux_mod[i], flux_packed[i], 1e-12) );

        // flux_jump_func via the function-pointer interface (this is what
        // wave_prop's F-wave path uses). Both sides of the interface are
        // the same in this single-cell test, so the jump is identically 0.
        double fj_packed[71], fj_mod[5];
        double amax_packed = gr_euler_packed->flux_jump(gr_euler_packed,
          q_local_packed, q_local_packed, fj_packed);
        double amax_mod    = gr_euler_mod->flux_jump(gr_euler_mod,
          q_local_mod_l, q_local_mod_l, fj_mod);
        TEST_CHECK( gkyl_compare(amax_mod, amax_packed, 1e-12) );
        for (int i = 0; i < 5; i++) {
          TEST_CHECK( gkyl_compare(fj_mod[i],     0.0, 1e-12) );
          TEST_CHECK( gkyl_compare(fj_packed[i],  0.0, 1e-12) );
        }

        // Locally-rotated spacetime products match packed q_local[5..66] for
        // both left and right scratch buffers (same input in single-cell run).
        check_rotated_spacetime(q_local_packed, grm->prodl_local);
        check_rotated_spacetime(q_local_packed, grm->prodr_local);

        // Derived blocks: spatial_det is invariant, INV_GIJ rotated should
        // satisfy INV * G = I in the rotated frame (internal consistency).
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

        // cons_to_riem cycle identity (currently identity on both paths).
        double w_packed[71], q_packed_cycle[71];
        gr_euler_packed->cons_to_riem(gr_euler_packed, q_local_packed, q_local_packed, w_packed);
        gr_euler_packed->riem_to_cons(gr_euler_packed, q_local_packed, w_packed, q_packed_cycle);
        for (int i = 0; i < 71; i++)
          TEST_CHECK( gkyl_compare(q_local_packed[i], q_packed_cycle[i], 1e-12) );

        double w_mod[5], q_mod_cycle[5];
        gr_euler_mod->cons_to_riem(gr_euler_mod, q_local_mod_l, q_local_mod_l, w_mod);
        gr_euler_mod->riem_to_cons(gr_euler_mod, q_local_mod_l, w_mod, q_mod_cycle);
        for (int i = 0; i < 5; i++)
          TEST_CHECK( gkyl_compare(q_local_mod_l[i], q_mod_cycle[i], 1e-12) );
      }
    }
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(gr_euler_mod);
  gkyl_wv_eqn_release(gr_euler_packed);
}

// Sanity that the constructor wires types, num_equations, num_waves, and the
// optional setters correctly.
static void
test_gr_euler_mod_construction()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  int lower[1] = { 0 };
  int upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_euler_mod_new(5.0/3.0, &conf_range, false);
  TEST_CHECK( eqn->type == GKYL_EQN_GR_EULER_MOD );
  TEST_CHECK( eqn->num_equations == 5 );
  TEST_CHECK( eqn->num_diag == 5 );
  TEST_CHECK( eqn->num_waves == 2 );      // HLL default
  TEST_CHECK( eqn->set_interface_idx_func != NULL );
  TEST_CHECK( eqn->set_cell_idx_func != NULL );
  TEST_CHECK( eqn->rotate_to_local_func != NULL );
  TEST_CHECK( eqn->rotate_to_global_func != NULL );
  TEST_CHECK( eqn->cons_to_diag != NULL );

  TEST_CHECK( gkyl_compare(gkyl_wv_gr_euler_mod_gas_gamma(eqn), 5.0/3.0, 1e-15) );

  // Different RP types yield different num_waves.
  struct gkyl_wv_eqn *eqn_lax = gkyl_wv_gr_euler_mod_inew(&(struct gkyl_wv_gr_euler_mod_inp){
    .gas_gamma = 1.4, .conf_range = conf_range,
    .rp_type = WV_GR_EULER_RP_LAX, .use_gpu = false,
  });
  TEST_CHECK( eqn_lax->num_waves == 2 );

  struct gkyl_wv_eqn *eqn_roe = gkyl_wv_gr_euler_mod_inew(&(struct gkyl_wv_gr_euler_mod_inp){
    .gas_gamma = 1.4, .conf_range = conf_range,
    .rp_type = WV_GR_EULER_RP_ROE, .use_gpu = false,
  });
  TEST_CHECK( eqn_roe->num_waves == 3 );

  gkyl_wv_eqn_release(eqn);
  gkyl_wv_eqn_release(eqn_lax);
  gkyl_wv_eqn_release(eqn_roe);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_rotation_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_rotation_equivalence(spacetime, 10, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_rotation_schwarzschild()
{
  // Same blackhole parameters as the packed test.
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_rotation_equivalence(spacetime, 10, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_rotation_kerr()
{
  // Same Kerr parameters as the packed test.
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_rotation_equivalence(spacetime, 10, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

// Riemann-solver equivalence: build a non-trivial 1D Riemann problem (left
// and right hydro states differ) sharing the spacetime at a single cell, run
// both packed and mod waves+qfluct at one face direction, and compare the
// hydro portion of each wave plus the fluctuations. Speeds must agree
// exactly (both equations compute the same max_abs_speed bound for Lax).
//
// For an N-wave solver, the packed `waves` array is stored num_equations
// stripes (71 per wave); the mod stores 5 per wave. We compare the first 5
// components of each stripe.
static void
run_riemann_equivalence(struct gkyl_gr_spacetime *spacetime,
  enum gkyl_wv_gr_euler_rp rp_type, int num_waves)
{
  double gas_gamma = 5.0 / 3.0;

  // Packed equation with the chosen Riemann solver.
  struct gkyl_wv_eqn *gr_euler_packed = gkyl_wv_gr_euler_inew(
    &(struct gkyl_wv_gr_euler_inp){
      .gas_gamma = gas_gamma,
      .spacetime_gauge = GKYL_STATIC_GAUGE,
      .reinit_freq = 0,
      .spacetime = spacetime,
      .rp_type = rp_type,
      .use_gpu = false,
    });

  // Single-cell range for the mod equation.
  int lower[1] = { 0 }, upper[1] = { 0 };
  struct gkyl_range conf_range;
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);

  struct gkyl_wv_eqn *gr_euler_mod = gkyl_wv_gr_euler_mod_inew(
    &(struct gkyl_wv_gr_euler_mod_inp){
      .gas_gamma = gas_gamma,
      .conf_range = conf_range,
      .rp_type = rp_type,
      .use_gpu = false,
    });
  gkyl_gr_euler_mod_set_auxfields(gr_euler_mod,
    (struct gkyl_wv_gr_euler_mod_auxfields){ .prods = prods });

  TEST_CHECK( gr_euler_mod->num_waves == num_waves );
  TEST_CHECK( gr_euler_packed->num_waves == num_waves );

  struct wv_gr_euler_mod *grm = container_of(gr_euler_mod,
    struct wv_gr_euler_mod, eqn);

  // Single direction (x). The rotation tests already exercise all three; this
  // test focuses on the solver math.
  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  // Sweep a small (x,y) grid; for each point use the same spacetime but
  // distinct left/right hydro states, mirroring how a Riemann problem looks
  // at an interface between two cells that happen to share spacetime.
  for (int x_ind = -5; x_ind < 6; x_ind++) {
    for (int y_ind = -5; y_ind < 6; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double ql_packed[71], qr_packed[71];
      double *prods_row = gkyl_array_fetch(prods, 0);
      fill_spacetime(spacetime, x, y, 0.0, ql_packed, prods_row);

      // ql and qr share the spacetime layout (single cell), so copy.
      for (int i = 5; i < 71; i++) qr_packed[i] = ql_packed[i];

      if (ql_packed[27] < 0.0) continue;

      // Different hydro states. Use modest velocities to stay subsonic and
      // inside the physically valid GR-Euler regime for all three spacetimes.
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

      // Rotate to local frame on both paths.
      double ql_packed_local[71], qr_packed_local[71];
      gr_euler_packed->rotate_to_local_func(gr_euler_packed, tau1, tau2, norm, ql_packed, ql_packed_local);
      gr_euler_packed->rotate_to_local_func(gr_euler_packed, tau1, tau2, norm, qr_packed, qr_packed_local);

      int idx_l[1] = { 0 }, idx_r[1] = { 0 };
      gr_euler_mod->set_interface_idx_func(gr_euler_mod, idx_l, idx_r);
      double ql_mod_local[5], qr_mod_local[5];
      gr_euler_mod->rotate_to_local_func(gr_euler_mod, tau1, tau2, norm, ql_mod, ql_mod_local);
      gr_euler_mod->rotate_to_local_func(gr_euler_mod, tau1, tau2, norm, qr_mod, qr_mod_local);

      // Delta needed by some Riemann solvers; Lax ignores it but pass anyway.
      double delta_packed[71], delta_mod[5];
      for (int i = 0; i < 71; i++) delta_packed[i] = qr_packed_local[i] - ql_packed_local[i];
      for (int i = 0; i <  5; i++) delta_mod[i]    = qr_mod_local[i]    - ql_mod_local[i];

      double waves_packed[3 /* max waves across Lax/HLL/Roe */ * 71];
      double speeds_packed[3 /* max waves across Lax/HLL/Roe */];
      double waves_mod[3 /* max waves across Lax/HLL/Roe */ * 5];
      double speeds_mod[3 /* max waves across Lax/HLL/Roe */];

      double maxs_packed = gr_euler_packed->waves_func(gr_euler_packed,
        GKYL_WV_HIGH_ORDER_FLUX, delta_packed,
        ql_packed_local, qr_packed_local, 1.0, 1.0,
        waves_packed, speeds_packed);
      double maxs_mod    = gr_euler_mod->waves_func(gr_euler_mod,
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

      // Fluctuations.
      double amdq_packed[71], apdq_packed[71];
      double amdq_mod[5],     apdq_mod[5];
      gr_euler_packed->qfluct_func(gr_euler_packed, GKYL_WV_HIGH_ORDER_FLUX,
        ql_packed_local, qr_packed_local, 1.0, 1.0,
        waves_packed, speeds_packed, amdq_packed, apdq_packed);
      gr_euler_mod->qfluct_func(gr_euler_mod, GKYL_WV_HIGH_ORDER_FLUX,
        ql_mod_local, qr_mod_local, 1.0, 1.0,
        waves_mod, speeds_mod, amdq_mod, apdq_mod);

      for (int i = 0; i < 5; i++) {
        TEST_CHECK( gkyl_compare(amdq_mod[i], amdq_packed[i], 1e-12) );
        TEST_CHECK( gkyl_compare(apdq_mod[i], apdq_packed[i], 1e-12) );
      }

      // Consistency: sum of waves equals the conserved-state jump (hydro
      // portion). True by construction for Lax (∑ w = ∆q) and HLL (qm − qL
      // + qR − qm = qR − qL). Not guaranteed for the GR Roe solver — its
      // eigenvector basis is the SR Roe linearization and does not span the
      // GR conserved jumps exactly, so we skip this check for Roe.
      if (rp_type != WV_GR_EULER_RP_ROE) {
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
  gkyl_wv_eqn_release(gr_euler_mod);
  gkyl_wv_eqn_release(gr_euler_packed);
}

void
test_gr_euler_mod_lax_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_LAX, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_lax_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_LAX, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_lax_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_LAX, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_hll_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_HLL, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_hll_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_HLL, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_hll_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_HLL, 2);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_roe_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_ROE, 3);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_roe_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_ROE, 3);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_mod_roe_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_equivalence(spacetime, WV_GR_EULER_RP_ROE, 3);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_euler_mod_construction",           test_gr_euler_mod_construction },
  { "gr_euler_mod_rotation_minkowski",     test_gr_euler_mod_rotation_minkowski },
  { "gr_euler_mod_rotation_schwarzschild", test_gr_euler_mod_rotation_schwarzschild },
  { "gr_euler_mod_rotation_kerr",          test_gr_euler_mod_rotation_kerr },
  { "gr_euler_mod_lax_minkowski",          test_gr_euler_mod_lax_minkowski },
  { "gr_euler_mod_lax_schwarzschild",      test_gr_euler_mod_lax_schwarzschild },
  { "gr_euler_mod_lax_kerr",               test_gr_euler_mod_lax_kerr },
  { "gr_euler_mod_hll_minkowski",          test_gr_euler_mod_hll_minkowski },
  { "gr_euler_mod_hll_schwarzschild",      test_gr_euler_mod_hll_schwarzschild },
  { "gr_euler_mod_hll_kerr",               test_gr_euler_mod_hll_kerr },
  { "gr_euler_mod_roe_minkowski",          test_gr_euler_mod_roe_minkowski },
  { "gr_euler_mod_roe_schwarzschild",      test_gr_euler_mod_roe_schwarzschild },
  { "gr_euler_mod_roe_kerr",               test_gr_euler_mod_roe_kerr },
  { NULL, NULL },
};
