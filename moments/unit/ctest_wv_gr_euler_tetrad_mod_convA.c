// Convention-A audit tests for wv_gr_euler_tetrad_mod (modular tetrad GR
// Euler) with q[1..3] = √γ·γ_ij·ρhW²·v^j (genuine covariant momentum).
//
// Three test groups:
//   1. Banyuls flux consistency. Build a state from primitives, run the
//      code's flux + flux_correction, compare to the analytic Banyuls
//      flux F^x = α√γ·(D·v̂^x, S_i·v̂^x + p·δ_i^x, τ·v̂^x + p·v^x) with
//      v̂^x = v^x − β^x/α.
//   2. Standard Riemann-solver properties for Lax and HLL: wave sum,
//      flux jump, fluctuation balance, eigenvalue ordering, max-abs-
//      speed dominance, trivial Riemann (qL = qR), L↔R sign symmetry.
//   3. Positivity sweep: curated 1D Riemann problems with hard near-
//      vacuum / near-luminal / cold-hot states; perform a single one-
//      sided Godunov update on each cell and check D > 0, s² > 0, τ > 0.
//
// All tests use Convention A consistently. There is no comparison to the
// packed tetrad path (which still uses Convention B).

#include <acutest.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_array.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler_prim_priv.h>
#include <gkyl_wv_gr_euler_tetrad_mod.h>
#include <gkyl_wv_gr_euler_tetrad_mod_priv.h>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Fill a single prods row from an analytic spacetime at (x, y, z).
static void
fill_prods_at(struct gkyl_gr_spacetime *spacetime,
  double x, double y, double z, double *prods)
{
  gkyl_moment_spacetime_coupling_fill_products_analytic(
    spacetime, 0.0, x, y, z, prods);
}

// Construct a Convention-A densitized state from primitives at the cell's
// spacetime products row. v_co is the contravariant 3-velocity v^i.
//
//   q[0]   = √γ · ρ · W
//   q[i+1] = √γ · γ_ij · ρhW² · v^j         (covariant momentum)
//   q[4]   = √γ · (ρhW² − p − ρW)
//
// W is computed self-consistently from |v|² = γ_ij v^i v^j.
static void
build_state_convA(double gas_gamma,
  double rho, const double v_co[3], double p,
  const double *prods, double q[5])
{
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);

  double v_lo[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v_lo[i] += prods[GKYL_GR_SP_GIJ + 3*i + j] * v_co[j];

  double vsq = 0.0;
  for (int i = 0; i < 3; i++) vsq += v_lo[i] * v_co[i];
  double W = 1.0 / sqrt(1.0 - vsq);
  double h = 1.0 + (p / rho) * (gas_gamma / (gas_gamma - 1.0));
  double rhohW2 = rho * h * W * W;

  q[0] = sqrt_det * rho * W;
  q[1] = sqrt_det * rhohW2 * v_lo[0];
  q[2] = sqrt_det * rhohW2 * v_lo[1];
  q[3] = sqrt_det * rhohW2 * v_lo[2];
  q[4] = sqrt_det * (rhohW2 - p - rho * W);
}

// Analytic Banyuls flux F^x[U] for the Convention-A state described by
// (rho, v_co, p) at the prods row. v̂^x = v^x − β^x/α.
//
//   F^x[D]   = α√γ · ρ W · v̂^x
//   F^x[S_i] = α√γ · (γ_ij·ρhW²·v^j · v̂^x + p · δ_i^x)
//   F^x[τ]   = α√γ · (τ · v̂^x + p · v^x)
static void
analytic_banyuls_flux(double gas_gamma,
  double rho, const double v_co[3], double p,
  const double *prods, double F[5])
{
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);
  double alpha    = prods[GKYL_GR_SP_LAPSE];
  double beta_x   = prods[GKYL_GR_SP_SHIFT + 0];

  double v_lo[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v_lo[i] += prods[GKYL_GR_SP_GIJ + 3*i + j] * v_co[j];

  double vsq = 0.0;
  for (int i = 0; i < 3; i++) vsq += v_lo[i] * v_co[i];
  double W = 1.0 / sqrt(1.0 - vsq);
  double h = 1.0 + (p / rho) * (gas_gamma / (gas_gamma - 1.0));
  double rhohW2 = rho * h * W * W;

  double v_tilde_x = v_co[0] - beta_x / alpha;
  double prefac = alpha * sqrt_det;

  F[0] = prefac * (rho * W) * v_tilde_x;
  F[1] = prefac * (rhohW2 * v_lo[0] * v_tilde_x + p);
  F[2] = prefac * (rhohW2 * v_lo[1] * v_tilde_x);
  F[3] = prefac * (rhohW2 * v_lo[2] * v_tilde_x);
  double tau = rhohW2 - p - rho * W;
  F[4] = prefac * (tau * v_tilde_x + p * v_co[0]);
}

// Set up a one-cell prods array bound to the equation. Returns the prods
// array (caller must release) and the prods row pointer.
static struct gkyl_array *
make_one_cell_prods(struct gkyl_wv_eqn *eqn, struct gkyl_range *conf_range,
  double *prods_row_copy)
{
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(conf_range, 1, lower, upper);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range->volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  (void)prods_row_copy;
  return prods;
}

static struct gkyl_wv_eqn *
make_eqn(double gas_gamma, struct gkyl_range conf_range,
  enum gkyl_wv_gr_euler_tetrad_rp rp)
{
  return gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = gas_gamma,
      .conf_range = conf_range,
      .rp_type = rp,
      .use_gpu = false,
    });
}

// ---------------------------------------------------------------------------
// 1. Banyuls flux consistency
// ---------------------------------------------------------------------------

static void
run_banyuls_flux_consistency(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z,
  enum gkyl_wv_gr_euler_tetrad_rp rp)
{
  double gas_gamma = 5.0 / 3.0;

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, x, y, z, prods_row);

  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }

  // A small grid of (rho, v_co, p) primitives.
  double rhos[]   = { 1.0, 0.5, 3.0 };
  double pres[]   = { 0.5, 1.5 };
  double vels[][3] = {
    { 0.10, 0.20, 0.30 },
    { 0.30, 0.10, 0.05 },
    { -0.20, 0.15, 0.40 },
    { 0.05, 0.05, 0.05 },
  };

  double max_diff = 0.0;
  for (size_t ir = 0; ir < sizeof(rhos)/sizeof(*rhos); ir++) {
    for (size_t ip = 0; ip < sizeof(pres)/sizeof(*pres); ip++) {
      for (size_t iv = 0; iv < sizeof(vels)/sizeof(*vels); iv++) {
        double rho = rhos[ir], p = pres[ip];
        const double *v = vels[iv];

        // Skip configurations that aren't subluminal at this metric.
        double v_lo[3] = { 0.0, 0.0, 0.0 };
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            v_lo[i] += prods_row[GKYL_GR_SP_GIJ + 3*i + j] * v[j];
        double vsq = 0.0;
        for (int i = 0; i < 3; i++) vsq += v_lo[i] * v[i];
        if (!(vsq < 1.0 - 1.0e-6)) continue;

        double q[5];
        build_state_convA(gas_gamma, rho, v, p, prods_row, q);

        double f_sr[5], f_gr[5];
        gkyl_gr_euler_tetrad_mod_flux(gas_gamma, q, prods_row, f_sr);
        gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, q, prods_row, f_sr, f_gr);

        double f_ref[5];
        analytic_banyuls_flux(gas_gamma, rho, v, p, prods_row, f_ref);

        for (int i = 0; i < 5; i++) {
          double d = fabs(f_gr[i] - f_ref[i]);
          if (d > max_diff) max_diff = d;
        }
      }
    }
  }

  TEST_CHECK_( max_diff < 1.0e-10,
    "[%s @ (%g,%g,%g)] Banyuls flux residual: max |F_code − F_analytic| = %.3e",
    label, x, y, z, max_diff );

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_banyuls_flux_consistency_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_banyuls_flux_consistency(st, "Mink-Lax", 0.3, 0.0, 0.0, WV_GR_EULER_TETRAD_RP_LAX);
  gkyl_gr_spacetime_release(st);
}
void test_banyuls_flux_consistency_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_banyuls_flux_consistency(st, "Schw-Lax", 0.3, 0.2, 0.0, WV_GR_EULER_TETRAD_RP_LAX);
  run_banyuls_flux_consistency(st, "Schw-Lax", 0.5, 0.0, 0.0, WV_GR_EULER_TETRAD_RP_LAX);
  run_banyuls_flux_consistency(st, "Schw-Lax", 0.4, 0.4, 0.0, WV_GR_EULER_TETRAD_RP_LAX);
  gkyl_gr_spacetime_release(st);
}
void test_banyuls_flux_consistency_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_banyuls_flux_consistency(st, "Kerr-Lax", 0.3, 0.2, 0.0, WV_GR_EULER_TETRAD_RP_LAX);
  run_banyuls_flux_consistency(st, "Kerr-Lax", 0.5, 0.0, 0.0, WV_GR_EULER_TETRAD_RP_LAX);
  run_banyuls_flux_consistency(st, "Kerr-Lax", 0.4, 0.4, 0.0, WV_GR_EULER_TETRAD_RP_LAX);
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 2. Standard Riemann-solver property tests for Lax / HLL
// ---------------------------------------------------------------------------

// Returns the curved-frame Banyuls flux ΔF for a given pair of states using
// the same (flux + flux_correction) machinery the wave construction uses.
static void
banyuls_delta_flux(double gas_gamma, struct wv_gr_euler_tetrad_mod *grm,
  const double qL[5], const double qR[5], double dF[5])
{
  double fL_sr[5], fR_sr[5], fL_gr[5], fR_gr[5];
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qL, grm->prodl_local, fL_sr);
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qR, grm->prodr_local, fR_sr);
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, qL, grm->prodl_local, fL_sr, fL_gr);
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, qR, grm->prodr_local, fR_sr, fR_gr);
  for (int i = 0; i < 5; i++) dF[i] = fR_gr[i] - fL_gr[i];
}

static void
run_riemann_properties_lax_hll(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z,
  enum gkyl_wv_gr_euler_tetrad_rp rp, int num_waves)
{
  double gas_gamma = 5.0 / 3.0;

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, x, y, z, prods_row);
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }

  TEST_CHECK( eqn->num_waves == num_waves );

  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  // Pick a single L/R primitive pair sufficient to exercise all properties.
  double rho_L = 1.0, p_L = 1.5;  double v_L[3] = { 0.10, 0.20, 0.30 };
  double rho_R = 0.5, p_R = 0.7;  double v_R[3] = { 0.05, 0.10, 0.15 };

  double qL_glob[5], qR_glob[5];
  build_state_convA(gas_gamma, rho_L, v_L, p_L, prods_row, qL_glob);
  build_state_convA(gas_gamma, rho_R, v_R, p_R, prods_row, qR_glob);

  int idx[1] = { 0 };
  eqn->set_interface_idx_func(eqn, idx, idx);

  double qL[5], qR[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL_glob, qL);
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR_glob, qR);

  double delta[5];
  for (int i = 0; i < 5; i++) delta[i] = qR[i] - qL[i];

  double waves[3 * 5], speeds[3];
  double maxs = eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    delta, qL, qR, 1.0, 1.0, waves, speeds);

  // (a) Wave sum: Σ w_k = Δq
  for (int i = 0; i < 5; i++) {
    double sum = 0.0;
    for (int k = 0; k < num_waves; k++) sum += waves[k * 5 + i];
    TEST_CHECK_( fabs(sum - delta[i]) < 1.0e-12,
      "[%s] wave sum: comp %d, |Σw − Δq| = %.3e", label, i, fabs(sum - delta[i]) );
  }

  // (b) Flux jump: Σ s_k · w_k = ΔF_GR (Banyuls)
  double dF[5];
  banyuls_delta_flux(gas_gamma, grm, qL, qR, dF);
  for (int i = 0; i < 5; i++) {
    double sw = 0.0;
    for (int k = 0; k < num_waves; k++) sw += speeds[k] * waves[k * 5 + i];
    TEST_CHECK_( fabs(sw - dF[i]) < 1.0e-10,
      "[%s] flux jump: comp %d, |Σs·w − ΔF| = %.3e", label, i, fabs(sw - dF[i]) );
  }

  // (c) Fluctuation balance: amdq + apdq = ΔF_GR
  double amdq[5], apdq[5];
  eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    qL, qR, 1.0, 1.0, waves, speeds, amdq, apdq);
  for (int i = 0; i < 5; i++) {
    TEST_CHECK_( fabs(amdq[i] + apdq[i] - dF[i]) < 1.0e-10,
      "[%s] fluct balance: comp %d, |amdq+apdq − ΔF| = %.3e",
      label, i, fabs(amdq[i] + apdq[i] - dF[i]) );
  }

  // (d) Eigenvalue ordering
  for (int k = 0; k + 1 < num_waves; k++) {
    TEST_CHECK_( speeds[k] <= speeds[k + 1] + 1.0e-14,
      "[%s] ordering: speeds[%d]=%g > speeds[%d]=%g",
      label, k, speeds[k], k + 1, speeds[k + 1] );
  }

  // (e) max-abs-speed dominates
  for (int k = 0; k < num_waves; k++) {
    TEST_CHECK_( maxs + 1.0e-12 >= fabs(speeds[k]),
      "[%s] maxs=%g < |speeds[%d]|=%g", label, maxs, k, fabs(speeds[k]) );
  }

  // (f) Trivial Riemann: qL = qR ⇒ all waves zero, fluctuations zero.
  double qE[5];
  build_state_convA(gas_gamma, 1.0, (double[]){0.1, 0.2, 0.05}, 1.0, prods_row, qE);
  double qE_loc[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qE, qE_loc);
  double dE[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  double wavesE[3 * 5], speedsE[3];
  double maxsE = eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    dE, qE_loc, qE_loc, 1.0, 1.0, wavesE, speedsE);
  for (int k = 0; k < num_waves; k++) {
    for (int i = 0; i < 5; i++) {
      TEST_CHECK_( fabs(wavesE[k * 5 + i]) < 1.0e-12,
        "[%s] trivial RP: wave[%d][%d] = %.3e", label, k, i, wavesE[k * 5 + i] );
    }
  }
  double amdqE[5], apdqE[5];
  eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    qE_loc, qE_loc, 1.0, 1.0, wavesE, speedsE, amdqE, apdqE);
  for (int i = 0; i < 5; i++) {
    TEST_CHECK_( fabs(amdqE[i]) < 1.0e-12 && fabs(apdqE[i]) < 1.0e-12,
      "[%s] trivial RP: |amdq[%d]|=%.3e |apdq[%d]|=%.3e",
      label, i, fabs(amdqE[i]), i, fabs(apdqE[i]) );
  }
  TEST_CHECK( isfinite(maxsE) );

  // (g) L↔R sign symmetry: swap states; amdq and apdq swap roles AND signs
  // for advection-style fluxes. The cleanest invariant is
  //     amdq(qL,qR) + apdq(qL,qR) = − [amdq(qR,qL) + apdq(qR,qL)]
  // since both sides equal ±ΔF and ΔF flips sign under L↔R swap.
  double waves_swap[3 * 5], speeds_swap[3];
  double delta_swap[5];
  for (int i = 0; i < 5; i++) delta_swap[i] = qL[i] - qR[i];
  eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    delta_swap, qR, qL, 1.0, 1.0, waves_swap, speeds_swap);
  double amdq_swap[5], apdq_swap[5];
  eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    qR, qL, 1.0, 1.0, waves_swap, speeds_swap, amdq_swap, apdq_swap);
  for (int i = 0; i < 5; i++) {
    double lhs = amdq[i] + apdq[i];
    double rhs = -(amdq_swap[i] + apdq_swap[i]);
    TEST_CHECK_( fabs(lhs - rhs) < 1.0e-10,
      "[%s] sign symmetry: comp %d, |ΔF + ΔF_swap| = %.3e",
      label, i, fabs(lhs - rhs) );
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_riemann_properties_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_riemann_properties_lax_hll(st, "Mink-Lax", 0.3, 0.0, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_properties_lax_hll(st, "Schw-Lax", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  run_riemann_properties_lax_hll(st, "Schw-Lax", 0.5, 0.0, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  run_riemann_properties_lax_hll(st, "Schw-Lax", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_properties_lax_hll(st, "Kerr-Lax", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  run_riemann_properties_lax_hll(st, "Kerr-Lax", 0.5, 0.0, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  run_riemann_properties_lax_hll(st, "Kerr-Lax", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(st);
}

void test_riemann_properties_hll_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_riemann_properties_lax_hll(st, "Mink-HLL", 0.3, 0.0, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_hll_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_properties_lax_hll(st, "Schw-HLL", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  run_riemann_properties_lax_hll(st, "Schw-HLL", 0.5, 0.0, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  run_riemann_properties_lax_hll(st, "Schw-HLL", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_hll_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_riemann_properties_lax_hll(st, "Kerr-HLL", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  run_riemann_properties_lax_hll(st, "Kerr-HLL", 0.5, 0.0, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  run_riemann_properties_lax_hll(st, "Kerr-HLL", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 3. Positivity preservation sweep
// ---------------------------------------------------------------------------

// One Riemann problem in primitive variables. v_co is the contravariant
// velocity; the rest are scalars.
struct rp_case {
  const char *name;
  double rho_L, p_L; double v_L[3];
  double rho_R, p_R; double v_R[3];
};

static void
record_admissibility(const double inv_g[3][3],
  const double q_undens[5],
  bool *D_ok, bool *S2_ok, bool *tau_ok)
{
  double D = q_undens[0], Sx = q_undens[1], Sy = q_undens[2];
  double Sz = q_undens[3], tau = q_undens[4];
  *D_ok   = (D > 0.0);
  *tau_ok = (tau >= 0.0);
  double mom_sq = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
  double s_sq = (D + tau) * (D + tau) - mom_sq;
  *S2_ok  = (s_sq > 0.0);
}

static void
run_positivity_for_rp(struct gkyl_wv_eqn *eqn,
  struct wv_gr_euler_tetrad_mod *grm,
  const double *prods,
  const struct rp_case *rc,
  int *D_violations, int *S2_violations, int *tau_violations)
{
  double gas_gamma = 5.0 / 3.0;
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g[3][3] = {
    { prods[GKYL_GR_SP_INV_GIJ + 0], prods[GKYL_GR_SP_INV_GIJ + 1], prods[GKYL_GR_SP_INV_GIJ + 2] },
    { prods[GKYL_GR_SP_INV_GIJ + 3], prods[GKYL_GR_SP_INV_GIJ + 4], prods[GKYL_GR_SP_INV_GIJ + 5] },
    { prods[GKYL_GR_SP_INV_GIJ + 6], prods[GKYL_GR_SP_INV_GIJ + 7], prods[GKYL_GR_SP_INV_GIJ + 8] },
  };

  // Build admissible L/R states.
  double qL[5], qR[5];
  build_state_convA(gas_gamma, rc->rho_L, rc->v_L, rc->p_L, prods, qL);
  build_state_convA(gas_gamma, rc->rho_R, rc->v_R, rc->p_R, prods, qR);

  // Sanity: both inputs admissible.
  double qL_und[5], qR_und[5];
  for (int i = 0; i < 5; i++) {
    qL_und[i] = qL[i] / sqrt_det;
    qR_und[i] = qR[i] / sqrt_det;
  }
  bool d_ok, s_ok, t_ok;
  record_admissibility(inv_g, qL_und, &d_ok, &s_ok, &t_ok);
  if (!(d_ok && s_ok && t_ok)) {
    fprintf(stderr, "  [positivity] WARNING: input qL inadmissible for case '%s' "
      "(D_ok=%d, S2_ok=%d, tau_ok=%d) — skipping\n", rc->name, d_ok, s_ok, t_ok);
    return;
  }
  record_admissibility(inv_g, qR_und, &d_ok, &s_ok, &t_ok);
  if (!(d_ok && s_ok && t_ok)) {
    fprintf(stderr, "  [positivity] WARNING: input qR inadmissible for case '%s' "
      "(D_ok=%d, S2_ok=%d, tau_ok=%d) — skipping\n", rc->name, d_ok, s_ok, t_ok);
    return;
  }

  // Compute waves and fluctuations at the interface. rotate_to_local
  // has the side-effect of filling prodl_local / prodr_local, which
  // wave_hll/wave_lax read. With norm = +x̂ this is a no-op rotation.
  int idx[1] = { 0 };
  eqn->set_interface_idx_func(eqn, idx, idx);
  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };
  double qL_loc[5], qR_loc[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL, qL_loc);
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR, qR_loc);

  double delta[5];
  for (int i = 0; i < 5; i++) delta[i] = qR_loc[i] - qL_loc[i];
  double waves[3 * 5], speeds[3];
  double maxs = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    delta, qL_loc, qR_loc, 1.0, 1.0, waves, speeds);
  double amdq[5], apdq[5];
  eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    qL_loc, qR_loc, 1.0, 1.0, waves, speeds, amdq, apdq);

  // CFL = 0.9 update against the interface.
  double cfl = 0.9;
  if (!(maxs > 0.0)) return;
  double dt_dx = cfl / maxs;

  double qL_new[5], qR_new[5];
  for (int i = 0; i < 5; i++) {
    qL_new[i] = qL_loc[i] - dt_dx * amdq[i];
    qR_new[i] = qR_loc[i] - dt_dx * apdq[i];
  }

  // Check both updated cells.
  for (int side = 0; side < 2; side++) {
    const double *qn = (side == 0) ? qL_new : qR_new;
    const char *who  = (side == 0) ? "L" : "R";
    double qn_und[5];
    for (int i = 0; i < 5; i++) qn_und[i] = qn[i] / sqrt_det;
    bool D_ok, S2_ok, tau_ok;
    record_admissibility(inv_g, qn_und, &D_ok, &S2_ok, &tau_ok);
    if (!D_ok)   { (*D_violations)++;   fprintf(stderr,
      "  [positivity] '%s' cell %s: D = %.3e\n", rc->name, who, qn_und[0]); }
    if (!S2_ok)  { (*S2_violations)++;  fprintf(stderr,
      "  [positivity] '%s' cell %s: s² = %.3e (D+τ = %.3e, |S|² = %.3e)\n",
      rc->name, who, (qn_und[0]+qn_und[4])*(qn_und[0]+qn_und[4])
                       - gkyl_gr_euler_mom_sq(inv_g, qn_und[1], qn_und[2], qn_und[3]),
      qn_und[0]+qn_und[4],
      gkyl_gr_euler_mom_sq(inv_g, qn_und[1], qn_und[2], qn_und[3])); }
    if (!tau_ok) { (*tau_violations)++; fprintf(stderr,
      "  [positivity] '%s' cell %s: τ = %.3e (D = %.3e)\n",
      rc->name, who, qn_und[4], qn_und[0]); }
  }
  (void)grm;
}

static const struct rp_case g_positivity_cases[] = {
  // Mild cases — should always preserve.
  { "uniform-mild",   1.0, 1.0, { 0.0, 0.0, 0.0 },  1.0, 1.0, { 0.0, 0.0, 0.0 } },
  { "smooth-shock",   2.0, 1.5, { 0.10, 0.05, 0.0 },  1.0, 0.7, { 0.0, 0.0, 0.0 } },
  { "smooth-rare",    1.0, 1.0, { -0.10, 0.0, 0.0 },  1.0, 1.0, {  0.10, 0.0, 0.0 } },

  // Strong shock: large rho/p ratio, head-on collision.
  { "strong-shock",   3.0, 5.0, { 0.30, 0.0, 0.0 },  0.5, 0.5, { -0.30, 0.0, 0.0 } },

  // Hot/cold contact.
  { "hot-cold",       1.0, 100.0, { 0.0, 0.0, 0.0 },  1.0, 0.001, { 0.0, 0.0, 0.0 } },

  // Near-vacuum: very low density on one side.
  { "near-vacuum-R",  1.0, 1.0,    { 0.0, 0.0, 0.0 },  1.0e-3, 1.0e-3, { 0.0, 0.0, 0.0 } },
  { "near-vacuum-L",  1.0e-3, 1.0e-3, { 0.0, 0.0, 0.0 },  1.0, 1.0,    { 0.0, 0.0, 0.0 } },

  // Strong rarefaction — pulls cells apart.
  { "strong-rare",    1.0, 1.0, { -0.50, 0.0, 0.0 },  1.0, 1.0, {  0.50, 0.0, 0.0 } },

  // Mildly relativistic shock.
  { "rel-shock",      5.0, 50.0, { 0.50, 0.0, 0.0 },  0.1, 0.01, { -0.50, 0.0, 0.0 } },
};

// Physical sample points for curved-spacetime positivity probes. Chosen
// to span on-axis, mid-radius, and corner geometry (off-diagonal γ from
// Kerr-Schild) without dipping into excision.
struct positivity_point { double x, y, z; };
static const struct positivity_point g_positivity_points[] = {
  { 0.0, 0.0, 0.0 },  // origin — Minkowski only; curved spacetimes filter
                       // via the in-excision-region check.
  { 0.3, 0.2, 0.0 },
  { 0.5, 0.0, 0.0 },
  { 0.4, 0.4, 0.0 },
};

static void
run_positivity_sweep(struct gkyl_gr_spacetime *spacetime,
  enum gkyl_wv_gr_euler_tetrad_rp rp,
  const char *label)
{
  double gas_gamma = 5.0 / 3.0;

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);

  int D_violations = 0, S2_violations = 0, tau_violations = 0;
  int total_cells = 0;

  int n_pts = sizeof(g_positivity_points) / sizeof(*g_positivity_points);
  int n_cases = sizeof(g_positivity_cases) / sizeof(*g_positivity_cases);

  double *prods_row = gkyl_array_fetch(prods, 0);

  for (int pi = 0; pi < n_pts; pi++) {
    const struct positivity_point *pt = &g_positivity_points[pi];
    fill_prods_at(spacetime, pt->x, pt->y, pt->z, prods_row);
    if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) continue;

    for (int c = 0; c < n_cases; c++) {
      run_positivity_for_rp(eqn, grm, prods_row, &g_positivity_cases[c],
        &D_violations, &S2_violations, &tau_violations);
      total_cells += 2;
    }
  }

  fprintf(stderr,
    "[positivity %s] across %d cell-updates: "
    "D<=0: %d (%.1f%%), s²<=0: %d (%.1f%%), τ<0: %d (%.1f%%)\n",
    label, total_cells,
    D_violations,   100.0 * D_violations   / total_cells,
    S2_violations,  100.0 * S2_violations  / total_cells,
    tau_violations, 100.0 * tau_violations / total_cells);

  // D and s² should provably be preserved by Lax/HLL — assert hard.
  TEST_CHECK_( D_violations  == 0, "D > 0 violated %d times for %s",
    D_violations, label );
  TEST_CHECK_( S2_violations == 0, "s² > 0 violated %d times for %s",
    S2_violations, label );
  // τ ≥ 0: report-only — Lax has no theoretical guarantee for τ.
  TEST_MSG( "τ ≥ 0 violated %d times for %s — info-only",
    tau_violations, label );

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_positivity_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_positivity_sweep(st, WV_GR_EULER_TETRAD_RP_LAX, "Lax-Mink");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_positivity_sweep(st, WV_GR_EULER_TETRAD_RP_LAX, "Lax-Schw");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_positivity_sweep(st, WV_GR_EULER_TETRAD_RP_LAX, "Lax-Kerr");
  gkyl_gr_spacetime_release(st);
}

void test_positivity_hll_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_positivity_sweep(st, WV_GR_EULER_TETRAD_RP_HLL, "HLL-Mink");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hll_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_positivity_sweep(st, WV_GR_EULER_TETRAD_RP_HLL, "HLL-Schw");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hll_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_positivity_sweep(st, WV_GR_EULER_TETRAD_RP_HLL, "HLL-Kerr");
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 4. Near-floor positivity sweep
//
// In the BHL run, the source step pushes cells to the τ-floor; the
// cascade-repair returns them at τ = GR_EULER_PRESSURE_FLOOR (currently
// 1e-14). Those cells then enter wave_prop with τ ≈ p_floor on one side
// and a typical state on the other. If the Lax/HLL fluctuation moves
// τ even slightly downward, the new τ goes below zero and wave_prop's
// post-positivity-sweep repair fires.
//
// This test sweeps tau_min ∈ {1e-14, …, 1e-6} on one side, with the
// other side a typical p ~ 1 state, and counts D / s² / τ violations
// per (solver, floor) pair. The output tells us the minimum tau_min
// at which Lax / HLL preserve admissibility.
// ---------------------------------------------------------------------------

struct near_floor_rp {
  const char *name;
  // Floored side (deliberately near the admissibility boundary).
  double rho_min, p_min;
  double v_min[3];
  // Typical-state side (drives the flux).
  double rho_typ, p_typ;
  double v_typ[3];
  bool floor_on_left;
};

static const struct near_floor_rp g_near_floor_cases[] = {
  // ρ_min held above the density floor (1e-8); only the pressure is
  // pushed near the τ-floor. v_min = 0 (the cell just got floored,
  // typically near-zero motion). The neighbor drives compression or
  // rarefaction.
  { "compress-from-right",  1.0e-4, 0.0, { 0.0, 0.0, 0.0 },
                            1.0,    1.0, { -0.30, 0.0, 0.0 },  true },
  { "compress-from-left",   1.0e-4, 0.0, { 0.0, 0.0, 0.0 },
                            1.0,    1.0, {  0.30, 0.0, 0.0 },  false },
  { "rarefy-to-right",      1.0e-4, 0.0, { 0.0, 0.0, 0.0 },
                            1.0,    1.0, {  0.30, 0.0, 0.0 },  true },
  { "rarefy-to-left",       1.0e-4, 0.0, { 0.0, 0.0, 0.0 },
                            1.0,    1.0, { -0.30, 0.0, 0.0 },  false },
  { "transverse-shear",     1.0e-4, 0.0, { 0.0, 0.0, 0.0 },
                            1.0,    1.0, {  0.10, 0.30, 0.0 }, true },
  { "hot-neighbor",         1.0e-4, 0.0, { 0.0, 0.0, 0.0 },
                            1.0,   50.0, {  0.0, 0.0, 0.0 },   true },
  // Floored cell has a slight outflow — tests whether the post-step
  // τ can be driven negative.
  { "small-outflow",        1.0e-4, 0.0, { 0.05, 0.0, 0.0 },
                            1.0,    1.0, {  0.30, 0.0, 0.0 },  true },
};

static void
run_near_floor_for_floor_value(struct gkyl_wv_eqn *eqn,
  struct wv_gr_euler_tetrad_mod *grm,
  const double *prods,
  double tau_min_target,
  int *D_v, int *S2_v, int *tau_v, int *total)
{
  double gas_gamma = 5.0 / 3.0;
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g[3][3] = {
    { prods[GKYL_GR_SP_INV_GIJ + 0], prods[GKYL_GR_SP_INV_GIJ + 1], prods[GKYL_GR_SP_INV_GIJ + 2] },
    { prods[GKYL_GR_SP_INV_GIJ + 3], prods[GKYL_GR_SP_INV_GIJ + 4], prods[GKYL_GR_SP_INV_GIJ + 5] },
    { prods[GKYL_GR_SP_INV_GIJ + 6], prods[GKYL_GR_SP_INV_GIJ + 7], prods[GKYL_GR_SP_INV_GIJ + 8] },
  };

  int n = sizeof(g_near_floor_cases) / sizeof(*g_near_floor_cases);
  for (int c = 0; c < n; c++) {
    const struct near_floor_rp *rc = &g_near_floor_cases[c];

    // Construct the "floored-side" state: target τ ≈ tau_min_target.
    // Cold-flow inversion (W≈1, h≈1) gives τ ≈ p/(γ-1), so set p_min
    // accordingly. Cap to ≥ p_min in the case rho_min was supplied.
    double p_min = (gas_gamma - 1.0) * tau_min_target;

    double qL[5], qR[5];
    if (rc->floor_on_left) {
      build_state_convA(gas_gamma, rc->rho_min, rc->v_min, p_min, prods, qL);
      build_state_convA(gas_gamma, rc->rho_typ, rc->v_typ, rc->p_typ, prods, qR);
    } else {
      build_state_convA(gas_gamma, rc->rho_typ, rc->v_typ, rc->p_typ, prods, qL);
      build_state_convA(gas_gamma, rc->rho_min, rc->v_min, p_min,    prods, qR);
    }

    // Skip if either side fails admissibility going in (shouldn't happen,
    // but defensive).
    {
      bool d_ok, s_ok, t_ok;
      double qu[5];
      for (int i = 0; i < 5; i++) qu[i] = qL[i] / sqrt_det;
      record_admissibility(inv_g, qu, &d_ok, &s_ok, &t_ok);
      if (!(d_ok && s_ok && t_ok)) continue;
      for (int i = 0; i < 5; i++) qu[i] = qR[i] / sqrt_det;
      record_admissibility(inv_g, qu, &d_ok, &s_ok, &t_ok);
      if (!(d_ok && s_ok && t_ok)) continue;
    }

    int idx[1] = { 0 };
    eqn->set_interface_idx_func(eqn, idx, idx);

    // rotate_to_local has the side-effect of filling prodl_local /
    // prodr_local on the equation, which wave_hll reads. Skipping this
    // call leaves uninitialized data and makes wave_hll mis-detect
    // excision. With norm = +x̂ this is a no-op rotation.
    double norm[3] = { 1.0, 0.0, 0.0 };
    double tau1[3] = { 0.0, 1.0, 0.0 };
    double tau2[3] = { 0.0, 0.0, 1.0 };
    double qL_loc[5], qR_loc[5];
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL, qL_loc);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR, qR_loc);

    double delta[5];
    for (int i = 0; i < 5; i++) delta[i] = qR_loc[i] - qL_loc[i];
    double waves[3 * 5], speeds[3];
    double maxs = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      delta, qL_loc, qR_loc, 1.0, 1.0, waves, speeds);
    double amdq[5], apdq[5];
    eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      qL_loc, qR_loc, 1.0, 1.0, waves, speeds, amdq, apdq);

    if (!(maxs > 0.0)) continue;
    double dt_dx = 0.95 / maxs;  // BHL uses cfl_frac = 0.95

    double qL_new[5], qR_new[5];
    for (int i = 0; i < 5; i++) {
      qL_new[i] = qL_loc[i] - dt_dx * amdq[i];
      qR_new[i] = qR_loc[i] - dt_dx * apdq[i];
    }

    for (int side = 0; side < 2; side++) {
      const double *qn = (side == 0) ? qL_new : qR_new;
      double qu[5];
      for (int i = 0; i < 5; i++) qu[i] = qn[i] / sqrt_det;
      bool D_ok, S2_ok, tau_ok;
      record_admissibility(inv_g, qu, &D_ok, &S2_ok, &tau_ok);
      if (!D_ok)   (*D_v)++;
      if (!S2_ok)  (*S2_v)++;
      if (!tau_ok) (*tau_v)++;
      (*total)++;
    }
  }
  (void)grm;
}

static void
run_near_floor_sweep(enum gkyl_wv_gr_euler_tetrad_rp rp,
  const char *rp_name)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, 0.0, 0.0, 0.0, prods_row);

  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);

  double tau_sweep[] = {
    1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10,
    1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14
  };

  fprintf(stderr, "[near-floor %s] tau_target  total  D<=0  s^2<=0  tau<0\n", rp_name);
  for (size_t i = 0; i < sizeof(tau_sweep)/sizeof(*tau_sweep); i++) {
    int D_v = 0, S2_v = 0, tau_v = 0, total = 0;
    run_near_floor_for_floor_value(eqn, grm, prods_row,
      tau_sweep[i], &D_v, &S2_v, &tau_v, &total);
    fprintf(stderr, "  %-8s    %.0e  %4d  %4d  %4d   %4d\n",
      rp_name, tau_sweep[i], total, D_v, S2_v, tau_v);
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
  gkyl_gr_spacetime_release(spacetime);
}

void test_near_floor_lax(void)
{
  run_near_floor_sweep(WV_GR_EULER_TETRAD_RP_LAX, "Lax");
}
void test_near_floor_hll(void)
{
  run_near_floor_sweep(WV_GR_EULER_TETRAD_RP_HLL, "HLL");
}

// ---------------------------------------------------------------------------
// 5. Audit tests for the contravariant-x triad correctness
// ---------------------------------------------------------------------------

// (a) Round-trip identity: q_curved → q_tet → q_curved should recover the
//     original (modulo √γ factor — q_tet is undensitized in scalars,
//     curved q is densitized). Catches sign/index errors in M, M_inv
//     that flux-jump wouldn't see if they cancel on Σ s·w.
static void
run_round_trip(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z)
{
  double gas_gamma = 5.0 / 3.0;
  double prods[GKYL_GR_SP_NCOMP_BASE];
  fill_prods_at(spacetime, x, y, z, prods);
  if (prods[GKYL_GR_SP_EXCISION] < 0.0) return;

  double rhos[] = { 1.0, 0.3 };
  double pres[] = { 0.5, 2.0 };
  double vels[][3] = {
    { 0.10, 0.20, 0.30 },
    { -0.20, 0.05, 0.0 },
    { 0.0, 0.5, 0.0 },         // pure transverse v
    { 0.0, 0.0, 0.0 },         // static
  };

  // Build triad.
  double g_ij[3][3], inv_g[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_ij[i][j]  = prods[GKYL_GR_SP_GIJ + 3*i + j];
      inv_g[i][j] = prods[GKYL_GR_SP_INV_GIJ + 3*i + j];
    }
  double M[3][3], M_inv[3][3];
  gkyl_gr_euler_tetrad_mod_build_triad_contravariant_x(g_ij, inv_g, M, M_inv);
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);

  // Sanity: M·M_inv = I, M^T·γ·M = I.
  double max_id_err = 0.0;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      double id_ab = M_inv[a][0]*M[0][b] + M_inv[a][1]*M[1][b] + M_inv[a][2]*M[2][b];
      double target = (a == b) ? 1.0 : 0.0;
      double e = fabs(id_ab - target);
      if (e > max_id_err) max_id_err = e;
    }
  }
  TEST_CHECK_( max_id_err < 1.0e-12,
    "[%s @ (%g,%g,%g)] M_inv·M ≠ I: max err = %.3e",
    label, x, y, z, max_id_err );

  // M^T γ M = I (orthonormality).
  double max_orth_err = 0.0;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      double inner = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          inner += g_ij[i][j] * M[i][a] * M[j][b];
      double target = (a == b) ? 1.0 : 0.0;
      double e = fabs(inner - target);
      if (e > max_orth_err) max_orth_err = e;
    }
  }
  TEST_CHECK_( max_orth_err < 1.0e-12,
    "[%s @ (%g,%g,%g)] M^T γ M ≠ I: max err = %.3e",
    label, x, y, z, max_orth_err );

  // Round trip on conserved state. The forward+back identity holds at
  // the WAVE level (w_tet round-tripped through back-transform gives
  // back the input w_tet, modulo √γ). For the STATE itself, the
  // forward map undensitizes (q_tet[i] = q_GR[i]/√γ for scalars; with
  // contravariant transform for momentum). Back-transform on a "wave"
  // with all 5 slots equal to q_tet then re-densitizes.
  double max_rt_err = 0.0;
  for (size_t ir = 0; ir < sizeof(rhos)/sizeof(*rhos); ir++) {
    for (size_t ip = 0; ip < sizeof(pres)/sizeof(*pres); ip++) {
      for (size_t iv = 0; iv < sizeof(vels)/sizeof(*vels); iv++) {
        double rho = rhos[ir], p = pres[ip];
        const double *v = vels[iv];

        // Skip superluminal.
        double v_lo[3] = {0,0,0};
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            v_lo[i] += g_ij[i][j] * v[j];
        double vsq = 0.0;
        for (int i = 0; i < 3; i++) vsq += v_lo[i] * v[i];
        if (!(vsq < 1.0 - 1.0e-6)) continue;

        double q[5], q_tet[5], q_back[5];
        build_state_convA(gas_gamma, rho, v, p, prods, q);
        gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(q, sqrt_det, inv_g, M_inv, q_tet);
        // The "wave_to_curved" maps tetrad-frame jumps to curved-frame
        // jumps. We use it on q_tet directly (which represents "the
        // full wave from 0 to q_tet"): the result should be the
        // densitized curved state.
        gkyl_gr_euler_tetrad_mod_wave_to_curved_contra(q_tet, sqrt_det, M_inv, q_back);
        for (int i = 0; i < 5; i++) {
          double e = fabs(q_back[i] - q[i]);
          if (e > max_rt_err) max_rt_err = e;
        }
      }
    }
  }
  TEST_CHECK_( max_rt_err < 1.0e-11,
    "[%s @ (%g,%g,%g)] state round-trip max err = %.3e",
    label, x, y, z, max_rt_err );
}

void test_round_trip_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_round_trip(st, "Mink", 0.3, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_round_trip_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_round_trip(st, "Schw", 0.3, 0.2, 0.0);
  run_round_trip(st, "Schw", 0.5, 0.0, 0.0);
  run_round_trip(st, "Schw", 0.4, 0.4, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_round_trip_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_round_trip(st, "Kerr", 0.3, 0.2, 0.0);
  run_round_trip(st, "Kerr", 0.5, 0.0, 0.0);
  run_round_trip(st, "Kerr", 0.4, 0.4, 0.0);
  gkyl_gr_spacetime_release(st);
}

// (b) Primitive-recovery consistency between curved-frame prim_vars and
//     SR Newton on q_tet. If both Newton solves are correct, they
//     should agree on (ρ, p, W) to machine precision. v_curved^i and
//     v_tet^a are related by the tetrad transform; v_tet^0 should be
//     v^x / √γ^{xx} by our construction.
static void
run_prim_consistency(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z)
{
  double gas_gamma = 5.0 / 3.0;
  double prods[GKYL_GR_SP_NCOMP_BASE];
  fill_prods_at(spacetime, x, y, z, prods);
  if (prods[GKYL_GR_SP_EXCISION] < 0.0) return;

  double g_ij[3][3], inv_g[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_ij[i][j]  = prods[GKYL_GR_SP_GIJ + 3*i + j];
      inv_g[i][j] = prods[GKYL_GR_SP_INV_GIJ + 3*i + j];
    }
  double M[3][3], M_inv[3][3];
  gkyl_gr_euler_tetrad_mod_build_triad_contravariant_x(g_ij, inv_g, M, M_inv);
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);

  double rho_in = 1.0, p_in = 1.5;
  double v_in[3] = { 0.12, 0.07, 0.03 };

  // Build state.
  double q[5];
  build_state_convA(gas_gamma, rho_in, v_in, p_in, prods, q);

  // Curved-frame recovery via shared helper.
  double D    = q[0] / sqrt_det;
  double momx = q[1] / sqrt_det;
  double momy = q[2] / sqrt_det;
  double momz = q[3] / sqrt_det;
  double tau  = q[4] / sqrt_det;
  struct gkyl_gr_euler_prim prim_curved;
  gkyl_gr_euler_recover_primitives(gas_gamma, D, momx, momy, momz, tau,
    inv_g, &prim_curved);

  // Tetrad-frame recovery: transform to tetrad, then run flat-space
  // Newton (replicate the body of sr_hll_minkowski).
  double q_tet[5];
  gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(q, sqrt_det, inv_g, M_inv, q_tet);

  double D_t = q_tet[0], Sx_t = q_tet[1], Sy_t = q_tet[2], Sz_t = q_tet[3], tau_t = q_tet[4];
  // Flat-space Banyuls Newton (inline).
  double s_sq_t = ((tau_t + D_t)*(tau_t + D_t))
                - (Sx_t*Sx_t + Sy_t*Sy_t + Sz_t*Sz_t);
  double C, C0;
  if (s_sq_t < pow(10.0, -10.0)) {
    C  = D_t / sqrt(pow(10.0, -10.0));
    C0 = (D_t + tau_t) / sqrt(pow(10.0, -10.0));
  } else {
    C  = D_t / sqrt(s_sq_t);
    C0 = (D_t + tau_t) / sqrt(s_sq_t);
  }
  double a0 = -1.0 / (gas_gamma*gas_gamma);
  double a1 = -2.0 * C * ((gas_gamma - 1.0)/(gas_gamma*gas_gamma));
  double a2 = ((gas_gamma - 2.0)/gas_gamma) * ((C0*C0) - 1.0) + 1.0
            - (C*C)*((gas_gamma-1.0)/gas_gamma)*((gas_gamma-1.0)/gas_gamma);
  double a4 = (C0*C0) - 1.0;
  double eta = 2.0 * C * ((gas_gamma - 1.0)/gas_gamma);
  double g = 1.0;
  for (int it = 0; it < 100; it++) {
    double poly = (a4*g*g*g)*(g-eta) + (a2*g*g) + (a1*g) + a0;
    double dpoly = a1 + 2.0*a2*g + 4.0*a4*g*g*g - 3.0*eta*a4*g*g;
    double gn = g - poly/dpoly;
    if (fabs(g - gn) < pow(10.0, -14.0)) { g = gn; break; }
    g = gn;
  }
  double W_t = 0.5 * C0 * g * (1.0 + sqrt(1.0
    + 4.0*((gas_gamma-1.0)/gas_gamma)*((1.0 - C*g)/(C0*C0*g*g))));
  double h_t = 1.0 / (C * g);
  double rho_t = D_t / W_t;
  double p_t  = (rho_t * h_t * W_t * W_t) - D_t - tau_t;

  // (ρ, p, W, h) should match between curved-frame and tetrad-frame
  // recovery to round-off.
  TEST_CHECK_( fabs(prim_curved.rho - rho_t) < 1.0e-12,
    "[%s] ρ mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.rho, rho_t, fabs(prim_curved.rho - rho_t) );
  TEST_CHECK_( fabs(prim_curved.p - p_t) < 1.0e-12,
    "[%s] p mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.p, p_t, fabs(prim_curved.p - p_t) );
  TEST_CHECK_( fabs(prim_curved.W - W_t) < 1.0e-12,
    "[%s] W mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.W, W_t, fabs(prim_curved.W - W_t) );
  TEST_CHECK_( fabs(prim_curved.h - h_t) < 1.0e-12,
    "[%s] h mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.h, h_t, fabs(prim_curved.h - h_t) );

  // The tetrad-frame normal velocity v_tet^x should equal v^x/√γ^{xx}.
  double vx_tet_expected = v_in[0] / sqrt(inv_g[0][0]);
  double vx_tet_computed = Sx_t / (rho_t * h_t * W_t * W_t);
  TEST_CHECK_( fabs(vx_tet_expected - vx_tet_computed) < 1.0e-12,
    "[%s] v_tet^x: expected v^x/√γ^{xx} = %g, got %g (diff %.3e)",
    label, vx_tet_expected, vx_tet_computed,
    fabs(vx_tet_expected - vx_tet_computed) );
}

void test_prim_consistency_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_prim_consistency(st, "Mink", 0.3, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_prim_consistency_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_prim_consistency(st, "Schw @ (0.3,0.2,0)", 0.3, 0.2, 0.0);
  run_prim_consistency(st, "Schw @ (0.5,0.0,0)", 0.5, 0.0, 0.0);
  run_prim_consistency(st, "Schw @ (0.4,0.4,0)", 0.4, 0.4, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_prim_consistency_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_prim_consistency(st, "Kerr @ (0.3,0.2,0)", 0.3, 0.2, 0.0);
  run_prim_consistency(st, "Kerr @ (0.5,0.0,0)", 0.5, 0.0, 0.0);
  run_prim_consistency(st, "Kerr @ (0.4,0.4,0)", 0.4, 0.4, 0.0);
  gkyl_gr_spacetime_release(st);
}

// (c) BHL-regime Riemann states: small τ, transverse-dominant velocity,
//     near-vacuum on one side, large γ off-diagonal. Tests whether the
//     standard Riemann properties (wave sum, flux jump, fluct balance,
//     etc.) hold on states the production solver actually sees.
static void
run_bhl_regime_states(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z,
  enum gkyl_wv_gr_euler_tetrad_rp rp, int num_waves)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, x, y, z, prods_row);
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };
  int idx[1] = { 0 };

  // BHL-like state pairs: small τ on one side, transverse v, etc.
  struct {
    const char *name;
    double rL, pL; double vL[3];
    double rR, pR; double vR[3];
  } cases[] = {
    // Small τ (near-floor) on left, typical on right.
    { "near-floor-τ-L",  0.01, 1.0e-10, { 0.05, 0.10, 0.0 },
                          1.0,  0.5,    { 0.10, 0.20, 0.0 } },
    // Transverse-dominant velocity on both sides.
    { "transverse-dom",   1.0,  0.5,    { 0.05, 0.40, 0.10 },
                          0.5,  0.3,    { 0.05, 0.45, 0.15 } },
    // Asymmetric strong shock with non-trivial transverse v.
    { "asym-strong",      3.0,  2.0,    { 0.20, 0.30, 0.05 },
                          0.3,  0.05,   { -0.10, 0.10, 0.0 } },
    // Near-equal states (small perturbation): linearization test.
    { "small-perturb",    1.0,  1.0,    { 0.10, 0.10, 0.10 },
                          1.001, 1.001, { 0.101, 0.100, 0.099 } },
  };

  for (size_t c = 0; c < sizeof(cases)/sizeof(*cases); c++) {
    // Build states.
    double v_co_l[3] = { cases[c].vL[0], cases[c].vL[1], cases[c].vL[2] };
    double v_co_r[3] = { cases[c].vR[0], cases[c].vR[1], cases[c].vR[2] };
    double qL_glob[5], qR_glob[5];
    build_state_convA(gas_gamma, cases[c].rL, v_co_l, cases[c].pL, prods_row, qL_glob);
    build_state_convA(gas_gamma, cases[c].rR, v_co_r, cases[c].pR, prods_row, qR_glob);

    eqn->set_interface_idx_func(eqn, idx, idx);
    double qL[5], qR[5];
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL_glob, qL);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR_glob, qR);

    double delta[5];
    for (int i = 0; i < 5; i++) delta[i] = qR[i] - qL[i];

    double waves[3 * 5], speeds[3];
    double maxs = eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      delta, qL, qR, 1.0, 1.0, waves, speeds);
    (void)maxs;

    // Wave sum.
    for (int i = 0; i < 5; i++) {
      double sum = 0.0;
      for (int k = 0; k < num_waves; k++) sum += waves[k * 5 + i];
      TEST_CHECK_( fabs(sum - delta[i]) < 1.0e-11,
        "[%s/%s] wave sum: comp %d, |Σw − Δq| = %.3e",
        label, cases[c].name, i, fabs(sum - delta[i]) );
    }

    // Flux jump.
    double dF[5];
    banyuls_delta_flux(gas_gamma, grm, qL, qR, dF);
    for (int i = 0; i < 5; i++) {
      double sw = 0.0;
      for (int k = 0; k < num_waves; k++) sw += speeds[k] * waves[k * 5 + i];
      // Near-floor τ cases incur catastrophic-cancellation precision
      // loss in the Banyuls Newton (Δ at 1e-10 squared minus terms at
      // 1e-8), so we loosen the tolerance to 1e-7 for those.
      double tol = strstr(cases[c].name, "near-floor") ? 1.0e-7 : 1.0e-9;
      TEST_CHECK_( fabs(sw - dF[i]) < tol,
        "[%s/%s] flux jump: comp %d, |Σs·w − ΔF| = %.3e",
        label, cases[c].name, i, fabs(sw - dF[i]) );
    }

    // Fluctuation balance.
    double amdq[5], apdq[5];
    eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      qL, qR, 1.0, 1.0, waves, speeds, amdq, apdq);
    for (int i = 0; i < 5; i++) {
      double tol_fluct = strstr(cases[c].name, "near-floor") ? 1.0e-7 : 1.0e-9;
      TEST_CHECK_( fabs(amdq[i] + apdq[i] - dF[i]) < tol_fluct,
        "[%s/%s] fluct: comp %d, |amdq+apdq − ΔF| = %.3e",
        label, cases[c].name, i, fabs(amdq[i] + apdq[i] - dF[i]) );
    }
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_bhl_regime_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_bhl_regime_states(st, "Schw-Lax", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  run_bhl_regime_states(st, "Schw-Lax", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(st);
}
void test_bhl_regime_hll_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_bhl_regime_states(st, "Schw-HLL", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  run_bhl_regime_states(st, "Schw-HLL", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(st);
}
void test_bhl_regime_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_bhl_regime_states(st, "Kerr-Lax", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  run_bhl_regime_states(st, "Kerr-Lax", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(st);
}
void test_bhl_regime_hll_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_bhl_regime_states(st, "Kerr-HLL", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  run_bhl_regime_states(st, "Kerr-HLL", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(st);
}

// (d) Precision sweep: how does flux-jump residual scale with τ on the
//     low side? Report-only — diagnostic of catastrophic cancellation in
//     the Banyuls Newton when τ is near the admissibility boundary. Tells
//     us what the "safe" repair floor is from a finite-precision
//     standpoint — below this value, fluxes degrade and feed back into
//     larger downstream errors.
static void
run_floor_precision_sweep(struct gkyl_gr_spacetime *spacetime,
  enum gkyl_wv_gr_euler_tetrad_rp rp, const char *rp_name)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, 0.3, 0.2, 0.0, prods_row);
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };
  int idx[1] = { 0 };

  // Sweep pL (which sets τ_L in the cold limit). Other side is a typical
  // state. Same per-side primitives otherwise.
  double pL_values[] = {
    1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10,
    1.0e-12, 1.0e-14, 1.0e-16
  };
  double rho_L = 0.01;
  double v_L[3] = { 0.05, 0.10, 0.0 };
  double rho_R = 1.0, p_R = 0.5;
  double v_R[3] = { 0.10, 0.20, 0.0 };

  fprintf(stderr,
    "[floor-precision %s] p_L      |Σs·w − ΔF|_max  |amdq+apdq − ΔF|_max\n",
    rp_name);
  for (size_t i = 0; i < sizeof(pL_values)/sizeof(*pL_values); i++) {
    double p_L = pL_values[i];
    double qL_glob[5], qR_glob[5];
    build_state_convA(gas_gamma, rho_L, v_L, p_L, prods_row, qL_glob);
    build_state_convA(gas_gamma, rho_R, v_R, p_R, prods_row, qR_glob);

    eqn->set_interface_idx_func(eqn, idx, idx);
    double qL[5], qR[5];
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL_glob, qL);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR_glob, qR);

    double delta[5];
    for (int k = 0; k < 5; k++) delta[k] = qR[k] - qL[k];
    double waves[3 * 5], speeds[3];
    eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      delta, qL, qR, 1.0, 1.0, waves, speeds);

    double dF[5];
    banyuls_delta_flux(gas_gamma, grm, qL, qR, dF);

    double max_fj = 0.0, max_fb = 0.0;
    for (int k = 0; k < 5; k++) {
      double sw = 0.0;
      for (int w = 0; w < 2; w++) sw += speeds[w] * waves[w * 5 + k];
      double r = fabs(sw - dF[k]);
      if (r > max_fj) max_fj = r;
    }
    double amdq[5], apdq[5];
    eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      qL, qR, 1.0, 1.0, waves, speeds, amdq, apdq);
    for (int k = 0; k < 5; k++) {
      double r = fabs(amdq[k] + apdq[k] - dF[k]);
      if (r > max_fb) max_fb = r;
    }

    fprintf(stderr, "  %-8s   %.0e  %.3e         %.3e\n",
      rp_name, p_L, max_fj, max_fb);
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_floor_precision_lax(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_floor_precision_sweep(st, WV_GR_EULER_TETRAD_RP_LAX, "Schw-Lax");
  gkyl_gr_spacetime_release(st);
}
void test_floor_precision_hll(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_floor_precision_sweep(st, WV_GR_EULER_TETRAD_RP_HLL, "Schw-HLL");
  gkyl_gr_spacetime_release(st);
}

// (e) τ/D << 1 positivity sweep — documents the regime where Lax (and
// HLL) provably do NOT preserve τ > 0 on SR/GR Euler. The Banyuls
// τ-flux F[τ] = (τ + p)·v^x can be much larger than τ when p ≫ τ,
// breaking the convex-bound proof for the τ component. This regime
// is exactly what the BHL bow shock exercises: τ/D ~ 1–5% with
// p/τ ~ 10–15. Mignone & Bodo (2005) document this as motivation
// for HLLC.
//
// This is a REPORT-ONLY test (TEST_MSG, not TEST_CHECK). The goal is to
// quantify the τ-violation as a function of τ/D, not to assert that
// Lax/HLL preserve τ in this regime (they don't).
static void
run_small_tau_over_D_sweep(struct gkyl_gr_spacetime *spacetime,
  enum gkyl_wv_gr_euler_tetrad_rp rp, const char *label)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, 0.0, 0.0, 0.0, prods_row);  // Minkowski only
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  double sqrt_det = sqrt(prods_row[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g[3][3] = {
    { prods_row[GKYL_GR_SP_INV_GIJ + 0], prods_row[GKYL_GR_SP_INV_GIJ + 1], prods_row[GKYL_GR_SP_INV_GIJ + 2] },
    { prods_row[GKYL_GR_SP_INV_GIJ + 3], prods_row[GKYL_GR_SP_INV_GIJ + 4], prods_row[GKYL_GR_SP_INV_GIJ + 5] },
    { prods_row[GKYL_GR_SP_INV_GIJ + 6], prods_row[GKYL_GR_SP_INV_GIJ + 7], prods_row[GKYL_GR_SP_INV_GIJ + 8] },
  };

  // Sweep p/ρ ratio. For cold (v=0) the relation is τ ≈ p/(γ−1) and
  // D = ρ, so τ/D = (p/ρ)/(γ−1). For non-zero v there's a (W²−1) term
  // that boosts τ — we want to be in the regime where p ≫ τ which
  // happens when v is not too small but p/ρ stays bounded.
  //
  // To mimic the BHL bow-shock cells exactly (D=2.98, τ=0.10 in the
  // dumped failure case), use ρ=1, p=1.5, v=(0.5, 0, 0). That gives
  // W=1.155 and τ ≈ 0.6 — still not the τ/D=0.04 regime. To reach that
  // regime we need stronger velocity and smaller p.
  struct {
    const char *name;
    double rho_L, p_L; double v_L[3];
    double rho_R, p_R; double v_R[3];
  } cases[] = {
    // Symmetric setup: same ρ both sides, different p, "head-on" v.
    // Sweep over decreasing p (decreasing τ/D).
    { "p=2.0  (τ/D~3)",    1.0, 2.0,    { 0.30, 0.0, 0.0 },  1.0, 2.0,   { -0.30, 0.0, 0.0 } },
    { "p=0.5  (τ/D~1)",    1.0, 0.5,    { 0.30, 0.0, 0.0 },  1.0, 0.5,   { -0.30, 0.0, 0.0 } },
    { "p=0.1  (τ/D~0.2)",  1.0, 0.1,    { 0.30, 0.0, 0.0 },  1.0, 0.1,   { -0.30, 0.0, 0.0 } },
    { "p=0.01 (τ/D~0.02)", 1.0, 0.01,   { 0.30, 0.0, 0.0 },  1.0, 0.01,  { -0.30, 0.0, 0.0 } },
    // BHL-bow-shock proxy: ρ moderate, p ≫ τ, strong relativistic v.
    { "BHL-proxy",         1.0, 1.0,    { 0.60, 0.10, 0.0 }, 1.0, 1.0,   { -0.60, 0.10, 0.0 } },
    // Asymmetric (one side hot, other cold) — the configuration that
    // routinely produces the bow-shock τ-flip in BHL.
    { "asym-hot-cold",     1.0, 2.0,    { 0.10, 0.0, 0.0 },  0.1, 0.001, { -0.10, 0.0, 0.0 } },
  };

  fprintf(stderr,
    "[small-τ %s]  case                    τ_L_in/D_L_in   τ_L_new   τ_R_in/D_R_in   τ_R_new\n",
    label);
  for (size_t c = 0; c < sizeof(cases)/sizeof(*cases); c++) {
    double qL[5], qR[5];
    build_state_convA(gas_gamma, cases[c].rho_L, cases[c].v_L, cases[c].p_L,
      prods_row, qL);
    build_state_convA(gas_gamma, cases[c].rho_R, cases[c].v_R, cases[c].p_R,
      prods_row, qR);

    // Sanity: inputs admissible.
    double qL_und[5], qR_und[5];
    for (int i = 0; i < 5; i++) { qL_und[i] = qL[i]/sqrt_det; qR_und[i] = qR[i]/sqrt_det; }
    bool d_ok, s_ok, t_ok;
    record_admissibility(inv_g, qL_und, &d_ok, &s_ok, &t_ok);
    if (!(d_ok && s_ok && t_ok)) { fprintf(stderr, "  [%s] L inadmissible — skip\n", cases[c].name); continue; }
    record_admissibility(inv_g, qR_und, &d_ok, &s_ok, &t_ok);
    if (!(d_ok && s_ok && t_ok)) { fprintf(stderr, "  [%s] R inadmissible — skip\n", cases[c].name); continue; }

    double tau_over_D_L = qL_und[4] / qL_und[0];
    double tau_over_D_R = qR_und[4] / qR_und[0];

    // Wave step.
    int idx[1] = { 0 };
    eqn->set_interface_idx_func(eqn, idx, idx);
    double norm[3] = { 1.0, 0.0, 0.0 };
    double tau1v[3] = { 0.0, 1.0, 0.0 };
    double tau2v[3] = { 0.0, 0.0, 1.0 };
    double qL_loc[5], qR_loc[5];
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qL, qL_loc);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qR, qR_loc);
    double delta[5];
    for (int i = 0; i < 5; i++) delta[i] = qR_loc[i] - qL_loc[i];
    double waves[3 * 5], speeds[3];
    double maxs = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      delta, qL_loc, qR_loc, 1.0, 1.0, waves, speeds);
    double amdq[5], apdq[5];
    eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      qL_loc, qR_loc, 1.0, 1.0, waves, speeds, amdq, apdq);

    if (!(maxs > 0.0)) { fprintf(stderr, "  [%s] degenerate — skip\n", cases[c].name); continue; }
    double dt_dx = 0.5 / maxs;  // 2D directional-splitting CFL

    // Undensitized new τ for each cell.
    double qL_new[5], qR_new[5];
    for (int i = 0; i < 5; i++) {
      qL_new[i] = qL_loc[i] - dt_dx * amdq[i];
      qR_new[i] = qR_loc[i] - dt_dx * apdq[i];
    }
    double tau_L_new = qL_new[4] / sqrt_det;
    double tau_R_new = qR_new[4] / sqrt_det;

    fprintf(stderr,
      "  [%s/%-20s] %.4f          %+.3e   %.4f          %+.3e\n",
      label, cases[c].name, tau_over_D_L, tau_L_new, tau_over_D_R, tau_R_new);
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_small_tau_over_D_lax(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_small_tau_over_D_sweep(st, WV_GR_EULER_TETRAD_RP_LAX, "Mink-Lax");
  gkyl_gr_spacetime_release(st);
}
void test_small_tau_over_D_hll(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_small_tau_over_D_sweep(st, WV_GR_EULER_TETRAD_RP_HLL, "Mink-HLL");
  gkyl_gr_spacetime_release(st);
}

// Direct-conservative variant: bypass `build_state_convA`'s primitive →
// conserved map and set (D, S_i, τ) directly. In Minkowski (flat), the
// constraint s² > 0 means |S| < D+τ, so for small τ/D we'd need
// vanishing |S|, contradicting the BHL bow-shock regime where Sx ~ D.
// The Schwarzschild γ_xx > 1 amplifies S relative to D+τ, so curved
// spacetime is the natural setting for τ/D << 1 with large Sx.
static void
run_direct_state_sweep(struct gkyl_gr_spacetime *spacetime,
  const char *spacetime_label,
  enum gkyl_wv_gr_euler_tetrad_rp rp, const char *rp_label)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  // Use (0.3, 0.2, 0) which is the same Schw sample point we've used
  // throughout the convA tests.
  fill_prods_at(spacetime, 0.3, 0.2, 0.0, prods_row);
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  double sqrt_det = sqrt(prods_row[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g[3][3] = {
    { prods_row[GKYL_GR_SP_INV_GIJ + 0], prods_row[GKYL_GR_SP_INV_GIJ + 1], prods_row[GKYL_GR_SP_INV_GIJ + 2] },
    { prods_row[GKYL_GR_SP_INV_GIJ + 3], prods_row[GKYL_GR_SP_INV_GIJ + 4], prods_row[GKYL_GR_SP_INV_GIJ + 5] },
    { prods_row[GKYL_GR_SP_INV_GIJ + 6], prods_row[GKYL_GR_SP_INV_GIJ + 7], prods_row[GKYL_GR_SP_INV_GIJ + 8] },
  };

  // Three states modeled directly after the dumped BHL repair fires
  // (see SESSION_NOTES_2.md §17). qin_C is the cell that fails;
  // qin_L, qin_R are its neighbors. Values are DENSITIZED (in code
  // form), as they come out of qin in wave_prop.c.
  struct {
    const char *name;
    double D_C, Sx_C, Sy_C, Sz_C, tau_C;
    double D_L, Sx_L, Sy_L, Sz_L, tau_L;
    double D_R, Sx_R, Sy_R, Sz_R, tau_R;
  } cases[] = {
    // BHL repair #0 dump verbatim.
    { "bhl-repair-#0",
      2.978, 2.499, 0.118, 0.000, 0.122,
      2.541, 2.146, 0.045, 0.000, 0.136,
      3.395, 2.827, 0.222, 0.000, 0.095 },
    // BHL repair #4 dump.
    { "bhl-repair-#4",
      3.407, 2.633, 0.209, 0.000, 0.021,
      3.123, 2.165, 0.119, 0.000, 0.025,
      3.009, 3.319, -0.766, 0.000, 0.526 },
    // BHL repair #19 dump — has very large negative Sx_R.
    { "bhl-repair-#19",
      3.505, 2.969, -1.664, 0.000, 0.193,
      3.872, 2.574, 0.285, 0.000, 0.025,
      2.251, 2.424, -0.887, 0.000, 0.392 },
  };

  fprintf(stderr,
    "[direct-state %s/%s]  case            τ_C_in/D_C_in    τ_C_new (sign+mag)\n",
    spacetime_label, rp_label);
  for (size_t c = 0; c < sizeof(cases)/sizeof(*cases); c++) {
    double qL_glob[5] = { cases[c].D_L, cases[c].Sx_L, cases[c].Sy_L, cases[c].Sz_L, cases[c].tau_L };
    double qC_glob[5] = { cases[c].D_C, cases[c].Sx_C, cases[c].Sy_C, cases[c].Sz_C, cases[c].tau_C };
    double qR_glob[5] = { cases[c].D_R, cases[c].Sx_R, cases[c].Sy_R, cases[c].Sz_R, cases[c].tau_R };

    // Check input admissibility on all three cells.
    bool all_ok = true;
    for (int side = 0; side < 3; side++) {
      const double *q = (side == 0) ? qL_glob : (side == 1) ? qC_glob : qR_glob;
      double qu[5];
      for (int i = 0; i < 5; i++) qu[i] = q[i] / sqrt_det;
      bool d_ok, s_ok, t_ok;
      record_admissibility(inv_g, qu, &d_ok, &s_ok, &t_ok);
      if (!(d_ok && s_ok && t_ok)) {
        fprintf(stderr, "  [%s] WARNING: cell %d inadmissible — skip\n", cases[c].name, side);
        all_ok = false; break;
      }
    }
    if (!all_ok) continue;

    double D_C_und  = qC_glob[0] / sqrt_det;
    double tau_C_und = qC_glob[4] / sqrt_det;
    double tau_over_D = tau_C_und / D_C_und;

    // Compute fluctuations at both interfaces (L-C and C-R), then
    // joint update on cell C.
    int idx[1] = { 0 };
    double norm[3] = { 1.0, 0.0, 0.0 };
    double tau1v[3] = { 0.0, 1.0, 0.0 };
    double tau2v[3] = { 0.0, 0.0, 1.0 };

    eqn->set_interface_idx_func(eqn, idx, idx);
    double qL_loc[5], qC_loc_LC[5], qC_loc_CR[5], qR_loc[5];
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qL_glob, qL_loc);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qC_glob, qC_loc_LC);
    double delta_LC[5];
    for (int i = 0; i < 5; i++) delta_LC[i] = qC_loc_LC[i] - qL_loc[i];
    double waves_LC[3 * 5], speeds_LC[3];
    double maxs_LC = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      delta_LC, qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC);
    double amdq_LC[5], apdq_LC[5];
    eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC, amdq_LC, apdq_LC);

    eqn->set_interface_idx_func(eqn, idx, idx);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qC_glob, qC_loc_CR);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qR_glob, qR_loc);
    double delta_CR[5];
    for (int i = 0; i < 5; i++) delta_CR[i] = qR_loc[i] - qC_loc_CR[i];
    double waves_CR[3 * 5], speeds_CR[3];
    double maxs_CR = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      delta_CR, qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR);
    double amdq_CR[5], apdq_CR[5];
    eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
      qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR, amdq_CR, apdq_CR);

    double maxs = fmax(maxs_LC, maxs_CR);
    if (!(maxs > 0.0)) { fprintf(stderr, "  [%s] degenerate\n", cases[c].name); continue; }
    double dt_dx = 0.5 / maxs;  // 2D directional-splitting CFL

    double qC_new[5];
    for (int i = 0; i < 5; i++) {
      qC_new[i] = qC_loc_CR[i] - dt_dx * (amdq_CR[i] + apdq_LC[i]);
    }
    double tau_C_new = qC_new[4] / sqrt_det;

    fprintf(stderr,
      "  [%s/%-16s] τ_C/D_C = %.4f    τ_C_new = %+.3e   %s\n",
      rp_label, cases[c].name, tau_over_D, tau_C_new,
      tau_C_new < 0.0 ? "** τ<0 VIOLATION **" : "");

    // Report-only TEST_MSG so the test PASSES while documenting the
    // pathology.
    TEST_MSG( "[%s/%s] τ/D = %.4f, τ_new = %+.3e",
      rp_label, cases[c].name, tau_over_D, tau_C_new );
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_direct_state_lax(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_LAX, "Lax");
  gkyl_gr_spacetime_release(st);
}
void test_direct_state_hll(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_HLL, "HLL");
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 6. Three-cell positivity test (mimics wave_prop's actual cell update)
// ---------------------------------------------------------------------------
// Real wave_prop does:
//   qM_new = qM - dt_dx * (amdq_R_interface + apdq_L_interface)
// where amdq_R is the LEFT-going fluctuation from M's right interface (M-R)
// and apdq_L is the RIGHT-going fluctuation from M's left interface (L-M).
// So the middle cell receives flux contributions from BOTH neighbors
// simultaneously. The 2-cell positivity test only exercised ONE interface
// per cell.
//
// This test sets up 3 cells (L, M, R) with M at the τ-floor and both
// neighbors on either side, then applies the joint flux update to M.
// The "outflow on both sides" pattern (v_L < 0, v_R > 0) is the
// pathological case Lax can fail on at the boundary of CFL.

static void
run_three_cell_positivity(struct gkyl_wv_eqn *eqn,
  struct wv_gr_euler_tetrad_mod *grm,
  const double *prods,
  double rho_M, double tau_M, double v_M[3],
  double rho_L, double tau_L, double v_L[3],
  double rho_R, double tau_R, double v_R[3],
  const char *case_name,
  int *D_v, int *S2_v, int *tau_v, int *total)
{
  double gas_gamma = 5.0 / 3.0;
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g[3][3] = {
    { prods[GKYL_GR_SP_INV_GIJ + 0], prods[GKYL_GR_SP_INV_GIJ + 1], prods[GKYL_GR_SP_INV_GIJ + 2] },
    { prods[GKYL_GR_SP_INV_GIJ + 3], prods[GKYL_GR_SP_INV_GIJ + 4], prods[GKYL_GR_SP_INV_GIJ + 5] },
    { prods[GKYL_GR_SP_INV_GIJ + 6], prods[GKYL_GR_SP_INV_GIJ + 7], prods[GKYL_GR_SP_INV_GIJ + 8] },
  };

  // Build the three states from primitives. Convert τ to p via the cold
  // limit (p = (γ−1)·τ). This isn't exact for hot flow but is the right
  // order-of-magnitude conversion for floor-level states.
  double p_M = (gas_gamma - 1.0) * tau_M;
  double p_L = (gas_gamma - 1.0) * tau_L;
  double p_R = (gas_gamma - 1.0) * tau_R;
  double qL[5], qM[5], qR[5];
  build_state_convA(gas_gamma, rho_L, v_L, p_L, prods, qL);
  build_state_convA(gas_gamma, rho_M, v_M, p_M, prods, qM);
  build_state_convA(gas_gamma, rho_R, v_R, p_R, prods, qR);

  // Verify all three input states admissible.
  bool d_ok, s_ok, t_ok;
  for (int side = 0; side < 3; side++) {
    const double *q = (side == 0) ? qL : (side == 1) ? qM : qR;
    double qu[5];
    for (int i = 0; i < 5; i++) qu[i] = q[i] / sqrt_det;
    record_admissibility(inv_g, qu, &d_ok, &s_ok, &t_ok);
    if (!(d_ok && s_ok && t_ok)) {
      fprintf(stderr, "  [3cell %s] WARNING: input cell %d inadmissible "
        "(D=%d, S²=%d, τ=%d) — skipping\n",
        case_name, side, d_ok, s_ok, t_ok);
      return;
    }
  }

  // Compute fluctuations at each interface via the equation.
  int idx[1] = { 0 };
  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  // L-M interface: ql = qL, qr = qM.
  eqn->set_interface_idx_func(eqn, idx, idx);
  double qL_loc[5], qM_loc[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL, qL_loc);
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qM, qM_loc);
  double delta_LM[5];
  for (int i = 0; i < 5; i++) delta_LM[i] = qM_loc[i] - qL_loc[i];
  double waves_LM[3 * 5], speeds_LM[3];
  double maxs_LM = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    delta_LM, qL_loc, qM_loc, 1.0, 1.0, waves_LM, speeds_LM);
  double amdq_LM[5], apdq_LM[5];
  eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    qL_loc, qM_loc, 1.0, 1.0, waves_LM, speeds_LM, amdq_LM, apdq_LM);

  // M-R interface: ql = qM, qr = qR.
  eqn->set_interface_idx_func(eqn, idx, idx);
  double qR_loc[5], qM_loc2[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qM, qM_loc2);
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR, qR_loc);
  double delta_MR[5];
  for (int i = 0; i < 5; i++) delta_MR[i] = qR_loc[i] - qM_loc2[i];
  double waves_MR[3 * 5], speeds_MR[3];
  double maxs_MR = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    delta_MR, qM_loc2, qR_loc, 1.0, 1.0, waves_MR, speeds_MR);
  double amdq_MR[5], apdq_MR[5];
  eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    qM_loc2, qR_loc, 1.0, 1.0, waves_MR, speeds_MR, amdq_MR, apdq_MR);

  // CFL = 0.5 with the maximum of the two interface speeds (same as
  // wave_prop's per-direction global CFL constraint).
  double maxs = fmax(maxs_LM, maxs_MR);
  if (!(maxs > 0.0)) return;
  double cfl = 0.5;
  double dt_dx = cfl / maxs;

  // Joint update: cell M receives apdq from L-M (right-going from L into M)
  // AND amdq from M-R (left-going from R into M).
  double qM_new[5];
  for (int i = 0; i < 5; i++) {
    qM_new[i] = qM_loc2[i] - dt_dx * (amdq_MR[i] + apdq_LM[i]);
  }

  // Check admissibility on the joint update.
  double qu[5];
  for (int i = 0; i < 5; i++) qu[i] = qM_new[i] / sqrt_det;
  bool D_ok, S2_ok, tau_ok;
  record_admissibility(inv_g, qu, &D_ok, &S2_ok, &tau_ok);
  if (!D_ok)   { (*D_v)++;   fprintf(stderr,
    "  [3cell %s] D = %.3e\n", case_name, qu[0]); }
  if (!S2_ok)  { (*S2_v)++;  fprintf(stderr,
    "  [3cell %s] s² = %.3e\n", case_name,
    (qu[0]+qu[4])*(qu[0]+qu[4]) - gkyl_gr_euler_mom_sq(inv_g, qu[1], qu[2], qu[3])); }
  if (!tau_ok) { (*tau_v)++; fprintf(stderr,
    "  [3cell %s] τ = %.3e (D = %.3e, was %.3e)\n",
    case_name, qu[4], qu[0], tau_M); }
  (*total)++;
}

// Curated 3-cell scenarios. The "outflow-both-sides" case is the
// pathological one: middle cell at floor with both neighbors flowing
// outward away from M.
struct three_cell_case {
  const char *name;
  double rho_M, tau_M; double v_M[3];
  double rho_L, tau_L; double v_L[3];
  double rho_R, tau_R; double v_R[3];
};

static const struct three_cell_case g_three_cell_cases[] = {
  // M at floor, both neighbors at floor, outflow on both sides — the
  // pathological "vacuum suction" pattern.
  { "outflow-both-floor",
      0.01, 1.0e-6, { 0.0,   0.0, 0.0 },
      0.01, 1.0e-6, { -0.5,  0.0, 0.0 },
      0.01, 1.0e-6, {  0.5,  0.0, 0.0 } },

  // M at floor, neighbors hot with outflow.
  { "outflow-hot-neighbors",
      0.01, 1.0e-6, { 0.0,  0.0, 0.0 },
      1.0,  1.0,    { -0.3, 0.0, 0.0 },
      1.0,  1.0,    {  0.3, 0.0, 0.0 } },

  // M at floor with INFLOW from both sides — should be safe.
  { "inflow-both-hot",
      0.01, 1.0e-6, { 0.0,  0.0, 0.0 },
      1.0,  1.0,    {  0.3, 0.0, 0.0 },
      1.0,  1.0,    { -0.3, 0.0, 0.0 } },

  // M at floor, asymmetric: neighbor on one side hot, other at floor.
  { "asym-hot-L-floor-R",
      0.01, 1.0e-6, { 0.0,  0.0, 0.0 },
      1.0,  1.0,    {  0.0, 0.0, 0.0 },
      0.01, 1.0e-6, {  0.5, 0.0, 0.0 } },

  // M at floor, transverse outflow (v_y != 0 on neighbors).
  { "transverse-outflow",
      0.01, 1.0e-6, { 0.0,  0.0,  0.0 },
      0.01, 1.0e-6, { -0.3, 0.4,  0.0 },
      0.01, 1.0e-6, {  0.3, 0.4,  0.0 } },
};

static void
run_three_cell_sweep(struct gkyl_gr_spacetime *spacetime,
  enum gkyl_wv_gr_euler_tetrad_rp rp,
  const char *label)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(gas_gamma, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, 0.0, 0.0, 0.0, prods_row);
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);

  int D_v = 0, S2_v = 0, tau_v = 0, total = 0;
  int n = sizeof(g_three_cell_cases)/sizeof(*g_three_cell_cases);
  for (int c = 0; c < n; c++) {
    const struct three_cell_case *tc = &g_three_cell_cases[c];
    double v_M[3] = { tc->v_M[0], tc->v_M[1], tc->v_M[2] };
    double v_L[3] = { tc->v_L[0], tc->v_L[1], tc->v_L[2] };
    double v_R[3] = { tc->v_R[0], tc->v_R[1], tc->v_R[2] };
    run_three_cell_positivity(eqn, grm, prods_row,
      tc->rho_M, tc->tau_M, v_M,
      tc->rho_L, tc->tau_L, v_L,
      tc->rho_R, tc->tau_R, v_R,
      tc->name, &D_v, &S2_v, &tau_v, &total);
  }

  fprintf(stderr,
    "[3cell %s] across %d updates: D≤0: %d, s²≤0: %d, τ<0: %d\n",
    label, total, D_v, S2_v, tau_v);

  TEST_CHECK_( D_v == 0,   "D ≤ 0 violated %d times in 3-cell test", D_v );
  TEST_CHECK_( S2_v == 0,  "s² ≤ 0 violated %d times in 3-cell test", S2_v );
  TEST_CHECK_( tau_v == 0, "τ < 0 violated %d times in 3-cell test", tau_v );

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_three_cell_lax(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_three_cell_sweep(st, WV_GR_EULER_TETRAD_RP_LAX, "Lax-Mink");
  gkyl_gr_spacetime_release(st);
}
void test_three_cell_hll(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_three_cell_sweep(st, WV_GR_EULER_TETRAD_RP_HLL, "HLL-Mink");
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------

TEST_LIST = {
  { "banyuls_flux_consistency_lax_minkowski",     test_banyuls_flux_consistency_lax_minkowski },
  { "banyuls_flux_consistency_lax_schwarzschild", test_banyuls_flux_consistency_lax_schwarzschild },
  { "banyuls_flux_consistency_lax_kerr",          test_banyuls_flux_consistency_lax_kerr },

  { "riemann_properties_lax_minkowski",     test_riemann_properties_lax_minkowski },
  { "riemann_properties_lax_schwarzschild", test_riemann_properties_lax_schwarzschild },
  { "riemann_properties_lax_kerr",          test_riemann_properties_lax_kerr },
  { "riemann_properties_hll_minkowski",     test_riemann_properties_hll_minkowski },
  { "riemann_properties_hll_schwarzschild", test_riemann_properties_hll_schwarzschild },
  { "riemann_properties_hll_kerr",          test_riemann_properties_hll_kerr },

  { "positivity_lax_minkowski",     test_positivity_lax_minkowski },
  { "positivity_lax_schwarzschild", test_positivity_lax_schwarzschild },
  { "positivity_lax_kerr",          test_positivity_lax_kerr },
  { "positivity_hll_minkowski",     test_positivity_hll_minkowski },
  { "positivity_hll_schwarzschild", test_positivity_hll_schwarzschild },
  { "positivity_hll_kerr",          test_positivity_hll_kerr },

  { "near_floor_lax", test_near_floor_lax },
  { "near_floor_hll", test_near_floor_hll },

  { "round_trip_minkowski",     test_round_trip_minkowski },
  { "round_trip_schwarzschild", test_round_trip_schwarzschild },
  { "round_trip_kerr",          test_round_trip_kerr },

  { "prim_consistency_minkowski",     test_prim_consistency_minkowski },
  { "prim_consistency_schwarzschild", test_prim_consistency_schwarzschild },
  { "prim_consistency_kerr",          test_prim_consistency_kerr },

  { "bhl_regime_lax_schwarzschild", test_bhl_regime_lax_schwarzschild },
  { "bhl_regime_hll_schwarzschild", test_bhl_regime_hll_schwarzschild },
  { "bhl_regime_lax_kerr",          test_bhl_regime_lax_kerr },
  { "bhl_regime_hll_kerr",          test_bhl_regime_hll_kerr },

  { "floor_precision_lax", test_floor_precision_lax },
  { "floor_precision_hll", test_floor_precision_hll },

  { "three_cell_lax", test_three_cell_lax },
  { "three_cell_hll", test_three_cell_hll },

  { "small_tau_over_D_lax", test_small_tau_over_D_lax },
  { "small_tau_over_D_hll", test_small_tau_over_D_hll },

  { "direct_state_lax", test_direct_state_lax },
  { "direct_state_hll", test_direct_state_hll },

  { NULL, NULL },
};
