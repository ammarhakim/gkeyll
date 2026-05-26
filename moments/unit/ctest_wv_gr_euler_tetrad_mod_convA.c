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

// EOS modes exercised by the parametrized runners below. The three
// production EOS configurations: IDEAL γ=5/3, APPROXIMATE_SYNGE with
// use_rcc=false (Mathews-Taub cubic closure), and APPROXIMATE_SYNGE
// with use_rcc=true (Ryu-Chattopadhyay Newton refinement on top of TM).
#define NUM_EOS_MODES 3
static const struct gkyl_gr_euler_eos eos_modes[NUM_EOS_MODES] = {
  { .type = GR_EULER_EOS_IDEAL,             .gas_gamma = 5.0/3.0, .use_rcc = false },
  { .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .gas_gamma = 0.0,     .use_rcc = false },
  { .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .gas_gamma = 0.0,     .use_rcc = true  },
};
static const char *
eos_label_for(struct gkyl_gr_euler_eos eos)
{
  if (eos.type == GR_EULER_EOS_IDEAL) return "IDEAL";
  return eos.use_rcc ? "RCC" : "TM";
}

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
// W is computed self-consistently from |v|² = γ_ij v^i v^j; h via the
// EOS-specific enthalpy formula so the seed conservative state matches
// what the code-side recovery would produce for the same primitives.
static void
build_state_convA(struct gkyl_gr_euler_eos eos,
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
  double h = gkyl_gr_euler_eos_enthalpy(eos, rho, p);
  double rhohW2 = rho * h * W * W;

  q[0] = sqrt_det * rho * W;
  q[1] = sqrt_det * rhohW2 * v_lo[0];
  q[2] = sqrt_det * rhohW2 * v_lo[1];
  q[3] = sqrt_det * rhohW2 * v_lo[2];
  q[4] = sqrt_det * (rhohW2 - p - rho * W);
}

// Analytic Banyuls flux F^x[U] for the Convention-A state described by
// (rho, v_co, p) at the prods row. v̂^x = v^x − β^x/α. Enthalpy h is
// computed via the EOS-specific formula so the analytic flux matches
// the code-path flux for the same EOS.
//
//   F^x[D]   = α√γ · ρ W · v̂^x
//   F^x[S_i] = α√γ · (γ_ij·ρhW²·v^j · v̂^x + p · δ_i^x)
//   F^x[τ]   = α√γ · (τ · v̂^x + p · v^x)
static void
analytic_banyuls_flux(struct gkyl_gr_euler_eos eos,
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
  double h = gkyl_gr_euler_eos_enthalpy(eos, rho, p);
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
make_eqn(struct gkyl_gr_euler_eos eos, struct gkyl_range conf_range,
  enum gkyl_wv_gr_euler_tetrad_rp rp)
{
  return gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp){
      .gas_gamma = eos.gas_gamma,  // legacy IDEAL shortcut (eos overrides)
      .eos = eos,
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
  struct gkyl_gr_euler_eos eos, enum gkyl_wv_gr_euler_tetrad_rp rp)
{

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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
        build_state_convA(eos, rho, v, p, prods_row, q);

        double f_sr[5], f_gr[5];
        gkyl_gr_euler_tetrad_mod_flux(eos, q, prods_row, NULL, f_sr);
        gkyl_gr_euler_tetrad_mod_flux_correction(eos, q, prods_row, NULL, f_sr, f_gr);

        double f_ref[5];
        analytic_banyuls_flux(eos, rho, v, p, prods_row, f_ref);

        for (int i = 0; i < 5; i++) {
          double d = fabs(f_gr[i] - f_ref[i]);
          if (d > max_diff) max_diff = d;
        }
      }
    }
  }

  TEST_CHECK_( max_diff < 1.0e-10,
    "[%s @ (%g,%g,%g) EOS=%s] Banyuls flux residual: max |F_code − F_analytic| = %.3e",
    label, x, y, z, eos_label_for(eos), max_diff );

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

// Each wrapper loops over [IDEAL, TM, RCC] so every curated state is
// exercised in all three production EOS configurations.
void test_banyuls_flux_consistency_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_banyuls_flux_consistency(st, "Mink-Lax", 0.3, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX);
  gkyl_gr_spacetime_release(st);
}
void test_banyuls_flux_consistency_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_banyuls_flux_consistency(st, "Schw-Lax", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX);
    run_banyuls_flux_consistency(st, "Schw-Lax", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX);
    run_banyuls_flux_consistency(st, "Schw-Lax", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX);
  }
  gkyl_gr_spacetime_release(st);
}
void test_banyuls_flux_consistency_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_banyuls_flux_consistency(st, "Kerr-Lax", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX);
    run_banyuls_flux_consistency(st, "Kerr-Lax", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX);
    run_banyuls_flux_consistency(st, "Kerr-Lax", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX);
  }
  gkyl_gr_spacetime_release(st);
}

// HLLC variants. Banyuls flux consistency depends only on
// flux+flux_correction, so the rp choice is just a factory parameter
// (the wave/qfluct callbacks aren't exercised). Still worth running
// against HLLC as a smoke test that the dispatch wires up cleanly.
void test_banyuls_flux_consistency_hllc_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_banyuls_flux_consistency(st, "Mink-HLLC", 0.3, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC);
  gkyl_gr_spacetime_release(st);
}
void test_banyuls_flux_consistency_hllc_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_banyuls_flux_consistency(st, "Schw-HLLC", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC);
    run_banyuls_flux_consistency(st, "Schw-HLLC", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC);
    run_banyuls_flux_consistency(st, "Schw-HLLC", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC);
  }
  gkyl_gr_spacetime_release(st);
}
void test_banyuls_flux_consistency_hllc_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_banyuls_flux_consistency(st, "Kerr-HLLC", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC);
    run_banyuls_flux_consistency(st, "Kerr-HLLC", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC);
    run_banyuls_flux_consistency(st, "Kerr-HLLC", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC);
  }
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 2. Standard Riemann-solver property tests for Lax / HLL
// ---------------------------------------------------------------------------

// Returns the curved-frame Banyuls flux ΔF for a given pair of states using
// the same (flux + flux_correction) machinery the wave construction uses.
static void
banyuls_delta_flux(struct gkyl_gr_euler_eos eos, struct wv_gr_euler_tetrad_mod *grm,
  const double qL[5], const double qR[5], double dF[5])
{
  double fL_sr[5], fR_sr[5], fL_gr[5], fR_gr[5];
  gkyl_gr_euler_tetrad_mod_flux(eos, qL, grm->prodl_local, NULL, fL_sr);
  gkyl_gr_euler_tetrad_mod_flux(eos, qR, grm->prodr_local, NULL, fR_sr);
  gkyl_gr_euler_tetrad_mod_flux_correction(eos, qL, grm->prodl_local, NULL, fL_sr, fL_gr);
  gkyl_gr_euler_tetrad_mod_flux_correction(eos, qR, grm->prodr_local, NULL, fR_sr, fR_gr);
  for (int i = 0; i < 5; i++) dF[i] = fR_gr[i] - fL_gr[i];
}

static void
run_riemann_properties_lax_hll(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z,
  struct gkyl_gr_euler_eos eos, enum gkyl_wv_gr_euler_tetrad_rp rp,
  int num_waves)
{
  double gas_gamma = eos.gas_gamma;

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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
  build_state_convA(eos, rho_L, v_L, p_L, prods_row, qL_glob);
  build_state_convA(eos, rho_R, v_R, p_R, prods_row, qR_glob);

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
  banyuls_delta_flux(eos, grm, qL, qR, dF);
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
  build_state_convA(eos, 1.0, (double[]){0.1, 0.2, 0.05}, 1.0, prods_row, qE);
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

// Each wrapper loops over [IDEAL, TM, RCC] for full EOS coverage.
void test_riemann_properties_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_riemann_properties_lax_hll(st, "Mink-Lax", 0.3, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_lax_hll(st, "Schw-Lax", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
    run_riemann_properties_lax_hll(st, "Schw-Lax", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
    run_riemann_properties_lax_hll(st, "Schw-Lax", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
  }
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_lax_hll(st, "Kerr-Lax", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
    run_riemann_properties_lax_hll(st, "Kerr-Lax", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
    run_riemann_properties_lax_hll(st, "Kerr-Lax", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
  }
  gkyl_gr_spacetime_release(st);
}

void test_riemann_properties_hll_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_riemann_properties_lax_hll(st, "Mink-HLL", 0.3, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_hll_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_lax_hll(st, "Schw-HLL", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
    run_riemann_properties_lax_hll(st, "Schw-HLL", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
    run_riemann_properties_lax_hll(st, "Schw-HLL", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
  }
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_hll_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_lax_hll(st, "Kerr-HLL", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
    run_riemann_properties_lax_hll(st, "Kerr-HLL", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
    run_riemann_properties_lax_hll(st, "Kerr-HLL", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
  }
  gkyl_gr_spacetime_release(st);
}

// HLLC variants. The runner is parameterized by num_waves so the same
// machinery works for HLLC (3 waves) as it does for Lax/HLL (2). The
// flux-jump check is the key one for finite-β verification: it confirms
// that the wave-based update reproduces the Banyuls flux exactly,
// including the moving-interface β·U term that lives in the speed shift
// (see SESSION_NOTES_BETA_FLUX_JUMP.md).
void test_riemann_properties_hllc_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_riemann_properties_lax_hll(st, "Mink-HLLC", 0.3, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_hllc_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_lax_hll(st, "Schw-HLLC", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
    run_riemann_properties_lax_hll(st, "Schw-HLLC", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
    run_riemann_properties_lax_hll(st, "Schw-HLLC", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
  }
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_hllc_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_lax_hll(st, "Kerr-HLLC", 0.3, 0.2, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
    run_riemann_properties_lax_hll(st, "Kerr-HLLC", 0.5, 0.0, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
    run_riemann_properties_lax_hll(st, "Kerr-HLLC", 0.4, 0.4, 0.0, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
  }
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 2b. Excision-boundary absorbing BC (Lax, HLL)
//
// Verifies that the new excision policy (zero-state Riemann with active-cell
// geometry — SESSION_NOTES_3.md §11.10–§11.13) gives correct one-sided upwind
// behavior at the BH horizon. Three invariants per case:
//   (i)   Wave-sum:   Σ w_k = q_R − q_L  (one side is zero)
//   (ii)  Flux-jump:  Σ s_k·w_k = F(q_R) − F(q_L)  with F(q_excised) ≡ 0
//   (iii) No spurious inflow into the active cell: under the wave_prop sign
//         convention, the active cell's update via the excision interface
//         (qL_new -= dt/dx · amdq if active on left, qR_new -= dt/dx · apdq
//         if active on right) must NOT increase the active cell's D when
//         matter is flowing into the BH.
//
// Setup: a 2-cell conf range with both cells filled from the same spacetime
// point — only EXCISION flag differs. This isolates the excision logic from
// any cross-cell metric variation.
// ---------------------------------------------------------------------------

static void
run_excision_absorbing_for_rp(struct gkyl_gr_spacetime *spacetime,
  const char *spacetime_label,
  enum gkyl_wv_gr_euler_tetrad_rp rp, int num_waves, const char *rp_label,
  double x, double y, double z)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);

  // 2-cell conf range so we have distinct slots for active + excised.
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  // Fill both cells at the same point. Active keeps the natural EXCISION
  // flag from the spacetime; mock cell flips it to -1.
  double *prods_active  = gkyl_array_fetch(prods, 0);
  double *prods_excised = gkyl_array_fetch(prods, 1);
  fill_prods_at(spacetime, x, y, z, prods_active);
  fill_prods_at(spacetime, x, y, z, prods_excised);
  if (prods_active[GKYL_GR_SP_EXCISION] < 0.0) {
    // Sample point is itself excised — can't test from here.
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }
  prods_excised[GKYL_GR_SP_EXCISION] = -1.0;

  double sqrt_det = sqrt(prods_active[GKYL_GR_SP_SPATIAL_DET]);
  (void)sqrt_det;

  // Active state — matter with rightward velocity (i.e. when on the LEFT
  // of an excised cell, this means matter flows INTO the BH).
  double rho = 1.0, p = 0.5;
  double v_into_excised[3] = { 0.4, 0.1, 0.0 };
  double q_active_glob[5];
  build_state_convA(eos, rho, v_into_excised, p, prods_active, q_active_glob);
  double q_excised_glob[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  int idx_active[1]  = { 0 };
  int idx_excised[1] = { 1 };
  double norm[3]  = { 1.0, 0.0, 0.0 };
  double tau1v[3] = { 0.0, 1.0, 0.0 };
  double tau2v[3] = { 0.0, 0.0, 1.0 };

  // ---- Subtest 1: excised on RIGHT (matter into BH) ----
  // L = active, R = excised. Matter flowing right (into BH).
  eqn->set_interface_idx_func(eqn, idx_active, idx_excised);
  double qL_loc[5], qR_loc[5];
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_active_glob,  qL_loc);  // → prodl_local (cell 0, active)
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_excised_glob, qR_loc);  // → prodr_local (cell 1, excised)

  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);

  double delta_R[5];
  for (int i = 0; i < 5; i++) delta_R[i] = qR_loc[i] - qL_loc[i];

  double waves_R[3 * 5], speeds_R[3];
  double maxs_R = eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    delta_R, qL_loc, qR_loc, 1.0, 1.0, waves_R, speeds_R);
  (void)maxs_R;

  // (i) Wave sum: Σ w_k = Δq.
  for (int i = 0; i < 5; i++) {
    double sum = 0.0;
    for (int k = 0; k < num_waves; k++) sum += waves_R[k * 5 + i];
    TEST_CHECK_( fabs(sum - delta_R[i]) < 1.0e-12,
      "[%s/%s/R-excised] wave sum: comp %d, |Σw − Δq| = %.3e",
      spacetime_label, rp_label, i, fabs(sum - delta_R[i]) );
  }

  // (ii) Flux jump: Σ s·w = F(qR) − F(qL) with F(qR_excised) = 0.
  //      Compute F(qL_active) via the production flux + flux_correction
  //      using the active-cell prods.
  double fL_sr[5], fL_gr[5];
  gkyl_gr_euler_tetrad_mod_flux(eos, qL_loc, grm->prodl_local, NULL, fL_sr);
  gkyl_gr_euler_tetrad_mod_flux_correction(eos, qL_loc, grm->prodl_local, NULL, fL_sr, fL_gr);
  double dF_R[5];
  for (int i = 0; i < 5; i++) dF_R[i] = -fL_gr[i];

  for (int i = 0; i < 5; i++) {
    double sw = 0.0;
    for (int k = 0; k < num_waves; k++) sw += speeds_R[k] * waves_R[k * 5 + i];
    TEST_CHECK_( fabs(sw - dF_R[i]) < 1.0e-10,
      "[%s/%s/R-excised] flux jump: comp %d, |Σs·w − ΔF| = %.3e",
      spacetime_label, rp_label, i, fabs(sw - dF_R[i]) );
  }

  // (iii) No spurious inflow into the active cell (on LEFT of excision).
  //       Cell L's update from this interface: qL_new -= dt/dx · amdq.
  //       For matter flowing INTO BH, D should not grow on L → amdq[D] ≥ 0.
  double amdq_R[5], apdq_R[5];
  eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    qL_loc, qR_loc, 1.0, 1.0, waves_R, speeds_R, amdq_R, apdq_R);
  TEST_CHECK_( amdq_R[0] >= -1.0e-12,
    "[%s/%s/R-excised] spurious D inflow: amdq[D]=%.3e (should be ≥ 0)",
    spacetime_label, rp_label, amdq_R[0] );

  // Conservation: amdq + apdq = ΔF.
  for (int i = 0; i < 5; i++) {
    TEST_CHECK_( fabs(amdq_R[i] + apdq_R[i] - dF_R[i]) < 1.0e-10,
      "[%s/%s/R-excised] fluct balance: comp %d, |amdq+apdq − ΔF| = %.3e",
      spacetime_label, rp_label, i, fabs(amdq_R[i] + apdq_R[i] - dF_R[i]) );
  }

  // ---- Subtest 2: excised on LEFT (active state on right) ----
  // L = excised, R = active. Active has v_x > 0 (now pointing away from BH).
  // Matter cannot emerge from the BH; the wave decomposition with q=0 on left
  // should not produce inflow into the active cell from the excised side.
  eqn->set_interface_idx_func(eqn, idx_excised, idx_active);
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_excised_glob, qL_loc); // → prodl_local (excised)
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_active_glob,  qR_loc); // → prodr_local (active)

  double delta_L[5];
  for (int i = 0; i < 5; i++) delta_L[i] = qR_loc[i] - qL_loc[i];

  double waves_L[3 * 5], speeds_L[3];
  double maxs_L = eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    delta_L, qL_loc, qR_loc, 1.0, 1.0, waves_L, speeds_L);
  (void)maxs_L;

  for (int i = 0; i < 5; i++) {
    double sum = 0.0;
    for (int k = 0; k < num_waves; k++) sum += waves_L[k * 5 + i];
    TEST_CHECK_( fabs(sum - delta_L[i]) < 1.0e-12,
      "[%s/%s/L-excised] wave sum: comp %d, |Σw − Δq| = %.3e",
      spacetime_label, rp_label, i, fabs(sum - delta_L[i]) );
  }

  double fR_sr[5], fR_gr[5];
  gkyl_gr_euler_tetrad_mod_flux(eos, qR_loc, grm->prodr_local, NULL, fR_sr);
  gkyl_gr_euler_tetrad_mod_flux_correction(eos, qR_loc, grm->prodr_local, NULL, fR_sr, fR_gr);
  double dF_L[5];
  for (int i = 0; i < 5; i++) dF_L[i] = fR_gr[i];  // F(qR_active) − F(qL_excised=0)

  for (int i = 0; i < 5; i++) {
    double sw = 0.0;
    for (int k = 0; k < num_waves; k++) sw += speeds_L[k] * waves_L[k * 5 + i];
    TEST_CHECK_( fabs(sw - dF_L[i]) < 1.0e-10,
      "[%s/%s/L-excised] flux jump: comp %d, |Σs·w − ΔF| = %.3e",
      spacetime_label, rp_label, i, fabs(sw - dF_L[i]) );
  }

  // No spurious inflow into active cell (on RIGHT). Cell R's update from
  // this interface: qR_new -= dt/dx · apdq. The "absorbing BC" principle
  // is that matter cannot emerge from the BH (vacuum cannot be a SOURCE
  // of matter for the active cell). So apdq[D] ≥ 0 is required — the
  // active cell either stays put (apdq[D] = 0) or loses matter to the
  // vacuum (apdq[D] > 0, e.g. via a left-going rarefaction when the
  // active flow is subsonic outflow into vacuum). What's forbidden is
  // apdq[D] < 0, which would have the BH spuriously supplying matter
  // to the active cell.
  double amdq_L[5], apdq_L[5];
  eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    qL_loc, qR_loc, 1.0, 1.0, waves_L, speeds_L, amdq_L, apdq_L);
  TEST_CHECK_( apdq_L[0] >= -1.0e-12,
    "[%s/%s/L-excised] spurious D from BH: apdq[D]=%.3e (should be ≥ 0)",
    spacetime_label, rp_label, apdq_L[0] );

  for (int i = 0; i < 5; i++) {
    TEST_CHECK_( fabs(amdq_L[i] + apdq_L[i] - dF_L[i]) < 1.0e-10,
      "[%s/%s/L-excised] fluct balance: comp %d, |amdq+apdq − ΔF| = %.3e",
      spacetime_label, rp_label, i, fabs(amdq_L[i] + apdq_L[i] - dF_L[i]) );
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_excision_absorbing_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_excision_absorbing_for_rp(st, "Mink", WV_GR_EULER_TETRAD_RP_LAX, 2, "Lax",
    0.3, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_excision_absorbing_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_excision_absorbing_for_rp(st, "Schw", WV_GR_EULER_TETRAD_RP_LAX, 2, "Lax",
    0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_excision_absorbing_hll_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_excision_absorbing_for_rp(st, "Mink", WV_GR_EULER_TETRAD_RP_HLL, 2, "HLL",
    0.3, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_excision_absorbing_hll_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_excision_absorbing_for_rp(st, "Schw", WV_GR_EULER_TETRAD_RP_HLL, 2, "HLL",
    0.3, 0.2, 0.0);
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

// Per-cell-metric positivity check. prods_L and prods_R are the
// spacetime products rows for the LEFT and RIGHT cells (production
// wave_prop pipes adjacent-cell metrics into prodl_local/prodr_local
// independently). Each cell's conservatives are densitized with its
// OWN √γ and admissibility-checked in its OWN inverse spatial metric,
// matching what wave_prop does post-POSITIVITY_SWEEP.
static void
run_positivity_for_rp(struct gkyl_wv_eqn *eqn,
  struct wv_gr_euler_tetrad_mod *grm,
  const double *prods_L, const double *prods_R,
  const struct rp_case *rc,
  int *D_violations, int *S2_violations, int *tau_violations)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  double sqrt_det_L = sqrt(prods_L[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_R = sqrt(prods_R[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g_L[3][3] = {
    { prods_L[GKYL_GR_SP_INV_GIJ + 0], prods_L[GKYL_GR_SP_INV_GIJ + 1], prods_L[GKYL_GR_SP_INV_GIJ + 2] },
    { prods_L[GKYL_GR_SP_INV_GIJ + 3], prods_L[GKYL_GR_SP_INV_GIJ + 4], prods_L[GKYL_GR_SP_INV_GIJ + 5] },
    { prods_L[GKYL_GR_SP_INV_GIJ + 6], prods_L[GKYL_GR_SP_INV_GIJ + 7], prods_L[GKYL_GR_SP_INV_GIJ + 8] },
  };
  double inv_g_R[3][3] = {
    { prods_R[GKYL_GR_SP_INV_GIJ + 0], prods_R[GKYL_GR_SP_INV_GIJ + 1], prods_R[GKYL_GR_SP_INV_GIJ + 2] },
    { prods_R[GKYL_GR_SP_INV_GIJ + 3], prods_R[GKYL_GR_SP_INV_GIJ + 4], prods_R[GKYL_GR_SP_INV_GIJ + 5] },
    { prods_R[GKYL_GR_SP_INV_GIJ + 6], prods_R[GKYL_GR_SP_INV_GIJ + 7], prods_R[GKYL_GR_SP_INV_GIJ + 8] },
  };

  // Build L/R states from primitives, each in its OWN metric. Same
  // primitives → different densitized conservatives because
  // build_state_convA uses √γ and γ_ij of the supplied prods row.
  double qL[5], qR[5];
  build_state_convA(eos, rc->rho_L, rc->v_L, rc->p_L, prods_L, qL);
  build_state_convA(eos, rc->rho_R, rc->v_R, rc->p_R, prods_R, qR);

  // Sanity: each cell admissible in its own metric.
  double qL_und[5], qR_und[5];
  for (int i = 0; i < 5; i++) {
    qL_und[i] = qL[i] / sqrt_det_L;
    qR_und[i] = qR[i] / sqrt_det_R;
  }
  bool d_ok, s_ok, t_ok;
  record_admissibility(inv_g_L, qL_und, &d_ok, &s_ok, &t_ok);
  if (!(d_ok && s_ok && t_ok)) {
    fprintf(stderr, "  [positivity] WARNING: input qL inadmissible for case '%s' "
      "(D_ok=%d, S2_ok=%d, tau_ok=%d) — skipping\n", rc->name, d_ok, s_ok, t_ok);
    return;
  }
  record_admissibility(inv_g_R, qR_und, &d_ok, &s_ok, &t_ok);
  if (!(d_ok && s_ok && t_ok)) {
    fprintf(stderr, "  [positivity] WARNING: input qR inadmissible for case '%s' "
      "(D_ok=%d, S2_ok=%d, tau_ok=%d) — skipping\n", rc->name, d_ok, s_ok, t_ok);
    return;
  }

  // Per-cell metric: idxl=0 → prods_L, idxr=1 → prods_R. rotate_to_local
  // populates prodl_local/prodr_local via the parity flip; with norm=+x̂
  // the rotation is a no-op on q components.
  int idx_L[1] = { 0 }, idx_R[1] = { 1 };
  eqn->set_interface_idx_func(eqn, idx_L, idx_R);
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

  // Check each updated cell in its OWN metric.
  for (int side = 0; side < 2; side++) {
    const double *qn = (side == 0) ? qL_new : qR_new;
    const char *who  = (side == 0) ? "L" : "R";
    double sqrt_det_s = (side == 0) ? sqrt_det_L : sqrt_det_R;
    const double (*inv_g_s)[3] = (side == 0) ? inv_g_L : inv_g_R;
    double qn_und[5];
    for (int i = 0; i < 5; i++) qn_und[i] = qn[i] / sqrt_det_s;
    bool D_ok, S2_ok, tau_ok;
    record_admissibility(inv_g_s, qn_und, &D_ok, &S2_ok, &tau_ok);
    if (!D_ok)   { (*D_violations)++;   fprintf(stderr,
      "  [positivity] '%s' cell %s: D = %.3e\n", rc->name, who, qn_und[0]); }
    if (!S2_ok)  { (*S2_violations)++;  fprintf(stderr,
      "  [positivity] '%s' cell %s: s² = %.3e (D+τ = %.3e, |S|² = %.3e)\n",
      rc->name, who, (qn_und[0]+qn_und[4])*(qn_und[0]+qn_und[4])
                       - gkyl_gr_euler_mom_sq(inv_g_s, qn_und[1], qn_und[2], qn_und[3]),
      qn_und[0]+qn_und[4],
      gkyl_gr_euler_mom_sq(inv_g_s, qn_und[1], qn_und[2], qn_und[3])); }
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

// Physical sample points for curved-spacetime positivity probes.
// Each point is the *interface centroid*: the test fills prods_L at
// (x - dx/2, y, z) and prods_R at (x + dx/2, y, z), mirroring how
// production wave_prop reads metrics from adjacent cells.
//
// Sample points span:
//   - origin (Minkowski only — curved spacetimes' excision filters)
//   - mid-radius axial / off-axis (γ_ij with off-diagonals from Kerr-Schild)
//   - BHL bow-shock front (just outside r_h for M=0.3 — sharp √γ jump)
//   - BHL inflow region (further out from horizon)
struct positivity_point { double x, y, z; };
static const struct positivity_point g_positivity_points[] = {
  { 0.0,   0.0,   0.0 },
  { 0.3,   0.2,   0.0 },
  { 0.5,   0.0,   0.0 },
  { 0.4,   0.4,   0.0 },
  // BHL-production-relevant: M=0.3 → r_h=0.6 in Schw/Kerr. These points
  // are valid for the M=0.3 case constructed below; Minkowski/Kerr-low-M
  // tests use the same coords but the excision check filters them.
  { 0.65,  0.0,   0.0 },   // axial, just outside r_h=0.6
  {-0.05, -0.6,   0.0 },   // bow-shock front, sharp √γ variation
  {-0.05,  0.6,   0.0 },   // mirror
  { 1.0,   1.0,   0.0 },   // off-axis bulk (γ_xy nonzero)
};

// Spacing between L and R cell-center metric samples. Matches production
// BHL grid (Lx=5, Nx=256 → dx≈0.0195). Small enough that √γ_L ≈ √γ_R
// for smooth metrics, but large enough at sharp-curvature points (near
// horizon) to expose the densitization mismatch hidden by single-cell
// tests.
#define GR_EULER_POSITIVITY_DX (0.01953)

static void
run_positivity_sweep(struct gkyl_gr_spacetime *spacetime,
  struct gkyl_gr_euler_eos eos,
  enum gkyl_wv_gr_euler_tetrad_rp rp,
  const char *label)
{
  double gas_gamma = eos.gas_gamma;

  // 2-cell conf_range so set_interface_idx(0, 1) routes prods_L and
  // prods_R to wave_lax_curved's prodl_local / prodr_local.
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);

  for (int pi = 0; pi < n_pts; pi++) {
    const struct positivity_point *pt = &g_positivity_points[pi];
    // L cell at (x - dx/2, y, z); R cell at (x + dx/2, y, z).
    double dx_half = 0.5 * GR_EULER_POSITIVITY_DX;
    fill_prods_at(spacetime, pt->x - dx_half, pt->y, pt->z, prods_L);
    fill_prods_at(spacetime, pt->x + dx_half, pt->y, pt->z, prods_R);
    if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
        prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;

    for (int c = 0; c < n_cases; c++) {
      run_positivity_for_rp(eqn, grm, prods_L, prods_R, &g_positivity_cases[c],
        &D_violations, &S2_violations, &tau_violations);
      total_cells += 2;
    }
  }

  fprintf(stderr,
    "[positivity %s] across %d cell-updates (per-cell-metric, dx=%.4f): "
    "D<=0: %d (%.1f%%), s²<=0: %d (%.1f%%), τ<0: %d (%.1f%%)\n",
    label, total_cells, GR_EULER_POSITIVITY_DX,
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

// Each positivity sweep loops over [IDEAL, TM, RCC] — admissibility
// preservation is an algorithm property that should hold across all
// production EOS configurations.
void test_positivity_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, "Lax-Mink");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, "Lax-Schw");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, "Lax-Kerr");
  gkyl_gr_spacetime_release(st);
}

void test_positivity_hll_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, "HLL-Mink");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hll_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, "HLL-Schw");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hll_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, "HLL-Kerr");
  gkyl_gr_spacetime_release(st);
}

// HLLC positivity sweep — the headline test: HLLC's star-state
// construction is supposed to give explicit τ-positivity from
// admissible inputs (MB05 §3.1.2). If HLLC genuinely preserves
// τ-positivity in regimes where Lax/HLL fail, we'll see zero
// τ-violations in the sweep where Lax/HLL produce some.
void test_positivity_hllc_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, "HLLC-Mink");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hllc_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, "HLLC-Schw");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hllc_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, "HLLC-Kerr");
  gkyl_gr_spacetime_release(st);
}

// BHL-production variants: M=0.3 (r_h=0.6) Schwarzschild matches the
// rt_gr_bhl_static_tetrad_mod.lua input. The sample-point list in
// g_positivity_points contains points just outside r_h=0.6 and on the
// bow-shock front, which only become physical under the M=0.3 metric.
// (Under M=0.1 those same points are deeper in the strong-field region
// but the excision check filters anything inside r_h=0.2.)
void test_positivity_lax_bhl(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, "Lax-BHL(M=0.3)");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hll_bhl(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, "HLL-BHL(M=0.3)");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_hllc_bhl(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, "HLLC-BHL(M=0.3)");
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
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
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
      build_state_convA(eos, rc->rho_min, rc->v_min, p_min, prods, qL);
      build_state_convA(eos, rc->rho_typ, rc->v_typ, rc->p_typ, prods, qR);
    } else {
      build_state_convA(eos, rc->rho_typ, rc->v_typ, rc->p_typ, prods, qL);
      build_state_convA(eos, rc->rho_min, rc->v_min, p_min,    prods, qR);
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
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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

void test_near_floor_hllc(void)
{
  run_near_floor_sweep(WV_GR_EULER_TETRAD_RP_HLLC, "HLLC");
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
  const char *label, double x, double y, double z,
  struct gkyl_gr_euler_eos eos)
{
  double gas_gamma = eos.gas_gamma;
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
        build_state_convA(eos, rho, v, p, prods, q);
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

// Each round-trip wrapper loops over [IDEAL, TM, RCC] — primitive
// recovery roundtrip is the most direct test of the EOS dispatch.
void test_round_trip_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_round_trip(st, "Mink", 0.3, 0.0, 0.0, eos_modes[ei]);
  gkyl_gr_spacetime_release(st);
}
void test_round_trip_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_round_trip(st, "Schw", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_round_trip(st, "Schw", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_round_trip(st, "Schw", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}
void test_round_trip_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_round_trip(st, "Kerr", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_round_trip(st, "Kerr", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_round_trip(st, "Kerr", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}

// (b) Primitive-recovery consistency between curved-frame prim_vars and
//     SR Newton on q_tet. If both Newton solves are correct, they
//     should agree on (ρ, p, W) to machine precision. v_curved^i and
//     v_tet^a are related by the tetrad transform; v_tet^0 should be
//     v^x / √γ^{xx} by our construction.
static void
run_prim_consistency(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z,
  struct gkyl_gr_euler_eos eos)
{
  double gas_gamma = eos.gas_gamma;
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
  build_state_convA(eos, rho_in, v_in, p_in, prods, q);

  // Curved-frame recovery via shared helper.
  double D    = q[0] / sqrt_det;
  double momx = q[1] / sqrt_det;
  double momy = q[2] / sqrt_det;
  double momz = q[3] / sqrt_det;
  double tau  = q[4] / sqrt_det;
  struct gkyl_gr_euler_prim prim_curved;
  gkyl_gr_euler_recover_primitives(eos,
    D, momx, momy, momz, tau, inv_g, NULL, &prim_curved);

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

// Each prim-consistency wrapper loops over [IDEAL, TM, RCC].
void test_prim_consistency_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_prim_consistency(st, "Mink", 0.3, 0.0, 0.0, eos_modes[ei]);
  gkyl_gr_spacetime_release(st);
}
void test_prim_consistency_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_prim_consistency(st, "Schw @ (0.3,0.2,0)", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Schw @ (0.5,0.0,0)", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Schw @ (0.4,0.4,0)", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}
void test_prim_consistency_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_prim_consistency(st, "Kerr @ (0.3,0.2,0)", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Kerr @ (0.5,0.0,0)", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Kerr @ (0.4,0.4,0)", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
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
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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
    build_state_convA(eos, cases[c].rL, v_co_l, cases[c].pL, prods_row, qL_glob);
    build_state_convA(eos, cases[c].rR, v_co_r, cases[c].pR, prods_row, qR_glob);

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
    banyuls_delta_flux(eos, grm, qL, qR, dF);
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
void test_bhl_regime_hllc_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_bhl_regime_states(st, "Schw-HLLC", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_HLLC, 3);
  run_bhl_regime_states(st, "Schw-HLLC", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_HLLC, 3);
  gkyl_gr_spacetime_release(st);
}
void test_bhl_regime_hllc_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_bhl_regime_states(st, "Kerr-HLLC", 0.3, 0.2, 0.0,
    WV_GR_EULER_TETRAD_RP_HLLC, 3);
  run_bhl_regime_states(st, "Kerr-HLLC", 0.4, 0.4, 0.0,
    WV_GR_EULER_TETRAD_RP_HLLC, 3);
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
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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
    build_state_convA(eos, rho_L, v_L, p_L, prods_row, qL_glob);
    build_state_convA(eos, rho_R, v_R, p_R, prods_row, qR_glob);

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
    banyuls_delta_flux(eos, grm, qL, qR, dF);

    double max_fj = 0.0, max_fb = 0.0;
    int nw = eqn->num_waves;
    for (int k = 0; k < 5; k++) {
      double sw = 0.0;
      for (int w = 0; w < nw; w++) sw += speeds[w] * waves[w * 5 + k];
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

void test_floor_precision_hllc(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_floor_precision_sweep(st, WV_GR_EULER_TETRAD_RP_HLLC, "HLLC");
  gkyl_gr_spacetime_release(st);
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
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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
    build_state_convA(eos, cases[c].rho_L, cases[c].v_L, cases[c].p_L,
      prods_row, qL);
    build_state_convA(eos, cases[c].rho_R, cases[c].v_R, cases[c].p_R,
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
// HLLC: this is the headline τ-positivity test. Lax and HLL produce
// τ-violations in the τ/D ≪ 1 regime (BHL bow-shock corner). HLLC's
// star-state construction is *supposed* to make τ ≥ 0 from admissible
// inputs; if it works here in Minkowski, that's strong evidence it will
// fix the wave_prop violations on the BHL run as well.
void test_small_tau_over_D_hllc(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_small_tau_over_D_sweep(st, WV_GR_EULER_TETRAD_RP_HLLC, "Mink-HLLC");
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
  enum gkyl_wv_gr_euler_tetrad_rp rp, const char *rp_label,
  enum gkyl_wv_flux_type flux_type,
  double sx, double sy, double sz)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  // Sample point passed in by caller. Default canonical convA point is
  // (0.3, 0.2, 0); use a near-horizon point (e.g. 0.08, 0.04, 0 for
  // M=0.1 Schw with r_h=0.2) to test admissibility under extreme β/γ.
  fill_prods_at(spacetime, sx, sy, sz, prods_row);
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
    // BHL repair #2 dump from rpType=hll, τ_floor=1e-8 run. Production
    // log shows qt[4]=+0.693 (τ>0) but s² ≈ −0.52 (failure mode is s²<0,
    // not τ<0). All three input cells are admissible in Minkowski
    // (s²_mink > 1 each). Use this case to isolate whether HLL HIGH
    // can generate s²<0 from admissible inputs at CFL ≤ 1.
    { "bhl-repair-s2-#2",
      3.764, 3.765, -1.130, 0.000, 0.418,
      3.245, 3.185, -1.073, 0.000, 0.332,
      3.510, 2.847,  0.308, 0.000, 0.127 },
    // Dump #10 from rpType=hll + curved Lax LOW (BHL bow-shock cell
    // i=125 near the BH). Pure s²-only failure: all 3 inputs are
    // admissible (s²_cart > 2.0 each), qt has τ=+1.86 (positive)
    // but |S|²_cart = 32.25 > (D+τ)² = 30.78, so s²_cart = −1.47.
    // This is the canonical case showing curved Lax DOES still leak
    // s²-positivity on supersonic bow-shock cells with sharp gradients,
    // even with admissible inputs.
    { "bhl-repair-curvedLax-s2-#10",
      4.313, 4.221,  0.779, 0.000, 1.084,
      4.357, 4.743, -1.245, 0.000, 1.187,
      1.778, 1.109,  0.346, 0.000, 0.067 },
    // ----- Production fires from rpType="hll" + curved Lax LOW with
    // full-3D amax (SESSION_NOTES_3 §17). 8749 wave_prop s²<0 fires
    // remain even with the wider envelope, so the dominant mechanism is
    // NOT off-diagonal-γ_xy with non-aligned momentum. Three new dumps
    // exercise the residual pathology:
    //
    // full3d-s2-#10  dir=0 i=125 — large-D x-aligned (S_y ≈ 0 in qC) but
    //   neighbors carry strong y-momentum (qin_L has S_y=-1.54). qt has
    //   τ_local=+1.07 but s²_cart = -0.31. Pure x-sweep on a cell whose
    //   neighbors have transverse momentum.
    // full3d-s2-#12  dir=0 i=126 — moderate-D mixed momentum, qC has
    //   strong x and modest y; qt has τ_local=+0.55 but s²_cart = -0.05.
    // full3d-s2-#17  dir=1 i=92 — y-SWEEP (norm=(0,1,0))! qC has strong
    //   y-momentum (S_y=-0.76) and even stronger x-momentum (S_x=2.6);
    //   qt has τ_local=+0.42 but s²_cart = -0.04. Inputs entered into
    //   the unit test in unrotated GLOBAL frame; the rotate_to_local in
    //   the test below uses norm=(1,0,0) so this case effectively maps
    //   a "y-sweep production state" onto an "x-sweep test state" —
    //   inconsistent with the actual dir=1 sweep in production, but
    //   still useful to check whether curved Lax fails on the raw state.
    { "bhl-repair-full3d-s2-#10",
      2.528, 2.695, -0.005, 0.000, 0.655,
      3.578, 3.814, -1.539, 0.000, 0.987,
      0.878, 0.670,  0.193, 0.000, 0.026 },
    { "bhl-repair-full3d-s2-#12",
      1.059, 0.935,  0.167, 0.000, 0.117,
      2.364, 3.174, -0.673, 0.000, 0.934,
      0.244, 0.182,  0.082, 0.000, 0.015 },
    { "bhl-repair-full3d-s2-#17",
      1.695, 2.594, -0.757, 0.000, 0.986,
      1.194, 1.231, -0.184, 0.000, 0.303,
      2.087, 2.902,  0.585, 0.000, 1.045 },
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
    double maxs_LC = eqn->waves_func(eqn, flux_type,
      delta_LC, qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC);
    double amdq_LC[5], apdq_LC[5];
    eqn->qfluct_func(eqn, flux_type,
      qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC, amdq_LC, apdq_LC);

    eqn->set_interface_idx_func(eqn, idx, idx);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qC_glob, qC_loc_CR);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qR_glob, qR_loc);
    double delta_CR[5];
    for (int i = 0; i < 5; i++) delta_CR[i] = qR_loc[i] - qC_loc_CR[i];
    double waves_CR[3 * 5], speeds_CR[3];
    double maxs_CR = eqn->waves_func(eqn, flux_type,
      delta_CR, qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR);
    double amdq_CR[5], apdq_CR[5];
    eqn->qfluct_func(eqn, flux_type,
      qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR, amdq_CR, apdq_CR);

    double maxs = fmax(maxs_LC, maxs_CR);
    if (!(maxs > 0.0)) { fprintf(stderr, "  [%s] degenerate\n", cases[c].name); continue; }
    double dt_dx = 0.5 / maxs;  // 2D directional-splitting CFL

    double qC_new[5];
    for (int i = 0; i < 5; i++) {
      qC_new[i] = qC_loc_CR[i] - dt_dx * (amdq_CR[i] + apdq_LC[i]);
    }
    double tau_C_new = qC_new[4] / sqrt_det;

    // Also check s² invariant: s² = (D+τ)² − γ^{ij}·S_i·S_j (undensitized).
    double qC_und[5];
    for (int i = 0; i < 5; i++) qC_und[i] = qC_new[i] / sqrt_det;
    double mom_sq_new = gkyl_gr_euler_mom_sq(inv_g, qC_und[1], qC_und[2], qC_und[3]);
    double s_sq_new   = (qC_und[0] + qC_und[4]) * (qC_und[0] + qC_und[4]) - mom_sq_new;

    fprintf(stderr,
      "  [%s/%-16s] τ_C/D_C = %.4f    τ_C_new = %+.3e %s   s²_new = %+.3e %s\n",
      rp_label, cases[c].name, tau_over_D, tau_C_new,
      tau_C_new < 0.0 ? "**τ<0**" : "      ",
      s_sq_new,
      s_sq_new <= 0.0 ? "**s²≤0**" : "       ");

    // Report-only TEST_MSG so the test PASSES while documenting the
    // pathology.
    TEST_MSG( "[%s/%s] τ/D = %.4f, τ_new = %+.3e, s²_new = %+.3e",
      rp_label, cases[c].name, tau_over_D, tau_C_new, s_sq_new );
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_direct_state_lax(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  // LOW_ORDER call. wave_lax_l ignores type, so this exercises Lax.
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_LAX, "Lax",
    GKYL_WV_LOW_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}

// HLLC fallback probe. Walks the same three BHL repair cases through
// gkyl_gr_euler_tetrad_mod_sr_hllc_minkowski directly with the
// did_fallback outparam so we can see, per L-C and C-R interface,
// whether HLLC's full star-state construction was used or the cold-gas
// fallback to HLL kicked in. Helps disambiguate τ-positivity failures:
// if a τ<0 case has fallback=1 on the violating interface, the failure
// is HLL's known weakness; if fallback=0, the failure is HLLC's own
// star-state construction interacting with the curved-frame back-
// transform.
//
// The setup mirrors run_direct_state_sweep step-for-step so the inputs
// are byte-identical: load a Schw spacetime at (0.3, 0.2, 0), build the
// interface metric (mean of L+C and C+R), forward-transform to the
// Gram-Schmidt-on-γ⁻¹ tetrad, and call sr_hllc_minkowski directly.
static void
run_direct_state_hllc_fallback_probe(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, WV_GR_EULER_TETRAD_RP_HLLC);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(st, 0.3, 0.2, 0.0, prods_row);
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    gkyl_gr_spacetime_release(st);
    return;
  }

  // Same three BHL cases as run_direct_state_sweep (verbatim wave_prop
  // dumps). qC is the cell that fails in the BHL run; qL, qR neighbors.
  // Values are DENSITIZED (in code form, as they came out of qin in
  // wave_prop.c).
  struct {
    const char *name;
    double D_C, Sx_C, Sy_C, Sz_C, tau_C;
    double D_L, Sx_L, Sy_L, Sz_L, tau_L;
    double D_R, Sx_R, Sy_R, Sz_R, tau_R;
  } cases[] = {
    { "bhl-repair-#0",
      2.978, 2.499, 0.118, 0.000, 0.122,
      2.541, 2.146, 0.045, 0.000, 0.136,
      3.395, 2.827, 0.222, 0.000, 0.095 },
    { "bhl-repair-#4",
      3.407, 2.633, 0.209, 0.000, 0.021,
      3.123, 2.165, 0.119, 0.000, 0.025,
      3.009, 3.319, -0.766, 0.000, 0.526 },
    { "bhl-repair-#19",
      3.505, 2.969, -1.664, 0.000, 0.193,
      3.872, 2.574, 0.285, 0.000, 0.025,
      2.251, 2.424, -0.887, 0.000, 0.392 },
  };

  // Build the interface metric and triad ONCE. All three cases share
  // the same prods row so γ_ij, α, β are uniform across the cells —
  // i.e. the Gram-Schmidt triad is the same on the L-C and C-R
  // interfaces and on every case. This is the same triad the HLLC
  // wave_hllc dispatch builds at runtime.
  double g_iface[3][3], inv_g_iface[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_iface[i][j]     = prods_row[GKYL_GR_SP_GIJ + 3*i + j];
      inv_g_iface[i][j] = prods_row[GKYL_GR_SP_INV_GIJ + 3*i + j];
    }
  double sqrt_det = sqrt(prods_row[GKYL_GR_SP_SPATIAL_DET]);
  double M[3][3], M_inv[3][3];
  gkyl_gr_euler_tetrad_mod_build_triad_contravariant_x(
    g_iface, inv_g_iface, M, M_inv);

  fprintf(stderr,
    "[hllc-fallback-probe Schw] case             L-C-fallback   C-R-fallback   notes\n");
  for (size_t c = 0; c < sizeof(cases)/sizeof(*cases); c++) {
    // Densitized states (as they live in code form). The state in the
    // tetrad-frame solver is undensitized — q_to_tetrad_contra divides
    // by √γ, so we pass densitized q as it expects.
    double qL[5] = { cases[c].D_L, cases[c].Sx_L, cases[c].Sy_L, cases[c].Sz_L, cases[c].tau_L };
    double qC[5] = { cases[c].D_C, cases[c].Sx_C, cases[c].Sy_C, cases[c].Sz_C, cases[c].tau_C };
    double qR[5] = { cases[c].D_R, cases[c].Sx_R, cases[c].Sy_R, cases[c].Sz_R, cases[c].tau_R };

    // Forward-transform to tetrad. (Same √γ for all sides since uniform
    // background; matches the production dispatch's per-side strip.)
    double qL_tet[5], qC_tet[5], qR_tet[5];
    gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(qL, sqrt_det, inv_g_iface, M_inv, qL_tet);
    gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(qC, sqrt_det, inv_g_iface, M_inv, qC_tet);
    gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(qR, sqrt_det, inv_g_iface, M_inv, qR_tet);

    // L-C interface.
    double waves_LC[3 * 5], speeds_LC[3];
    struct gkyl_gr_euler_tetrad_mod_hllc_diag diag_LC = {0};
    gkyl_gr_euler_tetrad_mod_sr_hllc_minkowski(
      eos, qL_tet, qC_tet,
      NULL, waves_LC, speeds_LC, &diag_LC);

    // C-R interface.
    double waves_CR[3 * 5], speeds_CR[3];
    struct gkyl_gr_euler_tetrad_mod_hllc_diag diag_CR = {0};
    gkyl_gr_euler_tetrad_mod_sr_hllc_minkowski(
      eos, qC_tet, qR_tet,
      NULL, waves_CR, speeds_CR, &diag_CR);

    // fallback_reason key (post tightening of the fallback policy):
    //   1 = lam_diff < 1e-14 (degenerate Davis bracket — λ_L ≈ λ_R)
    //   2 = λ* not finite (sqrt of negative discriminant; should be
    //       caught by the disc<0 → 0 clamp inside, so rare)
    //   3 = |λ_L − λ*| < tol (would blow up 1/(λ_L − λ*) in U_L*)
    //   4 = |λ_R − λ*| < tol (would blow up 1/(λ_R − λ*) in U_R*)
    // λ* outside [λ_L, λ_R] (e.g. supersonic flow) is NOT a fallback
    // trigger — the wave decomposition is still conservative there.
    static const char *reason_str[] = {
      "—",                  // 0
      "lam_diff~0",         // 1
      "λ* not finite",      // 2
      "|λ_L−λ*| < tol",     // 3
      "|λ_R−λ*| < tol",     // 4
    };
    fprintf(stderr,
      "  [%-16s] τ_C/D_C=%.4f\n", cases[c].name, cases[c].tau_C / cases[c].D_C);
    fprintf(stderr,
      "    L-C: fb=%d (%s)  λ_L=%+.4f  λ_R=%+.4f  λ*=%+.4f\n",
      diag_LC.did_fallback, reason_str[diag_LC.fallback_reason],
      diag_LC.lambda_L, diag_LC.lambda_R, diag_LC.lambda_star);
    fprintf(stderr,
      "    C-R: fb=%d (%s)  λ_L=%+.4f  λ_R=%+.4f  λ*=%+.4f\n",
      diag_CR.did_fallback, reason_str[diag_CR.fallback_reason],
      diag_CR.lambda_L, diag_CR.lambda_R, diag_CR.lambda_star);
    TEST_MSG( "[%s] L-C fb=%d (reason %d), C-R fb=%d (reason %d)",
      cases[c].name,
      diag_LC.did_fallback, diag_LC.fallback_reason,
      diag_CR.did_fallback, diag_CR.fallback_reason );
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
  gkyl_gr_spacetime_release(st);
}

void test_direct_state_hllc_fallback_probe(void)
{
  run_direct_state_hllc_fallback_probe();
}
// HLLC LOW_ORDER. Post the wave_hllc_l fix (HIGH→HLLC, LOW→Lax), this
// invokes Lax (the POSITIVITY_SWEEP cleanup branch). So it should give
// identical τ_C_new to test_direct_state_lax — sanity-check that the
// dispatch is wired up correctly.
void test_direct_state_hllc(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_HLLC, "HLLC-LOW",
    GKYL_WV_LOW_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}
// HLLC HIGH_ORDER. The production first-sweep behavior. wave_prop runs
// HIGH_ORDER at the FIRST_SWEEP stage; if a cell fails admissibility,
// POSITIVITY_SWEEP retries with LOW_ORDER (Lax). The composite
// production behavior is the union of this test's pass-set with
// test_direct_state_lax's pass-set: a cell passes overall iff at least
// one of (HLLC HIGH, Lax LOW) gives an admissible result.
void test_direct_state_hllc_high_order(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_HLLC, "HLLC-HIGH",
    GKYL_WV_HIGH_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}

// HLL HIGH_ORDER: this is what wave_prop runs in FIRST_SWEEP for rpType="hll".
// In particular this should show what's responsible for the production
// s²<0 wave_prop fires (when run in Mink, the bhl-repair-s2-#2 case mimics
// dump #2 from the τ_floor=1e-8 production run almost exactly).
void test_direct_state_hll_high_order(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_HLL, "HLL-HIGH",
    GKYL_WV_HIGH_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}

// Same in Minkowski — pure flat SR baseline for the s²<0 question.
// If HLL HIGH in Mink reproduces s²<0 on admissible inputs, that's a
// real Mignone-Bodo invariant-domain violation. If not, the production
// s²<0 fires likely come from MUSCL reconstruction or curved-metric
// wrap, not the underlying SR HLL.
void test_direct_state_hll_high_order_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_direct_state_sweep(st, "Mink", WV_GR_EULER_TETRAD_RP_HLL, "HLL-HIGH",
    GKYL_WV_HIGH_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_direct_state_lax_high_order_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_direct_state_sweep(st, "Mink", WV_GR_EULER_TETRAD_RP_LAX, "Lax-HIGH",
    GKYL_WV_HIGH_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}

// HIGH_ORDER Lax on Schwarzschild = tetrad-first Lax (the path that
// produced the production s²<0 dumps). Pairs with test_direct_state_lax
// (which is LOW_ORDER → curved Lax) so the two together give an A/B
// comparison: HIGH_ORDER row should show s²<0 on bhl-repair-s2-#2,
// LOW_ORDER row should show s²>0 (the curved-frame Lax preserving A_γ).
void test_direct_state_lax_high_order_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_LAX, "Lax-HIGH",
    GKYL_WV_HIGH_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}
// LOW_ORDER Lax on Schwarzschild = wave_lax_curved (the new direct
// curved-frame Lax). After the wiring change, test_direct_state_lax
// already exercises this — but make the contrast explicit with a
// matching wrapper that has the new name baked in.
void test_direct_state_lax_curved_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, 0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}

// Near-horizon sweep at multiple radii outside r_h = 2M = 0.2 for M=0.1
// Schwarzschild. Tests whether wave_lax_curved leaks s²-positivity at
// extreme β^x / inflated γ_xx regimes — the regime the BHL production
// run actually has at the bow-shock cells immediately outside the BH
// horizon.
// ---------------------------------------------------------------------------
// Per-cell-metric variant of run_direct_state_sweep
// ---------------------------------------------------------------------------
// run_direct_state_sweep uses a SINGLE metric for L, C, R cells (the
// equation's conf_range has volume 1 and idxl=idxr=0). Production wave_prop
// uses adjacent cells with DIFFERENT metrics: prodl_local comes from
// idx_L=(i-1, j) and prodr_local from idx_R=(i, j). The curved-Banyuls
// flux is F = α·√γ · (...), so F_L uses (α_L, √γ_L) and F_R uses
// (α_R, √γ_R). The Lax envelope wraps Δq = q_R - q_L, where q_packed
// is densitized with the cell's own √γ — so Δq is NOT just a state jump,
// it includes a metric variation. This may break the Lax convex-
// combination argument when √γ varies across the interface (the
// "densitization tax" — SESSION_NOTES_3 §16.6).
//
// This wrapper sets up a 3-cell conf_range, fills prods at three
// different sample points (centered at (sx,sy,sz) with x-offsets
// ±dx_metric), and computes fluctuations with set_interface_idx(0,1)
// and (1,2) so wave_lax_curved sees per-side metrics. If the same
// dumps that PASS in run_direct_state_sweep (shared metric) FAIL here
// (per-cell metric), the densitization tax is the production s²<0
// mechanism.
static void
run_per_cell_metric_sweep(struct gkyl_gr_spacetime *spacetime,
  const char *spacetime_label,
  enum gkyl_wv_gr_euler_tetrad_rp rp, const char *rp_label,
  enum gkyl_wv_flux_type flux_type,
  double sx, double sy, double sz, double dx_metric)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 2 };  // 3 cells: L=0, C=1, R=2
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  // Fill prods at three sample points along x: L at (sx-dx, sy, sz),
  // C at (sx, sy, sz), R at (sx+dx, sy, sz).
  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_C = gkyl_array_fetch(prods, 1);
  double *prods_R = gkyl_array_fetch(prods, 2);
  fill_prods_at(spacetime, sx - dx_metric, sy, sz, prods_L);
  fill_prods_at(spacetime, sx,              sy, sz, prods_C);
  fill_prods_at(spacetime, sx + dx_metric, sy, sz, prods_R);

  // Skip if any cell is excised — the test needs all three active.
  if (prods_L[GKYL_GR_SP_EXCISION] < 0.0
   || prods_C[GKYL_GR_SP_EXCISION] < 0.0
   || prods_R[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  double sqrt_det_L = sqrt(prods_L[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_C = sqrt(prods_C[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_R = sqrt(prods_R[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g_C[3][3] = {
    { prods_C[GKYL_GR_SP_INV_GIJ + 0], prods_C[GKYL_GR_SP_INV_GIJ + 1], prods_C[GKYL_GR_SP_INV_GIJ + 2] },
    { prods_C[GKYL_GR_SP_INV_GIJ + 3], prods_C[GKYL_GR_SP_INV_GIJ + 4], prods_C[GKYL_GR_SP_INV_GIJ + 5] },
    { prods_C[GKYL_GR_SP_INV_GIJ + 6], prods_C[GKYL_GR_SP_INV_GIJ + 7], prods_C[GKYL_GR_SP_INV_GIJ + 8] },
  };

  // Same dump cases as run_direct_state_sweep, copied for clarity.
  struct {
    const char *name;
    double D_C, Sx_C, Sy_C, Sz_C, tau_C;
    double D_L, Sx_L, Sy_L, Sz_L, tau_L;
    double D_R, Sx_R, Sy_R, Sz_R, tau_R;
  } cases[] = {
    { "full3d-s2-#10",
      2.528, 2.695, -0.005, 0.000, 0.655,
      3.578, 3.814, -1.539, 0.000, 0.987,
      0.878, 0.670,  0.193, 0.000, 0.026 },
    { "full3d-s2-#12",
      1.059, 0.935,  0.167, 0.000, 0.117,
      2.364, 3.174, -0.673, 0.000, 0.934,
      0.244, 0.182,  0.082, 0.000, 0.015 },
    { "full3d-s2-#17",
      1.695, 2.594, -0.757, 0.000, 0.986,
      1.194, 1.231, -0.184, 0.000, 0.303,
      2.087, 2.902,  0.585, 0.000, 1.045 },
    { "curvedLax-s2-#10",
      4.313, 4.221,  0.779, 0.000, 1.084,
      4.357, 4.743, -1.245, 0.000, 1.187,
      1.778, 1.109,  0.346, 0.000, 0.067 },
  };

  fprintf(stderr,
    "[per-cell-metric %s/%s dx_metric=%.4f]  case        τ_C_new (sign+mag)\n",
    spacetime_label, rp_label, dx_metric);
  for (size_t c = 0; c < sizeof(cases)/sizeof(*cases); c++) {
    // qin in GLOBAL frame, packed with EACH cell's own √γ. The dump qin
    // values are already densitized in production using the source cell's
    // √γ, so we apply the same per-cell densitization here.
    double qL_glob[5] = { cases[c].D_L, cases[c].Sx_L, cases[c].Sy_L, cases[c].Sz_L, cases[c].tau_L };
    double qC_glob[5] = { cases[c].D_C, cases[c].Sx_C, cases[c].Sy_C, cases[c].Sz_C, cases[c].tau_C };
    double qR_glob[5] = { cases[c].D_R, cases[c].Sx_R, cases[c].Sy_R, cases[c].Sz_R, cases[c].tau_R };

    int idx_L[1] = { 0 }, idx_C[1] = { 1 }, idx_R[1] = { 2 };
    double norm[3] = { 1.0, 0.0, 0.0 };
    double tau1v[3] = { 0.0, 1.0, 0.0 };
    double tau2v[3] = { 0.0, 0.0, 1.0 };

    // L-C interface: prodl from cell L, prodr from cell C
    eqn->set_interface_idx_func(eqn, idx_L, idx_C);
    double qL_loc[5], qC_loc_LC[5];
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qL_glob, qL_loc);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qC_glob, qC_loc_LC);
    double delta_LC[5];
    for (int i = 0; i < 5; i++) delta_LC[i] = qC_loc_LC[i] - qL_loc[i];
    double waves_LC[3 * 5], speeds_LC[3];
    double maxs_LC = eqn->waves_func(eqn, flux_type,
      delta_LC, qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC);
    double amdq_LC[5], apdq_LC[5];
    eqn->qfluct_func(eqn, flux_type,
      qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC, amdq_LC, apdq_LC);

    // C-R interface: prodl from cell C, prodr from cell R
    eqn->set_interface_idx_func(eqn, idx_C, idx_R);
    double qC_loc_CR[5], qR_loc[5];
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qC_glob, qC_loc_CR);
    eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qR_glob, qR_loc);
    double delta_CR[5];
    for (int i = 0; i < 5; i++) delta_CR[i] = qR_loc[i] - qC_loc_CR[i];
    double waves_CR[3 * 5], speeds_CR[3];
    double maxs_CR = eqn->waves_func(eqn, flux_type,
      delta_CR, qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR);
    double amdq_CR[5], apdq_CR[5];
    eqn->qfluct_func(eqn, flux_type,
      qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR, amdq_CR, apdq_CR);

    double maxs = fmax(maxs_LC, maxs_CR);
    if (!(maxs > 0.0)) { fprintf(stderr, "  [%s] degenerate\n", cases[c].name); continue; }
    double dt_dx = 0.5 / maxs;

    double qC_new[5];
    for (int i = 0; i < 5; i++) {
      qC_new[i] = qC_loc_CR[i] - dt_dx * (amdq_CR[i] + apdq_LC[i]);
    }
    // Undensitize qC_new with CENTER cell's √γ (POSITIVITY_SWEEP checks
    // invariants in the cell's own undensitized coords).
    double qC_und[5];
    for (int i = 0; i < 5; i++) qC_und[i] = qC_new[i] / sqrt_det_C;
    double tau_C_new = qC_und[4];
    double mom_sq_new = gkyl_gr_euler_mom_sq(inv_g_C, qC_und[1], qC_und[2], qC_und[3]);
    double s_sq_new   = (qC_und[0] + qC_und[4]) * (qC_und[0] + qC_und[4]) - mom_sq_new;

    fprintf(stderr,
      "  [%s/%-16s] √γ_L=%.4f √γ_C=%.4f √γ_R=%.4f   τ_C_new = %+.3e %s   s²_new = %+.3e %s\n",
      rp_label, cases[c].name, sqrt_det_L, sqrt_det_C, sqrt_det_R,
      tau_C_new,
      tau_C_new < 0.0 ? "**τ<0**" : "      ",
      s_sq_new,
      s_sq_new <= 0.0 ? "**s²≤0**" : "       ");
  }
  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

// ---------------------------------------------------------------------------
// Faithful production reproducer: takes verbatim wave_prop dump info
// (per-cell coords + per-cell q values + dir) and reruns the curved-Lax
// POSITIVITY_SWEEP step on cell C using the EXACT same per-cell metrics
// that production used. Reports the resulting (qC_new) admissibility.
//
// The wave_prop dump format (after the §17 instrumentation upgrade) gives
// us xcL, xcC, xcR and the dir of the sweep — everything we need to call
// fill_prods_at at the exact production cell-centers and run the same
// flux/fluctuation pipeline. If THIS reproducer fails to produce s²<0,
// the production fire is NOT coming from wave_lax_curved on the dumped
// qin inputs — it must be from a tier or transformation we haven't
// instrumented yet (MUSCL HIGH inadmissible reconstruction, sweep
// directional interaction, etc.).
//
// dir=0 → x-sweep (norm = +x̂)
// dir=1 → y-sweep (norm = +ŷ, tau1 = -x̂, tau2 = +ẑ) — matches wave_geom
//        convention for axis-aligned grids.
struct dump_repro {
  const char *name;
  int dir;
  double xL, yL, zL;
  double xC, yC, zC;
  double xR, yR, zR;
  double qC[5];
  double qL[5];
  double qR[5];
  // Production-dumped dt/dx — when nonzero, the reproducer uses this
  // exact dt_dx instead of 0.5/maxs_local. Lets us match production
  // exactly for cells where the global-vs-local CFL difference matters.
  double dt;   // 0 means "use 0.5/maxs"
  double dx;
};

static void
run_dump_reproducer(struct gkyl_gr_spacetime *spacetime,
  const char *spacetime_label,
  enum gkyl_wv_gr_euler_tetrad_rp rp, const char *rp_label,
  enum gkyl_wv_flux_type flux_type,
  const struct dump_repro *d)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 2 };  // 3 cells: L=0, C=1, R=2
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_mod_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_mod_auxfields){ .prods = prods });

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_C = gkyl_array_fetch(prods, 1);
  double *prods_R = gkyl_array_fetch(prods, 2);
  fill_prods_at(spacetime, d->xL, d->yL, d->zL, prods_L);
  fill_prods_at(spacetime, d->xC, d->yC, d->zC, prods_C);
  fill_prods_at(spacetime, d->xR, d->yR, d->zR, prods_R);

  double sqrt_det_C = sqrt(prods_C[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g_C[3][3] = {
    { prods_C[GKYL_GR_SP_INV_GIJ + 0], prods_C[GKYL_GR_SP_INV_GIJ + 1], prods_C[GKYL_GR_SP_INV_GIJ + 2] },
    { prods_C[GKYL_GR_SP_INV_GIJ + 3], prods_C[GKYL_GR_SP_INV_GIJ + 4], prods_C[GKYL_GR_SP_INV_GIJ + 5] },
    { prods_C[GKYL_GR_SP_INV_GIJ + 6], prods_C[GKYL_GR_SP_INV_GIJ + 7], prods_C[GKYL_GR_SP_INV_GIJ + 8] },
  };

  // Rotation triad keyed off dir. MUST match wave_geom (gkyl_wave_geom_priv.h
  // calc_geom_2d_from_nodes) for production-faithful reproduction:
  //   dir=0 (x-normal face): norm=+x̂, tau1=+ŷ, tau2=+ẑ
  //   dir=1 (y-normal face): norm=+ŷ, tau1=+x̂, tau2=-ẑ  ← cross-product
  //     ordering picks tau2=-ẑ so that norm = tau1 × tau2 = ŷ.
  // (The basis only changes the local-frame components, not the
  // physical answer, but matching production lets us cross-check the
  // rotation pipeline as well.)
  double norm[3], tau1v[3], tau2v[3];
  switch (d->dir) {
    case 0:
      norm[0]=1; norm[1]=0; norm[2]=0;
      tau1v[0]=0; tau1v[1]=1; tau1v[2]=0;
      tau2v[0]=0; tau2v[1]=0; tau2v[2]=1;
      break;
    case 1:
      norm[0]=0; norm[1]=1; norm[2]=0;
      tau1v[0]=1; tau1v[1]=0; tau1v[2]=0;
      tau2v[0]=0; tau2v[1]=0; tau2v[2]=-1;
      break;
    default:
      norm[0]=0; norm[1]=0; norm[2]=1;
      tau1v[0]=1; tau1v[1]=0; tau1v[2]=0;
      tau2v[0]=0; tau2v[1]=1; tau2v[2]=0;
      break;
  }

  double qL_glob[5], qC_glob[5], qR_glob[5];
  for (int i = 0; i < 5; i++) {
    qL_glob[i] = d->qL[i];
    qC_glob[i] = d->qC[i];
    qR_glob[i] = d->qR[i];
  }

  int idx_L[1] = { 0 }, idx_C[1] = { 1 }, idx_R[1] = { 2 };

  // L-C interface: prodl from cell L, prodr from cell C
  eqn->set_interface_idx_func(eqn, idx_L, idx_C);
  double qL_loc[5], qC_loc_LC[5];
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qL_glob, qL_loc);
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qC_glob, qC_loc_LC);
  double delta_LC[5];
  for (int i = 0; i < 5; i++) delta_LC[i] = qC_loc_LC[i] - qL_loc[i];
  double waves_LC[3 * 5], speeds_LC[3];
  double maxs_LC = eqn->waves_func(eqn, flux_type,
    delta_LC, qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC);
  double amdq_LC[5], apdq_LC[5];
  eqn->qfluct_func(eqn, flux_type,
    qL_loc, qC_loc_LC, 1.0, 1.0, waves_LC, speeds_LC, amdq_LC, apdq_LC);

  // C-R interface: prodl from cell C, prodr from cell R
  eqn->set_interface_idx_func(eqn, idx_C, idx_R);
  double qC_loc_CR[5], qR_loc[5];
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qC_glob, qC_loc_CR);
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, qR_glob, qR_loc);
  double delta_CR[5];
  for (int i = 0; i < 5; i++) delta_CR[i] = qR_loc[i] - qC_loc_CR[i];
  double waves_CR[3 * 5], speeds_CR[3];
  double maxs_CR = eqn->waves_func(eqn, flux_type,
    delta_CR, qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR);
  double amdq_CR[5], apdq_CR[5];
  eqn->qfluct_func(eqn, flux_type,
    qC_loc_CR, qR_loc, 1.0, 1.0, waves_CR, speeds_CR, amdq_CR, apdq_CR);

  double maxs = fmax(maxs_LC, maxs_CR);
  if (!(maxs > 0.0)) {
    fprintf(stderr, "  [%s] degenerate (maxs<=0)\n", d->name);
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }
  // If the dump captured production's exact dt/dx, use it. Otherwise
  // fall back to local CFL=0.5.
  double dt_dx = (d->dt > 0.0 && d->dx > 0.0) ? d->dt / d->dx : 0.5 / maxs;

  double qC_new[5];
  for (int i = 0; i < 5; i++) {
    qC_new[i] = qC_loc_CR[i] - dt_dx * (amdq_CR[i] + apdq_LC[i]);
  }
  // qC_new is in the LOCAL-rotated frame. The wave_prop check_inv test
  // also looks at qt in the local frame (per the dump comment), so we
  // check admissibility on qC_new directly without rotating back.
  double qC_und[5];
  for (int i = 0; i < 5; i++) qC_und[i] = qC_new[i] / sqrt_det_C;

  // For per-cell metric in rotated frame: when norm=ŷ (dir=1), the
  // rotated inv_g has γ^{yy} in the [0][0] slot. The cleanest way is
  // to check in the ORIGINAL (unrotated) frame: rotate qC_new back to
  // global, divide by √γ_C, check with the unrotated inv_g.
  double qC_glob_new[5];
  eqn->rotate_to_global_func(eqn, tau1v, tau2v, norm, qC_new, qC_glob_new);
  double qC_glob_und[5];
  for (int i = 0; i < 5; i++) qC_glob_und[i] = qC_glob_new[i] / sqrt_det_C;
  double tau_new = qC_glob_und[4];
  double mom_sq_new = gkyl_gr_euler_mom_sq(inv_g_C,
    qC_glob_und[1], qC_glob_und[2], qC_glob_und[3]);
  double s_sq_new = (qC_glob_und[0] + qC_glob_und[4]) * (qC_glob_und[0] + qC_glob_und[4]) - mom_sq_new;

  fprintf(stderr,
    "  [%s] %s/%s dir=%d  "
    "xC=(%.3f,%.3f) dt/dx=%.4f √γ_C=%.4f   "
    "qC_new_local=[%.6e %.6e %.6e %.6e %.6e]   "
    "qC_new_glob=[%.6e %.6e %.6e %.6e %.6e]   "
    "s²_new = %+.3e %s\n",
    d->name, spacetime_label, rp_label, d->dir,
    d->xC, d->yC, dt_dx, sqrt_det_C,
    qC_new[0], qC_new[1], qC_new[2], qC_new[3], qC_new[4],
    qC_glob_new[0], qC_glob_new[1], qC_glob_new[2], qC_glob_new[3], qC_glob_new[4],
    s_sq_new,
    s_sq_new <= 0.0 ? "**s²≤0**" : "       ");
  (void)qC_und;

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

// Production dumps captured by the §17-instrumented wave_prop POSITIVITY
// _SWEEP. Each entry is a verbatim reproducer: feeding these into
// run_dump_reproducer should produce s²<0 in cell C if curved Lax is
// genuinely the source of the production failure. If they all PASS in
// isolation, the production fire is caused by something other than
// wave_lax_curved on these inputs.
//
// Entries TBD — fill from the bhl_dump_with_coords.log after the run
// completes.
static const struct dump_repro g_production_dumps[] = {
  // Verbatim dumps from BHL run with §17 instrumented wave_prop (M=0.3,
  // r_h=0.6, BH at (2.5,2.5), dx=0.01953). Coords are in the GLOBAL
  // grid frame; run_dump_reproducer shifts by (-2.5,-2.5,0) to feed
  // into BH-relative spacetime API.
  // dir=0 x-sweep s2-only fires:
  { "repair#10-s2", 0,
    /*L*/2.412109, 1.865234, 0.0,
    /*C*/2.431641, 1.865234, 0.0,
    /*R*/2.451172, 1.865234, 0.0,
    /*qC*/ { 2.528, 2.695, -0.005315, 0.0, 0.6550 },
    /*qL*/ { 3.578, 3.814, -1.539,    0.0, 0.9866 },
    /*qR*/ { 0.8776, 0.6700, 0.1928,  0.0, 0.02617 } },
  { "repair#13-s2", 0,
    2.529297, 1.865234, 0.0,
    2.548828, 1.865234, 0.0,
    2.568359, 1.865234, 0.0,
    { 0.1612, 0.2580,  0.06944, 0.0, 0.1378 },
    { 0.3345, 0.6247,  0.06834, 0.0, 0.3298 },
    { 0.03990, 0.04751, 0.04237, 0.0, 0.03538 } },
  { "repair#16-s2", 0,
    2.548828, 3.173828, 0.0,
    2.568359, 3.173828, 0.0,
    2.587891, 3.173828, 0.0,
    { 0.9301, 1.117, -0.2134, 0.0, 0.3867 },
    { 1.466,  1.809, -0.09216, 0.0, 0.4412 },
    { 0.3511, 0.4189, -0.04190, 0.0, 0.1490 } },
  // dir=1 y-sweep s2-only fires:
  { "repair#17-s2", 1,
    /*L*/2.529297, 1.767578, 0.0,
    /*C*/2.529297, 1.787109, 0.0,
    /*R*/2.529297, 1.806641, 0.0,
    /*qC*/ { 1.695, 2.594, -0.7574, 0.0, 0.9862 },
    /*qL*/ { 1.194, 1.231, -0.1838, 0.0, 0.3025 },
    /*qR*/ { 2.087, 2.902,  0.5848, 0.0, 1.045 } },
  { "repair#19-s2", 1,
    2.548828, 1.767578, 0.0,
    2.548828, 1.787109, 0.0,
    2.548828, 1.806641, 0.0,
    { 0.7084, 1.316, -0.2676, 0.0, 0.6415 },
    { 0.2359, 0.4153, -0.02211, 0.0, 0.2145 },
    { 1.204,  1.855,  0.01382, 0.0, 0.7738 } },
  // dir=1 tau+s2 fires:
  { "repair#21-taus2", 1,
    2.626953, 1.689453, 0.0,
    2.626953, 1.708984, 0.0,
    2.626953, 1.728516, 0.0,
    { 0.08045, 0.06253, -0.02530, 0.0, 0.003853 },
    { 0.02293, 0.02563,  0.03113, 0.0, 0.02994 },
    { 0.3296, 0.3375,   -0.1184,  0.0, 0.05345 } },
  // === Residual dumps from margin=1e-6 + 1x amax (Stage A winner). ===
  // 12 wave_prop s²<0 remain. Of the 10 dumped (5 mirror pairs), only
  // these distinct physical events are added; mirrors omitted.
  // m6-#10 (dir=1): ratio_cart=0.96 — near-boundary case, margin fix
  //   couldn't fully protect.
  { "m6-repair#10-s2", 1,
    2.509766, 1.806641, 0.0,
    2.509766, 1.826172, 0.0,
    2.509766, 1.845703, 0.0,
    { 0.9599, 1.999, -0.4613, 0.0, 1.128 },
    { 0.5331, 0.4811,  0.1294, 0.0, 0.09473 },
    { 1.802,  2.330,   0.5787, 0.0, 0.7766 } },
  // m6-#12 (dir=1): ratio_cart=1.01 — beyond Cartesian, admissible only
  //   in curved metric. Different physical event from #10 (later step).
  { "m6-repair#12-s2", 1,
    2.509766, 1.806641, 0.0,
    2.509766, 1.826172, 0.0,
    2.509766, 1.845703, 0.0,
    { 1.158, 2.045, -0.5523, 0.0, 0.9492 },
    { 0.6707, 0.6638, -0.02188, 0.0, 0.1875 },
    { 1.649, 2.391,   0.6050,  0.0, 0.9514 } },
  // m6-#14 (dir=0): high-precision dump. Production qt: D+τ=0.5444,
  //   |S|²_cart=0.2974 → s²_cart=-0.0010 (real failure of magnitude 1e-3).
  //   Input qin_C is well inside Cartesian boundary (ratio_cart=0.77)
  //   but the sharp L/C gradient (qin_L much larger) drives the Lax
  //   update toward the boundary in the curved frame.
  { "m6-repair#14-s2", 0,
    2.5488281250, 1.8457031250, 0.0,
    2.5683593750, 1.8457031250, 0.0,
    2.5878906250, 1.8457031250, 0.0,
    { 2.210722125432e-01, 3.147104175096e-01, 3.089140027060e-02, 0.0, 1.421972464263e-01 },
    { 5.639505503463e-01, 1.143098909165e+00, -5.529981844672e-02, 0.0, 5.887044319329e-01 },
    { 6.785050901076e-02, 9.288799320179e-02, 3.444705579548e-02, 0.0, 6.077562973962e-02 },
    6.403732456302e-03, 1.953125000000e-02 },
  // m6-#15 (dir=0): high-precision dump. Production qt: D+τ=0.1693,
  //   |S|²_cart=0.0290 → s²_cart=-0.0004. Very low absolute scale —
  //   D=0.067, τ=0.061.
  { "m6-repair#15-s2", 0,
    2.5878906250, 1.8652343750, 0.0,
    2.6074218750, 1.8652343750, 0.0,
    2.6269531250, 1.8652343750, 0.0,
    { 6.747470649607e-02, 9.495031964782e-02, 5.204785154879e-02, 0.0, 6.099348705102e-02 },
    { 2.262458644598e-01, 2.933156770815e-01, -2.095106241374e-03, 0.0, 8.806511095333e-02 },
    { 3.391654654402e-02, 4.068520153585e-02, 4.308902146394e-02, 0.0, 3.676802280504e-02 },
    6.403732456302e-03, 1.953125000000e-02 },
  // m6-#18 (dir=1): ratio_cart=1.02 — beyond Cartesian.
  { "m6-repair#18-s2", 1,
    2.568359, 1.826172, 0.0,
    2.568359, 1.845703, 0.0,
    2.568359, 1.865234, 0.0,
    { 0.2595, 0.5448, 0.01668, 0.0, 0.2849 },
    { 0.1358, 0.3683, -0.01390, 0.0, 0.2644 },
    { 0.4606, 0.7802, 0.1149,  0.0, 0.3782 } },
};

void test_direct_state_lax_curved_production_reproducer(void)
{
  // Match BHL production: M=0.3, spin=0, BH center at (2.5, 2.5).
  // The dumps log cell-center coords in the GLOBAL grid frame (origin
  // at the grid lower-left), so we must shift to BH-centered coords
  // before calling fill_prods_at.
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  const double bh_pos_x = 2.5, bh_pos_y = 2.5, bh_pos_z = 0.0;

  fprintf(stderr, "\n[production reproducer — LAX (LOW): what wave_lax_curved would produce]\n");
  size_t n_dumps = sizeof(g_production_dumps) / sizeof(*g_production_dumps);
  for (size_t i = 0; i < n_dumps; i++) {
    struct dump_repro d = g_production_dumps[i];
    if (strcmp(d.name, "PLACEHOLDER") == 0) continue;
    d.xL -= bh_pos_x; d.yL -= bh_pos_y; d.zL -= bh_pos_z;
    d.xC -= bh_pos_x; d.yC -= bh_pos_y; d.zC -= bh_pos_z;
    d.xR -= bh_pos_x; d.yR -= bh_pos_y; d.zR -= bh_pos_z;
    run_dump_reproducer(st, "Schw(M=0.3)", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
      GKYL_WV_LOW_ORDER_FLUX, &d);
  }

  // Production s2-only dumps with redo_fluct=0/0 actually use HIGH HLL
  // (first-order, no MUSCL) — POSITIVITY_SWEEP wholesale-resets qout=qin
  // and the cell update uses stored FIRST_SWEEP amdq/apdq which are
  // HIGH HLL first-order. The HIGH MUSCL second-order correction that
  // kept the cell admissible in FIRST_SWEEP is gone. So to reproduce
  // production qt for these cells, call with rp=HLL and flux=HIGH_ORDER.
  fprintf(stderr, "\n[production reproducer — HLL HIGH first-order (no MUSCL):\n"
                  "  reproduces production qt for cells with redo_fluct=0/0,\n"
                  "  i.e., cells admissible after FIRST_SWEEP that POSITIVITY_SWEEP\n"
                  "  re-updated with stored HIGH HLL fluxes without MUSCL]\n");
  for (size_t i = 0; i < n_dumps; i++) {
    struct dump_repro d = g_production_dumps[i];
    if (strcmp(d.name, "PLACEHOLDER") == 0) continue;
    d.xL -= bh_pos_x; d.yL -= bh_pos_y; d.zL -= bh_pos_z;
    d.xC -= bh_pos_x; d.yC -= bh_pos_y; d.zC -= bh_pos_z;
    d.xR -= bh_pos_x; d.yR -= bh_pos_y; d.zR -= bh_pos_z;
    run_dump_reproducer(st, "Schw(M=0.3)", WV_GR_EULER_TETRAD_RP_HLL, "HLL-HIGH",
      GKYL_WV_HIGH_ORDER_FLUX, &d);
  }
  gkyl_gr_spacetime_release(st);
}

void test_direct_state_lax_curved_per_cell_metric(void)
{
  // Match BHL production: M=0.3, r_h=0.6, dx_production = 5/256 ≈ 0.01953.
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  // Test at a few sample points with PRODUCTION-MATCHING dx_metric:
  // (1) bow-shock front (-0.05, -0.6, 0): √γ varies sharply across an
  //     x-interface here.
  fprintf(stderr, "\n[per-cell-metric Schw/Lax-CURVED at (-0.05,-0.6,0), dx_metric=0.0195]\n");
  run_per_cell_metric_sweep(st, "Schw-bow", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, -0.05, -0.6, 0.0, 0.01953);
  // (2) Just outside r_h=0.6 on the x-axis
  fprintf(stderr, "\n[per-cell-metric Schw/Lax-CURVED at (0.65,0,0), dx_metric=0.0195]\n");
  run_per_cell_metric_sweep(st, "Schw-r0.65", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, 0.65, 0.0, 0.0, 0.01953);
  // (3) Off-axis (1, 1, 0) — already failed with shared metric before §17
  //     fix; check it now passes with per-cell metric too.
  fprintf(stderr, "\n[per-cell-metric Schw/Lax-CURVED at (1,1,0), dx_metric=0.0195]\n");
  run_per_cell_metric_sweep(st, "Schw-(1,1,0)", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, 1.0, 1.0, 0.0, 0.01953);
  // (4) Exaggerated dx_metric: 10x production to amplify the metric jump
  //     between L and R. If densitization tax is real, this should fail
  //     more dramatically.
  fprintf(stderr, "\n[per-cell-metric Schw/Lax-CURVED at (-0.05,-0.6,0), dx_metric=0.2 (10x exaggerated)]\n");
  run_per_cell_metric_sweep(st, "Schw-bow-EXAGG", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, -0.05, -0.6, 0.0, 0.2);
  gkyl_gr_spacetime_release(st);
}

void test_direct_state_lax_curved_near_horizon(void)
{
  // Match BHL production input (rt_gr_bhl_static_tetrad_mod.lua): M=0.3,
  // spin=0 → r_h = 2M = 0.6 in Schwarzschild Kerr-Schild coords. Production
  // setup is centered at (pos_x=2.5, pos_y=2.5); the spacetime API uses
  // BH-relative coords so points (sx, sy, sz) here correspond to global
  // cell positions (2.5+sx, 2.5+sy, sz).
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  // (1) On the x-axis just outside the horizon: R=0.65, sy=sz=0. Mostly β^x,
  //     diagonal γ. Tests "pure radial outflow at horizon front".
  fprintf(stderr, "\n[direct-state Schw/Lax-CURVED at (sx,sy,sz)=(0.65, 0, 0)  R=0.65]\n");
  run_direct_state_sweep(st, "Schw-(0.65,0,0)", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, 0.65, 0.0, 0.0);
  // (2) Production bow-shock front: x-sweep dump cell i=125 has x≈2.45,
  //     so sx≈-0.05. The bow-shock sits at sqrt(sx²+sy²) ≈ r_h=0.6, which
  //     puts sy ≈ ±0.6. γ_xy is nonzero here, and β has both x- and y-
  //     components — a much harsher test than (0.65,0,0).
  fprintf(stderr, "\n[direct-state Schw/Lax-CURVED at (sx,sy,sz)=(-0.05, -0.6, 0)  bow-shock front]\n");
  run_direct_state_sweep(st, "Schw-bow(-,-)", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, -0.05, -0.6, 0.0);
  fprintf(stderr, "\n[direct-state Schw/Lax-CURVED at (sx,sy,sz)=(-0.05, +0.6, 0)  bow-shock front (mirror)]\n");
  run_direct_state_sweep(st, "Schw-bow(-,+)", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, -0.05, 0.6, 0.0);
  // (3) Closer to horizon (just outside r=0.6 with both sx and sy contributions)
  fprintf(stderr, "\n[direct-state Schw/Lax-CURVED at (sx,sy,sz)=(-0.05, -0.605, 0)  hugging horizon]\n");
  run_direct_state_sweep(st, "Schw-hug(-,-)", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, -0.05, -0.605, 0.0);
  // (4) Control: far from horizon
  fprintf(stderr, "\n[direct-state Schw/Lax-CURVED at (sx,sy,sz)=(1.0, 1.0, 0)  control (far)]\n");
  run_direct_state_sweep(st, "Schw-(1,1,0)", WV_GR_EULER_TETRAD_RP_LAX, "Lax-CURVED",
    GKYL_WV_LOW_ORDER_FLUX, 1.0, 1.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_direct_state_hll(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_direct_state_sweep(st, "Schw", WV_GR_EULER_TETRAD_RP_HLL, "HLL-LOW",
    GKYL_WV_LOW_ORDER_FLUX, 0.3, 0.2, 0.0);
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
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
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
  build_state_convA(eos, rho_L, v_L, p_L, prods, qL);
  build_state_convA(eos, rho_M, v_M, p_M, prods, qM);
  build_state_convA(eos, rho_R, v_R, p_R, prods, qR);

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
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
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

void test_three_cell_hllc(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_three_cell_sweep(st, WV_GR_EULER_TETRAD_RP_HLLC, "HLLC");
  gkyl_gr_spacetime_release(st);
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
// EOS-accuracy validation vs the true Synge gas.
//
// The Synge gas (single-component perfect relativistic gas) has the
// closed-form specific enthalpy
//   h_Synge(θ) = K3(1/θ) / K2(1/θ),     θ = p/ρ.
// where K_n is the modified Bessel function of the second kind. Both
// MATHEWS_TAUB (TM) and RYU_CHATTOPADHYAY (RCC) closures are designed
// to approximate this curve in the trans-relativistic regime
// θ ∈ [O(0.1), O(10)] with stated maximum relative errors:
//   TM (Mignone+ 2005)         ~ 2%   on h
//   RCC (Ryu+ 2006)            ~ 0.8% on h
//
// This test sweeps θ ∈ [0.1, 10], runs a primitive recovery roundtrip
// at each point under each EOS, then compares the EOS's enthalpy at
// the recovered primitives against the Synge value. We also assert
// the roundtrip recovery is accurate (the recovery should exactly
// invert the EOS-specific build_state_convA seeding to machine
// precision for the chosen EOS).
//
// Coverage:
//   TM no-flow (W = 1): tests the closed-form TM cubic recovery in
//                       the static-fluid M = 0 short-circuit.
//   RCC relativistic   (|v| = 0.9, W ≈ 2.29): tests the TM cubic
//                       warm-start → RC Newton refinement → paper-
//                       grounded (a)/(b)/(c) physicality check path.
// ---------------------------------------------------------------------------

// Modified Bessel functions K_0, K_1 via the Abramowitz & Stegun
// polynomial fits (9.8.1-9.8.8). Accurate to <2e-7 over all x > 0
// — far better than our test tolerances of 2% (TM) and 0.8% (RCC).
// libm provides these on Linux/glibc as XSI extensions but Apple's
// libSystem does NOT, so we implement them locally for portability.
static double bessel_i0(double x)
{
  if (x < 0.0) x = -x;
  if (x < 3.75) {
    double t = (x / 3.75) * (x / 3.75);
    return 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492
      + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
  }
  // x ≥ 3.75 — asymptotic form. Not needed for our K_0 path (we only
  // call I_0 from K_0's "small x" branch where x ≤ 2 < 3.75), but kept
  // for completeness.
  double t = 3.75 / x;
  return (exp(x) / sqrt(x)) * (0.39894228 + t * (0.01328592
    + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281
    + t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633
    + t * 0.00392377))))))));
}
static double bessel_i1(double x)
{
  double ax = (x < 0.0) ? -x : x;
  double res;
  if (ax < 3.75) {
    double t = (ax / 3.75) * (ax / 3.75);
    res = ax * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934
      + t * (0.02658733 + t * (0.00301532 + t * 0.00032411))))));
  } else {
    double t = 3.75 / ax;
    res = (exp(ax) / sqrt(ax)) * (0.39894228 + t * (-0.03988024
      + t * (-0.00362018 + t * (0.00163801 + t * (-0.01031555
      + t * (0.02282967 + t * (-0.02895312 + t * (0.01787654
      + t * (-0.00420059)))))))));
  }
  return (x < 0.0) ? -res : res;
}
static double bessel_k0(double x)
{
  if (x <= 2.0) {
    double t = (x / 2.0) * (x / 2.0);
    return (-log(x / 2.0) * bessel_i0(x)) + (-0.57721566
      + t * (0.42278420 + t * (0.23069756 + t * (0.03488590
      + t * (0.00262698 + t * (0.00010750 + t * 0.00000740))))));
  }
  double t = 2.0 / x;
  return (exp(-x) / sqrt(x)) * (1.25331414 + t * (-0.07832358
    + t * (0.02189568 + t * (-0.01062446 + t * (0.00587872
    + t * (-0.00251540 + t * 0.00053208))))));
}
static double bessel_k1(double x)
{
  if (x <= 2.0) {
    double t = (x / 2.0) * (x / 2.0);
    return (log(x / 2.0) * bessel_i1(x)) + (1.0 / x) * (1.0
      + t * (0.15443144 + t * (-0.67278579 + t * (-0.18156897
      + t * (-0.01919402 + t * (-0.00110404 + t * (-0.00004686))))))) ;
  }
  double t = 2.0 / x;
  return (exp(-x) / sqrt(x)) * (1.25331414 + t * (0.23498619
    + t * (-0.03655620 + t * (0.01504268 + t * (-0.00780353
    + t * (0.00325614 + t * (-0.00068245)))))));
}

// K_2, K_3 via forward recurrence K_{n+1}(x) = K_{n-1}(x) + (2n/x)·K_n(x).
// Stable in the forward direction for K_n (unlike I_n).
//
// Numerical range used here: θ ∈ [0.1, 10] ⇒ x = 1/θ ∈ [0.1, 10].
//   - At x = 10: K_n(10) is small (~1e-5) but well above underflow;
//     the ratio K3/K2 is well-conditioned.
//   - At x = 0.1: K_n(0.1) is large (K0 ~ 2.4, K3 ~ 8e3); still finite,
//     ratio K3/K2 → 4/x = 4θ (ultra-relativistic limit).
// Both extremes evaluate cleanly in double precision.
static double bessel_k2_arg(double x)
{
  return bessel_k0(x) + (2.0 / x) * bessel_k1(x);
}
static double bessel_k3_arg(double x)
{
  return bessel_k1(x) + (4.0 / x) * bessel_k2_arg(x);
}
static double synge_enthalpy(double theta)
{
  double x = 1.0 / theta;
  return bessel_k3_arg(x) / bessel_k2_arg(x);
}

// Sweep θ ∈ [0.1, 10] (16 logspaced points) and validate that:
//   (a) primitive recovery roundtrip preserves (ρ, v, p) to round_trip_tol
//   (b) the EOS-specific enthalpy at the recovered (ρ, p) matches the
//       true Synge gas enthalpy to within max_rel_err.
static void
run_eos_synge_validation(struct gkyl_gr_euler_eos eos,
  const double v_co[3], double round_trip_tol, double max_rel_err,
  const char *label)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  double prods[GKYL_GR_SP_NCOMP_BASE];
  fill_prods_at(st, 0.3, 0.0, 0.0, prods);

  // Sanity-check velocity is subluminal at this metric (Minkowski → just
  // need |v|² < 1).
  double v_lo[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v_lo[i] += prods[GKYL_GR_SP_GIJ + 3*i + j] * v_co[j];
  double vsq = 0.0;
  for (int i = 0; i < 3; i++) vsq += v_lo[i] * v_co[i];
  TEST_CHECK_( vsq < 1.0,
    "[%s] velocity is superluminal (γ_ij v^i v^j = %g)", label, vsq );

  // Logspaced θ scan over [0.1, 10] (16 points).
  static const double thetas[] = {
    0.10, 0.135, 0.18, 0.24, 0.32, 0.42, 0.56, 0.75,
    1.00, 1.33,  1.78, 2.37, 3.16, 4.22, 5.62, 10.0
  };
  const int n_theta = sizeof(thetas) / sizeof(*thetas);

  double max_rt_dr = 0.0, max_rt_dp = 0.0, max_rt_dv = 0.0;
  double max_h_err = 0.0;
  double max_h_err_theta = 0.0;

  for (int it = 0; it < n_theta; it++) {
    double theta = thetas[it];
    double rho = 1.0;
    double p = theta;  // → p/ρ = θ

    // Seed the conservative state using the EOS-specific enthalpy.
    double q[5];
    build_state_convA(eos, rho, v_co, p, prods, q);

    // Run the production recovery dispatch.
    double sd  = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);
    double D   = q[0] / sd;
    double Sx  = q[1] / sd;
    double Sy  = q[2] / sd;
    double Sz  = q[3] / sd;
    double tau = q[4] / sd;
    const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
    double inv_g[3][3] = {
      { ig[0], ig[1], ig[2] },
      { ig[3], ig[4], ig[5] },
      { ig[6], ig[7], ig[8] },
    };
    struct gkyl_gr_euler_prim prim;
    gkyl_gr_euler_recover_primitives(eos, D, Sx, Sy, Sz, tau, inv_g,
      NULL, &prim);

    // (a) Roundtrip: recovered primitives should match the inputs.
    double dr = fabs(prim.rho - rho);
    double dp = fabs(prim.p - p);
    double dv = 0.0;
    for (int i = 0; i < 3; i++) {
      double d = fabs(prim.v[i] - v_co[i]);
      if (d > dv) dv = d;
    }
    if (dr > max_rt_dr) max_rt_dr = dr;
    if (dp > max_rt_dp) max_rt_dp = dp;
    if (dv > max_rt_dv) max_rt_dv = dv;

    // (b) EOS enthalpy at recovered primitives vs Synge enthalpy at
    // the recovered θ. Using recovered values keeps the comparison
    // self-consistent if (a) has any FP slack.
    double theta_r = prim.p / prim.rho;
    double h_eos   = gkyl_gr_euler_eos_enthalpy(eos, prim.rho, prim.p);
    double h_synge = synge_enthalpy(theta_r);
    double rel_err = fabs(h_eos - h_synge) / h_synge;
    if (rel_err > max_h_err) {
      max_h_err = rel_err;
      max_h_err_theta = theta;
    }
  }

  TEST_CHECK_( max_rt_dr < round_trip_tol,
    "[%s] roundtrip ρ residual %.3e exceeds %.3e", label, max_rt_dr, round_trip_tol );
  TEST_CHECK_( max_rt_dp < round_trip_tol,
    "[%s] roundtrip p residual %.3e exceeds %.3e", label, max_rt_dp, round_trip_tol );
  TEST_CHECK_( max_rt_dv < round_trip_tol,
    "[%s] roundtrip v residual %.3e exceeds %.3e", label, max_rt_dv, round_trip_tol );
  TEST_CHECK_( max_h_err < max_rel_err,
    "[%s] h_EOS vs h_Synge max rel err %.3e exceeds %.3e (worst at θ=%.3g)",
    label, max_h_err, max_rel_err, max_h_err_theta );

  // Informational dump of the worst-case errors (visible only on fail
  // or via TEST_MSG always-shown channel for postmortem).
  TEST_MSG( "[%s] roundtrip residuals max: ρ=%.3e p=%.3e v=%.3e | "
            "max |h_EOS − h_Synge|/h_Synge = %.3e at θ=%.3g",
    label, max_rt_dr, max_rt_dp, max_rt_dv, max_h_err, max_h_err_theta );

  gkyl_gr_spacetime_release(st);
}

// TM (use_rcc=false) with no flow (v = 0, W = 1). Exercises the
// closed-form TM cubic recovery's M = 0 static-fluid short-circuit.
// Tolerance 2.1% rather than the literature 2%: a direct numerical
// sweep |h_TM − h_Synge|/h_Synge over θ ∈ [0.1, 10] peaks at ~2.011%
// near θ ≈ 0.52. Mignone+ 2005 quotes "≲ 2%" as a rounded-down
// description of this same peak; the extra 0.1% margin in the test
// captures the actual worst-case the bound is approximating.
void test_eos_tm_synge_no_flow(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = false };
  double v[3] = { 0.0, 0.0, 0.0 };
  run_eos_synge_validation(eos, v, 1.0e-12, 2.1e-2, "TM no-flow");
}

// RCC (use_rcc=true) with relativistic flow (|v|² = 0.81, W ≈ 2.29).
// Exercises the TM cubic warm-start → RC Newton refinement → paper-
// grounded physicality check path. Tolerance 0.8% per Ryu+ 2006.
void test_eos_rcc_synge_relativistic(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = true };
  // Components chosen so |v|² = (0.9)² with v isotropic across x/y/z.
  double v_each = 0.9 / sqrt(3.0);
  double v[3] = { v_each, v_each, v_each };
  run_eos_synge_validation(eos, v, 1.0e-10, 8.0e-3, "RCC relativistic");
}

// RCC at higher Lorentz factor (|v| = 0.99 → W ≈ 7.09). Stresses the
// pressure recovery at large W where the (E²−M²)/E² ~ 1/W² ≈ 0.02
// conditioning measure starts to bite — historically the regime where
// TM cubic precision-loss appeared. RC Newton on the squared
// polynomial is supposed to converge robustly here without the TM
// precision-loss feedback loop. Round-trip and enthalpy tolerances
// stay tight (RCC's Synge accuracy is an enthalpy-formula property,
// independent of W; precision-loss would show up first in the
// roundtrip ρ/p residual).
void test_eos_rcc_synge_high_W(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = true };
  double v_each = 0.99 / sqrt(3.0);
  double v[3] = { v_each, v_each, v_each };
  run_eos_synge_validation(eos, v, 1.0e-10, 8.0e-3, "RCC high-W");
}

// RCC at extreme Lorentz factor (|v| = 0.9999 → W ≈ 70.71). The
// (E²−M²)/E² ~ 1/W² ≈ 2e-4 conditioning measure is now squarely in
// the regime where the TM cubic's coefficients lose ~3 digits to
// cancellation. SESSION_NOTES_EOS_IMPROVEMENTS §2.2 catalogs cases
// where TM gave W = 5×10⁵ and RC corrected it to W = 50; this test
// ensures RC Newton's refinement + paper-grounded physicality check
// still locks onto the right root at extreme W.
//
// Round-trip tolerance loosened to 1e-6 to reflect realistic FP
// precision at extreme W:
//   - ρ residual scales as ULP × (D/ρ ratio) ~ 1e-12 × W ~ 7e-11
//     compounded with the W-solve precision; measured ≈ 7e-8.
//   - p residual is worse because p = ρhW² − τ − ρW is a catastrophic
//     cancellation: |ρhW²| ≈ |τ| ≈ ρ·4θ·W² ~ 2×10⁴ at the worst θ,
//     while p ≈ θ ~ 1, so p inherits relative precision degraded by
//     ρhW²/p ≈ 2×10⁴; measured ≈ 5e-7.
// What this test is really checking: RC Newton at extreme W does NOT
// diverge or fall back to TM with W = 5×10⁵ (a documented failure
// mode at this regime). Enthalpy accuracy stays at 0.8% — that's an
// algebraic property of the recovered θ and unaffected by W.
void test_eos_rcc_synge_extreme_W(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = true };
  double v_each = 0.9999 / sqrt(3.0);
  double v[3] = { v_each, v_each, v_each };
  run_eos_synge_validation(eos, v, 1.0e-6, 8.0e-3, "RCC extreme-W");
}

// ---------------------------------------------------------------------------

TEST_LIST = {
  { "banyuls_flux_consistency_lax_minkowski",      test_banyuls_flux_consistency_lax_minkowski },
  { "banyuls_flux_consistency_lax_schwarzschild",  test_banyuls_flux_consistency_lax_schwarzschild },
  { "banyuls_flux_consistency_lax_kerr",           test_banyuls_flux_consistency_lax_kerr },
  { "banyuls_flux_consistency_hllc_minkowski",     test_banyuls_flux_consistency_hllc_minkowski },
  { "banyuls_flux_consistency_hllc_schwarzschild", test_banyuls_flux_consistency_hllc_schwarzschild },
  { "banyuls_flux_consistency_hllc_kerr",          test_banyuls_flux_consistency_hllc_kerr },

  { "riemann_properties_lax_minkowski",     test_riemann_properties_lax_minkowski },
  { "riemann_properties_lax_schwarzschild", test_riemann_properties_lax_schwarzschild },
  { "riemann_properties_lax_kerr",          test_riemann_properties_lax_kerr },
  { "riemann_properties_hll_minkowski",     test_riemann_properties_hll_minkowski },
  { "riemann_properties_hll_schwarzschild", test_riemann_properties_hll_schwarzschild },
  { "riemann_properties_hll_kerr",          test_riemann_properties_hll_kerr },
  { "riemann_properties_hllc_minkowski",     test_riemann_properties_hllc_minkowski },
  { "riemann_properties_hllc_schwarzschild", test_riemann_properties_hllc_schwarzschild },
  { "riemann_properties_hllc_kerr",          test_riemann_properties_hllc_kerr },

  { "excision_absorbing_lax_minkowski",     test_excision_absorbing_lax_minkowski },
  { "excision_absorbing_lax_schwarzschild", test_excision_absorbing_lax_schwarzschild },
  { "excision_absorbing_hll_minkowski",     test_excision_absorbing_hll_minkowski },
  { "excision_absorbing_hll_schwarzschild", test_excision_absorbing_hll_schwarzschild },

  { "positivity_lax_minkowski",     test_positivity_lax_minkowski },
  { "positivity_lax_schwarzschild", test_positivity_lax_schwarzschild },
  { "positivity_lax_kerr",          test_positivity_lax_kerr },
  { "positivity_hll_minkowski",     test_positivity_hll_minkowski },
  { "positivity_hll_schwarzschild", test_positivity_hll_schwarzschild },
  { "positivity_hll_kerr",          test_positivity_hll_kerr },
  { "positivity_hllc_minkowski",     test_positivity_hllc_minkowski },
  { "positivity_hllc_schwarzschild", test_positivity_hllc_schwarzschild },
  { "positivity_hllc_kerr",          test_positivity_hllc_kerr },
  { "positivity_lax_bhl",            test_positivity_lax_bhl },
  { "positivity_hll_bhl",            test_positivity_hll_bhl },
  { "positivity_hllc_bhl",           test_positivity_hllc_bhl },

  { "near_floor_lax",  test_near_floor_lax },
  { "near_floor_hll",  test_near_floor_hll },
  { "near_floor_hllc", test_near_floor_hllc },

  { "round_trip_minkowski",     test_round_trip_minkowski },
  { "round_trip_schwarzschild", test_round_trip_schwarzschild },
  { "round_trip_kerr",          test_round_trip_kerr },

  { "prim_consistency_minkowski",     test_prim_consistency_minkowski },
  { "prim_consistency_schwarzschild", test_prim_consistency_schwarzschild },
  { "prim_consistency_kerr",          test_prim_consistency_kerr },

  { "bhl_regime_lax_schwarzschild",  test_bhl_regime_lax_schwarzschild },
  { "bhl_regime_hll_schwarzschild",  test_bhl_regime_hll_schwarzschild },
  { "bhl_regime_lax_kerr",           test_bhl_regime_lax_kerr },
  { "bhl_regime_hll_kerr",           test_bhl_regime_hll_kerr },
  { "bhl_regime_hllc_schwarzschild", test_bhl_regime_hllc_schwarzschild },
  { "bhl_regime_hllc_kerr",          test_bhl_regime_hllc_kerr },

  { "floor_precision_lax",  test_floor_precision_lax },
  { "floor_precision_hll",  test_floor_precision_hll },
  { "floor_precision_hllc", test_floor_precision_hllc },

  { "three_cell_lax",  test_three_cell_lax },
  { "three_cell_hll",  test_three_cell_hll },
  { "three_cell_hllc", test_three_cell_hllc },

  { "small_tau_over_D_lax",  test_small_tau_over_D_lax },
  { "small_tau_over_D_hll",  test_small_tau_over_D_hll },
  { "small_tau_over_D_hllc", test_small_tau_over_D_hllc },

  { "direct_state_lax",  test_direct_state_lax },
  { "direct_state_hll",  test_direct_state_hll },
  { "direct_state_hllc", test_direct_state_hllc },
  { "direct_state_hllc_high_order", test_direct_state_hllc_high_order },
  { "direct_state_hll_high_order",  test_direct_state_hll_high_order },
  { "direct_state_hll_high_order_minkowski",  test_direct_state_hll_high_order_minkowski },
  { "direct_state_lax_high_order_minkowski",  test_direct_state_lax_high_order_minkowski },
  { "direct_state_lax_high_order_schwarzschild",  test_direct_state_lax_high_order_schwarzschild },
  { "direct_state_lax_curved_schwarzschild",  test_direct_state_lax_curved_schwarzschild },
  { "direct_state_lax_curved_near_horizon",  test_direct_state_lax_curved_near_horizon },
  { "direct_state_lax_curved_per_cell_metric",  test_direct_state_lax_curved_per_cell_metric },
  { "direct_state_lax_curved_production_reproducer",  test_direct_state_lax_curved_production_reproducer },
  { "direct_state_hllc_fallback_probe", test_direct_state_hllc_fallback_probe },

  // EOS accuracy vs the true Synge gas (Bessel-function reference).
  { "eos_tm_synge_no_flow",          test_eos_tm_synge_no_flow },
  { "eos_rcc_synge_relativistic",    test_eos_rcc_synge_relativistic },
  { "eos_rcc_synge_high_W",          test_eos_rcc_synge_high_W },
  { "eos_rcc_synge_extreme_W",       test_eos_rcc_synge_extreme_W },

  { NULL, NULL },
};
