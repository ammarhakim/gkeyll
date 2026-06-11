// Positivity / invariant-domain tests for wv_gr_euler_tetrad (Convention A).
// Split out of ctest_wv_gr_euler_tetrad.c; see CTEST_SPLIT_PLAN.md.
//
//   1. positivity_curved_lax_* — LOW_ORDER sweep over the shared state
//      table (CFL = 0.9, hard D/s² asserts)
//   2. near_floor_*           — τ-floor-value sweep through the fallback
//   3. small_tau_over_D_*     — report-only τ/D ≪ 1 regime (HIGH_ORDER,
//      per-rp: the HLLC-vs-HLL τ-positivity question)
//   4. production reproducer  — flagged 3-cell debugging harness
//      (replicates wave_prop verbatim; the one deliberate 3-cell exception)
//   5. positivity_registry_*  — the per-(condition × flux) failure matrix,
//      extended with the coarse-dx and full margin axes that absorbed the
//      retired lax_convexity_tetrad_vs_curved_* and lax_margin_* tests
//      (measured history preserved in the registry header comment)

#include <acutest.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler_prim_priv.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>

#include "ctest_gr_euler_pos_states.h"
#include "ctest_gr_euler_test_util.h"

// ---------------------------------------------------------------------------
// 3. Positivity preservation sweep
// ---------------------------------------------------------------------------

// The Riemann-problem state table (struct gr_euler_pos_rp_case, in
// primitives with contravariant velocity) is shared with the mod-variant
// ctest via ctest_gr_euler_pos_states.h.

// Per-cell-metric positivity check. prods_L and prods_R are the
// spacetime products rows for the LEFT and RIGHT cells (production
// wave_prop pipes adjacent-cell metrics into prodl_local/prodr_local
// independently). Each cell's conservatives are densitized with its
// OWN √γ and admissibility-checked in its OWN inverse spatial metric,
// matching what wave_prop does post-POSITIVITY_SWEEP.
//
// `ftype` selects which Lax this probes: GKYL_WV_LOW_ORDER_FLUX routes
// to the direct curved-frame Lax fallback (wave_lax_curved);
// GKYL_WV_HIGH_ORDER_FLUX routes to the variant's tetrad-first solver
// (for rp_type=LAX, the Gorard+25 transform + SR Lax kernel).
//
// Returns the number of cell-updates actually performed (2, or 0 when
// the case is skipped because it is unbuildable / inadmissible at this
// pair of metric samples), so callers can report honest percentages.
static int
run_positivity_for_rp(struct gkyl_wv_eqn *eqn,
  struct wv_gr_euler_tetrad *grm,
  const double *prods_L, const double *prods_R,
  const struct gr_euler_pos_rp_case *rc,
  enum gkyl_wv_flux_type ftype, double cfl,
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

  // Silent skip when either side's velocity is superluminal in its own
  // cell metric. The shared state table deliberately carries near-luminal
  // flat-space entries (v = 0.85–0.95) that exceed c once contracted with
  // a near-horizon γ_ij — those combinations are unbuildable, not failed.
  {
    double vsq_L = 0.0, vsq_R = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        vsq_L += prods_L[GKYL_GR_SP_GIJ + 3*i + j] * rc->v_L[i] * rc->v_L[j];
        vsq_R += prods_R[GKYL_GR_SP_GIJ + 3*i + j] * rc->v_R[i] * rc->v_R[j];
      }
    if (!(vsq_L < 1.0 - 1.0e-6) || !(vsq_R < 1.0 - 1.0e-6)) return 0;
  }

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
    return 0;
  }
  record_admissibility(inv_g_R, qR_und, &d_ok, &s_ok, &t_ok);
  if (!(d_ok && s_ok && t_ok)) {
    fprintf(stderr, "  [positivity] WARNING: input qR inadmissible for case '%s' "
      "(D_ok=%d, S2_ok=%d, tau_ok=%d) — skipping\n", rc->name, d_ok, s_ok, t_ok);
    return 0;
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
  double maxs = eqn->waves_func(eqn, ftype,
    delta, qL_loc, qR_loc, 1.0, 1.0, waves, speeds);
  double amdq[5], apdq[5];
  eqn->qfluct_func(eqn, ftype,
    qL_loc, qR_loc, 1.0, 1.0, waves, speeds, amdq, apdq);

  // CFL-scaled update against the interface (callers pass 0.9 for the
  // production-like probe, 1.0 for the extremal convexity probe).
  if (!(maxs > 0.0)) return 0;
  double dt_dx = 1.0 / maxs;  // CFL = 1: the extremal convexity probe

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
  return 2;
}

// LOW_ORDER curved-Lax positivity sweep. NOTE: this sweep is
// rp-independent BY CONSTRUCTION — the LOW_ORDER branch of
// wave_tetrad_dispatch routes to wave_lax_curved / qfluct_lax_curved for
// every rp_type, without reading sr_kernel / num_waves / excision_policy.
// (Historical: this runner used to be invoked three times as
// positivity_{lax,hll,hllc}_* — three names for the identical code path;
// in particular the old "HLLC headline" comment promised a star-state
// τ-positivity probe this sweep never performed. The per-scheme
// HIGH_ORDER positivity story lives in lax_convexity_tetrad_vs_curved_*
// below; the HLLC-vs-HLL τ-positivity question lives in
// small_tau_over_D_*, which genuinely runs HIGH_ORDER.)
static void
run_positivity_sweep(struct gkyl_gr_spacetime *spacetime,
  struct gkyl_gr_euler_eos eos,
  const char *label)
{
  double gas_gamma = eos.gas_gamma;

  // 2-cell conf_range so set_interface_idx(0, 1) routes prods_L and
  // prods_R to wave_lax_curved's prodl_local / prodr_local.
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range,
    WV_GR_EULER_TETRAD_RP_LAX);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);

  int D_violations = 0, S2_violations = 0, tau_violations = 0;
  int total_cells = 0;

  int n_pts = sizeof(g_positivity_points) / sizeof(*g_positivity_points);

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

    for (int c = 0; c < GR_EULER_POS_NUM_CASES; c++) {
      total_cells += run_positivity_for_rp(eqn, grm, prods_L, prods_R,
        &gr_euler_pos_cases[c], GKYL_WV_LOW_ORDER_FLUX, 0.9,
        &D_violations, &S2_violations, &tau_violations);
    }
  }

  fprintf(stderr,
    "[positivity %s] across %d cell-updates (per-cell-metric, dx=%.4f): "
    "D<=0: %d (%.1f%%), s²<=0: %d (%.1f%%), τ<0: %d (%.1f%%)\n",
    label, total_cells, GR_EULER_POSITIVITY_DX,
    D_violations,   100.0 * D_violations   / total_cells,
    S2_violations,  100.0 * S2_violations  / total_cells,
    tau_violations, 100.0 * tau_violations / total_cells);

  // D and s² preservation: the convex-combination argument for curved
  // Lax assumes one flux function on both sides, which holds only up to
  // O(Δgeom) — exact at production dx (this sweep), degrading at coarse
  // near-horizon Δγ (measured in lax_convexity_tetrad_vs_curved_*).
  // Assert hard at production dx.
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
void test_positivity_curved_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-Mink");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_curved_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-Schw");
  gkyl_gr_spacetime_release(st);
}
void test_positivity_curved_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-Kerr");
  gkyl_gr_spacetime_release(st);
}

// BHL-production variant: M=0.3 (r_h=0.6) Schwarzschild matches the
// rt_gr_bhl_static_tetrad.lua input. The sample-point list in
// g_positivity_points contains points just outside r_h=0.6 and on the
// bow-shock front, which only become physical under the M=0.3 metric.
// (Under M=0.1 those same points are deeper in the strong-field region
// but the excision check filters anything inside r_h=0.2.)
void test_positivity_curved_lax_bhl(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-BHL(M=0.3)");
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
  const double *prods_L, const double *prods_R,
  double tau_min_target,
  int *D_v, int *S2_v, int *tau_v, int *total)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);

  double sqrt_det_L = sqrt(prods_L[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_R = sqrt(prods_R[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g_L[3][3], inv_g_R[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      inv_g_L[i][j] = prods_L[GKYL_GR_SP_INV_GIJ + 3*i + j];
      inv_g_R[i][j] = prods_R[GKYL_GR_SP_INV_GIJ + 3*i + j];
    }

  int n = sizeof(g_near_floor_cases) / sizeof(*g_near_floor_cases);
  for (int c = 0; c < n; c++) {
    const struct near_floor_rp *rc = &g_near_floor_cases[c];

    // Construct the "floored-side" state: target τ ≈ tau_min_target.
    // Cold-flow inversion (W≈1, h≈1) gives τ ≈ p/(γ-1), so set p_min
    // accordingly.
    double p_min = (gas_gamma - 1.0) * tau_min_target;

    // Each side is built with its OWN cell metric (per-cell-metrics
    // policy — production wave_prop reads adjacent cells' prods).
    double qL[5], qR[5];
    if (rc->floor_on_left) {
      build_state_convA(eos, rc->rho_min, rc->v_min, p_min, prods_L, qL);
      build_state_convA(eos, rc->rho_typ, rc->v_typ, rc->p_typ, prods_R, qR);
    } else {
      build_state_convA(eos, rc->rho_typ, rc->v_typ, rc->p_typ, prods_L, qL);
      build_state_convA(eos, rc->rho_min, rc->v_min, p_min,    prods_R, qR);
    }

    // Skip if either side fails admissibility going in (defensive).
    {
      bool d_ok, s_ok, t_ok;
      double qu[5];
      for (int i = 0; i < 5; i++) qu[i] = qL[i] / sqrt_det_L;
      record_admissibility(inv_g_L, qu, &d_ok, &s_ok, &t_ok);
      if (!(d_ok && s_ok && t_ok)) continue;
      for (int i = 0; i < 5; i++) qu[i] = qR[i] / sqrt_det_R;
      record_admissibility(inv_g_R, qu, &d_ok, &s_ok, &t_ok);
      if (!(d_ok && s_ok && t_ok)) continue;
    }

    int idxl[1] = { 0 }, idxr[1] = { 1 };
    eqn->set_interface_idx_func(eqn, idxl, idxr);

    // rotate_to_local fills prodl_local / prodr_local on the equation
    // (per-cell rows). With norm = +x̂ this is a no-op on q components.
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

    // Per-side single-interface updates, each checked in its OWN metric.
    for (int side = 0; side < 2; side++) {
      double qn[5], qu[5];
      double sd = (side == 0) ? sqrt_det_L : sqrt_det_R;
      const double (*ig)[3] = (side == 0) ? inv_g_L : inv_g_R;
      for (int i = 0; i < 5; i++)
        qn[i] = (side == 0) ? qL_loc[i] - dt_dx * amdq[i]
                            : qR_loc[i] - dt_dx * apdq[i];
      for (int i = 0; i < 5; i++) qu[i] = qn[i] / sd;
      bool D_ok, S2_ok, tau_ok;
      record_admissibility(ig, qu, &D_ok, &S2_ok, &tau_ok);
      if (!D_ok)   (*D_v)++;
      if (!S2_ok)  (*S2_v)++;
      if (!tau_ok) (*tau_v)++;
      (*total)++;
    }
  }
}

// Floor-value sweep through the LOW_ORDER curved-Lax fallback — the path
// that processes post-repair near-floor cells in production. LOW_ORDER is
// rp-independent (wave_tetrad_dispatch routes every rp_type to
// wave_lax_curved), so this runs ONE equation; the old per-rp trio
// (near_floor_{lax,hll,hllc}) ran identical math three times. Cells L/R
// carry per-cell metrics at the production grid spacing.
static void
run_near_floor_sweep(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range,
    WV_GR_EULER_TETRAD_RP_LAX);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);
  fill_prods_at(spacetime, x - 0.5*GR_EULER_POSITIVITY_DX, y, z, prods_L);
  fill_prods_at(spacetime, x + 0.5*GR_EULER_POSITIVITY_DX, y, z, prods_R);
  if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
      prods_R[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }

  double tau_sweep[] = {
    1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10,
    1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14
  };

  fprintf(stderr, "[near-floor %s] tau_target  total  D<=0  s^2<=0  tau<0\n", label);
  for (size_t i = 0; i < sizeof(tau_sweep)/sizeof(*tau_sweep); i++) {
    int D_v = 0, S2_v = 0, tau_v = 0, total = 0;
    run_near_floor_for_floor_value(eqn, prods_L, prods_R,
      tau_sweep[i], &D_v, &S2_v, &tau_v, &total);
    fprintf(stderr, "  %-12s  %.0e  %4d  %4d  %4d   %4d\n",
      label, tau_sweep[i], total, D_v, S2_v, tau_v);
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_near_floor_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_near_floor_sweep(st, "Mink", 0.0, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
// Production τ-floor cells live near the horizon — probe the curved
// near-horizon regime the old single-cell Minkowski version never saw.
void test_near_floor_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_near_floor_sweep(st, "Schw-near", 0.23, 0.0, 0.0);
  run_near_floor_sweep(st, "Schw-mid",  0.33, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_near_floor_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_near_floor_sweep(st, "Kerr-near",    0.17, 0.0,  0.0);
  run_near_floor_sweep(st, "Kerr-offaxis", 0.23, 0.10, 0.0);
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
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
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

// ---------------------------------------------------------------------------
// DELIBERATE 3-CELL EXCEPTION: this is a debugging HARNESS, not scheme
// coverage — it replicates the production wave_prop update verbatim
// (exact dump coords, per-cell metrics, sweep dir) so future production
// fires can be reproduced. All scheme-coverage tests are 2-cell with
// per-cell metrics (per-side CFL = 1 probes dominate the two-interface
// form by the convex half-split argument).
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
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

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


// ---------------------------------------------------------------------------
// 9. Unified positivity registry
//
// One sweep over ALL the positivity condition kinds this suite has
// accumulated — generic extreme primitives (shared table), τ-floor
// neighborhoods (the near_floor family's states), small-τ/D ladder
// states, metric-aware margin states, and captured production states
// (direct_state dumps) — run through MULTIPLE Riemann solves of the same
// equation family:
//
//   tetrad-Lax  (rp=LAX,  HIGH_ORDER)
//   tetrad-HLL  (rp=HLL,  HIGH_ORDER)
//   curved-Lax  (LOW_ORDER — rp-independent fallback)
//
// The registry RETAINS which conditions fail for which Riemann solve:
// every (condition × flux) cell of the matrix is counted independently
// and the full matrix is printed at the end of the sweep (grep
// "REGISTRY" in the test output). The printed matrix — not a hardcoded
// expectation table — is the retained knowledge: per-cell pass/fail at
// near-boundary margins is FP-fragile under -ffast-math, so the binding
// asserts are the structural invariants only:
//
//   - Minkowski control: every flux, every condition, zero violations.
//   - D > 0 holds for every flux on every condition (curved too).
//   - Fine-dx shared-table clean (assert_shared_clean callers): zero s²
//     for tetrad-Lax and curved-Lax on the shared PRIM block — the
//     production-resolution contract inherited from the retired
//     lax_convexity_tetrad_vs_curved_* tests.
//   - Coarse-dx gross ceiling (callers passing a coarse ladder): s²
//     violations ≤ 10% of updates per Lax flux — exact counts are
//     FP-fragile; the ceiling guards against gross regression only.
//
// Probes run at CFL = 1 (extremal).
//
// AXES ABSORBED DURING THE CTEST SPLIT (CTEST_SPLIT_PLAN.md):
//
// dx axis (from lax_convexity_tetrad_vs_curved_*, retired): the shared
// PRIM table re-swept at GR_EULER_POS_COARSE_DX = 0.15 over near-horizon
// coarse ladders ("coarse <case>" rows). Measured history (2026-06-10,
// CFL = 1): production dx — both Lax schemes clean everywhere; coarse dx
// — the ordering REVERSES: Schw s² fails 1 (tetrad) vs 3 (curved), Kerr
// a=0.9 2 (tetrad) vs 7 (curved), all high-momentum states, all
// downstream cells. Mechanism: MODE A (metric tilt — interface inverse
// metric contracts cell-built momentum, tilting the cone O(Δγ); hurts
// tetrad-Lax only at near-floor margins) vs MODE B (cross-geometry flux
// content in curved Lax's per-cell ΔF, not bounded by a_max·Δq; grows
// with |S|² and Δγ — see SESSION_NOTES_POSITIVITY_UNIFICATION.md §4,
// §4.1). "Curved Lax is the superior fallback" is a REGIME statement —
// true at production resolution, reversed at coarse Δγ.
//
// Full margin axes (from lax_margin_*, retired): W ∈ {10, 100, 500} ×
// dirs {x, x+y, x+z} × thermo {cold p=1e-6, warm p=1}, velocity
// metric-normalized PER CELL so γ_ij·v^i·v^j = 1 − 1/W² (relative s²
// margin ~1/W², the production post-repair scale). Measured history
// (2026-06-10, Schw/Kerr ladders, production dx, CFL = 1):
//
//            tetrad s²-fails    curved s²-fails
//   W=10        5/84,  6/144      27/84, 49/144
//   W=100      64/84, 84/144      42/84, 71/144
//   W=500      82/84, 84/144      42/84, 71/144
//
// (Curved-Lax columns recorded under the pre-adoption full-3D amax; the
// x-only normal-direction penalization adopted later that day — §9 of
// the session notes — lowers them slightly, e.g. Schw 25/40/40, Kerr
// 40/54/55 at the first post-split run. Tetrad columns are amax-
// independent and reproduce exactly. The live REGISTRY-MARGIN print is
// the current record.)
//
// Neither scheme preserves the invariant domain at CFL = 1 once the
// margin 1/W² approaches O(Δγ); Mode B dominates at moderate margins and
// PLATEAUS, Mode A keeps growing and saturates. Production sees
// post-repair states at W ≲ 3 (margins ≳ 0.1, where the table is clean)
// at CFL < 1 — the repair cascade is structural; no flux construction
// makes near-floor states safe at the CFL limit. Minkowski is the
// control: Δγ = 0, so a violation there at ANY margin or thermo is a
// genuine solver bug, not Mode A/B.
// ---------------------------------------------------------------------------

// Bracket provenance: the whole SR kernel family uses the exact
// tangential-aware eigenvalue estimate (MB05 eqs 21–23) since
// 2026-06-11; the historical 1D velocity-addition form was DELETED
// after the bracket A/B study (HLLC_AUDIT_PLAN.md — registry columns
// for both brackets, plus production BHL static/spinning runs, showed
// the exact form stable with slightly lower repair activity; tetrad-Lax
// margin rows are near-identical between brackets, tetrad-HLL W=10
// rows shift, e.g. Schw 18→25/84). Measured-history tables in this
// header predating that date were recorded under the addition bracket.
enum gr_pos_flux {
  GRPF_TET_LAX = 0, GRPF_TET_HLL = 1, GRPF_TET_HLLC = 2, GRPF_CUR_LAX = 3,
  GRPF_N = 4
};
static const char *gr_pos_flux_name[GRPF_N] =
  { "tetrad-Lax", "tetrad-HLL", "tetrad-HLLC", "curved-Lax" };

#define GRPR_MAX_CONDS 128
// Coarse-Δγ axis spacing (near-horizon production-coarse cells; the
// regime where Mode B reverses the Lax ordering — see header comment).
#define GR_EULER_POS_COARSE_DX (0.15)
struct gr_pos_registry_row {
  const char *name;
  int n[GRPF_N], D[GRPF_N], S2[GRPF_N], tau[GRPF_N];
};

// Probe one (qL, qR) pair of DENSITIZED conservatives through one flux
// path at CFL = 1; counts like run_positivity_for_rp's core. Inputs are
// guard-checked for admissibility in their own cell metrics (skip = 0).
static int
gr_pos_probe_qpair(struct gkyl_wv_eqn *eqn,
  const double *prods_L, const double *prods_R,
  const double qL[5], const double qR[5],
  enum gkyl_wv_flux_type ftype, const char *cond_name,
  int *D_v, int *S2_v, int *tau_v)
{
  double sqrt_det_L = sqrt(prods_L[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_R = sqrt(prods_R[GKYL_GR_SP_SPATIAL_DET]);
  double inv_g_L[3][3], inv_g_R[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      inv_g_L[i][j] = prods_L[GKYL_GR_SP_INV_GIJ + 3*i + j];
      inv_g_R[i][j] = prods_R[GKYL_GR_SP_INV_GIJ + 3*i + j];
    }

  double qL_und[5], qR_und[5];
  for (int i = 0; i < 5; i++) {
    qL_und[i] = qL[i] / sqrt_det_L;
    qR_und[i] = qR[i] / sqrt_det_R;
  }
  bool d_ok, s_ok, t_ok;
  record_admissibility(inv_g_L, qL_und, &d_ok, &s_ok, &t_ok);
  if (!(d_ok && s_ok && t_ok)) return 0;
  record_admissibility(inv_g_R, qR_und, &d_ok, &s_ok, &t_ok);
  if (!(d_ok && s_ok && t_ok)) return 0;

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
  double maxs = eqn->waves_func(eqn, ftype,
    delta, qL_loc, qR_loc, 1.0, 1.0, waves, speeds);
  if (!(maxs > 0.0)) return 0;
  double amdq[5], apdq[5];
  eqn->qfluct_func(eqn, ftype,
    qL_loc, qR_loc, 1.0, 1.0, waves, speeds, amdq, apdq);

  double dt_dx = 1.0 / maxs;  // CFL = 1: the extremal convexity probe
  for (int side = 0; side < 2; side++) {
    double qn[5];
    for (int i = 0; i < 5; i++)
      qn[i] = (side == 0) ? qL_loc[i] - dt_dx * amdq[i]
                          : qR_loc[i] - dt_dx * apdq[i];
    double sd = (side == 0) ? sqrt_det_L : sqrt_det_R;
    const double (*ig)[3] = (side == 0) ? inv_g_L : inv_g_R;
    double qn_und[5];
    for (int i = 0; i < 5; i++) qn_und[i] = qn[i] / sd;
    bool D_ok, S2_ok, tau_ok;
    record_admissibility(ig, qn_und, &D_ok, &S2_ok, &tau_ok);
    if (!D_ok)   (*D_v)++;
    if (!S2_ok)  (*S2_v)++;
    if (!tau_ok) (*tau_v)++;
  }
  (void)cond_name;
  return 2;
}

// Small-τ/D ladder rows (provenance: the small_tau_over_D sweep) and the
// captured direct-state pairs are folded into the registry alongside the
// shared PRIM table, the near_floor table, and margin-generated states.
static const struct gr_euler_pos_rp_case gr_pos_small_tau_rows[] = {
  { "sTau p=2.0",    1.0, 2.0,  { 0.30, 0.0, 0.0 }, 1.0, 2.0,  { -0.30, 0.0, 0.0 } },
  { "sTau p=0.5",    1.0, 0.5,  { 0.30, 0.0, 0.0 }, 1.0, 0.5,  { -0.30, 0.0, 0.0 } },
  { "sTau p=0.1",    1.0, 0.1,  { 0.30, 0.0, 0.0 }, 1.0, 0.1,  { -0.30, 0.0, 0.0 } },
  { "sTau p=0.01",   1.0, 0.01, { 0.30, 0.0, 0.0 }, 1.0, 0.01, { -0.30, 0.0, 0.0 } },
  { "sTau BHL-proxy",1.0, 1.0,  { 0.60, 0.10, 0.0 },1.0, 1.0,  { -0.60, 0.10, 0.0 } },
  { "sTau asym-hot-cold", 1.0, 2.0, { 0.10, 0.0, 0.0 }, 0.1, 0.001, { -0.10, 0.0, 0.0 } },
};

// Captured production dumps (densitized conservatives; provenance:
// direct_state family, BHL repair fires, SESSION_NOTES_2 §17 onward).
// Each dump gives (C; L; R); the registry probes the L–C and C–R
// interfaces at the canonical near-horizon sample point (0.3, 0.2, 0)
// on Schw M=0.1, matching the single-metric direct_state tests.
struct gr_pos_direct_case {
  const char *name;
  double qC[5], qL[5], qR[5];
};
static const struct gr_pos_direct_case gr_pos_direct_cases[] = {
  { "direct bhl-repair-#0",
    { 2.978, 2.499, 0.118, 0.000, 0.122 },
    { 2.541, 2.146, 0.045, 0.000, 0.136 },
    { 3.395, 2.827, 0.222, 0.000, 0.095 } },
  { "direct bhl-repair-s2-#2",
    { 3.764, 3.765, -1.130, 0.000, 0.418 },
    { 3.245, 3.185, -1.073, 0.000, 0.332 },
    { 3.510, 2.847,  0.308, 0.000, 0.127 } },
  { "direct curvedLax-s2-#10",
    { 4.313, 4.221,  0.779, 0.000, 1.084 },
    { 4.357, 4.743, -1.245, 0.000, 1.187 },
    { 1.778, 1.109,  0.346, 0.000, 0.067 } },
  { "direct bhl-repair-#4",
    { 3.407, 2.633,  0.209, 0.000, 0.021 },
    { 3.123, 2.165,  0.119, 0.000, 0.025 },
    { 3.009, 3.319, -0.766, 0.000, 0.526 } },
  { "direct bhl-repair-#19",
    { 3.505, 2.969, -1.664, 0.000, 0.193 },
    { 3.872, 2.574,  0.285, 0.000, 0.025 },
    { 2.251, 2.424, -0.887, 0.000, 0.392 } },
  { "direct full3d-s2-#10",
    { 2.528, 2.695, -0.005, 0.000, 0.655 },
    { 3.578, 3.814, -1.539, 0.000, 0.987 },
    { 0.878, 0.670,  0.193, 0.000, 0.026 } },
  { "direct full3d-s2-#12",
    { 1.059, 0.935,  0.167, 0.000, 0.117 },
    { 2.364, 3.174, -0.673, 0.000, 0.934 },
    { 0.244, 0.182,  0.082, 0.000, 0.015 } },
  { "direct full3d-s2-#17",
    { 1.695, 2.594, -0.757, 0.000, 0.986 },
    { 1.194, 1.231, -0.184, 0.000, 0.303 },
    { 2.087, 2.902,  0.585, 0.000, 1.045 } },
};

// coarse_pts/n_coarse: near-horizon ladder for the coarse-dx axis (pass
// NULL/0 to skip — Minkowski has no Δγ, BHL has no measured coarse
// contract). assert_shared_clean: enforce the fine-dx production
// contract (zero s² for both Lax schemes on the shared PRIM block);
// false for BHL, where the convexity contract was never measured at
// CFL = 1 (positivity_curved_lax_bhl covers curved Lax at CFL = 0.9).
static void
run_positivity_registry(struct gkyl_gr_spacetime *spacetime,
  const char *label,
  const struct positivity_point *pts, int n_pts,
  const struct positivity_point *coarse_pts, int n_coarse,
  bool is_minkowski_control, bool assert_shared_clean)
{
  struct gkyl_gr_euler_eos eos = eos_modes[0];

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  // One equation per HIGH_ORDER rp; curved-Lax runs through the LAX
  // equation's LOW_ORDER branch (rp-independent).
  struct gkyl_wv_eqn *eqn_lax = make_eqn(eos, conf_range,
    WV_GR_EULER_TETRAD_RP_LAX);
  struct gkyl_wv_eqn *eqn_hll = make_eqn(eos, conf_range,
    WV_GR_EULER_TETRAD_RP_HLL);
  struct gkyl_wv_eqn *eqn_hllc = make_eqn(eos, conf_range,
    WV_GR_EULER_TETRAD_RP_HLLC);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn_lax,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
  gkyl_gr_euler_tetrad_set_auxfields(eqn_hll,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
  gkyl_gr_euler_tetrad_set_auxfields(eqn_hllc,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
  struct wv_gr_euler_tetrad *grm_lax = container_of(eqn_lax,
    struct wv_gr_euler_tetrad, eqn);
  struct wv_gr_euler_tetrad *grm_hll = container_of(eqn_hll,
    struct wv_gr_euler_tetrad, eqn);
  struct wv_gr_euler_tetrad *grm_hllc = container_of(eqn_hllc,
    struct wv_gr_euler_tetrad, eqn);

  struct gkyl_wv_eqn *flux_eqn[GRPF_N] = { eqn_lax, eqn_hll, eqn_hllc, eqn_lax };
  struct wv_gr_euler_tetrad *flux_grm[GRPF_N] =
    { grm_lax, grm_hll, grm_hllc, grm_lax };
  enum gkyl_wv_flux_type flux_ft[GRPF_N] = {
    GKYL_WV_HIGH_ORDER_FLUX, GKYL_WV_HIGH_ORDER_FLUX, GKYL_WV_HIGH_ORDER_FLUX,
    GKYL_WV_LOW_ORDER_FLUX };

  struct gr_pos_registry_row reg[GRPR_MAX_CONDS];
  int n_conds = 0;
  memset(reg, 0, sizeof(reg));

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);

  // --- PRIM-like conditions: shared table + τ-floor rows + small-τ/D
  //     ladder rows, swept over the point ladder.
  struct prim_block { const struct gr_euler_pos_rp_case *rows; int n; };
  // τ-floor rows materialized from g_near_floor_cases (floored side as
  // given; p_min = 0 builds cleanly with h = 1).
  struct gr_euler_pos_rp_case nf_rows[16];
  int n_nf = (int)(sizeof(g_near_floor_cases) / sizeof(*g_near_floor_cases));
  for (int i = 0; i < n_nf; i++) {
    const struct near_floor_rp *nf = &g_near_floor_cases[i];
    struct gr_euler_pos_rp_case *r = &nf_rows[i];
    r->name = nf->name;
    if (nf->floor_on_left) {
      r->rho_L = nf->rho_min; r->p_L = nf->p_min;
      r->rho_R = nf->rho_typ; r->p_R = nf->p_typ;
      for (int k = 0; k < 3; k++) { r->v_L[k] = nf->v_min[k]; r->v_R[k] = nf->v_typ[k]; }
    } else {
      r->rho_L = nf->rho_typ; r->p_L = nf->p_typ;
      r->rho_R = nf->rho_min; r->p_R = nf->p_min;
      for (int k = 0; k < 3; k++) { r->v_L[k] = nf->v_typ[k]; r->v_R[k] = nf->v_min[k]; }
    }
  }
  struct prim_block blocks[] = {
    { gr_euler_pos_cases, GR_EULER_POS_NUM_CASES },
    { nf_rows, n_nf },
    { gr_pos_small_tau_rows,
      (int)(sizeof(gr_pos_small_tau_rows) / sizeof(*gr_pos_small_tau_rows)) },
  };

  for (size_t b = 0; b < sizeof(blocks)/sizeof(*blocks); b++) {
    for (int c = 0; c < blocks[b].n; c++) {
      const struct gr_euler_pos_rp_case *rc = &blocks[b].rows[c];
      struct gr_pos_registry_row *row = &reg[n_conds++];
      row->name = rc->name;
      for (int pi = 0; pi < n_pts; pi++) {
        double dxh = 0.5 * GR_EULER_POSITIVITY_DX;
        fill_prods_at(spacetime, pts[pi].x - dxh, pts[pi].y, pts[pi].z, prods_L);
        fill_prods_at(spacetime, pts[pi].x + dxh, pts[pi].y, pts[pi].z, prods_R);
        if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
            prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;
        for (int fl = 0; fl < GRPF_N; fl++) {
          row->n[fl] += run_positivity_for_rp(flux_eqn[fl], flux_grm[fl],
            prods_L, prods_R, rc, flux_ft[fl], 1.0,
            &row->D[fl], &row->S2[fl], &row->tau[fl]);
        }
      }
    }
  }

  // --- MARGIN conditions: metric-aware W targets, full axes (absorbed
  //     from the retired lax_margin_* — see header comment): 3 dirs
  //     (face-normal, in-plane diagonal activating γ_xy, out-of-plane
  //     diagonal activating γ_xz) × {cold, warm} thermo. Velocity is
  //     metric-normalized PER CELL so both cells carry the prescribed W
  //     in their own frames.
  static const double Wt[] = { 10.0, 100.0, 500.0 };
  static const double mdirs[][3] = {
    { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 1.0, 0.0, 1.0 },
  };
  static const struct { double rho, p; const char *tag; } mthermo[] = {
    { 1.0, 1.0e-6, "cold" }, { 1.0, 1.0, "warm" },
  };
  enum { N_WT = 3, N_MDIRS = 3, N_MTH = 2 };
  static char margin_names[N_WT * N_MDIRS * N_MTH][40];
  int margin_S2[N_WT][GRPF_N], margin_n[N_WT][GRPF_N];
  memset(margin_S2, 0, sizeof(margin_S2));
  memset(margin_n, 0, sizeof(margin_n));
  for (size_t wi = 0; wi < N_WT; wi++) {
    for (size_t di = 0; di < N_MDIRS; di++) {
      for (size_t ti = 0; ti < N_MTH; ti++) {
        struct gr_pos_registry_row *row = &reg[n_conds++];
        char *nm = margin_names[(wi * N_MDIRS + di) * N_MTH + ti];
        snprintf(nm, 40, "margin W=%g dir=%zu %s", Wt[wi], di, mthermo[ti].tag);
        row->name = nm;
        double vsq_t = 1.0 - 1.0/(Wt[wi]*Wt[wi]);
        for (int pi = 0; pi < n_pts; pi++) {
          double dxh = 0.5 * GR_EULER_POSITIVITY_DX;
          fill_prods_at(spacetime, pts[pi].x - dxh, pts[pi].y, pts[pi].z, prods_L);
          fill_prods_at(spacetime, pts[pi].x + dxh, pts[pi].y, pts[pi].z, prods_R);
          if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
              prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;
          struct gr_euler_pos_rp_case rc = {
            .name = row->name,
            .rho_L = mthermo[ti].rho, .p_L = mthermo[ti].p,
            .rho_R = mthermo[ti].rho, .p_R = mthermo[ti].p };
          double nsq_L = 0.0, nsq_R = 0.0;
          for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
              nsq_L += prods_L[GKYL_GR_SP_GIJ + 3*i + j] * mdirs[di][i] * mdirs[di][j];
              nsq_R += prods_R[GKYL_GR_SP_GIJ + 3*i + j] * mdirs[di][i] * mdirs[di][j];
            }
          for (int i = 0; i < 3; i++) {
            rc.v_L[i] = mdirs[di][i] * sqrt(vsq_t / nsq_L);
            rc.v_R[i] = mdirs[di][i] * sqrt(vsq_t / nsq_R);
          }
          for (int fl = 0; fl < GRPF_N; fl++) {
            row->n[fl] += run_positivity_for_rp(flux_eqn[fl], flux_grm[fl],
              prods_L, prods_R, &rc, flux_ft[fl], 1.0,
              &row->D[fl], &row->S2[fl], &row->tau[fl]);
          }
        }
        for (int fl = 0; fl < GRPF_N; fl++) {
          margin_S2[wi][fl] += row->S2[fl];
          margin_n[wi][fl] += row->n[fl];
        }
      }
    }
  }

  // --- COARSE-DX conditions (absorbed from the retired
  //     lax_convexity_tetrad_vs_curved_* — see header comment): the
  //     shared PRIM table re-swept at GR_EULER_POS_COARSE_DX over the
  //     near-horizon coarse ladder. Mode-B regime: s² violations are
  //     EXPECTED, predominantly for curved Lax; bounded only by the
  //     gross-regression ceiling asserted below.
  int first_coarse_cond = n_conds;
  static char coarse_names[GR_EULER_POS_NUM_CASES][56];
  for (int c = 0; c < (n_coarse > 0 ? GR_EULER_POS_NUM_CASES : 0); c++) {
    const struct gr_euler_pos_rp_case *rc = &gr_euler_pos_cases[c];
    struct gr_pos_registry_row *row = &reg[n_conds++];
    snprintf(coarse_names[c], 56, "coarse-dx %s", rc->name);
    row->name = coarse_names[c];
    for (int pi = 0; pi < n_coarse; pi++) {
      double dxh = 0.5 * GR_EULER_POS_COARSE_DX;
      fill_prods_at(spacetime, coarse_pts[pi].x - dxh,
        coarse_pts[pi].y, coarse_pts[pi].z, prods_L);
      fill_prods_at(spacetime, coarse_pts[pi].x + dxh,
        coarse_pts[pi].y, coarse_pts[pi].z, prods_R);
      if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
          prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;
      for (int fl = 0; fl < GRPF_N; fl++) {
        row->n[fl] += run_positivity_for_rp(flux_eqn[fl], flux_grm[fl],
          prods_L, prods_R, rc, flux_ft[fl], 1.0,
          &row->D[fl], &row->S2[fl], &row->tau[fl]);
      }
    }
  }

  // --- DIRECT_Q conditions: captured production dumps around the
  //     canonical near-horizon point, with PER-CELL metrics (L at
  //     x − dx/2, R at x + dx/2 — production wave_prop reads adjacent
  //     cells' metrics, never one shared metric). Probe the L–C and C–R
  //     interfaces. These are fixed conservatives, so their
  //     admissibility margins do NOT adapt to the local metric — they
  //     sit near the cone boundary in any frame (they were dumped
  //     because they produced failures). The Minkowski-control assert
  //     below therefore excludes them.
  int first_direct_cond = n_conds;
  for (size_t c = 0; c < sizeof(gr_pos_direct_cases)/sizeof(*gr_pos_direct_cases); c++) {
    const struct gr_pos_direct_case *dc = &gr_pos_direct_cases[c];
    struct gr_pos_registry_row *row = &reg[n_conds++];
    row->name = dc->name;
    fill_prods_at(spacetime, pts[0].x - 0.5*GR_EULER_POSITIVITY_DX,
      pts[0].y, pts[0].z, prods_L);
    fill_prods_at(spacetime, pts[0].x + 0.5*GR_EULER_POSITIVITY_DX,
      pts[0].y, pts[0].z, prods_R);
    if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
        prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;

    // Print the dump's input relative s² margins under THIS metric, so
    // the matrix below self-documents: fixed conservatives carry their
    // margins with them, and a ~10% flat-space margin failing at CFL = 1
    // is the state being near-boundary, not a solver bug.
    {
      double sd = sqrt(prods_L[GKYL_GR_SP_SPATIAL_DET]);
      double ig[3][3];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          ig[i][j] = prods_L[GKYL_GR_SP_INV_GIJ + 3*i + j];
      const double *qs[3] = { dc->qC, dc->qL, dc->qR };
      const char *who[3] = { "C", "L", "R" };
      fprintf(stderr, "  REGISTRY-DIRECT %-24s input rel s² margins:",
        dc->name);
      for (int s = 0; s < 3; s++) {
        double q_und[5];
        for (int i = 0; i < 5; i++) q_und[i] = qs[s][i] / sd;
        double dt2 = (q_und[0] + q_und[4]) * (q_und[0] + q_und[4]);
        double msq = gkyl_gr_euler_mom_sq(ig, q_und[1], q_und[2], q_und[3]);
        fprintf(stderr, "  %s=%.3f", who[s], (dt2 - msq) / dt2);
      }
      fprintf(stderr, "\n");
    }

    for (int fl = 0; fl < GRPF_N; fl++) {
      row->n[fl] += gr_pos_probe_qpair(flux_eqn[fl], prods_L, prods_R,
        dc->qL, dc->qC, flux_ft[fl], dc->name,
        &row->D[fl], &row->S2[fl], &row->tau[fl]);
      row->n[fl] += gr_pos_probe_qpair(flux_eqn[fl], prods_L, prods_R,
        dc->qC, dc->qR, flux_ft[fl], dc->name,
        &row->D[fl], &row->S2[fl], &row->tau[fl]);
    }
  }

  // --- Print the registry matrix: the retained "what fails for which
  //     Riemann solve" record. Only conditions with at least one
  //     violation are listed (plus a global summary line).
  int tot_S2[GRPF_N] = {0}, tot_D[GRPF_N] = {0}, tot_tau[GRPF_N] = {0}, tot_n[GRPF_N] = {0};
  int ctl_S2[GRPF_N] = {0}, ctl_tau[GRPF_N] = {0};  // control sums exclude DIRECT_Q
  fprintf(stderr, "[REGISTRY %s] conditions with violations "
    "(s²/τ per flux; n = updates):\n", label);
  for (int c = 0; c < n_conds; c++) {
    bool any = false;
    for (int fl = 0; fl < GRPF_N; fl++) {
      tot_S2[fl] += reg[c].S2[fl]; tot_D[fl] += reg[c].D[fl];
      tot_tau[fl] += reg[c].tau[fl]; tot_n[fl] += reg[c].n[fl];
      if (c < first_direct_cond) {
        ctl_S2[fl] += reg[c].S2[fl]; ctl_tau[fl] += reg[c].tau[fl];
      }
      if (reg[c].S2[fl] || reg[c].tau[fl] || reg[c].D[fl]) any = true;
    }
    if (any) {
      fprintf(stderr, "  REGISTRY %-24s", reg[c].name);
      for (int fl = 0; fl < GRPF_N; fl++)
        fprintf(stderr, "  %s: s2=%d τ=%d (n=%d)", gr_pos_flux_name[fl],
          reg[c].S2[fl], reg[c].tau[fl], reg[c].n[fl]);
      fprintf(stderr, "\n");
    }
  }
  for (int fl = 0; fl < GRPF_N; fl++)
    fprintf(stderr, "  REGISTRY-TOTAL %-12s s2=%d τ=%d D=%d across %d updates\n",
      gr_pos_flux_name[fl], tot_S2[fl], tot_tau[fl], tot_D[fl], tot_n[fl]);

  // Per-W margin aggregate (the retired lax_margin_* table, now per flux).
  for (int wi = 0; wi < N_WT; wi++) {
    fprintf(stderr, "  REGISTRY-MARGIN W=%-5g (margin %.0e):", Wt[wi],
      1.0/(Wt[wi]*Wt[wi]));
    for (int fl = 0; fl < GRPF_N; fl++)
      fprintf(stderr, "  %s %d/%d", gr_pos_flux_name[fl],
        margin_S2[wi][fl], margin_n[wi][fl]);
    fprintf(stderr, "\n");
  }

  // Structural asserts (see header comment).
  for (int fl = 0; fl < GRPF_N; fl++) {
    TEST_CHECK_( tot_D[fl] == 0, "%s: D > 0 violated %d times for %s",
      gr_pos_flux_name[fl], tot_D[fl], label );
    if (is_minkowski_control) {
      // Metric-relative conditions (PRIM/MARGIN) only — captured dumps
      // are near-boundary in any frame and are exempt (see above).
      TEST_CHECK_( ctl_S2[fl] == 0 && ctl_tau[fl] == 0,
        "%s: Minkowski control violated (s² %d, τ %d) for %s",
        gr_pos_flux_name[fl], ctl_S2[fl], ctl_tau[fl], label );
    }
  }

  // Fine-dx production contract on the shared PRIM block (inherited from
  // the retired lax_convexity_* clean asserts; both Lax schemes).
  if (assert_shared_clean) {
    static const int lax_fluxes[] = { GRPF_TET_LAX, GRPF_CUR_LAX };
    for (size_t k = 0; k < sizeof(lax_fluxes)/sizeof(*lax_fluxes); k++) {
      int fl = lax_fluxes[k], s2 = 0;
      for (int c = 0; c < GR_EULER_POS_NUM_CASES; c++) s2 += reg[c].S2[fl];
      TEST_CHECK_( s2 == 0,
        "%s: fine-dx shared-table s² violated %d times for %s",
        gr_pos_flux_name[fl], s2, label );
    }
  }

  // Coarse-dx gross-regression ceiling (inherited from the retired
  // lax_convexity_* coarse asserts; exact counts are FP-fragile, the
  // ceiling guards order-of-magnitude regressions only).
  if (n_coarse > 0) {
    static const int lax_fluxes[] = { GRPF_TET_LAX, GRPF_CUR_LAX };
    for (size_t k = 0; k < sizeof(lax_fluxes)/sizeof(*lax_fluxes); k++) {
      int fl = lax_fluxes[k], s2 = 0, n = 0;
      for (int c = first_coarse_cond; c < first_direct_cond; c++) {
        s2 += reg[c].S2[fl]; n += reg[c].n[fl];
      }
      TEST_CHECK_( s2 <= n / 10,
        "%s: coarse-dx s² violations (%d of %d) exceed gross ceiling for %s",
        gr_pos_flux_name[fl], s2, n, label );
    }
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn_lax);
  gkyl_wv_eqn_release(eqn_hll);
  gkyl_wv_eqn_release(eqn_hllc);
}

void test_positivity_registry_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  // No coarse ladder: Minkowski has no Δγ, so the coarse-dx axis is
  // contentless (identical metrics at any spacing).
  run_positivity_registry(st, "Mink-registry",
    g_positivity_points,
    (int)(sizeof(g_positivity_points) / sizeof(*g_positivity_points)),
    NULL, 0,
    /*is_minkowski_control=*/true, /*assert_shared_clean=*/true);
  gkyl_gr_spacetime_release(st);
}

// BHL M = 0.3 (production spacetime): PRIM/MARGIN conditions sweep the
// bow-shock points; DIRECT_Q dumps are probed at the radial-front point
// (their capture geometry). This is the 2-cell per-cell-metric successor
// of the deleted 3-cell near-horizon/per-cell-metric direct-state tests
// — by convexity, per-side probes at CFL = 1 dominate the two-interface
// center-cell update at CFL = 0.5.
void test_positivity_registry_bhl(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  static const struct positivity_point pts[] = {
    {  0.65,  0.0,   0.0 },  // radial front (DIRECT_Q anchor)
    { -0.05, -0.6,   0.0 },  // bow-shock front
    { -0.05,  0.6,   0.0 },  // mirror
    { -0.05, -0.605, 0.0 },  // hugging horizon
    {  1.0,   1.0,   0.0 },  // control, far
  };
  // No coarse ladder and no fine-clean assert: the convexity contract
  // was never measured on BHL M=0.3 at CFL = 1 (positivity_curved_lax_bhl
  // covers curved Lax at CFL = 0.9); the matrix print is the record.
  run_positivity_registry(st, "BHL(M=0.3)-registry",
    pts, (int)(sizeof(pts) / sizeof(*pts)),
    NULL, 0,
    /*is_minkowski_control=*/false, /*assert_shared_clean=*/false);
  gkyl_gr_spacetime_release(st);
}

void test_positivity_registry_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  // Coarse ladder (from the retired convexity test): centroids chosen so
  // xL = centroid − dx/2 stays outside r_+ = 0.2 at dx = 0.15.
  static const struct positivity_point coarse[] = {
    { 0.28, 0.0, 0.0 }, { 0.33, 0.0, 0.0 }, { 0.41, 0.0, 0.0 },
    { 0.61, 0.0, 0.0 },
  };
  run_positivity_registry(st, "Schw-registry",
    g_ladder_schw_m01, GR_EULER_N_PTS(g_ladder_schw_m01),
    coarse, GR_EULER_N_PTS(coarse),
    /*is_minkowski_control=*/false, /*assert_shared_clean=*/true);
  gkyl_gr_spacetime_release(st);
}

void test_positivity_registry_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  // Coarse ladder (from the retired convexity test): xL stays outside
  // r_+ ≈ 0.144 at dx = 0.15; includes off-axis frame-dragging rows.
  static const struct positivity_point coarse[] = {
    { 0.23, 0.0, 0.0 },  { 0.31, 0.0, 0.0 },  { 0.46, 0.0, 0.0 },
    { 0.23, 0.10, 0.0 }, { 0.31, 0.10, 0.0 },
  };
  run_positivity_registry(st, "Kerr-registry",
    g_ladder_kerr_a09, GR_EULER_N_PTS(g_ladder_kerr_a09),
    coarse, GR_EULER_N_PTS(coarse),
    /*is_minkowski_control=*/false, /*assert_shared_clean=*/true);
  gkyl_gr_spacetime_release(st);
}

TEST_LIST = {
  { "positivity_curved_lax_minkowski",     test_positivity_curved_lax_minkowski },
  { "positivity_curved_lax_schwarzschild", test_positivity_curved_lax_schwarzschild },
  { "positivity_curved_lax_kerr",          test_positivity_curved_lax_kerr },
  { "positivity_curved_lax_bhl",           test_positivity_curved_lax_bhl },

  { "positivity_registry_minkowski",     test_positivity_registry_minkowski },
  { "positivity_registry_schwarzschild", test_positivity_registry_schwarzschild },
  { "positivity_registry_kerr",          test_positivity_registry_kerr },
  { "positivity_registry_bhl",           test_positivity_registry_bhl },

  { "near_floor_minkowski",     test_near_floor_minkowski },
  { "near_floor_schwarzschild", test_near_floor_schwarzschild },
  { "near_floor_kerr",          test_near_floor_kerr },

  { "small_tau_over_D_lax",  test_small_tau_over_D_lax },
  { "small_tau_over_D_hll",  test_small_tau_over_D_hll },
  { "small_tau_over_D_hllc", test_small_tau_over_D_hllc },

  { "direct_state_lax_curved_production_reproducer",
    test_direct_state_lax_curved_production_reproducer },

  { NULL, NULL },
};
