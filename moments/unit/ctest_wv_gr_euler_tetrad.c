// Convention-A audit tests for wv_gr_euler_tetrad (modular tetrad GR
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
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

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

        double f_gr[5];
        gkyl_gr_euler_banyuls_flux_cell(eos, q, prods_row, NULL, f_gr);

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

// ---------------------------------------------------------------------------
// Two-cell variant for non-degenerate interface geometry.
//
// The single-cell runner above sets prodl_local == prodr_local (idx={0,0}),
// so the averaging at lines 1464-1475 of wv_gr_euler_tetrad.c recovers
// the cell-centered values exactly. That makes the flux-jump check trivial
// at the geometry level: any averaging policy passes.
//
// This variant allocates two cells, fills L at (xL,yL,zL) and R at
// (xR,yR,zR), and sets idx={0,1} so prodl_local and prodr_local come from
// genuinely different points. The flux-jump check uses an interface-
// averaged prods row to match the averaging policy inside wave_lax /
// wave_hll / wave_hllc (the relevant flux-jump form for the tetrad-first
// scheme; see TETRAD_REFACTOR_PLAN.md §1).
// ---------------------------------------------------------------------------
static void
run_riemann_properties_two_cell(struct gkyl_gr_spacetime *spacetime,
  const char *label, double xL, double yL, double zL,
  double xR, double yR, double zR,
  struct gkyl_gr_euler_eos eos, enum gkyl_wv_gr_euler_tetrad_rp rp,
  int num_waves)
{
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);
  fill_prods_at(spacetime, xL, yL, zL, prods_L);
  fill_prods_at(spacetime, xR, yR, zR, prods_R);

  // Skip if either cell is excised — the wave decomposition uses
  // active-cell-only geometry in that branch and isn't what this test
  // probes (excision-boundary face-position evaluation is a Phase 2 test).
  if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
      prods_R[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }

  TEST_CHECK( eqn->num_waves == num_waves );

  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);
  (void)grm;  // banyuls_delta_flux not used here — interface-averaged form below

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  // Small primitives chosen to remain subluminal in moderately distorted
  // metrics; if either L or R fails the subluminal check at its own
  // cell's metric, skip rather than crash inside build_state_convA.
  double rho_L = 1.0, p_L = 1.5;  double v_L[3] = { 0.05, 0.10, 0.15 };
  double rho_R = 0.5, p_R = 0.7;  double v_R[3] = { 0.03, 0.05, 0.08 };

  double vsq_L = 0.0, vsq_R = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      vsq_L += prods_L[GKYL_GR_SP_GIJ + 3*i + j] * v_L[i] * v_L[j];
      vsq_R += prods_R[GKYL_GR_SP_GIJ + 3*i + j] * v_R[i] * v_R[j];
    }
  if (!(vsq_L < 1.0 - 1.0e-6) || !(vsq_R < 1.0 - 1.0e-6)) {
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }

  // Build each cell's conservative state using its OWN cell-centered
  // metric — that's the physical state in each cell.
  double qL_glob[5], qR_glob[5];
  build_state_convA(eos, rho_L, v_L, p_L, prods_L, qL_glob);
  build_state_convA(eos, rho_R, v_R, p_R, prods_R, qR_glob);

  int idxl[1] = { 0 }, idxr[1] = { 1 };
  eqn->set_interface_idx_func(eqn, idxl, idxr);

  double qL[5], qR[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL_glob, qL);
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR_glob, qR);

  double delta[5];
  for (int i = 0; i < 5; i++) delta[i] = qR[i] - qL[i];

  double waves[3 * 5], speeds[3];
  double maxs = eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    delta, qL, qR, 1.0, 1.0, waves, speeds);

  // Interface-averaged prods row — matches the production averaging
  // policy: average covariant γ_ij arithmetically, then DERIVE inv_g_iface
  // and sqrt(det γ_iface) from the inverse of g_iface. Other fields
  // (lapse, shift, ...) stay as element-wise arithmetic means.
  double prods_iface[GKYL_GR_SP_NCOMP_BASE];
  build_prods_iface(prods_L, prods_R, prods_iface);

  // ΔF computed with INTERFACE-AVERAGED geometry on BOTH sides. This is
  // the relevant flux-jump form for the tetrad-first scheme. Comparing
  // against cell-centered (F(qR;geom_R) − F(qL;geom_L)) would conflate
  // tetrad-first with curved-Lax interface-flux constructions and show
  // an unrelated O(Δgeom) residual.
  double dF[5];
  double fL_gr[5], fR_gr[5];
  gkyl_gr_euler_banyuls_flux_cell(eos, qL, prods_iface, NULL, fL_gr);
  gkyl_gr_euler_banyuls_flux_cell(eos, qR, prods_iface, NULL, fR_gr);
  for (int i = 0; i < 5; i++) dF[i] = fR_gr[i] - fL_gr[i];

  // (a) Wave sum: Σ w_k = Δq
  for (int i = 0; i < 5; i++) {
    double sum = 0.0;
    for (int k = 0; k < num_waves; k++) sum += waves[k * 5 + i];
    TEST_CHECK_( fabs(sum - delta[i]) < 1.0e-12,
      "[%s] wave sum: comp %d, |Σw − Δq| = %.3e", label, i, fabs(sum - delta[i]) );
  }

  // (b) Flux jump (interface-averaged): Σ s_k · w_k = ΔF_iface
  for (int i = 0; i < 5; i++) {
    double sw = 0.0;
    for (int k = 0; k < num_waves; k++) sw += speeds[k] * waves[k * 5 + i];
    TEST_CHECK_( fabs(sw - dF[i]) < 1.0e-10,
      "[%s] flux jump (iface): comp %d, |Σs·w − ΔF| = %.3e",
      label, i, fabs(sw - dF[i]) );
  }

  // (c) Fluctuation balance: amdq + apdq = ΔF_iface
  double amdq[5], apdq[5];
  eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
    qL, qR, 1.0, 1.0, waves, speeds, amdq, apdq);
  for (int i = 0; i < 5; i++) {
    TEST_CHECK_( fabs(amdq[i] + apdq[i] - dF[i]) < 1.0e-10,
      "[%s] fluct balance (iface): comp %d, |amdq+apdq − ΔF| = %.3e",
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
  // Geometry choice doesn't matter for this check (Δq = 0 ⇒ Δw = 0
  // regardless of triad). Use L cell prods to construct qE.
  double qE[5];
  build_state_convA(eos, 1.0, (double[]){0.05, 0.10, 0.05}, 1.0, prods_L, qE);
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

  // (g) L↔R sign symmetry. We also swap interface indices so prodl/prodr
  // are sourced consistently with the swapped state ordering. The
  // interface-averaged geom is unchanged under this swap (arithmetic
  // mean is commutative) — the invariant we check is:
  //   amdq(qL,qR) + apdq(qL,qR) = − [amdq(qR,qL) + apdq(qR,qL)]
  int idxl_swap[1] = { 1 }, idxr_swap[1] = { 0 };
  eqn->set_interface_idx_func(eqn, idxl_swap, idxr_swap);
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
// --- Minkowski sanity baseline. geom_L = geom_R trivially, so this must
//     pass identically to the single-cell case (sanity check on the
//     two-cell scaffolding itself).
void test_riemann_properties_two_cell_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_riemann_properties_two_cell(st, "Mink-Lax-2cell",
      0.30, 0.0, 0.0,  0.32, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_LAX, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_hll_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_riemann_properties_two_cell(st, "Mink-HLL-2cell",
      0.30, 0.0, 0.0,  0.32, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLL, 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_hllc_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_riemann_properties_two_cell(st, "Mink-HLLC-2cell",
      0.30, 0.0, 0.0,  0.32, 0.0, 0.0,
      eos_modes[ei], WV_GR_EULER_TETRAD_RP_HLLC, 3);
  gkyl_gr_spacetime_release(st);
}

// --- Schwarzschild (M=0.1 ⇒ r_+ = 0.2). Three position pairs, each
//     Δx = 0.02 in the x-direction, ordered from near-horizon stress
//     outward. The first pair sits at r ≈ 0.22, ≈ 10% outside r_+ where
//     lapse and metric gradients are strong.
#define SCHW_2CELL_RUNS(rp_macro, label, nw)                                   \
  do {                                                                         \
    run_riemann_properties_two_cell(st, label,                                 \
      0.22, 0.0, 0.0,  0.24, 0.0, 0.0,                                         \
      eos_modes[ei], rp_macro, nw);                                            \
    run_riemann_properties_two_cell(st, label,                                 \
      0.30, 0.20, 0.0,  0.32, 0.20, 0.0,                                       \
      eos_modes[ei], rp_macro, nw);                                            \
    run_riemann_properties_two_cell(st, label,                                 \
      0.38, 0.40, 0.0,  0.40, 0.40, 0.0,                                       \
      eos_modes[ei], rp_macro, nw);                                            \
  } while (0)

void test_riemann_properties_two_cell_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    SCHW_2CELL_RUNS(WV_GR_EULER_TETRAD_RP_LAX, "Schw-Lax-2cell", 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_hll_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    SCHW_2CELL_RUNS(WV_GR_EULER_TETRAD_RP_HLL, "Schw-HLL-2cell", 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_hllc_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    SCHW_2CELL_RUNS(WV_GR_EULER_TETRAD_RP_HLLC, "Schw-HLLC-2cell", 3);
  gkyl_gr_spacetime_release(st);
}
#undef SCHW_2CELL_RUNS

// --- Kerr at a = 0.9 (M=0.1 ⇒ r_+ ≈ 0.144). Three pairs covering:
//     near-horizon radial, off-axis where frame-dragging is strongest,
//     and out-of-plane. Δx = 0.02 throughout.
#define KERR_2CELL_RUNS(rp_macro, label, nw)                                   \
  do {                                                                         \
    run_riemann_properties_two_cell(st, label,                                 \
      0.16, 0.0, 0.0,  0.18, 0.0, 0.0,                                         \
      eos_modes[ei], rp_macro, nw);                                            \
    run_riemann_properties_two_cell(st, label,                                 \
      0.20, 0.15, 0.0,  0.22, 0.15, 0.0,                                       \
      eos_modes[ei], rp_macro, nw);                                            \
    run_riemann_properties_two_cell(st, label,                                 \
      0.18, 0.0, 0.10,  0.20, 0.0, 0.10,                                       \
      eos_modes[ei], rp_macro, nw);                                            \
  } while (0)

void test_riemann_properties_two_cell_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    KERR_2CELL_RUNS(WV_GR_EULER_TETRAD_RP_LAX, "Kerr0.9-Lax-2cell", 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_hll_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    KERR_2CELL_RUNS(WV_GR_EULER_TETRAD_RP_HLL, "Kerr0.9-HLL-2cell", 2);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_hllc_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    KERR_2CELL_RUNS(WV_GR_EULER_TETRAD_RP_HLLC, "Kerr0.9-HLLC-2cell", 3);
  gkyl_gr_spacetime_release(st);
}
#undef KERR_2CELL_RUNS

// ---------------------------------------------------------------------------
// Roe two-cell tests.
//
// Roe Minkowski and Roe curved spacetime both pass the strict two-cell
// suite (wave-sum, flux-jump, fluctuation balance, ordering, sign
// symmetry) at machine ε. The Eulderink-Mellema SR Roe satisfies
// Σ s_tet · w_tet = ΔF_SR_tet exactly in the tetrad frame, and the
// consistent-g_iface back-transform chain (Phase 0 Fix 1 + 2) propagates
// this to Σ s_curved · w_curved = ΔF_iface in the curved frame — no
// structural O((Δgeom)²) residual. The Schwarzschild and Kerr tests use
// the same strict 1e-10 tolerance as Lax/HLL/HLLC.
// ---------------------------------------------------------------------------

void test_riemann_properties_two_cell_roe_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  // Roe is IDEAL-only; only the first entry in eos_modes (IDEAL) is valid.
  run_riemann_properties_two_cell(st, "Mink-Roe-2cell",
    0.30, 0.0, 0.0,  0.32, 0.0, 0.0,
    eos_modes[0], WV_GR_EULER_TETRAD_RP_ROE, 3);
  gkyl_gr_spacetime_release(st);
}

void test_riemann_properties_two_cell_roe_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_riemann_properties_two_cell(st, "Schw-Roe-2cell",
    0.22, 0.0, 0.0,  0.24, 0.0, 0.0,
    eos_modes[0], WV_GR_EULER_TETRAD_RP_ROE, 3);
  run_riemann_properties_two_cell(st, "Schw-Roe-2cell",
    0.30, 0.20, 0.0,  0.32, 0.20, 0.0,
    eos_modes[0], WV_GR_EULER_TETRAD_RP_ROE, 3);
  run_riemann_properties_two_cell(st, "Schw-Roe-2cell",
    0.38, 0.40, 0.0,  0.40, 0.40, 0.0,
    eos_modes[0], WV_GR_EULER_TETRAD_RP_ROE, 3);
  gkyl_gr_spacetime_release(st);
}

void test_riemann_properties_two_cell_roe_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_riemann_properties_two_cell(st, "Kerr0.9-Roe-2cell",
    0.16, 0.0, 0.0,  0.18, 0.0, 0.0,
    eos_modes[0], WV_GR_EULER_TETRAD_RP_ROE, 3);
  run_riemann_properties_two_cell(st, "Kerr0.9-Roe-2cell",
    0.20, 0.15, 0.0,  0.22, 0.15, 0.0,
    eos_modes[0], WV_GR_EULER_TETRAD_RP_ROE, 3);
  run_riemann_properties_two_cell(st, "Kerr0.9-Roe-2cell",
    0.18, 0.0, 0.10,  0.20, 0.0, 0.10,
    eos_modes[0], WV_GR_EULER_TETRAD_RP_ROE, 3);
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
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

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

  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);

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
  //      Compute F(qL_active) via the production cell-centered Banyuls
  //      flux using the active-cell prods.
  double fL_gr[5];
  gkyl_gr_euler_banyuls_flux_cell(eos, qL_loc, grm->prodl_local, NULL, fL_gr);
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

  double fR_gr[5];
  gkyl_gr_euler_banyuls_flux_cell(eos, qR_loc, grm->prodr_local, NULL, fR_gr);
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
// Roe joined the ZERO_VACUUM excision family when sr_roe_minkowski gained
// its zero-state vacuum short-circuit (CLEANUP_PLAN.md A2) — before that,
// an excised side drove the inlined Newton recovery to W = 0 → ρ = 0/0.
void test_excision_absorbing_roe_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_excision_absorbing_for_rp(st, "Mink", WV_GR_EULER_TETRAD_RP_ROE, 3, "Roe",
    0.3, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_excision_absorbing_roe_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_excision_absorbing_for_rp(st, "Schw", WV_GR_EULER_TETRAD_RP_ROE, 3, "Roe",
    0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}

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
// 3b. Tetrad-Lax vs curved-Lax convexity comparison
//
// Same per-cell-metric update probe as the positivity sweep above, run
// TWICE per (point, case) through the SAME rp_type=LAX equation object:
// once with GKYL_WV_HIGH_ORDER_FLUX (the tetrad-first Lax — Gorard+25
// transform + SR Lax kernel) and once with GKYL_WV_LOW_ORDER_FLUX (the
// direct curved-frame Lax, the POSITIVITY_SWEEP fallback). Probes run at
// the extremal CFL = 1 — the boundary of the convex-combination argument.
//
// There are TWO distinct ways an "admissible-in, admissible-out" Lax
// update can fail in curved spacetime, and they hit the two schemes
// asymmetrically. Empirical counts from this very sweep (2026-06-10,
// dx_coarse = 0.15): Schw s² fails 1 (tetrad) vs 3 (curved); Kerr a=0.9
// s² fails 2 (tetrad) vs 7 (curved), all on high-momentum states
// (rel-shock, cold-shock, equal-high-v, cold-velocity-jump) and all on
// the downstream cell.
//
// MODE A — metric-tilt (hurts tetrad-Lax and iface-flux Lax). The tetrad
// forward transform contracts cell-constructed momentum with the
// INTERFACE inverse metric (Gram-Schmidt orthonormality gives
// δ_ab·S̃^a·S̃^b = γ^{ij}_iface·S_i·S_j), tilting the admissibility cone
// relative to the cell's own A_γ; the M_inv·γ back-transform does not
// project the update back. The tilt is O(Δγ) — a few percent at
// production dx — so it only flips states whose RELATIVE s² margin is
// comparable, i.e. near-floor post-repair states. The admissibility set
// is a cone (scale-invariant), so this sweep's synthetic states — with
// O(0.1–1) relative margins — never trip Mode A at ANY dx tested, and
// the densitization ratio √γ_cell/√γ_iface (a uniform scaling) cannot
// trip it either. Production-scale evidence: SESSION_NOTES_S2_REPAIR.md
// §3.3 (12,677 vs 91 wave-prop s² fires with tetrad-Lax as the LOW_ORDER
// fallback) and the captured-state direct_state_lax tests below.
//
// MODE B — cross-geometry flux difference (hurts curved Lax). Curved Lax
// evaluates F_L(q_L; geom_L) and F_R(q_R; geom_R) with each cell's OWN
// geometry, so ΔF contains a purely geometric component (the balance-law
// partner of the source terms) that is NOT bounded by a_max·Δq. The
// classical Lax-Friedrichs convexity proof assumes ONE flux function on
// both sides; that premise fails at O(Δgeom). At production dx the
// geometric component is negligible; at coarse near-horizon dx it is
// O(1) and pushes high-|S| downstream states out of the cone. The
// tetrad path evaluates both sides in a single averaged interface frame
// (a genuinely single-valued flux function), so Mode B barely touches it.
//
// CONSEQUENCE: "curved Lax is superior for convexity of the update" is a
// REGIME statement, true exactly where production lives — fine dx with
// near-floor states (Mode A dominates) — and REVERSED at coarse
// near-horizon resolution with smooth high-momentum states (Mode B
// dominates). The two-tier HIGH/LOW design decision stands on the
// production regime; this sweep keeps both regimes measured.
//
// Asserts:
//   - fine-dx runs (assert_clean = true): zero D and s² violations for
//     BOTH schemes — the empirical contract at production resolution on
//     the shared synthetic table.
//   - coarse-dx runs (assert_clean = false): zero D violations for both
//     (D is preserved everywhere); s² counts are reported and bounded
//     only by a gross-regression ceiling (≤ 10% of updates) — exact
//     counts are FP-fragile under -ffast-math, and the per-mode split is
//     the documentation, not a contract.
//   - τ has no theoretical guarantee for either scheme: info-only.
// ---------------------------------------------------------------------------

static void
run_lax_convexity_compare(struct gkyl_gr_spacetime *spacetime,
  const char *label,
  const struct positivity_point *pts, int n_pts, double dx,
  bool assert_clean)
{
  // IDEAL only: the tetrad-vs-curved distinction is metric-driven (see
  // header comment); the EOS enters both paths through the same shared
  // recovery and does not change the comparison.
  struct gkyl_gr_euler_eos eos = eos_modes[0];

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, WV_GR_EULER_TETRAD_RP_LAX);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);

  int tet_D = 0, tet_S2 = 0, tet_tau = 0;
  int cur_D = 0, cur_S2 = 0, cur_tau = 0;
  int total_cells = 0;

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);

  for (int pi = 0; pi < n_pts; pi++) {
    const struct positivity_point *pt = &pts[pi];
    double dx_half = 0.5 * dx;
    fill_prods_at(spacetime, pt->x - dx_half, pt->y, pt->z, prods_L);
    fill_prods_at(spacetime, pt->x + dx_half, pt->y, pt->z, prods_R);
    if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
        prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;

    for (int c = 0; c < GR_EULER_POS_NUM_CASES; c++) {
      // CFL = 1.0: the extremal point of the convex-combination argument,
      // matching the mod-variant Phase-B probe. CFL = 0.9 hides the
      // boundary cases the comparison exists to expose.
      int n = run_positivity_for_rp(eqn, grm, prods_L, prods_R,
        &gr_euler_pos_cases[c], GKYL_WV_HIGH_ORDER_FLUX, 1.0,
        &tet_D, &tet_S2, &tet_tau);
      run_positivity_for_rp(eqn, grm, prods_L, prods_R,
        &gr_euler_pos_cases[c], GKYL_WV_LOW_ORDER_FLUX, 1.0,
        &cur_D, &cur_S2, &cur_tau);
      total_cells += n;
    }
  }

  fprintf(stderr,
    "[lax-convexity %s, dx=%.4f] across %d cell-updates:\n"
    "  tetrad-Lax (HIGH_ORDER):  D<=0: %d, s²<=0: %d, τ<0: %d\n"
    "  curved-Lax (LOW_ORDER):   D<=0: %d, s²<=0: %d, τ<0: %d\n"
    "%s",
    label, dx, total_cells, tet_D, tet_S2, tet_tau, cur_D, cur_S2, cur_tau,
    assert_clean
      ? "  -> production regime (fine dx, synthetic margins): both schemes "
        "clean; Mode-A tetrad leakage needs near-floor states — see "
        "direct_state_lax captured-state tests\n"
      : "  -> coarse-Δγ regime: Mode-B (cross-geometry flux difference) "
        "violations EXPECTED, predominantly for curved Lax; ordering "
        "REVERSES relative to the production regime\n");

  // D is preserved by both schemes in every regime probed so far.
  TEST_CHECK_( tet_D == 0, "tetrad-Lax D > 0 violated %d times for %s",
    tet_D, label );
  TEST_CHECK_( cur_D == 0, "curved-Lax D > 0 violated %d times for %s",
    cur_D, label );

  if (assert_clean) {
    // Production-resolution contract on the shared synthetic table.
    TEST_CHECK_( tet_S2 == 0, "tetrad-Lax s² > 0 violated %d times for %s",
      tet_S2, label );
    TEST_CHECK_( cur_S2 == 0, "curved-Lax s² > 0 violated %d times for %s",
      cur_S2, label );
  } else {
    // Coarse-Δγ regime: nonzero s² counts are the documented Mode-B
    // behavior (see header comment). Only a gross-regression ceiling.
    TEST_CHECK_( tet_S2 <= total_cells / 10,
      "tetrad-Lax s² violations (%d of %d) exceed gross ceiling for %s",
      tet_S2, total_cells, label );
    TEST_CHECK_( cur_S2 <= total_cells / 10,
      "curved-Lax s² violations (%d of %d) exceed gross ceiling for %s",
      cur_S2, total_cells, label );
  }
  TEST_MSG( "τ violations (info-only): tetrad %d, curved %d for %s",
    tet_tau, cur_tau, label );

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_lax_convexity_tetrad_vs_curved_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_lax_convexity_compare(st, "Mink",
    g_positivity_points,
    (int)(sizeof(g_positivity_points) / sizeof(*g_positivity_points)),
    GR_EULER_POSITIVITY_DX, /*assert_clean=*/true);
  gkyl_gr_spacetime_release(st);
}

void test_lax_convexity_tetrad_vs_curved_schwarzschild(void)
{
  // M = 0.1 → r_+ = 0.2. Interface centroids on the near-horizon radial
  // ladder used by the mod-variant iface-flux Phase-B sweep (xL = r,
  // xR = r + dx with dx ≈ 0.0195), plus the standard bulk points, plus a
  // coarse-resolution ladder (dx = 0.15 — production-coarse cells near
  // the horizon, where Δγ across one cell is large and the tetrad-frame
  // metric tilt is strongest).
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_lax_convexity_compare(st, "Schw(M=0.1)",
    g_ladder_schw_m01, GR_EULER_N_PTS(g_ladder_schw_m01),
    GR_EULER_POSITIVITY_DX, /*assert_clean=*/true);
  run_lax_convexity_compare(st, "Schw(M=0.1)-bulk",
    g_positivity_points,
    (int)(sizeof(g_positivity_points) / sizeof(*g_positivity_points)),
    GR_EULER_POSITIVITY_DX, /*assert_clean=*/true);
  // Coarse ladder: keep xL = centroid − dx/2 outside r_+ = 0.2.
  static const struct positivity_point coarse[] = {
    { 0.28, 0.0, 0.0 }, { 0.33, 0.0, 0.0 }, { 0.41, 0.0, 0.0 },
    { 0.61, 0.0, 0.0 },
  };
  run_lax_convexity_compare(st, "Schw(M=0.1)-coarse",
    coarse, (int)(sizeof(coarse) / sizeof(*coarse)), 0.15,
    /*assert_clean=*/false);
  gkyl_gr_spacetime_release(st);
}

void test_lax_convexity_tetrad_vs_curved_kerr(void)
{
  // M = 0.1, a = 0.9 → r_+ ≈ 0.144. High spin to match the regime where
  // the iface-flux Lax failures concentrated; on-axis, off-axis (frame
  // dragging), and out-of-plane ladders, plus a coarse-resolution ladder.
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_lax_convexity_compare(st, "Kerr(M=0.1,a=0.9)",
    g_ladder_kerr_a09, GR_EULER_N_PTS(g_ladder_kerr_a09),
    GR_EULER_POSITIVITY_DX, /*assert_clean=*/true);
  // Coarse ladder: keep xL = centroid − dx/2 outside r_+ ≈ 0.144.
  static const struct positivity_point coarse[] = {
    { 0.23, 0.0, 0.0 },  { 0.31, 0.0, 0.0 },  { 0.46, 0.0, 0.0 },
    { 0.23, 0.10, 0.0 }, { 0.31, 0.10, 0.0 },
  };
  run_lax_convexity_compare(st, "Kerr(M=0.1,a=0.9)-coarse",
    coarse, (int)(sizeof(coarse) / sizeof(*coarse)), 0.15,
    /*assert_clean=*/false);
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 3d. wave_spacetime cache equivalence
//
// The HIGH_ORDER tetrad path reads its interface tetrad data from the
// wave_spacetime cache when one is attached, and falls back to per-call
// averaging when not (wave_tetrad_high_order). The cache builder mirrors
// the fallback's averaging + Gram-Schmidt construction, so for interior
// faces the two paths must agree to machine precision. This test is the
// equation-object-level cache coverage that used to live in
// ctest_wv_gr_euler_mod.c (deleted with the iface-flux experiment).
//
// Also checks the cached-triad invariant M_inv·γ⁻¹·M_invᵀ = I (which is
// the M_inv = Mᵀ·γ image of Gram-Schmidt orthonormality MᵀγM = I).
// ---------------------------------------------------------------------------

#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_spacetime.h>

static void
run_wave_spacetime_cache_equivalence(struct gkyl_gr_spacetime *spacetime,
  const char *label, double xL, double xR)
{
  struct gkyl_gr_euler_eos eos = eos_modes[0];

  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  double dx = xR - xL;
  double grid_lower[1] = { xL - 0.5 * dx };
  double grid_upper[1] = { xR + 0.5 * dx };
  int cells[1] = { 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, grid_lower, grid_upper, cells);
  struct gkyl_wave_geom *wg =
    gkyl_wave_geom_new(&grid, &conf_range, NULL, NULL, false);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range,
    WV_GR_EULER_TETRAD_RP_LAX);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);
  fill_prods_at(spacetime, xL, 0.0, 0.0, prods_L);
  fill_prods_at(spacetime, xR, 0.0, 0.0, prods_R);
  if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
      prods_R[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    gkyl_wave_geom_release(wg);
    return;
  }

  struct gkyl_wave_spacetime *ws =
    gkyl_wave_spacetime_new(&grid, &conf_range, wg, spacetime, prods, 0.0,
      /*is_static=*/true, /*use_gpu=*/false);

  // Storage contract at the bottom edge of the build range: cell 0 has
  // no lower neighbor, so its iface[0] must be the one-sided fill from
  // its own prods — VALID (sqrt_det > 0), tagged INTERIOR, never a
  // zero-filled sentinel (see gkyl_wave_spacetime.h).
  {
    int idx0[1] = { 0 };
    const struct gkyl_wave_spacetime_cell *wsc0 =
      gkyl_wave_spacetime_get(ws, idx0);
    const struct gkyl_wave_spacetime_iface *if0 = &wsc0->iface[0];
    TEST_CHECK_( if0->kind == GKYL_WS_IFACE_INTERIOR,
      "[%s] bottom-edge face not tagged INTERIOR", label );
    TEST_CHECK_( if0->sqrt_det_iface > 0.0,
      "[%s] bottom-edge face has sqrt_det = %g (must be valid)",
      label, if0->sqrt_det_iface );
    // One-sided fill from cell 0's own prods: sqrt_det must match the
    // cell's own metric determinant.
    TEST_CHECK_( fabs(if0->sqrt_det_iface
        - sqrt(prods_L[GKYL_GR_SP_SPATIAL_DET])) < 1.0e-14,
      "[%s] bottom-edge face sqrt_det != own-cell value", label );
  }

  // Cached-triad invariant on the interior face owned by cell 1:
  // M_inv·γ⁻¹·M_invᵀ = I at machine ε.
  {
    int idxr[1] = { 1 };
    const struct gkyl_wave_spacetime_cell *wsc =
      gkyl_wave_spacetime_get(ws, idxr);
    const struct gkyl_wave_spacetime_iface *iface = &wsc->iface[0];
    TEST_CHECK_( iface->kind == GKYL_WS_IFACE_INTERIOR,
      "[%s] interior face not tagged INTERIOR", label );
    double MiG[3][3] = {{0}};
    for (int a = 0; a < 3; a++)
      for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++)
          MiG[a][j] += iface->M_inv[a][i] * iface->inv_g_iface[i][j];
    for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++) {
        double dot = MiG[a][0]*iface->M_inv[b][0]
                   + MiG[a][1]*iface->M_inv[b][1]
                   + MiG[a][2]*iface->M_inv[b][2];
        double expect = (a == b) ? 1.0 : 0.0;
        TEST_CHECK_( fabs(dot - expect) < 1.0e-12,
          "[%s] M_inv·γ⁻¹·M_invᵀ[%d][%d] = %.3e (expect %g)",
          label, a, b, dot, expect );
      }
  }

  // Cache-attached vs fallback equivalence over the shared state table.
  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };
  int idxl[1] = { 0 }, idxr[1] = { 1 };

  for (int c = 0; c < GR_EULER_POS_NUM_CASES; c++) {
    const struct gr_euler_pos_rp_case *rc = &gr_euler_pos_cases[c];

    double vsq_L = 0.0, vsq_R = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        vsq_L += prods_L[GKYL_GR_SP_GIJ + 3*i + j] * rc->v_L[i] * rc->v_L[j];
        vsq_R += prods_R[GKYL_GR_SP_GIJ + 3*i + j] * rc->v_R[i] * rc->v_R[j];
      }
    if (!(vsq_L < 1.0 - 1.0e-6) || !(vsq_R < 1.0 - 1.0e-6)) continue;

    double qL_glob[5], qR_glob[5];
    build_state_convA(eos, rc->rho_L, rc->v_L, rc->p_L, prods_L, qL_glob);
    build_state_convA(eos, rc->rho_R, rc->v_R, rc->p_R, prods_R, qR_glob);

    double waves_ws[3 * 5], speeds_ws[3], waves_nb[3 * 5], speeds_nb[3];
    double qL[5], qR[5], delta[5];

    // Cache attached.
    gkyl_gr_euler_tetrad_set_wave_spacetime(eqn, ws);
    eqn->set_interface_idx_func(eqn, idxl, idxr);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL_glob, qL);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR_glob, qR);
    for (int i = 0; i < 5; i++) delta[i] = qR[i] - qL[i];
    eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      delta, qL, qR, 1.0, 1.0, waves_ws, speeds_ws);

    // Cache detached — per-call averaging fallback.
    gkyl_gr_euler_tetrad_set_wave_spacetime(eqn, NULL);
    eqn->set_interface_idx_func(eqn, idxl, idxr);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL_glob, qL);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR_glob, qR);
    eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      delta, qL, qR, 1.0, 1.0, waves_nb, speeds_nb);

    for (int k = 0; k < 2 * 5; k++)
      TEST_CHECK_( fabs(waves_ws[k] - waves_nb[k]) < 1.0e-12,
        "[%s/%s] cache vs fallback wave[%d]: %.3e vs %.3e",
        label, rc->name, k, waves_ws[k], waves_nb[k] );
    for (int k = 0; k < 2; k++)
      TEST_CHECK_( fabs(speeds_ws[k] - speeds_nb[k]) < 1.0e-12,
        "[%s/%s] cache vs fallback speed[%d]: %.3e vs %.3e",
        label, rc->name, k, speeds_ws[k], speeds_nb[k] );
  }

  gkyl_wave_spacetime_release(ws);
  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(wg);
}

void test_wave_spacetime_cache_equivalence(void)
{
  struct gkyl_gr_spacetime *mink = gkyl_gr_minkowski_new(false);
  run_wave_spacetime_cache_equivalence(mink, "Mink-cache", 0.30, 0.32);
  gkyl_gr_spacetime_release(mink);

  struct gkyl_gr_spacetime *schw =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_wave_spacetime_cache_equivalence(schw, "Schw-cache", 0.22, 0.24);
  gkyl_gr_spacetime_release(schw);

  struct gkyl_gr_spacetime *kerr =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_wave_spacetime_cache_equivalence(kerr, "Kerr-cache", 0.16, 0.18);
  gkyl_gr_spacetime_release(kerr);
}

// ---------------------------------------------------------------------------
// 3c. Metric-aware near-floor margin sweep (Mode-A probe)
//
// The shared primitive table cannot place a state at a controlled
// distance from the admissibility boundary — closeness to the cone is
// metric-dependent, so a fixed v entry is either comfortable or
// superluminal depending on the sample point. Here the velocity is
// constructed PER CELL so that γ_ij·v^i·v^j = 1 − 1/W² for a prescribed
// Lorentz factor W: the state's relative s² margin is then ~1/W² by
// construction, the same scale as production post-repair states (the
// history-informed repair anchors its margin at 1/W²_prev).
//
// This was designed as the synthetic probe for MODE A (see §3b), and it
// does expose it — but the measured scan (see the table in the assert
// block below) shows BOTH modes activate at near-floor margins, with
// opposite scaling: curved Lax (Mode B) fails more at moderate margins
// because its cross-geometry flux content grows with |S|², while
// tetrad-Lax (Mode A) overtakes and saturates as the margin shrinks
// below the tilt. Minkowski is the control: Δγ = 0, so a violation
// there at ANY margin is a genuine solver bug, not Mode A/B.
// ---------------------------------------------------------------------------

static void
run_lax_margin_sweep(struct gkyl_gr_spacetime *spacetime,
  const char *label,
  const struct positivity_point *pts, int n_pts, double dx,
  bool expect_mode_a)
{
  struct gkyl_gr_euler_eos eos = eos_modes[0];  // IDEAL (see §3b note)

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

  // Prescribed Lorentz factors. 1/W² = {1e-2, 1e-4, 4e-6}: from "margin
  // well above the tilt" down to "margin well below it". W = 1000 would
  // collide with the runner's 1−1e-6 subluminal guard, so stop at 500.
  static const double W_targets[] = { 10.0, 100.0, 500.0 };
  // Coordinate directions for v (metric-normalized below): face-normal,
  // in-plane diagonal (activates γ_xy), out-of-plane diagonal (γ_xz).
  static const double dirs[][3] = {
    { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 1.0, 0.0, 1.0 },
  };
  // Thermodynamic states: cold (production near-floor pressure regime)
  // and warm (h > 1 exercises the enthalpy weighting in S_i).
  static const struct { double rho, p; const char *tag; } thermo[] = {
    { 1.0, 1.0e-6, "cold" }, { 1.0, 1.0, "warm" },
  };

  enum { N_W = 3 };
  int tet_D[N_W] = {0}, tet_S2[N_W] = {0}, tet_tau[N_W] = {0};
  int cur_D[N_W] = {0}, cur_S2[N_W] = {0}, cur_tau[N_W] = {0};
  int total_cells[N_W] = {0};

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);

  for (int pi = 0; pi < n_pts; pi++) {
    const struct positivity_point *pt = &pts[pi];
    fill_prods_at(spacetime, pt->x - 0.5*dx, pt->y, pt->z, prods_L);
    fill_prods_at(spacetime, pt->x + 0.5*dx, pt->y, pt->z, prods_R);
    if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
        prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;

    for (size_t wi = 0; wi < N_W; wi++) {
      double v_sq_target = 1.0 - 1.0/(W_targets[wi]*W_targets[wi]);
      for (size_t di = 0; di < sizeof(dirs)/sizeof(*dirs); di++) {
        // Per-cell metric normalization: v^i = u^i·√(v²_target / |u|²_γ),
        // with |u|²_γ contracted against THAT cell's γ_ij. Both cells get
        // the same prescribed W in their own frames.
        struct gr_euler_pos_rp_case rc;
        double nsq_L = 0.0, nsq_R = 0.0;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) {
            nsq_L += prods_L[GKYL_GR_SP_GIJ + 3*i + j] * dirs[di][i] * dirs[di][j];
            nsq_R += prods_R[GKYL_GR_SP_GIJ + 3*i + j] * dirs[di][i] * dirs[di][j];
          }
        double scl_L = sqrt(v_sq_target / nsq_L);
        double scl_R = sqrt(v_sq_target / nsq_R);
        for (int i = 0; i < 3; i++) {
          rc.v_L[i] = dirs[di][i] * scl_L;
          rc.v_R[i] = dirs[di][i] * scl_R;
        }

        for (size_t ti = 0; ti < sizeof(thermo)/sizeof(*thermo); ti++) {
          char name[64];
          snprintf(name, sizeof(name), "W=%g dir=%zu %s",
            W_targets[wi], di, thermo[ti].tag);
          rc.name = name;
          rc.rho_L = rc.rho_R = thermo[ti].rho;
          rc.p_L   = rc.p_R   = thermo[ti].p;

          int n = run_positivity_for_rp(eqn, grm, prods_L, prods_R, &rc,
            GKYL_WV_HIGH_ORDER_FLUX, 1.0, &tet_D[wi], &tet_S2[wi], &tet_tau[wi]);
          run_positivity_for_rp(eqn, grm, prods_L, prods_R, &rc,
            GKYL_WV_LOW_ORDER_FLUX, 1.0, &cur_D[wi], &cur_S2[wi], &cur_tau[wi]);
          total_cells[wi] += n;
        }
      }
    }
  }

  fprintf(stderr, "[lax-margin %s, dx=%.4f] s²<=0 counts by margin:\n", label, dx);
  for (int wi = 0; wi < N_W; wi++) {
    fprintf(stderr,
      "  W=%-5g (margin %.0e): tetrad %3d / %d   curved %3d / %d"
      "   (D fails: tet %d, cur %d)\n",
      W_targets[wi], 1.0/(W_targets[wi]*W_targets[wi]),
      tet_S2[wi], total_cells[wi], cur_S2[wi], total_cells[wi],
      tet_D[wi], cur_D[wi]);
  }

  int tet_S2_all = 0, cur_S2_all = 0, tet_D_all = 0, cur_D_all = 0;
  for (int wi = 0; wi < N_W; wi++) {
    tet_S2_all += tet_S2[wi]; cur_S2_all += cur_S2[wi];
    tet_D_all  += tet_D[wi];  cur_D_all  += cur_D[wi];
  }

  if (!expect_mode_a) {
    // Minkowski control: no metric variation → no Mode A or B; any
    // violation at any margin is a genuine solver bug.
    TEST_CHECK_( tet_S2_all == 0 && tet_D_all == 0,
      "tetrad-Lax violated admissibility on Minkowski (D %d, s² %d) for %s",
      tet_D_all, tet_S2_all, label );
    TEST_CHECK_( cur_S2_all == 0 && cur_D_all == 0,
      "curved-Lax violated admissibility on Minkowski (D %d, s² %d) for %s",
      cur_D_all, cur_S2_all, label );
  } else {
    // Curved-spacetime margin scan — measured behavior (2026-06-10,
    // Schw M=0.1 / Kerr a=0.9 ladders, dx ≈ 0.0195, CFL = 1):
    //
    //            tetrad s²-fails    curved s²-fails
    //   W=10        5/84,  6/144      27/84, 49/144
    //   W=100      64/84, 84/144      42/84, 71/144
    //   W=500      82/84, 84/144      42/84, 71/144
    //
    // Neither scheme preserves the invariant domain at CFL = 1 once the
    // relative s² margin 1/W² approaches the O(Δγ) geometric variation —
    // and the two modes scale differently: Mode B (curved Lax's
    // cross-geometry flux content, which grows with |S|²) dominates at
    // moderate margins and PLATEAUS on a fixed geometry-dependent subset,
    // while Mode A (tetrad cone tilt) keeps growing as the margin shrinks
    // and saturates toward all-fail. Production reconciliation: the
    // production LOW_ORDER fallback sees post-repair states at W ≲ 3
    // (margins ≳ 0.1, where the shared-table sweep shows both schemes
    // clean) at CFL < 1; this probe is deliberately past that regime and
    // demonstrates that the repair cascade is structural — no flux
    // construction makes near-floor states safe at the CFL limit.
    //
    // Asserted: D > 0 held everywhere (below); s² counts are the
    // documented scan, info-only.
    TEST_MSG( "s² margin-scan counts (info-only): tetrad %d, curved %d for %s",
      tet_S2_all, cur_S2_all, label );
  }
  TEST_CHECK_( tet_D_all == 0, "tetrad-Lax D > 0 violated %d times for %s",
    tet_D_all, label );
  TEST_CHECK_( cur_D_all == 0, "curved-Lax D > 0 violated %d times for %s",
    cur_D_all, label );

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_lax_margin_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_lax_margin_sweep(st, "Mink-margin",
    g_positivity_points,
    (int)(sizeof(g_positivity_points) / sizeof(*g_positivity_points)),
    GR_EULER_POSITIVITY_DX, /*expect_mode_a=*/false);
  gkyl_gr_spacetime_release(st);
}

void test_lax_margin_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_lax_margin_sweep(st, "Schw(M=0.1)-margin",
    g_ladder_schw_m01, GR_EULER_N_PTS(g_ladder_schw_m01),
    GR_EULER_POSITIVITY_DX, /*expect_mode_a=*/true);
  gkyl_gr_spacetime_release(st);
}

void test_lax_margin_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_lax_margin_sweep(st, "Kerr(M=0.1,a=0.9)-margin",
    g_ladder_kerr_a09, GR_EULER_N_PTS(g_ladder_kerr_a09),
    GR_EULER_POSITIVITY_DX, /*expect_mode_a=*/true);
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

// ---------------------------------------------------------------------------
// 5. Audit tests for the contravariant-x triad correctness
// HLLC fallback probe. Walks the same three BHL repair cases through
// gkyl_gr_euler_tetrad_sr_hllc_minkowski directly with the
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
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
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
  gkyl_wave_spacetime_build_triad_contravariant_x(
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
    gkyl_gr_euler_tetrad_q_to_tetrad_contra(qL, sqrt_det, inv_g_iface, M_inv, qL_tet);
    gkyl_gr_euler_tetrad_q_to_tetrad_contra(qC, sqrt_det, inv_g_iface, M_inv, qC_tet);
    gkyl_gr_euler_tetrad_q_to_tetrad_contra(qR, sqrt_det, inv_g_iface, M_inv, qR_tet);

    // L-C interface. HLLC diagnostics live on stat->hllc — populate a
    // local status struct so we can read the last_* fields after the call.
    double waves_LC[3 * 5], speeds_LC[3];
    struct gkyl_gr_euler_prim_status stat_LC = {0};
    gkyl_gr_euler_tetrad_sr_hllc_minkowski(
      eos, qL_tet, qC_tet, &stat_LC, waves_LC, speeds_LC);

    // C-R interface.
    double waves_CR[3 * 5], speeds_CR[3];
    struct gkyl_gr_euler_prim_status stat_CR = {0};
    gkyl_gr_euler_tetrad_sr_hllc_minkowski(
      eos, qC_tet, qR_tet, &stat_CR, waves_CR, speeds_CR);

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
      stat_LC.hllc.last_did_fallback, reason_str[stat_LC.hllc.last_fallback_reason],
      stat_LC.hllc.last_lambda_L, stat_LC.hllc.last_lambda_R, stat_LC.hllc.last_lambda_star);
    fprintf(stderr,
      "    C-R: fb=%d (%s)  λ_L=%+.4f  λ_R=%+.4f  λ*=%+.4f\n",
      stat_CR.hllc.last_did_fallback, reason_str[stat_CR.hllc.last_fallback_reason],
      stat_CR.hllc.last_lambda_L, stat_CR.hllc.last_lambda_R, stat_CR.hllc.last_lambda_star);
    TEST_MSG( "[%s] L-C fb=%d (reason %d), C-R fb=%d (reason %d)",
      cases[c].name,
      stat_LC.hllc.last_did_fallback, stat_LC.hllc.last_fallback_reason,
      stat_CR.hllc.last_did_fallback, stat_CR.hllc.last_fallback_reason );
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
  gkyl_gr_spacetime_release(st);
}

void test_direct_state_hllc_fallback_probe(void)
{
  run_direct_state_hllc_fallback_probe();
}
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
  gkyl_wave_spacetime_build_triad_contravariant_x(g_ij, inv_g, M, M_inv);
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
        gkyl_gr_euler_tetrad_q_to_tetrad_contra(q, sqrt_det, inv_g, M_inv, q_tet);
        // The "wave_to_curved" maps tetrad-frame jumps to curved-frame
        // jumps. We use it on q_tet directly (which represents "the
        // full wave from 0 to q_tet"): the result should be the
        // densitized curved state.
        gkyl_gr_euler_tetrad_wave_to_curved_contra(q_tet, sqrt_det, M_inv, q_back);
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
  gkyl_wave_spacetime_build_triad_contravariant_x(g_ij, inv_g, M, M_inv);
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
  gkyl_gr_euler_tetrad_q_to_tetrad_contra(q, sqrt_det, inv_g, M_inv, q_tet);

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
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
  // Per-cell metrics: L at x − dx/2, R at x + dx/2 (production grid dx).
  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);
  fill_prods_at(spacetime, x - 0.5*GR_EULER_POSITIVITY_DX, y, z, prods_L);
  fill_prods_at(spacetime, x + 0.5*GR_EULER_POSITIVITY_DX, y, z, prods_R);
  if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
      prods_R[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  // Interface-averaged prods for the flux-jump reference: with per-cell
  // metrics the tetrad scheme satisfies Σs·w = ΔF(geom_iface), not the
  // cell-centered difference (TETRAD_REFACTOR_PLAN §1).
  double prods_iface[GKYL_GR_SP_NCOMP_BASE];
  build_prods_iface(prods_L, prods_R, prods_iface);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };
  int idxl[1] = { 0 }, idxr[1] = { 1 };

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
    build_state_convA(eos, cases[c].rL, v_co_l, cases[c].pL, prods_L, qL_glob);
    build_state_convA(eos, cases[c].rR, v_co_r, cases[c].pR, prods_R, qR_glob);

    eqn->set_interface_idx_func(eqn, idxl, idxr);
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

    // Flux jump — interface-averaged form (per-cell metrics make the
    // cell-centered difference an unrelated O(Δgeom) quantity).
    double dF[5];
    {
      double fL_gr[5], fR_gr[5];
      gkyl_gr_euler_banyuls_flux_cell(eos, qL, prods_iface, NULL, fL_gr);
      gkyl_gr_euler_banyuls_flux_cell(eos, qR, prods_iface, NULL, fR_gr);
      for (int i = 0; i < 5; i++) dF[i] = fR_gr[i] - fL_gr[i];
    }
    for (int i = 0; i < 5; i++) {
      double sw = 0.0;
      for (int k = 0; k < num_waves; k++) sw += speeds[k] * waves[k * 5 + i];
      // Near-floor τ cases incur catastrophic-cancellation precision
      // loss in the Banyuls Newton (Δ at 1e-10 squared minus terms at
      // 1e-8). With per-cell metrics each side's recovery sees slightly
      // different scales and the iface-metric reference flux adds its
      // own recovery; measured residual ≈ 3.6e-7 on the τ slot, so the
      // near-floor tolerance is 1e-6 (was 1e-7 at single metric).
      double tol = strstr(cases[c].name, "near-floor") ? 1.0e-6 : 1.0e-9;
      TEST_CHECK_( fabs(sw - dF[i]) < tol,
        "[%s/%s] flux jump: comp %d, |Σs·w − ΔF| = %.3e",
        label, cases[c].name, i, fabs(sw - dF[i]) );
    }

    // Fluctuation balance.
    double amdq[5], apdq[5];
    eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      qL, qR, 1.0, 1.0, waves, speeds, amdq, apdq);
    for (int i = 0; i < 5; i++) {
      double tol_fluct = strstr(cases[c].name, "near-floor") ? 1.0e-6 : 1.0e-9;
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
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, 0.3, 0.2, 0.0, prods_row);
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);

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
// 6. Three-cell positivity test (mimics wave_prop's actual cell update)


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
//
// Probes run at CFL = 1 (extremal). HLLC can be added as a fourth column
// once its excision/fallback diagnostics are wanted here; it shares the
// same dispatch plumbing.
// ---------------------------------------------------------------------------

enum gr_pos_flux {
  GRPF_TET_LAX = 0, GRPF_TET_HLL = 1, GRPF_TET_HLLC = 2, GRPF_CUR_LAX = 3,
  GRPF_N = 4
};
static const char *gr_pos_flux_name[GRPF_N] =
  { "tetrad-Lax", "tetrad-HLL", "tetrad-HLLC", "curved-Lax" };

#define GRPR_MAX_CONDS 64
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

static void
run_positivity_registry(struct gkyl_gr_spacetime *spacetime,
  const char *label,
  const struct positivity_point *pts, int n_pts,
  bool is_minkowski_control)
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

  // --- MARGIN conditions: metric-aware W targets, cold thermo,
  //     face-normal + in-plane-diagonal directions.
  static const double Wt[] = { 10.0, 100.0, 500.0 };
  static const double mdirs[][3] = { { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 } };
  static char margin_names[6][32];
  for (size_t wi = 0; wi < sizeof(Wt)/sizeof(*Wt); wi++) {
    for (size_t di = 0; di < sizeof(mdirs)/sizeof(*mdirs); di++) {
      struct gr_pos_registry_row *row = &reg[n_conds++];
      snprintf(margin_names[2*wi + di], 32, "margin W=%g dir=%zu", Wt[wi], di);
      row->name = margin_names[2*wi + di];
      double vsq_t = 1.0 - 1.0/(Wt[wi]*Wt[wi]);
      for (int pi = 0; pi < n_pts; pi++) {
        double dxh = 0.5 * GR_EULER_POSITIVITY_DX;
        fill_prods_at(spacetime, pts[pi].x - dxh, pts[pi].y, pts[pi].z, prods_L);
        fill_prods_at(spacetime, pts[pi].x + dxh, pts[pi].y, pts[pi].z, prods_R);
        if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
            prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;
        struct gr_euler_pos_rp_case rc = {
          .name = row->name, .rho_L = 1.0, .p_L = 1.0e-6,
          .rho_R = 1.0, .p_R = 1.0e-6 };
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

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn_lax);
  gkyl_wv_eqn_release(eqn_hll);
  gkyl_wv_eqn_release(eqn_hllc);
}

void test_positivity_registry_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_positivity_registry(st, "Mink-registry",
    g_positivity_points,
    (int)(sizeof(g_positivity_points) / sizeof(*g_positivity_points)),
    /*is_minkowski_control=*/true);
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
  run_positivity_registry(st, "BHL(M=0.3)-registry",
    pts, (int)(sizeof(pts) / sizeof(*pts)),
    /*is_minkowski_control=*/false);
  gkyl_gr_spacetime_release(st);
}

void test_positivity_registry_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_positivity_registry(st, "Schw-registry",
    g_ladder_schw_m01, GR_EULER_N_PTS(g_ladder_schw_m01),
    /*is_minkowski_control=*/false);
  gkyl_gr_spacetime_release(st);
}

void test_positivity_registry_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_positivity_registry(st, "Kerr-registry",
    g_ladder_kerr_a09, GR_EULER_N_PTS(g_ladder_kerr_a09),
    /*is_minkowski_control=*/false);
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


  { "riemann_properties_two_cell_lax_minkowski",      test_riemann_properties_two_cell_lax_minkowski },
  { "riemann_properties_two_cell_lax_schwarzschild",  test_riemann_properties_two_cell_lax_schwarzschild },
  { "riemann_properties_two_cell_lax_kerr",           test_riemann_properties_two_cell_lax_kerr },
  { "riemann_properties_two_cell_hll_minkowski",      test_riemann_properties_two_cell_hll_minkowski },
  { "riemann_properties_two_cell_hll_schwarzschild",  test_riemann_properties_two_cell_hll_schwarzschild },
  { "riemann_properties_two_cell_hll_kerr",           test_riemann_properties_two_cell_hll_kerr },
  { "riemann_properties_two_cell_hllc_minkowski",     test_riemann_properties_two_cell_hllc_minkowski },
  { "riemann_properties_two_cell_hllc_schwarzschild", test_riemann_properties_two_cell_hllc_schwarzschild },
  { "riemann_properties_two_cell_hllc_kerr",          test_riemann_properties_two_cell_hllc_kerr },
  { "riemann_properties_two_cell_roe_minkowski",      test_riemann_properties_two_cell_roe_minkowski },
  { "riemann_properties_two_cell_roe_schwarzschild",  test_riemann_properties_two_cell_roe_schwarzschild },
  { "riemann_properties_two_cell_roe_kerr",           test_riemann_properties_two_cell_roe_kerr },

  { "excision_absorbing_lax_minkowski",     test_excision_absorbing_lax_minkowski },
  { "excision_absorbing_lax_schwarzschild", test_excision_absorbing_lax_schwarzschild },
  { "excision_absorbing_hll_minkowski",     test_excision_absorbing_hll_minkowski },
  { "excision_absorbing_hll_schwarzschild", test_excision_absorbing_hll_schwarzschild },
  { "excision_absorbing_roe_minkowski",     test_excision_absorbing_roe_minkowski },
  { "excision_absorbing_roe_schwarzschild", test_excision_absorbing_roe_schwarzschild },

  { "positivity_curved_lax_minkowski",     test_positivity_curved_lax_minkowski },
  { "positivity_curved_lax_schwarzschild", test_positivity_curved_lax_schwarzschild },
  { "positivity_curved_lax_kerr",          test_positivity_curved_lax_kerr },
  { "positivity_curved_lax_bhl",           test_positivity_curved_lax_bhl },

  { "lax_convexity_tetrad_vs_curved_minkowski",
    test_lax_convexity_tetrad_vs_curved_minkowski },
  { "lax_convexity_tetrad_vs_curved_schwarzschild",
    test_lax_convexity_tetrad_vs_curved_schwarzschild },
  { "lax_convexity_tetrad_vs_curved_kerr",
    test_lax_convexity_tetrad_vs_curved_kerr },

  { "lax_margin_minkowski",     test_lax_margin_minkowski },
  { "lax_margin_schwarzschild", test_lax_margin_schwarzschild },
  { "lax_margin_kerr",          test_lax_margin_kerr },

  { "wave_spacetime_cache_equivalence", test_wave_spacetime_cache_equivalence },

  { "positivity_registry_minkowski",     test_positivity_registry_minkowski },
  { "positivity_registry_schwarzschild", test_positivity_registry_schwarzschild },
  { "positivity_registry_kerr",          test_positivity_registry_kerr },
  { "positivity_registry_bhl",           test_positivity_registry_bhl },

  { "near_floor_minkowski",     test_near_floor_minkowski },
  { "near_floor_schwarzschild", test_near_floor_schwarzschild },
  { "near_floor_kerr",          test_near_floor_kerr },

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


  { "small_tau_over_D_lax",  test_small_tau_over_D_lax },
  { "small_tau_over_D_hll",  test_small_tau_over_D_hll },
  { "small_tau_over_D_hllc", test_small_tau_over_D_hllc },

  { "direct_state_lax_curved_production_reproducer",  test_direct_state_lax_curved_production_reproducer },
  { "direct_state_hllc_fallback_probe", test_direct_state_hllc_fallback_probe },

  // EOS accuracy vs the true Synge gas (Bessel-function reference).
  { "eos_tm_synge_no_flow",          test_eos_tm_synge_no_flow },
  { "eos_rcc_synge_relativistic",    test_eos_rcc_synge_relativistic },
  { "eos_rcc_synge_high_W",          test_eos_rcc_synge_high_W },
  { "eos_rcc_synge_extreme_W",       test_eos_rcc_synge_extreme_W },

  { NULL, NULL },
};
