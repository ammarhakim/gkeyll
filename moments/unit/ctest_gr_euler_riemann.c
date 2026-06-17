// Riemann-solver property tests for wv_gr_euler_tetrad (Convention A) across
// the five solvers: curved Lax (LOW_ORDER fallback), tetrad Lax, HLL, HLLC,
// Roe. Split out of ctest_wv_gr_euler_tetrad.c; see CTEST_SPLIT_PLAN.md.
//
//   1.  Banyuls flux consistency — code flux vs analytic F^x at a point
//   2.  Two-cell Riemann properties, per-cell metrics: wave sum, flux jump
//       (interface-averaged ΔF for the tetrad HIGH_ORDER solvers), fluct
//       balance, ordering, max-speed dominance, trivial RP, sign symmetry
//   2a. Curved-Lax (LOW_ORDER) two-cell properties against its PER-CELL
//       flux-jump identity Σs·w = F(q_R;γ_R) − F(q_L;γ_L)
//   2b. Excision-boundary absorbing BC (all four solvers; HLLC routes
//       vacuum sides to its HLL fallback, reason 5)
//   3.  BHL-regime states — properties on production-relevant inputs
//   4.  wave_spacetime builder correctness
//   5.  SR HLLC kernel fallback probe (kernel-level diagnostics)

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
  struct test_ws tws = test_ws_new(&conf_range, prods, eqn);

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

  test_ws_release(&tws);
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
// The single-cell runner above uses idx={0,0}, so both sides see the
// same cell's geometry and the averaging recovers the cell-centered
// values exactly — the flux-jump check is trivial at the geometry
// level: any averaging policy passes.
//
// This variant allocates two cells, fills L at (xL,yL,zL) and R at
// (xR,yR,zR), and sets idx={0,1} so the cached per-cell rows come from
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

  struct test_ws tws = test_ws_new(&conf_range, prods, eqn);

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
    test_ws_release(&tws);
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

  test_ws_release(&tws);
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
// 2a. Curved-Lax (LOW_ORDER) two-cell Riemann properties — the fifth solver.
//
// The LOW_ORDER dispatch routes every rp_type to wave_lax_curved /
// qfluct_lax_curved, which work in curved-frame conservatives throughout
// and satisfy the PER-CELL flux-jump identity
//   Σ s_k·w_k = F(q_R; γ_R) − F(q_L; γ_L)
// (cell-centered Banyuls flux per side) — NOT the interface-averaged form
// the tetrad HIGH_ORDER solvers satisfy. This identity is what makes the
// positivity fallback a convex-combination update in A_γ
// (SESSION_NOTES_POSITIVITY_UNIFICATION.md §5) and is enforced nowhere
// else in the suite.
//
// rp_type is irrelevant on this path (the dispatch never reads sr_kernel /
// num_waves / excision_policy for LOW_ORDER), so a single RP_LAX
// instantiation covers all variants. wave_lax_curved always fills the
// ±amax pair: 2 waves regardless of the variant's HIGH_ORDER num_waves.
// ---------------------------------------------------------------------------
static void
run_riemann_properties_two_cell_curved_lax(struct gkyl_gr_spacetime *spacetime,
  const char *label, double xL, double yL, double zL,
  double xR, double yR, double zR, struct gkyl_gr_euler_eos eos)
{
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
  fill_prods_at(spacetime, xL, yL, zL, prods_L);
  fill_prods_at(spacetime, xR, yR, zR, prods_R);

  if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
      prods_R[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }

  struct test_ws tws = test_ws_new(&conf_range, prods, eqn);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  // Same per-side primitives as the tetrad two-cell runner.
  double rho_L = 1.0, p_L = 1.5;  double v_L[3] = { 0.05, 0.10, 0.15 };
  double rho_R = 0.5, p_R = 0.7;  double v_R[3] = { 0.03, 0.05, 0.08 };

  double vsq_L = 0.0, vsq_R = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      vsq_L += prods_L[GKYL_GR_SP_GIJ + 3*i + j] * v_L[i] * v_L[j];
      vsq_R += prods_R[GKYL_GR_SP_GIJ + 3*i + j] * v_R[i] * v_R[j];
    }
  if (!(vsq_L < 1.0 - 1.0e-6) || !(vsq_R < 1.0 - 1.0e-6)) {
    test_ws_release(&tws);
    gkyl_array_release(prods);
    gkyl_wv_eqn_release(eqn);
    return;
  }

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

  enum { NW_LAX = 2 };
  double waves[NW_LAX * 5], speeds[NW_LAX];
  double maxs = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    delta, qL, qR, 1.0, 1.0, waves, speeds);

  // PER-CELL ΔF: each side's Banyuls flux at its OWN cell-centered prods
  // (x-normal identity frame: the global rows are already face-local).
  double dF[5];
  banyuls_delta_flux(eos, prods_L, prods_R, qL, qR, dF);

  // (a) Wave sum: Σ w_k = Δq
  for (int i = 0; i < 5; i++) {
    double sum = 0.0;
    for (int k = 0; k < NW_LAX; k++) sum += waves[k * 5 + i];
    TEST_CHECK_( fabs(sum - delta[i]) < 1.0e-12,
      "[%s] wave sum: comp %d, |Σw − Δq| = %.3e", label, i, fabs(sum - delta[i]) );
  }

  // (b) Flux jump (PER-CELL): Σ s_k·w_k = F(q_R;γ_R) − F(q_L;γ_L)
  for (int i = 0; i < 5; i++) {
    double sw = 0.0;
    for (int k = 0; k < NW_LAX; k++) sw += speeds[k] * waves[k * 5 + i];
    TEST_CHECK_( fabs(sw - dF[i]) < 1.0e-10,
      "[%s] flux jump (per-cell): comp %d, |Σs·w − ΔF| = %.3e",
      label, i, fabs(sw - dF[i]) );
  }

  // (c) Fluctuation balance: amdq + apdq = ΔF_percell
  double amdq[5], apdq[5];
  eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    qL, qR, 1.0, 1.0, waves, speeds, amdq, apdq);
  for (int i = 0; i < 5; i++) {
    TEST_CHECK_( fabs(amdq[i] + apdq[i] - dF[i]) < 1.0e-10,
      "[%s] fluct balance (per-cell): comp %d, |amdq+apdq − ΔF| = %.3e",
      label, i, fabs(amdq[i] + apdq[i] - dF[i]) );
  }

  // (d) ±amax pair structure: s = {−amax, +amax}, maxs = amax > 0.
  TEST_CHECK_( speeds[0] <= 0.0 && speeds[1] >= 0.0,
    "[%s] speed signs: s0=%g s1=%g", label, speeds[0], speeds[1] );
  TEST_CHECK_( fabs(speeds[0] + speeds[1]) < 1.0e-14,
    "[%s] speeds not symmetric: s0=%g s1=%g", label, speeds[0], speeds[1] );
  TEST_CHECK_( maxs > 0.0 && fabs(maxs - speeds[1]) < 1.0e-14,
    "[%s] maxs=%g != +amax=%g", label, maxs, speeds[1] );

  // (e) Trivial RP — same state AND same cell. For the per-cell scheme,
  // zero waves require the geometry to match too: with q_L = q_R across
  // different metrics the wave content is exactly the Mode-B
  // cross-geometry flux term — see probe (f).
  //
  // Geometry is selected by set_interface_idx alone (cache fetch by
  // index — rotation is state-only and call-order-free since the
  // parity-contract retirement).
  int idx0[1] = { 0 };
  eqn->set_interface_idx_func(eqn, idx0, idx0);
  double qE[5];
  build_state_convA(eos, 1.0, (double[]){0.05, 0.10, 0.05}, 1.0, prods_L, qE);
  double qE_loc[5];
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qE, qE_loc);
  double dE[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  double wavesE[NW_LAX * 5], speedsE[NW_LAX];
  double maxsE = eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    dE, qE_loc, qE_loc, 1.0, 1.0, wavesE, speedsE);
  for (int k = 0; k < NW_LAX; k++) {
    for (int i = 0; i < 5; i++) {
      TEST_CHECK_( fabs(wavesE[k * 5 + i]) < 1.0e-12,
        "[%s] trivial RP: wave[%d][%d] = %.3e", label, k, i, wavesE[k * 5 + i] );
    }
  }
  double amdqE[5], apdqE[5];
  eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    qE_loc, qE_loc, 1.0, 1.0, wavesE, speedsE, amdqE, apdqE);
  for (int i = 0; i < 5; i++) {
    TEST_CHECK_( fabs(amdqE[i]) < 1.0e-12 && fabs(apdqE[i]) < 1.0e-12,
      "[%s] trivial RP: |amdq[%d]|=%.3e |apdq[%d]|=%.3e",
      label, i, fabs(amdqE[i]), i, fabs(apdqE[i]) );
  }
  TEST_CHECK( isfinite(maxsE) );

  // (f) Mode-B identity probe: same state, DIFFERENT cells. Δq = 0 but
  // ΔF = F(q;γ_R) − F(q;γ_L) ≠ 0 — the cross-geometry flux content that
  // drives curved-Lax repair near horizons. The per-cell identity must
  // hold exactly on it (the waves are genuinely nonzero; that is the
  // scheme, not a bug).
  eqn->set_interface_idx_func(eqn, idxl, idxr);
  double dB[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  double wavesB[NW_LAX * 5], speedsB[NW_LAX];
  eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    dB, qE_loc, qE_loc, 1.0, 1.0, wavesB, speedsB);
  double dFB[5];
  banyuls_delta_flux(eos, prods_L, prods_R, qE_loc, qE_loc, dFB);
  for (int i = 0; i < 5; i++) {
    double sumB = wavesB[0 * 5 + i] + wavesB[1 * 5 + i];
    double swB = speedsB[0] * wavesB[0 * 5 + i] + speedsB[1] * wavesB[1 * 5 + i];
    TEST_CHECK_( fabs(sumB) < 1.0e-12,
      "[%s] Mode-B probe: wave sum comp %d = %.3e (Δq = 0)", label, i, sumB );
    TEST_CHECK_( fabs(swB - dFB[i]) < 1.0e-10,
      "[%s] Mode-B probe: comp %d, |Σs·w − ΔF_percell| = %.3e",
      label, i, fabs(swB - dFB[i]) );
  }

  // (g) L↔R sign symmetry with swapped interface indices: ΔF_percell is
  // antisymmetric under the (states, cells) swap, so
  //   amdq(qL,qR) + apdq(qL,qR) = −[amdq(qR,qL) + apdq(qR,qL)].
  int idxl_swap[1] = { 1 }, idxr_swap[1] = { 0 };
  eqn->set_interface_idx_func(eqn, idxl_swap, idxr_swap);
  double delta_swap[5];
  for (int i = 0; i < 5; i++) delta_swap[i] = qL[i] - qR[i];
  double waves_swap[NW_LAX * 5], speeds_swap[NW_LAX];
  eqn->waves_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    delta_swap, qR, qL, 1.0, 1.0, waves_swap, speeds_swap);
  double amdq_swap[5], apdq_swap[5];
  eqn->qfluct_func(eqn, GKYL_WV_LOW_ORDER_FLUX,
    qR, qL, 1.0, 1.0, waves_swap, speeds_swap, amdq_swap, apdq_swap);
  for (int i = 0; i < 5; i++) {
    double lhs = amdq[i] + apdq[i];
    double rhs = -(amdq_swap[i] + apdq_swap[i]);
    TEST_CHECK_( fabs(lhs - rhs) < 1.0e-10,
      "[%s] sign symmetry: comp %d, |ΔF + ΔF_swap| = %.3e",
      label, i, fabs(lhs - rhs) );
  }

  test_ws_release(&tws);
  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

// Same position pairs as the tetrad two-cell tests; full EOS coverage.
void test_riemann_properties_two_cell_curved_lax_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_riemann_properties_two_cell_curved_lax(st, "Mink-CurvedLax-2cell",
      0.30, 0.0, 0.0,  0.32, 0.0, 0.0, eos_modes[ei]);
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_curved_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_two_cell_curved_lax(st, "Schw-CurvedLax-2cell",
      0.22, 0.0, 0.0,  0.24, 0.0, 0.0, eos_modes[ei]);
    run_riemann_properties_two_cell_curved_lax(st, "Schw-CurvedLax-2cell",
      0.30, 0.20, 0.0,  0.32, 0.20, 0.0, eos_modes[ei]);
    run_riemann_properties_two_cell_curved_lax(st, "Schw-CurvedLax-2cell",
      0.38, 0.40, 0.0,  0.40, 0.40, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}
void test_riemann_properties_two_cell_curved_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_riemann_properties_two_cell_curved_lax(st, "Kerr0.9-CurvedLax-2cell",
      0.16, 0.0, 0.0,  0.18, 0.0, 0.0, eos_modes[ei]);
    run_riemann_properties_two_cell_curved_lax(st, "Kerr0.9-CurvedLax-2cell",
      0.20, 0.15, 0.0,  0.22, 0.15, 0.0, eos_modes[ei]);
    run_riemann_properties_two_cell_curved_lax(st, "Kerr0.9-CurvedLax-2cell",
      0.18, 0.0, 0.10,  0.20, 0.0, 0.10, eos_modes[ei]);
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

  // Cache built AFTER the excision-flag mock so the cached rows carry it.
  struct test_ws tws = test_ws_new(&conf_range, prods, eqn);

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
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_active_glob,  qL_loc);
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_excised_glob, qR_loc);

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
  gkyl_gr_euler_banyuls_flux_cell(eos, qL_loc, prods_active, NULL, fL_gr);
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
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_excised_glob, qL_loc);
  eqn->rotate_to_local_func(eqn, tau1v, tau2v, norm, q_active_glob,  qR_loc);

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
  gkyl_gr_euler_banyuls_flux_cell(eos, qR_loc, prods_active, NULL, fR_gr);
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

  test_ws_release(&tws);
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

// HLLC routes vacuum sides to its in-kernel HLL fallback (reason 5) and
// runs under ZERO_VACUUM like the others since the absorbing-BC port
// (HLLC_AUDIT_PLAN.md Phase 1).
void test_excision_absorbing_hllc_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_excision_absorbing_for_rp(st, "Mink", WV_GR_EULER_TETRAD_RP_HLLC, 3, "HLLC",
    0.3, 0.0, 0.0);
  gkyl_gr_spacetime_release(st);
}
void test_excision_absorbing_hllc_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_excision_absorbing_for_rp(st, "Schw", WV_GR_EULER_TETRAD_RP_HLLC, 3, "HLLC",
    0.3, 0.2, 0.0);
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 3d. wave_spacetime builder correctness
//
// The cache is the ONLY geometry path (the per-call fallback retired
// with the rotation-parity contract — WAVE_SPACETIME_PARITY_PLAN.md;
// cache == fallback was validated by the historical version of this
// test before removal). Checks: the bottom-edge storage contract, the
// cached-triad invariant M_inv·γ⁻¹·M_invᵀ = I, the per-(cell, dir)
// face-local rows against a hand rotation, and the full-row rotation
// helper on a synthetic non-identity frame (1D wave_geom frames are
// identity, which would not exercise the tensor blocks).
// ---------------------------------------------------------------------------

#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_spacetime.h>

static void
run_wave_spacetime_builder(struct gkyl_gr_spacetime *spacetime,
  const char *label, double xL, double xR)
{
  struct gkyl_gr_euler_eos eos = eos_modes[0];

  // Range {1,2} so gkyl's cell-center convention puts cell 1 at xL and
  // cell 2 at xR — the excision face-eval check below needs grid
  // positions aligned with the prods sample points.
  struct gkyl_range conf_range;
  int lower[1] = { 1 }, upper[1] = { 2 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  double dx = xR - xL;
  double grid_lower[1] = { xL - 0.5 * dx };
  double grid_upper[1] = { xR + 0.5 * dx };
  int cells[1] = { 2 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, grid_lower, grid_upper, cells);
  struct gkyl_wave_geom *wg =
    gkyl_wave_geom_new(&grid, &conf_range, NULL, NULL, false);

  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);
  fill_prods_at(spacetime, xL, 0.0, 0.0, prods_L);
  fill_prods_at(spacetime, xR, 0.0, 0.0, prods_R);
  if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
      prods_R[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods);
    gkyl_wave_geom_release(wg);
    return;
  }
  (void)eos;

  struct gkyl_wave_spacetime *ws =
    gkyl_wave_spacetime_new(&grid, &conf_range, wg, spacetime, prods, 0.0,
      /*use_gpu=*/false);

  // Storage contract at the bottom edge of the build range: cell 0 has
  // no lower neighbor, so its iface[0] must be the one-sided fill from
  // its own prods — VALID (sqrt_det > 0), tagged INTERIOR, never a
  // zero-filled sentinel (see gkyl_wave_spacetime.h).
  {
    int idx0[1] = { 1 };
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
    int idxr[1] = { 2 };
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

  // Per-(cell, dir) face-local rows: cached row == hand rotation of the
  // cell's global row with the same wave_geom frames.
  for (int ci = 0; ci < 2; ci++) {
    int idx[1] = { ci + 1 };
    const struct gkyl_wave_spacetime_cell *wsc =
      gkyl_wave_spacetime_get(ws, idx);
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idx);
    const double *row = (ci == 0) ? prods_L : prods_R;
    double expect[GKYL_GR_SP_NCOMP_BASE];
    gkyl_wave_spacetime_rotate_prods_row(cg->tau1[0], cg->tau2[0],
      cg->norm[0], row, expect);
    for (int k = 0; k < GKYL_GR_SP_NCOMP_BASE; k++)
      TEST_CHECK_( wsc->cell_prods_local[0][k] == expect[k],
        "[%s] cell %d face-local row comp %d: cached %.17e vs %.17e",
        label, ci, k, wsc->cell_prods_local[0][k], expect[k] );
  }

  // Full-row rotation helper on a synthetic y-normal frame (pure
  // permutation — spot-check one component per tensor rank).
  {
    double nrm[3] = { 0.0, 1.0, 0.0 };
    double t1v[3] = { -1.0, 0.0, 0.0 };
    double t2v[3] = { 0.0, 0.0, 1.0 };
    double rot[GKYL_GR_SP_NCOMP_BASE];
    gkyl_wave_spacetime_rotate_prods_row(t1v, t2v, nrm, prods_L, rot);
    TEST_CHECK_( fabs(rot[GKYL_GR_SP_SHIFT + 0]
        - prods_L[GKYL_GR_SP_SHIFT + 1]) < 1.0e-15,
      "[%s] y-frame rank-1: shift_n != β_y", label );
    TEST_CHECK_( fabs(rot[GKYL_GR_SP_GIJ + 0]
        - prods_L[GKYL_GR_SP_GIJ + 4]) < 1.0e-15,
      "[%s] y-frame rank-2: γ_nn != γ_yy", label );
    TEST_CHECK_( fabs(rot[GKYL_GR_SP_DGIJ + 0]
        - prods_L[GKYL_GR_SP_DGIJ + 9*1 + 3*1 + 1]) < 1.0e-15,
      "[%s] y-frame rank-3: (∂γ)_nnn != ∂_y γ_yy", label );
  }

  // Excision-face geometry: with the analytic spacetime available, the
  // builder evaluates the metric AT THE FACE CENTROID. Mock cell 1 (at
  // xL) excised; the face between the cells sits at (xL+xR)/2 and its
  // iface entry must match the analytic values there, not cell 2's
  // center values.
  {
    double exc_save = prods_L[GKYL_GR_SP_EXCISION];
    prods_L[GKYL_GR_SP_EXCISION] = -1.0;
    struct gkyl_wave_spacetime *ws_ex =
      gkyl_wave_spacetime_new(&grid, &conf_range, wg, spacetime, prods, 0.0,
        /*use_gpu=*/false);
    int idxr[1] = { 2 };
    const struct gkyl_wave_spacetime_iface *ifx =
      &gkyl_wave_spacetime_get(ws_ex, idxr)->iface[0];
    TEST_CHECK_( ifx->kind == GKYL_WS_IFACE_EXCISION,
      "[%s] excision face not tagged EXCISION", label );
    double pF[GKYL_GR_SP_NCOMP_BASE];
    fill_prods_at(spacetime, 0.5 * (xL + xR), 0.0, 0.0, pF);
    TEST_CHECK_( fabs(ifx->alpha - pF[GKYL_GR_SP_LAPSE]) < 1.0e-14,
      "[%s] excision face α = %.17e != analytic %.17e at face centroid",
      label, ifx->alpha, pF[GKYL_GR_SP_LAPSE] );
    TEST_CHECK_( fabs(ifx->shift_n - pF[GKYL_GR_SP_SHIFT + 0]) < 1.0e-14,
      "[%s] excision face β_n != analytic at face centroid", label );
    TEST_CHECK_( fabs(ifx->sqrt_det_iface
        - sqrt(pF[GKYL_GR_SP_SPATIAL_DET])) < 1.0e-12,
      "[%s] excision face √det = %.17e != analytic %.17e",
      label, ifx->sqrt_det_iface, sqrt(pF[GKYL_GR_SP_SPATIAL_DET]) );
    gkyl_wave_spacetime_release(ws_ex);
    prods_L[GKYL_GR_SP_EXCISION] = exc_save;
  }

  gkyl_wave_spacetime_release(ws);
  gkyl_array_release(prods);
  gkyl_wave_geom_release(wg);
}

void test_wave_spacetime_builder(void)
{
  struct gkyl_gr_spacetime *mink = gkyl_gr_minkowski_new(false);
  run_wave_spacetime_builder(mink, "Mink-cache", 0.30, 0.32);
  gkyl_gr_spacetime_release(mink);

  struct gkyl_gr_spacetime *schw =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_wave_spacetime_builder(schw, "Schw-cache", 0.22, 0.24);
  gkyl_gr_spacetime_release(schw);

  struct gkyl_gr_spacetime *kerr =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_wave_spacetime_builder(kerr, "Kerr-cache", 0.16, 0.18);
  gkyl_gr_spacetime_release(kerr);
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
  struct test_ws tws = test_ws_new(&conf_range, prods, eqn);

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

  test_ws_release(&tws);
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

// ---------------------------------------------------------------------------
// 5. SR HLLC kernel fallback probe
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

    // fallback_reason key (see gkyl_wv_gr_euler_prim_priv.h):
    //   1 = λ* not finite (sqrt of negative discriminant; should be
    //       caught by the disc<0 → 0 clamp inside, so rare)
    //   2 = |λ_L − λ*| < tol (would blow up 1/(λ_L − λ*) in U_L*)
    //   3 = |λ_R − λ*| < tol (would blow up 1/(λ_R − λ*) in U_R*)
    //   4 = vacuum side, 5 = star-state admissibility guard
    // λ* outside [λ_L, λ_R] (e.g. supersonic flow) is NOT a fallback
    // trigger — the wave decomposition is still conservative there.
    static const char *reason_str[] = {
      "—",                  // 0
      "λ* not finite",      // 1
      "|λ_L−λ*| < tol",     // 2
      "|λ_R−λ*| < tol",     // 3
      "vacuum side",        // 4 (excision absorbing BC)
      "star inadmissible",  // 5 (star-state guard)
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

  // The fifth solver: LOW_ORDER curved Lax against its per-cell identity.
  { "riemann_properties_two_cell_curved_lax_minkowski",
    test_riemann_properties_two_cell_curved_lax_minkowski },
  { "riemann_properties_two_cell_curved_lax_schwarzschild",
    test_riemann_properties_two_cell_curved_lax_schwarzschild },
  { "riemann_properties_two_cell_curved_lax_kerr",
    test_riemann_properties_two_cell_curved_lax_kerr },

  { "excision_absorbing_lax_minkowski",     test_excision_absorbing_lax_minkowski },
  { "excision_absorbing_lax_schwarzschild", test_excision_absorbing_lax_schwarzschild },
  { "excision_absorbing_hll_minkowski",     test_excision_absorbing_hll_minkowski },
  { "excision_absorbing_hll_schwarzschild", test_excision_absorbing_hll_schwarzschild },
  { "excision_absorbing_roe_minkowski",     test_excision_absorbing_roe_minkowski },
  { "excision_absorbing_roe_schwarzschild", test_excision_absorbing_roe_schwarzschild },
  { "excision_absorbing_hllc_minkowski",     test_excision_absorbing_hllc_minkowski },
  { "excision_absorbing_hllc_schwarzschild", test_excision_absorbing_hllc_schwarzschild },

  { "bhl_regime_lax_schwarzschild",  test_bhl_regime_lax_schwarzschild },
  { "bhl_regime_hll_schwarzschild",  test_bhl_regime_hll_schwarzschild },
  { "bhl_regime_lax_kerr",           test_bhl_regime_lax_kerr },
  { "bhl_regime_hll_kerr",           test_bhl_regime_hll_kerr },
  { "bhl_regime_hllc_schwarzschild", test_bhl_regime_hllc_schwarzschild },
  { "bhl_regime_hllc_kerr",          test_bhl_regime_hllc_kerr },

  { "wave_spacetime_builder", test_wave_spacetime_builder },

  { "direct_state_hllc_fallback_probe", test_direct_state_hllc_fallback_probe },

  { NULL, NULL },
};
