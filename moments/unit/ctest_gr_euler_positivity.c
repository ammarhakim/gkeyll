// Invariant-domain tests for wv_gr_euler_tetrad (Convention A).
//
// The invariant domain is the p≥0 cone (SESSION_NOTES_P_GE_0_CONE.md):
//   D > 0,  τ ≥ 0,  s² = (D+τ)² − γ^{ij}S_iS_j > D²
// — the existence set for a p ≥ 0 primitive inversion. Production's
// update is the COMPOSITE solver step → admissibility check → cascade
// repair, and that composite is what these tests pin:
//
//   - D > 0 after the raw step: hard zero, every metric, every flux.
//   - Minkowski: the raw step preserves the cone EXACTLY (Δγ = 0 makes
//     the convex-combination argument exact — the cone is a norm cone,
//     hence convex). Any crossing on flat space is a solver bug.
//   - Curved metrics: the raw step may graze the cone (O(Δγ) flux
//     mismatch on cold states — the production attractor regime), and
//     the contract is the composite one, asserted per crossing:
//       (a) the cascade repair lands the state strictly inside, and
//       (b) the repair removes no more momentum than the state's actual
//           overshoot plus the landing margin: clip ≤ δ + margin, which
//           is an algebraic identity of the rescale — a violation means
//           the repair AMPLIFIES (the W_prev-anchor bug class).
//     Crossing counts and depths are printed as diagnostics (they read
//     directly against production clip histograms — same cone-budget
//     normalization) but are not contracts: crossing frequency is a
//     property of how hard the flow rides the cone, not of correctness.
//
// Test families:
//   1. positivity_curved_lax_* — LOW_ORDER curved-Lax sweep over the
//      shared PRIM table + the cone-marginal table (production dx,
//      CFL = 0.9), per EOS mode.
//   2. positivity_registry_*  — the (condition × flux) matrix over
//      tetrad-Lax / tetrad-HLL / tetrad-HLLC / curved-Lax at CFL = 1
//      (extremal), production dx + a coarse-dx near-horizon ladder.
//      The printed matrix is the retained record of which conditions
//      cross for which flux; the asserts are the composite contract.

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

#include "ctest_gr_euler_test_util.h"
#include "ctest_gr_euler_pos_states.h"

// Production landing margin for the cascade repair's S-rescale — must
// match the value repair_state passes in wv_gr_euler_tetrad.c. The
// composite contract below pins repair behavior at this margin.
#define GR_EULER_POS_REPAIR_MARGIN (1.0e-6)

// ---------------------------------------------------------------------------
// Composite-update contract: apply the production cascade repair to a
// cone-violating post-step state and assert
//   (a) the repaired state is strictly admissible, and
//   (b) on the pure S-rescale path, the momentum clip does not exceed
//       the state's overshoot plus the landing margin. Algebraically
//       clip = (δ + m)/(1 + δ) ≤ δ + m, so any excess means the repair
//       target itself is wrong (e.g. the retired W_prev-anchored margin
//       cut |S|² by 1 − 1/W_prev² regardless of δ).
// δ and the budget are measured pre-repair; on the rescale-only path D
// and τ are untouched, so the bound is exact up to roundoff.
// ---------------------------------------------------------------------------
static void
check_composite_repair(const double inv_g[3][3], const double qn_und[5],
  const char *cond_name, const char *who)
{
  double D = qn_und[0], Sx = qn_und[1], Sy = qn_und[2], Sz = qn_und[3];
  double tau = qn_und[4];
  double msq_pre = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
  double budget_pre = (D + tau) * (D + tau) - D * D;

  unsigned int fixed = gkyl_gr_euler_repair_admissibility_cascade(
    inv_g, GR_EULER_POS_REPAIR_MARGIN, &D, &Sx, &Sy, &Sz, &tau);

  TEST_CHECK_( gkyl_gr_euler_check_admissibility(D, Sx, Sy, Sz, tau, inv_g)
      == GR_EULER_ADM_OK,
    "'%s' cell %s: cascade repair left the state inadmissible",
    cond_name, who );

  if (fixed == GR_EULER_REPAIR_S2 && budget_pre > 0.0 && msq_pre > 0.0) {
    double msq_post = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
    double clip = 1.0 - msq_post / msq_pre;
    double overshoot = (msq_pre - budget_pre) / budget_pre;
    TEST_CHECK_( clip <= overshoot + GR_EULER_POS_REPAIR_MARGIN + 1.0e-12,
      "'%s' cell %s: repair clip %.3e exceeds overshoot %.3e + margin — "
      "repair amplifies", cond_name, who, clip, overshoot );
  }
}

// ---------------------------------------------------------------------------
// Per-cell-metric positivity probe. prods_L and prods_R are the
// spacetime products rows for the LEFT and RIGHT cells (production
// wave_prop exposes adjacent cells' metrics independently through the
// wave_spacetime cache). Each cell's conservatives are densitized with
// its OWN √γ and admissibility-checked in its OWN inverse spatial
// metric, matching what wave_prop does post-POSITIVITY_SWEEP.
//
// `ftype` selects which path this probes: GKYL_WV_LOW_ORDER_FLUX routes
// to the direct curved-frame Lax fallback (wave_lax_curved);
// GKYL_WV_HIGH_ORDER_FLUX routes to the equation's tetrad-first solver.
//
// Every cone crossing runs the composite-repair contract above. Counts
// and the max crossing depth (cone-budget units, same normalization as
// the production clip stats) are accumulated for the caller's
// diagnostics and for the Minkowski exactness assert.
//
// Returns the number of cell-updates actually performed (2, or 0 when
// the case is skipped because it is unbuildable at this pair of metric
// samples), so callers can report honest percentages.
// ---------------------------------------------------------------------------
static int
run_positivity_for_rp(struct gkyl_wv_eqn *eqn,
  const double *prods_L, const double *prods_R,
  const struct gr_euler_pos_rp_case *rc,
  enum gkyl_wv_flux_type ftype, double cfl,
  int *D_violations, int *S2_violations, int *tau_violations,
  double *max_cone_depth)
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

  // Sanity: each cell admissible in its own metric. Primitive-built
  // states with p > 0 are strictly cone-interior by the sign identity
  // sign(s² − D²) = sign(p); a skip here means an FP-degenerate build
  // (e.g. θ below double precision at this metric), not coverage loss.
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

  // Per-cell metric: idxl=0 → prods_L, idxr=1 → prods_R, fetched from
  // the wave_spacetime cache by (index, dir); with norm=+x̂ the state
  // rotation is a no-op on q components.
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

  // CFL-scaled update against the interface speed (cfl = 1 is the
  // extremal convexity probe; the sweep passes the production-like 0.9).
  if (!(maxs > 0.0)) return 0;
  double dt_dx = cfl / maxs;

  double qL_new[5], qR_new[5];
  for (int i = 0; i < 5; i++) {
    qL_new[i] = qL_loc[i] - dt_dx * amdq[i];
    qR_new[i] = qR_loc[i] - dt_dx * apdq[i];
  }

  // Check each updated cell in its OWN metric; run the composite
  // contract on every cone crossing.
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
    if (!S2_ok)  {
      (*S2_violations)++;
      double ms = gkyl_gr_euler_mom_sq(inv_g_s, qn_und[1], qn_und[2], qn_und[3]);
      double s2 = (qn_und[0]+qn_und[4])*(qn_und[0]+qn_und[4]) - ms;
      double D2 = qn_und[0]*qn_und[0];
      double budget = (qn_und[0]+qn_und[4])*(qn_und[0]+qn_und[4]) - D2;
      double depth = (budget > 0.0) ? (D2 - s2) / budget : 1.0;
      if (max_cone_depth && depth > *max_cone_depth)
        *max_cone_depth = depth;
      check_composite_repair(inv_g_s, qn_und, rc->name, who);
    }
    if (!tau_ok) (*tau_violations)++;
  }
  return 2;
}

// ---------------------------------------------------------------------------
// 1. LOW_ORDER curved-Lax positivity sweep. This sweep is rp-independent
// BY CONSTRUCTION — the LOW_ORDER branch of wave_tetrad_dispatch routes
// to wave_lax_curved / qfluct_lax_curved for every rp_type.
// ---------------------------------------------------------------------------
static void
run_positivity_sweep(struct gkyl_gr_spacetime *spacetime,
  struct gkyl_gr_euler_eos eos,
  const char *label, bool expect_exact)
{
  // 2-cell conf_range so set_interface_idx(0, 1) selects cells 0/1's
  // cached face-local rows for wave_lax_curved.
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 1 };
  gkyl_range_init(&conf_range, 1, lower, upper);

  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range,
    WV_GR_EULER_TETRAD_RP_LAX);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });

  struct test_ws tws = test_ws_new(&conf_range, prods, eqn);

  int D_violations = 0, S2_violations = 0, tau_violations = 0;
  int total_cells = 0;
  double max_cone_depth = 0.0;

  int n_pts = sizeof(g_positivity_points) / sizeof(*g_positivity_points);

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);

  struct cond_block { const struct gr_euler_pos_rp_case *rows; int n; };
  struct cond_block blocks[] = {
    { gr_euler_pos_cases, GR_EULER_POS_NUM_CASES },
    { gr_euler_pos_cone_marginal_cases, GR_EULER_POS_NUM_CONE_MARGINAL_CASES },
  };

  for (int pi = 0; pi < n_pts; pi++) {
    const struct positivity_point *pt = &g_positivity_points[pi];
    // L cell at (x - dx/2, y, z); R cell at (x + dx/2, y, z).
    double dx_half = 0.5 * GR_EULER_POSITIVITY_DX;
    fill_prods_at(spacetime, pt->x - dx_half, pt->y, pt->z, prods_L);
    fill_prods_at(spacetime, pt->x + dx_half, pt->y, pt->z, prods_R);
    if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
        prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;
    test_ws_refresh(&tws, prods);

    for (size_t b = 0; b < sizeof(blocks)/sizeof(*blocks); b++)
      for (int c = 0; c < blocks[b].n; c++)
        total_cells += run_positivity_for_rp(eqn, prods_L, prods_R,
          &blocks[b].rows[c], GKYL_WV_LOW_ORDER_FLUX, 0.9,
          &D_violations, &S2_violations, &tau_violations, &max_cone_depth);
  }

  fprintf(stderr,
    "[positivity %s] across %d cell-updates (per-cell-metric, dx=%.4f): "
    "D<=0: %d (%.1f%%), cone crossings: %d (%.1f%%, max depth %.2e), "
    "τ<0: %d (%.1f%%)\n",
    label, total_cells, GR_EULER_POSITIVITY_DX,
    D_violations,   100.0 * D_violations   / total_cells,
    S2_violations,  100.0 * S2_violations  / total_cells, max_cone_depth,
    tau_violations, 100.0 * tau_violations / total_cells);

  TEST_CHECK_( D_violations  == 0, "D > 0 violated %d times for %s",
    D_violations, label );
  // Minkowski: Δγ = 0 makes the convex-combination argument exact, so
  // the raw step must preserve the cone with NO crossings. On curved
  // metrics crossings are allowed (the composite contract was asserted
  // per crossing inside the probe).
  if (expect_exact)
    TEST_CHECK_( S2_violations == 0,
      "raw step crossed the cone %d times on flat space for %s",
      S2_violations, label );
  TEST_MSG( "τ ≥ 0 violated %d times for %s — info-only",
    tau_violations, label );

  test_ws_release(&tws);
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
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-Mink",
      /*expect_exact=*/true);
  gkyl_gr_spacetime_release(st);
}
void test_positivity_curved_lax_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-Schw",
      /*expect_exact=*/false);
  gkyl_gr_spacetime_release(st);
}
void test_positivity_curved_lax_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-Kerr",
      /*expect_exact=*/false);
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
    run_positivity_sweep(st, eos_modes[ei], "curvedLax-BHL(M=0.3)",
      /*expect_exact=*/false);
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// 2. Unified positivity registry: the (condition × flux) matrix.
//
// Conditions = shared PRIM table + cone-marginal table, swept at
// production dx over the metric point ladder, plus the same conditions
// re-swept at GR_EULER_POS_COARSE_DX over a near-horizon coarse ladder
// (coarse near-horizon grids are a resolution users can run; crossings
// there are deeper — the large-Δγ regime — and the composite contract
// still binds them).
//
// Fluxes: tetrad-Lax / tetrad-HLL / tetrad-HLLC (HIGH_ORDER) and
// curved-Lax (LOW_ORDER), all at CFL = 1 (extremal).
//
// The printed matrix is the retained "which conditions cross for which
// flux" record — counts at the boundary are FP-fragile under
// -ffast-math, so the binding asserts are structural:
//   - D > 0 for every flux on every condition.
//   - Minkowski control: zero crossings, every flux, every condition.
//   - The composite repair contract on every crossing (in-probe).
// ---------------------------------------------------------------------------
enum gr_pos_flux {
  GRPF_TET_LAX = 0, GRPF_TET_HLL = 1, GRPF_TET_HLLC = 2, GRPF_CUR_LAX = 3,
  GRPF_N = 4
};
static const char *gr_pos_flux_name[GRPF_N] =
  { "tetrad-Lax", "tetrad-HLL", "tetrad-HLLC", "curved-Lax" };

#define GRPR_MAX_CONDS 128
// Coarse-Δγ axis spacing (near-horizon production-coarse cells).
#define GR_EULER_POS_COARSE_DX (0.15)
struct gr_pos_registry_row {
  const char *name;
  int n[GRPF_N], D[GRPF_N], S2[GRPF_N], tau[GRPF_N];
  double depth[GRPF_N];  // max cone-crossing depth (cone-budget units)
};

static void
run_positivity_registry(struct gkyl_gr_spacetime *spacetime,
  const char *label,
  const struct positivity_point *pts, int n_pts,
  const struct positivity_point *coarse_pts, int n_coarse,
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
  // One cache shared by all three equations (same prods, same range).
  struct test_ws tws = test_ws_new(&conf_range, prods, eqn_lax);
  gkyl_gr_euler_tetrad_set_wave_spacetime(eqn_hll, tws.ws);
  gkyl_gr_euler_tetrad_set_wave_spacetime(eqn_hllc, tws.ws);

  struct gkyl_wv_eqn *flux_eqn[GRPF_N] = { eqn_lax, eqn_hll, eqn_hllc, eqn_lax };
  enum gkyl_wv_flux_type flux_ft[GRPF_N] = {
    GKYL_WV_HIGH_ORDER_FLUX, GKYL_WV_HIGH_ORDER_FLUX, GKYL_WV_HIGH_ORDER_FLUX,
    GKYL_WV_LOW_ORDER_FLUX };

  struct gr_pos_registry_row reg[GRPR_MAX_CONDS];
  int n_conds = 0;
  memset(reg, 0, sizeof(reg));

  double *prods_L = gkyl_array_fetch(prods, 0);
  double *prods_R = gkyl_array_fetch(prods, 1);

  struct cond_block { const struct gr_euler_pos_rp_case *rows; int n; };
  struct cond_block blocks[] = {
    { gr_euler_pos_cases, GR_EULER_POS_NUM_CASES },
    { gr_euler_pos_cone_marginal_cases, GR_EULER_POS_NUM_CONE_MARGINAL_CASES },
  };
  int n_block_rows = GR_EULER_POS_NUM_CASES
                   + GR_EULER_POS_NUM_CONE_MARGINAL_CASES;

  // Production-dx axis over the point ladder.
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
        test_ws_refresh(&tws, prods);
        for (int fl = 0; fl < GRPF_N; fl++) {
          row->n[fl] += run_positivity_for_rp(flux_eqn[fl],
            prods_L, prods_R, rc, flux_ft[fl], 1.0,
            &row->D[fl], &row->S2[fl], &row->tau[fl], &row->depth[fl]);
        }
      }
    }
  }

  // Coarse-dx axis: same conditions over the near-horizon coarse ladder.
  static char coarse_names[64][56];
  if (n_coarse > 0) {
    int ci = 0;
    for (size_t b = 0; b < sizeof(blocks)/sizeof(*blocks); b++) {
      for (int c = 0; c < blocks[b].n; c++, ci++) {
        const struct gr_euler_pos_rp_case *rc = &blocks[b].rows[c];
        struct gr_pos_registry_row *row = &reg[n_conds++];
        snprintf(coarse_names[ci], 56, "coarse-dx %s", rc->name);
        row->name = coarse_names[ci];
        for (int pi = 0; pi < n_coarse; pi++) {
          double dxh = 0.5 * GR_EULER_POS_COARSE_DX;
          fill_prods_at(spacetime, coarse_pts[pi].x - dxh,
            coarse_pts[pi].y, coarse_pts[pi].z, prods_L);
          fill_prods_at(spacetime, coarse_pts[pi].x + dxh,
            coarse_pts[pi].y, coarse_pts[pi].z, prods_R);
          if (prods_L[GKYL_GR_SP_EXCISION] < 0.0 ||
              prods_R[GKYL_GR_SP_EXCISION] < 0.0) continue;
          test_ws_refresh(&tws, prods);
          for (int fl = 0; fl < GRPF_N; fl++) {
            row->n[fl] += run_positivity_for_rp(flux_eqn[fl],
              prods_L, prods_R, rc, flux_ft[fl], 1.0,
              &row->D[fl], &row->S2[fl], &row->tau[fl], &row->depth[fl]);
          }
        }
      }
    }
  }
  TEST_ASSERT( n_conds <= GRPR_MAX_CONDS );
  TEST_ASSERT( n_block_rows <= 64 );

  // Print the registry matrix: only conditions with at least one
  // crossing are listed (plus a global summary line).
  int tot_S2[GRPF_N] = {0}, tot_D[GRPF_N] = {0}, tot_tau[GRPF_N] = {0},
      tot_n[GRPF_N] = {0};
  fprintf(stderr, "[REGISTRY %s] conditions with crossings "
    "(s²/τ per flux; n = updates):\n", label);
  for (int c = 0; c < n_conds; c++) {
    bool any = false;
    for (int fl = 0; fl < GRPF_N; fl++) {
      tot_S2[fl] += reg[c].S2[fl]; tot_D[fl] += reg[c].D[fl];
      tot_tau[fl] += reg[c].tau[fl]; tot_n[fl] += reg[c].n[fl];
      if (reg[c].S2[fl] || reg[c].tau[fl] || reg[c].D[fl]) any = true;
    }
    if (any) {
      fprintf(stderr, "  REGISTRY %-24s", reg[c].name);
      for (int fl = 0; fl < GRPF_N; fl++)
        fprintf(stderr, "  %s: s2=%d depth=%.1e (n=%d)",
          gr_pos_flux_name[fl], reg[c].S2[fl], reg[c].depth[fl], reg[c].n[fl]);
      fprintf(stderr, "\n");
    }
  }
  for (int fl = 0; fl < GRPF_N; fl++)
    fprintf(stderr, "  REGISTRY-TOTAL %-12s s2=%d τ=%d D=%d across %d updates\n",
      gr_pos_flux_name[fl], tot_S2[fl], tot_tau[fl], tot_D[fl], tot_n[fl]);

  // Structural asserts (see file header).
  for (int fl = 0; fl < GRPF_N; fl++) {
    TEST_CHECK_( tot_D[fl] == 0, "%s: D > 0 violated %d times for %s",
      gr_pos_flux_name[fl], tot_D[fl], label );
    if (is_minkowski_control) {
      TEST_CHECK_( tot_S2[fl] == 0 && tot_tau[fl] == 0,
        "%s: Minkowski control violated (s² %d, τ %d) for %s",
        gr_pos_flux_name[fl], tot_S2[fl], tot_tau[fl], label );
    }
  }

  gkyl_array_release(prods);
  test_ws_release(&tws);
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
    /*is_minkowski_control=*/true);
  gkyl_gr_spacetime_release(st);
}

// BHL M = 0.3 (production spacetime): conditions sweep the bow-shock
// points at production dx.
void test_positivity_registry_bhl(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.3, 0.0, 0.0, 0.0, 0.0);
  static const struct positivity_point pts[] = {
    {  0.65,  0.0,   0.0 },  // radial front
    { -0.05, -0.6,   0.0 },  // bow-shock front
    { -0.05,  0.6,   0.0 },  // mirror
    { -0.05, -0.605, 0.0 },  // hugging horizon
    {  1.0,   1.0,   0.0 },  // control, far
  };
  run_positivity_registry(st, "BHL(M=0.3)-registry",
    pts, (int)(sizeof(pts) / sizeof(*pts)),
    NULL, 0,
    /*is_minkowski_control=*/false);
  gkyl_gr_spacetime_release(st);
}

void test_positivity_registry_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  // Coarse ladder: centroids chosen so xL = centroid − dx/2 stays
  // outside r_+ = 0.2 at dx = 0.15.
  static const struct positivity_point coarse[] = {
    { 0.28, 0.0, 0.0 }, { 0.33, 0.0, 0.0 }, { 0.41, 0.0, 0.0 },
    { 0.61, 0.0, 0.0 },
  };
  run_positivity_registry(st, "Schw-registry",
    g_ladder_schw_m01, GR_EULER_N_PTS(g_ladder_schw_m01),
    coarse, GR_EULER_N_PTS(coarse),
    /*is_minkowski_control=*/false);
  gkyl_gr_spacetime_release(st);
}

void test_positivity_registry_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  // Coarse ladder: xL stays outside r_+ ≈ 0.144 at dx = 0.15; includes
  // off-axis frame-dragging rows.
  static const struct positivity_point coarse[] = {
    { 0.23, 0.0, 0.0 },  { 0.31, 0.0, 0.0 },  { 0.46, 0.0, 0.0 },
    { 0.23, 0.10, 0.0 }, { 0.31, 0.10, 0.0 },
  };
  run_positivity_registry(st, "Kerr-registry",
    g_ladder_kerr_a09, GR_EULER_N_PTS(g_ladder_kerr_a09),
    coarse, GR_EULER_N_PTS(coarse),
    /*is_minkowski_control=*/false);
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

  { NULL, NULL },
};
