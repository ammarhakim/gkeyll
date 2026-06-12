// Shared helpers for the GR-Euler tetrad unit-test suites
// (ctest_gr_euler_prim / ctest_gr_euler_riemann / ctest_gr_euler_positivity).
// Split out of the former ctest_wv_gr_euler_tetrad.c monolith; see
// CTEST_SPLIT_PLAN.md. All helpers are static — each translation unit gets
// its own copy, which is fine for test code.
//
// acutest.h is NOT included here: TEST_LIST macros live per test binary.
#pragma once

#include <math.h>

#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_wave_spacetime.h>
#include <gkyl_wv_gr_euler_prim_priv.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>

// EOS modes exercised by the parametrized runners. The three production
// EOS configurations: IDEAL γ=5/3, APPROXIMATE_SYNGE with use_rcc=false
// (Mathews-Taub cubic closure), and APPROXIMATE_SYNGE with use_rcc=true
// (Ryu-Chattopadhyay Newton refinement on top of TM).
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

static struct gkyl_wv_eqn *
make_eqn(struct gkyl_gr_euler_eos eos, struct gkyl_range conf_range,
  enum gkyl_wv_gr_euler_tetrad_rp rp)
{
  return gkyl_wv_gr_euler_tetrad_inew(
    &(struct gkyl_wv_gr_euler_tetrad_inp){
      .gas_gamma = eos.gas_gamma,  // legacy IDEAL shortcut (eos overrides)
      .eos = eos,
      .conf_range = conf_range,
      .rp_type = rp,
      .use_gpu = false,
    });
}

// Returns the curved-frame Banyuls flux ΔF for a given pair of states
// using the same cell-centered Banyuls flux helper the wave construction
// uses. prods_l/prods_r are the per-side FACE-LOCAL rows; for the
// x-normal identity frame every runner here uses, the global-frame rows
// are already face-local.
static void
banyuls_delta_flux(struct gkyl_gr_euler_eos eos,
  const double *prods_l, const double *prods_r,
  const double qL[5], const double qR[5], double dF[5])
{
  double fL_gr[5], fR_gr[5];
  gkyl_gr_euler_banyuls_flux_cell(eos, qL, prods_l, NULL, fL_gr);
  gkyl_gr_euler_banyuls_flux_cell(eos, qR, prods_r, NULL, fR_gr);
  for (int i = 0; i < 5; i++) dF[i] = fR_gr[i] - fL_gr[i];
}

// ---------------------------------------------------------------------------
// wave_spacetime scaffolding. The tetrad equation REQUIRES the cache
// (the rotation-parity contract is retired — WAVE_SPACETIME_PARITY_PLAN
// .md), so every runner that drives waves_func/qfluct/flux_jump attaches
// one. Build AFTER the prods array exists; call test_ws_refresh after
// every prods refill so the cached face-local rows track the new metric
// samples. The grid extents are arbitrary — geometry enters only through
// the prods rows; wave_geom supplies the axis-aligned face frames.
// ---------------------------------------------------------------------------
struct test_ws {
  struct gkyl_rect_grid grid;
  struct gkyl_wave_geom *wg;
  struct gkyl_wave_spacetime *ws;
};

static struct test_ws
test_ws_new(const struct gkyl_range *conf_range,
  const struct gkyl_array *prods, struct gkyl_wv_eqn *eqn)
{
  struct test_ws t;
  int ndim = conf_range->ndim;
  double lo[GKYL_MAX_CDIM] = { 0.0 }, up[GKYL_MAX_CDIM];
  int cells[GKYL_MAX_CDIM];
  for (int d = 0; d < ndim; d++) {
    cells[d] = conf_range->upper[d] - conf_range->lower[d] + 1;
    lo[d] = 0.0;
    up[d] = (double)cells[d];
  }
  gkyl_rect_grid_init(&t.grid, ndim, lo, up, cells);
  t.wg = gkyl_wave_geom_new(&t.grid, (struct gkyl_range *)conf_range,
    0, 0, false);
  t.ws = gkyl_wave_spacetime_new(&t.grid, conf_range, t.wg, NULL, prods,
    0.0, /*use_gpu=*/false);
  gkyl_gr_euler_tetrad_set_wave_spacetime(eqn, t.ws);
  return t;
}

static void
test_ws_refresh(struct test_ws *t, const struct gkyl_array *prods)
{
  gkyl_wave_spacetime_refresh(t->ws, &t->grid, t->wg, NULL, prods, 0.0);
}

static void
test_ws_release(struct test_ws *t)
{
  gkyl_wave_spacetime_release(t->ws);
  gkyl_wave_geom_release(t->wg);
}

// Interface-averaged prods row — the production averaging policy:
// element-wise arithmetic mean of all components, then DERIVE inv_g and
// sqrt(det γ) from the inverse of the averaged γ_ij (Phase 0 Fix 2:
// a consistent pair from one inversion).
static void
build_prods_iface(const double *prods_L, const double *prods_R,
  double *prods_iface)
{
  for (int k = 0; k < GKYL_GR_SP_NCOMP_BASE; k++)
    prods_iface[k] = 0.5 * (prods_L[k] + prods_R[k]);
  double g_m[3][3], inv_g_m[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g_m[i][j] = prods_iface[GKYL_GR_SP_GIJ + 3*i + j];
  double det = gkyl_wave_spacetime_invert_metric_3x3(g_m, inv_g_m);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods_iface[GKYL_GR_SP_INV_GIJ + 3*i + j] = inv_g_m[i][j];
  prods_iface[GKYL_GR_SP_SPATIAL_DET] = det;
}

// Admissibility of an undensitized Convention-A state in the cell's own
// inverse spatial metric: D > 0, τ ≥ 0, s² = (D+τ)² − S·S > 0.
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

// Near-horizon radial ladders (interface centroids; cells at x ∓ dx/2),
// shared by the coarse-dx registry axis and the metric-aware margin sweep.
// Schwarzschild M = 0.1 → r_+ = 0.2; Kerr M = 0.1, a = 0.9 → r_+ ≈ 0.144
// (on-axis, off-axis frame-dragging, and out-of-plane rows).
static const struct positivity_point g_ladder_schw_m01[] = {
  { 0.23, 0.0, 0.0 }, { 0.25, 0.0, 0.0 }, { 0.28, 0.0, 0.0 },
  { 0.33, 0.0, 0.0 }, { 0.41, 0.0, 0.0 }, { 0.61, 0.0, 0.0 },
  { 1.01, 0.0, 0.0 },
};
static const struct positivity_point g_ladder_kerr_a09[] = {
  { 0.17, 0.0,  0.0 },  { 0.19, 0.0,  0.0 },  { 0.23, 0.0,  0.0 },
  { 0.31, 0.0,  0.0 },  { 0.46, 0.0,  0.0 },  { 0.81, 0.0,  0.0 },
  { 0.17, 0.10, 0.0 },  { 0.23, 0.10, 0.0 },  { 0.31, 0.10, 0.0 },
  { 0.17, 0.0,  0.05 }, { 0.23, 0.0,  0.05 }, { 0.31, 0.0,  0.05 },
};
#define GR_EULER_N_PTS(arr) ((int)(sizeof(arr) / sizeof(*(arr))))
