#include <assert.h>
#include <math.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_spacetime_coupling.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_gr_euler_prim_priv.h>

// Single-call cascade: each repair_state invocation fixes all three
// admissibility constraints (D > 0, τ ≥ 0, s² > 0) in safe order, so the
// post-call state is always admissible. No iteration needed.

// ---------------------------------------------------------------------------
// Internal helpers.
// ---------------------------------------------------------------------------

// Compute the inverse of a 3x3 symmetric matrix using the same closed-form
// expression used by gkyl_gr_euler_inv_spatial_metric in wv_gr_euler.c. This
// is duplicated here (rather than calling the packed helper) so the
// spacetime-coupling object stays decoupled from the packed equation; once
// the packed wv_gr_euler is removed this remains the canonical path.
static inline void
compute_inv_spatial_metric(const double g[3][3], double inv_g[3][3])
{
  double spatial_det =
    g[0][0] * (g[1][1] * g[2][2] - g[2][1] * g[1][2])
  - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
  + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);

  double trace = g[0][0] + g[1][1] + g[2][2];

  double g_sq[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_sq[i][j] = 0.0;
      for (int k = 0; k < 3; k++) g_sq[i][j] += g[i][k] * g[k][j];
    }

  double sq_trace = g_sq[0][0] + g_sq[1][1] + g_sq[2][2];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double euclid = (i == j) ? 1.0 : 0.0;
      inv_g[i][j] = (1.0 / spatial_det) *
        ((0.5 * (trace * trace - sq_trace) * euclid)
         - (trace * g[i][j])
         + g_sq[i][j]);
    }
  }
}

void
gkyl_moment_spacetime_coupling_fill_products_analytic(
  struct gkyl_gr_spacetime *spacetime,
  double t_curr, double x, double y, double z, double *prods)
{
  double spatial_det, lapse;
  bool in_excision_region;

  // Allocate the heap arrays the gr_spacetime callbacks fill. Packed code
  // does the same allocation pattern; matching it keeps memory ownership
  // stories identical, and the cost is per-cell-per-step (cheap for typical
  // grid sizes).
  double *shift = gkyl_malloc(sizeof(double[3]));
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i]    = gkyl_malloc(sizeof(double[3]));
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

  // Match the finite-difference epsilon used by gr_euler_impose_gauge in
  // level_set.c so derive_products numerically agrees with what the packed
  // path reinitialises into q[5..66] each sweep.
  const double eps = pow(10.0, -8.0);

  spacetime->lapse_function_func(spacetime, t_curr, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, t_curr, x, y, z, &shift);
  spacetime->spatial_metric_det_func(spacetime, t_curr, x, y, z, &spatial_det);
  spacetime->excision_region_func(spacetime, t_curr, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, t_curr, x, y, z, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, t_curr, x, y, z, eps, eps, eps, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, t_curr, x, y, z, eps, eps, eps, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, t_curr, x, y, z, eps, eps, eps, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, t_curr, x, y, z, eps, eps, eps, &spatial_metric_der);

  prods[GKYL_GR_SP_LAPSE]    = lapse;
  prods[GKYL_GR_SP_SHIFT + 0] = shift[0];
  prods[GKYL_GR_SP_SHIFT + 1] = shift[1];
  prods[GKYL_GR_SP_SHIFT + 2] = shift[2];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_GIJ + 3*i + j] = spatial_metric[i][j];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_KIJ + 3*i + j] = extrinsic_curvature[i][j];

  prods[GKYL_GR_SP_EXCISION] = in_excision_region ? -1.0 : 1.0;

  for (int i = 0; i < 3; i++)
    prods[GKYL_GR_SP_DALPHA + i] = lapse_der[i];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_DBETA + 3*i + j] = shift_der[i][j];

  for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        prods[GKYL_GR_SP_DGIJ + 9*k + 3*i + j] = spatial_metric_der[k][i][j];

  // Derived blocks: inverse spatial metric and determinant. Computing them
  // once here saves recomputation in every source-update / wave call.
  double g[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g[i][j] = spatial_metric[i][j];

  double inv_g[3][3];
  compute_inv_spatial_metric(g, inv_g);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_INV_GIJ + 3*i + j] = inv_g[i][j];

  prods[GKYL_GR_SP_SPATIAL_DET] = spatial_det;

  prods[GKYL_GR_SP_CELLCTR + 0] = x;
  prods[GKYL_GR_SP_CELLCTR + 1] = y;
  prods[GKYL_GR_SP_CELLCTR + 2] = z;

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]);
    gkyl_free(shift_der[i]);
    for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]);
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(extrinsic_curvature);
  gkyl_free(shift_der);
  gkyl_free(spatial_metric_der);
  gkyl_free(shift);
  gkyl_free(lapse_der);
}

// Compute the geometric source rate vector S(q) for the modular GR Euler
// equation. dq/dt = S(q) for the Banyuls source terms; integrators wrap
// this to take a time step. Excised cells return zero.
//
// Returns the rate of change of the densitized conservative state
// d(√γ·U)/dt = α√γ·S(w), so callers get the same scaling the original
// inlined source step used.
//
// The eos bundle controls the primitive-variable recovery closure (IDEAL
// or MATHEWS_TAUB). Everything downstream of the recovery uses h directly
// from the recovered primitives — no further EOS dispatch needed.
// gamma_eff_cell (optional, may be NULL): per-cell γ_eff cache slot used by
// the recovery's Picard iteration as a warm-start initial guess for TM.
// When non-NULL, the converged γ_eff is written back after the recovery so
// the next call from this cell starts from the right neighborhood. NULL
// preserves the cold-flow default initial guess of γ=5/3.
static void
compute_source_rate(struct gkyl_gr_euler_eos eos, const double *prods,
  const double q[5], double *gamma_eff_cell, double S_rate[5])
{
  bool in_excision_region = prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (in_excision_region) {
    for (int i = 0; i < 5; i++) S_rate[i] = 0.0;
    return;
  }

  double lapse   = prods[GKYL_GR_SP_LAPSE];
  double shift_x = prods[GKYL_GR_SP_SHIFT + 0];
  double shift_y = prods[GKYL_GR_SP_SHIFT + 1];
  double shift_z = prods[GKYL_GR_SP_SHIFT + 2];
  double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];

  // Convention A primitive recovery via the shared helper.
  double D    = q[0] / sqrt(spatial_det);
  double momx = q[1] / sqrt(spatial_det);
  double momy = q[2] / sqrt(spatial_det);
  double momz = q[3] / sqrt(spatial_det);
  double Etot = q[4] / sqrt(spatial_det);

  const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
  double inv_g[3][3] = {
    { ig[0], ig[1], ig[2] },
    { ig[3], ig[4], ig[5] },
    { ig[6], ig[7], ig[8] },
  };

  struct gkyl_gr_euler_prim prim;
  gkyl_gr_euler_set_recovery_context(GR_EULER_CTX_SOURCE);
  gkyl_gr_euler_recover_primitives(eos,
    D, momx, momy, momz, Etot, inv_g, gamma_eff_cell, &prim);
  gkyl_gr_euler_set_recovery_context(GR_EULER_CTX_UNKNOWN);

  double rho = prim.rho;
  double vx  = prim.v[0];
  double vy  = prim.v[1];
  double vz  = prim.v[2];
  double p   = prim.p;
  double W   = prim.W;
  double h   = prim.h;

  // Spacetime 4-velocity and contravariant 4-metric.
  double u4[4];
  u4[0] = W / lapse;
  u4[1] = (W * vx) - (shift_x * (W / lapse));
  u4[2] = (W * vy) - (shift_y * (W / lapse));
  u4[3] = (W * vz) - (shift_z * (W / lapse));

  double shift[3] = { shift_x, shift_y, shift_z };

  double inv_g4[4][4];
  inv_g4[0][0] = - (1.0 / (lapse * lapse));
  for (int i = 0; i < 3; i++) {
    inv_g4[0][i + 1] = (1.0 / (lapse * lapse)) * shift[i];
    inv_g4[i + 1][0] = (1.0 / (lapse * lapse)) * shift[i];
  }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      inv_g4[i + 1][j + 1] = prods[GKYL_GR_SP_INV_GIJ + 3*i + j]
        - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);

  // Perfect-fluid stress-energy tensor.
  double T[4][4];
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      T[i][j] = (rho * h * u4[i] * u4[j]) + (p * inv_g4[i][j]);

  // Spacetime-derivative lookups from products.
  double lapse_der[3] = {
    prods[GKYL_GR_SP_DALPHA + 0],
    prods[GKYL_GR_SP_DALPHA + 1],
    prods[GKYL_GR_SP_DALPHA + 2],
  };

  double shift_der[3][3];
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
      shift_der[j][i] = prods[GKYL_GR_SP_DBETA + 3*j + i];

  double spatial_metric_der[3][3][3];
  for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        spatial_metric_der[k][i][j] = prods[GKYL_GR_SP_DGIJ + 9*k + 3*i + j];

  double extrinsic_curvature[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      extrinsic_curvature[i][j] = prods[GKYL_GR_SP_KIJ + 3*i + j];

  const double *g_ij = &prods[GKYL_GR_SP_GIJ];

  double shift_lower[3] = { 0.0, 0.0, 0.0 };
  for (int m = 0; m < 3; m++)
    for (int k = 0; k < 3; k++)
      shift_lower[m] += g_ij[3*m + k] * shift[k];

  for (int i = 0; i < 5; i++) S_rate[i] = 0.0;

  double prefac = sqrt(spatial_det) * lapse;  // = α√γ = √(-g)

  // Energy density source — see source_euler for derivation.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      S_rate[4] += prefac * (T[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
      S_rate[4] += prefac * (2.0 * T[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
      S_rate[4] += prefac * (T[i + 1][j + 1] * extrinsic_curvature[i][j]);
    }
    S_rate[4] -= prefac * (T[0][0] * shift[i] * lapse_der[i]);
    S_rate[4] -= prefac * (T[0][i + 1] * lapse_der[i]);
  }

  // Momentum density sources — full Banyuls (1/2)·T^{μν}·∂_j g_{μν}
  // expansion. See source_euler for the per-term key.
  for (int j = 0; j < 3; j++) {
    S_rate[1 + j] -= prefac * (T[0][0] * lapse * lapse_der[j]);

    for (int k = 0; k < 3; k++) {
      for (int l = 0; l < 3; l++) {
        S_rate[1 + j] += prefac * (0.5 * T[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
        S_rate[1 + j] += prefac * (0.5 * T[k + 1][l + 1] * spatial_metric_der[j][k][l]);
      }
      S_rate[1 + j] += prefac * (T[0][0] * shift_lower[k] * shift_der[j][k]);

      for (int i = 0; i < 3; i++) {
        S_rate[1 + j] += prefac * (T[0][i + 1] * g_ij[3*i + k] * shift_der[j][k]);
        S_rate[1 + j] += prefac * (T[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
      }
    }
  }
}

// Per-cell forward-Euler source step for the modular GR Euler equation
// with TAU + S² POSITIVITY LIMITERS.
//
// Naive forward Euler: q_new = q_old + dt·S(q_old).
//
// The limiters compute the largest α ∈ [0, 1] such that
//   τ_new = τ + α·dt·S_τ ≥ 0
// AND
//   s²(q_new) = (D+τ + α·δDτ)² − γ^{ij}(S_i + α δS_i)(S_j + α δS_j)
//             ≥ margin · (D+τ + α·δDτ)²
// (where δ = dt·S). Both constraints reduce to picking a scalar α; the
// final α = min(α_τ, α_s²) is applied uniformly to ALL components.
// The cell gets a partial, physically-directed source contribution
// rather than a hard clamp.
//
// SSP-RK3 compatibility: the outer integrator applies convex
// combinations of (q_old, q_stage1, q_stage2). A_γ is convex, so a
// convex combination of admissible states is admissible. Per-substage
// limiting is safe.
//
// Counters: separate tallies for τ and s² limiter activations.
static uint64_t s_tau_limiter_fires = 0;
static uint64_t s_s2_limiter_fires = 0;

uint64_t
gkyl_moment_spacetime_coupling_tau_limiter_fires(void)
{
  return s_tau_limiter_fires;
}

uint64_t
gkyl_moment_spacetime_coupling_s2_limiter_fires(void)
{
  return s_s2_limiter_fires;
}

// ---------------------------------------------------------------------------
// Recovery-iteration instrumentation storage. Defined here (rather than in a
// header) so the histograms can be referenced from the static-inline
// gkyl_gr_euler_recover_primitives via extern declarations in
// gkyl_wv_gr_euler_prim_priv.h. Persists for the lifetime of the process;
// gkyl_gr_euler_print_recovery_stats dumps the tallies.
// ---------------------------------------------------------------------------

const int gkyl_gr_euler_newton_bin_edges[GR_EULER_NEWTON_HIST_NBINS] = {
  4, 8, 16, 32, 64, 99, 100
};
const int gkyl_gr_euler_picard_bin_edges[GR_EULER_PICARD_HIST_NBINS] = {
  1, 2, 3, 4, 9, 14, 19, 29, 30
};

static uint64_t s_newton_hist[GR_EULER_NEWTON_HIST_NBINS] = {0};
static uint64_t s_picard_hist[GR_EULER_PICARD_HIST_NBINS] = {0};
static uint64_t s_newton_total_iters = 0;
static uint64_t s_picard_total_iters = 0;
static uint64_t s_newton_calls = 0;
static uint64_t s_picard_calls = 0;
static int      s_newton_max = 0;
static int      s_picard_max = 0;

static inline int
hist_bucket(int n, const int *edges, int nbins)
{
  for (int b = 0; b < nbins; b++) {
    if (n <= edges[b]) return b;
  }
  return nbins - 1;  // overflow goes into the last bin
}

void
gkyl_gr_euler_record_newton_iters(int n)
{
  s_newton_hist[hist_bucket(n, gkyl_gr_euler_newton_bin_edges,
    GR_EULER_NEWTON_HIST_NBINS)] += 1;
  s_newton_total_iters += (uint64_t)n;
  s_newton_calls += 1;
  if (n > s_newton_max) s_newton_max = n;
}

void
gkyl_gr_euler_record_picard_iters(int n)
{
  s_picard_hist[hist_bucket(n, gkyl_gr_euler_picard_bin_edges,
    GR_EULER_PICARD_HIST_NBINS)] += 1;
  s_picard_total_iters += (uint64_t)n;
  s_picard_calls += 1;
  if (n > s_picard_max) s_picard_max = n;
}

// Bounded ring buffer for EM-Newton cap-hit diagnostic events. Caps total
// recorded events so a pathological run doesn't flood memory; once the
// cap is reached we still increment the total counter (so the dump can
// report how many were lost) but stop recording payloads. Per-run reset
// not exposed; tests that exercise extreme states (stringent prim_vars)
// rarely fire this path, so cross-test pollution is bounded.
#define GR_EULER_CAPFAIL_RING_SIZE 32

static struct gkyl_gr_euler_newton_capfail s_capfail_ring[GR_EULER_CAPFAIL_RING_SIZE];
static uint64_t s_capfail_recorded = 0;  // calls to record_newton_capfail
static int      s_capfail_written  = 0;  // count actually written into ring (≤ size)

void
gkyl_gr_euler_record_newton_capfail(
  const struct gkyl_gr_euler_newton_capfail *event)
{
  s_capfail_recorded += 1;
  if (s_capfail_written < GR_EULER_CAPFAIL_RING_SIZE) {
    s_capfail_ring[s_capfail_written] = *event;
    s_capfail_written += 1;
  }
}

// RC EOS comparison diagnostic counters. See header docstring for
// gkyl_gr_euler_record_rc_compare.
static uint64_t s_rc_attempts        = 0;
static uint64_t s_rc_accepted        = 0;  // all acceptances (standard + escape)
static uint64_t s_rc_accepted_escape = 0;  // subset accepted via conditioning escape
static double   s_rc_sum_dW_acc      = 0.0;  // sum of |ΔΓ|/Γ_tm over accepted cells
static double   s_rc_max_dW_acc      = 0.0;  // max  |ΔΓ|/Γ_tm over accepted cells
static double   s_rc_max_dW_rej      = 0.0;  // max  |ΔΓ|/Γ_tm over rejected cells

void
gkyl_gr_euler_record_rc_compare(
  bool accepted, bool conditioning_escape, double dW_rel)
{
  s_rc_attempts += 1;
  if (accepted) {
    s_rc_accepted += 1;
    if (conditioning_escape) s_rc_accepted_escape += 1;
    s_rc_sum_dW_acc += dW_rel;
    if (dW_rel > s_rc_max_dW_acc) s_rc_max_dW_acc = dW_rel;
  } else {
    if (dW_rel > s_rc_max_dW_rej) s_rc_max_dW_rej = dW_rel;
  }
}

// RC rejection diagnostic ring buffer. Stores the first 32 rejected
// cells so we can dump their (W_tm, W_rc, conservatives, conditioning)
// to see whether rejection is caused by spurious-root Newton failures
// (rejection correct, RC is wrong) or by TM coefficient cancellation at
// high W (rejection wrong, RC is actually better and we're discarding
// a real refinement).
#define GR_EULER_RC_REJECT_RING_SIZE 32
static struct gkyl_gr_euler_rc_reject s_rc_reject_ring[GR_EULER_RC_REJECT_RING_SIZE];
static int s_rc_reject_written = 0;

void
gkyl_gr_euler_record_rc_reject(const struct gkyl_gr_euler_rc_reject *event)
{
  if (s_rc_reject_written < GR_EULER_RC_REJECT_RING_SIZE) {
    s_rc_reject_ring[s_rc_reject_written] = *event;
    s_rc_reject_written += 1;
  }
}

// Floor-hit counters. Track how often the outer recover_primitives
// clamp had to apply the ρ or p floor — answers "is the floor doing
// work or is it irrelevant?" for floor-value sweep experiments.
static uint64_t s_floor_calls       = 0;
static uint64_t s_rho_floor_hits    = 0;
static uint64_t s_p_floor_hits      = 0;
static uint64_t s_both_floor_hits   = 0;  // cells where ρ AND p both floored

// Per-context floor counters. s_recovery_context is set by callers via
// gkyl_gr_euler_set_recovery_context before each recover_primitives call;
// the floor diagnostic reads it to bin the hit. Lets us distinguish
// transient HLL-rejected-by-positivity hits from real source/Lax hits.
static int s_recovery_context = GR_EULER_CTX_UNKNOWN;
static uint64_t s_floor_calls_ctx[GR_EULER_CTX_COUNT]      = {0};
static uint64_t s_p_floor_hits_ctx[GR_EULER_CTX_COUNT]     = {0};
static uint64_t s_rho_floor_hits_ctx[GR_EULER_CTX_COUNT]   = {0};

void
gkyl_gr_euler_set_recovery_context(int ctx)
{
  if (ctx >= 0 && ctx < GR_EULER_CTX_COUNT) s_recovery_context = ctx;
  else                                       s_recovery_context = GR_EULER_CTX_UNKNOWN;
}

// D-magnitude histogram for floored cells. Answers "are floored cells
// all near-vacuum (D << 1) or do moderate/normal-density cells also
// get floored?" Decade bins from D < 1e-10 to D ≥ 1.
#define GR_EULER_FLOOR_D_BINS 11
static uint64_t s_floor_D_hist[GR_EULER_FLOOR_D_BINS] = {0};
// Bin edges: <1e-10, 1e-10..1e-8, 1e-8..1e-6, 1e-6..1e-4, 1e-4..1e-2,
//            1e-2..1e-1, 1e-1..1, 1..10, 10..100, 100..1000, ≥1000
// (we put the "≥1000" overflow in the last slot)

static int
floor_D_bin(double D)
{
  if (D < 1.0e-10)  return 0;
  if (D < 1.0e-8)   return 1;
  if (D < 1.0e-6)   return 2;
  if (D < 1.0e-4)   return 3;
  if (D < 1.0e-2)   return 4;
  if (D < 1.0e-1)   return 5;
  if (D < 1.0)      return 6;
  if (D < 10.0)     return 7;
  if (D < 100.0)    return 8;
  if (D < 1000.0)   return 9;
  return 10;
}

// (E²−M²)/E² histogram — conditioning measure (~1/W² for floor-typical
// cells). Tells us whether floored cells are near-luminal or moderate W.
#define GR_EULER_FLOOR_S2_BINS 8
static uint64_t s_floor_s2_hist[GR_EULER_FLOOR_S2_BINS] = {0};
// Bin edges: <1e-12, 1e-12..1e-10, 1e-10..1e-8, 1e-8..1e-6,
//            1e-6..1e-4, 1e-4..1e-2, 1e-2..1, ≥1

static int
floor_s2_bin(double s_sq_over_E2)
{
  if (s_sq_over_E2 < 1.0e-12) return 0;
  if (s_sq_over_E2 < 1.0e-10) return 1;
  if (s_sq_over_E2 < 1.0e-8)  return 2;
  if (s_sq_over_E2 < 1.0e-6)  return 3;
  if (s_sq_over_E2 < 1.0e-4)  return 4;
  if (s_sq_over_E2 < 1.0e-2)  return 5;
  if (s_sq_over_E2 < 1.0)     return 6;
  return 7;
}

// Ring buffer of 64 floor-hit snapshots — captures conservatives +
// raw recovered (ρ, p) so we can characterize what these cells "look
// like" physically without needing position info threaded through the
// recovery API.
struct gkyl_gr_euler_floor_event {
  double D, tau, mom_sq;
  double rho_raw, p_raw;
  double s_sq_over_E2;
};
#define GR_EULER_FLOOR_RING_SIZE 64
static struct gkyl_gr_euler_floor_event s_floor_event_ring[GR_EULER_FLOOR_RING_SIZE];
static int s_floor_event_written = 0;

void
gkyl_gr_euler_record_floor_hit(
  bool rho_floored, bool p_floored,
  double D, double tau, double mom_sq,
  double rho_raw, double p_raw)
{
  s_floor_calls += 1;
  if (rho_floored) s_rho_floor_hits += 1;
  if (p_floored)   s_p_floor_hits   += 1;
  if (rho_floored && p_floored) s_both_floor_hits += 1;

  // Per-context counters: tag the call site that triggered this floor.
  int ctx = (s_recovery_context >= 0 && s_recovery_context < GR_EULER_CTX_COUNT)
    ? s_recovery_context : GR_EULER_CTX_UNKNOWN;
  s_floor_calls_ctx[ctx] += 1;
  if (rho_floored) s_rho_floor_hits_ctx[ctx] += 1;
  if (p_floored)   s_p_floor_hits_ctx[ctx]   += 1;

  // Only build histograms / ring entries for cells that actually
  // triggered a floor — keeps the histograms a property of FLOORED
  // cells, not all cells.
  if (rho_floored || p_floored) {
    s_floor_D_hist[floor_D_bin(D)] += 1;

    double E_lab = tau + D;
    double s_sq = (D + tau) * (D + tau) - mom_sq;
    double s_sq_over_E2 = (E_lab * E_lab > 1.0e-300)
      ? s_sq / (E_lab * E_lab) : 0.0;
    if (s_sq_over_E2 < 0.0) s_sq_over_E2 = 0.0;
    s_floor_s2_hist[floor_s2_bin(s_sq_over_E2)] += 1;

    if (s_floor_event_written < GR_EULER_FLOOR_RING_SIZE) {
      s_floor_event_ring[s_floor_event_written] = (struct gkyl_gr_euler_floor_event){
        .D = D, .tau = tau, .mom_sq = mom_sq,
        .rho_raw = rho_raw, .p_raw = p_raw,
        .s_sq_over_E2 = s_sq_over_E2,
      };
      s_floor_event_written += 1;
    }
  }
}

void
gkyl_gr_euler_print_recovery_stats(FILE *fp)
{
  if (s_newton_calls == 0 && s_picard_calls == 0) return;

  fprintf(fp,
    "[gr_euler_recovery] inner-Newton stats: %llu calls, total %llu iters, "
    "avg %.2f iters/call, max %d\n",
    (unsigned long long)s_newton_calls,
    (unsigned long long)s_newton_total_iters,
    s_newton_calls > 0
      ? (double)s_newton_total_iters / (double)s_newton_calls : 0.0,
    s_newton_max);
  const char *newton_labels[] = {
    "  ≤4   ", "  5–8  ", "  9–16 ", "  17–32", "  33–64", "  65–99", "  100  "
  };
  for (int b = 0; b < GR_EULER_NEWTON_HIST_NBINS; b++) {
    fprintf(fp, "    %s : %llu\n", newton_labels[b],
      (unsigned long long)s_newton_hist[b]);
  }

  fprintf(fp,
    "[gr_euler_recovery] outer-Picard stats: %llu calls, total %llu iters, "
    "avg %.2f iters/call, max %d\n",
    (unsigned long long)s_picard_calls,
    (unsigned long long)s_picard_total_iters,
    s_picard_calls > 0
      ? (double)s_picard_total_iters / (double)s_picard_calls : 0.0,
    s_picard_max);
  const char *picard_labels[] = {
    "  1    ", "  2    ", "  3    ", "  4    ",
    "  5–9  ", "  10–14", "  15–19", "  20–29", "  30+  "
  };
  for (int b = 0; b < GR_EULER_PICARD_HIST_NBINS; b++) {
    fprintf(fp, "    %s : %llu\n", picard_labels[b],
      (unsigned long long)s_picard_hist[b]);
  }

  // Floor-hit stats — applies to all EOSs. Tracks how often the outer
  // recover_primitives clamp had to fire. Crucial for floor-value sweep
  // experiments (1e-8 / 1e-10 / 1e-12) to distinguish "lower floor is
  // doing work" from "lower floor irrelevant, no cells reach it".
  if (s_floor_calls > 0) {
    fprintf(fp,
      "[gr_euler_recovery] Floor hits over %llu recovery calls:\n"
      "    ρ floored: %llu (%.4f%%)\n"
      "    p floored: %llu (%.4f%%)\n"
      "    both:      %llu (%.4f%%)\n",
      (unsigned long long)s_floor_calls,
      (unsigned long long)s_rho_floor_hits,
      100.0 * (double)s_rho_floor_hits / (double)s_floor_calls,
      (unsigned long long)s_p_floor_hits,
      100.0 * (double)s_p_floor_hits / (double)s_floor_calls,
      (unsigned long long)s_both_floor_hits,
      100.0 * (double)s_both_floor_hits / (double)s_floor_calls);

    // Per-context breakdown — distinguishes transient HLL-rejected-by-
    // positivity-sweep hits from real source / Lax (committed) hits.
    const char *ctx_labels[GR_EULER_CTX_COUNT] = {
      "UNKNOWN", "SOURCE ", "PRIMS  ", "HLL    ", "LAX    ", "HLLC   "
    };
    fprintf(fp, "[gr_euler_recovery] Per-context floor breakdown:\n");
    fprintf(fp, "    context  | total calls | p floored | %% of ctx | %% of all p floors\n");
    for (int c = 0; c < GR_EULER_CTX_COUNT; c++) {
      if (s_floor_calls_ctx[c] == 0) continue;
      double pct_of_ctx = 100.0 * (double)s_p_floor_hits_ctx[c] / (double)s_floor_calls_ctx[c];
      double pct_of_all_pfloor = (s_p_floor_hits > 0)
        ? 100.0 * (double)s_p_floor_hits_ctx[c] / (double)s_p_floor_hits : 0.0;
      fprintf(fp, "    %s  | %11llu | %9llu | %6.2f%% | %6.2f%%\n",
        ctx_labels[c],
        (unsigned long long)s_floor_calls_ctx[c],
        (unsigned long long)s_p_floor_hits_ctx[c],
        pct_of_ctx, pct_of_all_pfloor);
    }

    // D-magnitude histogram of floored cells.
    uint64_t total_floor_hits = s_rho_floor_hits + s_p_floor_hits - s_both_floor_hits;
    if (total_floor_hits > 0) {
      const char *D_labels[GR_EULER_FLOOR_D_BINS] = {
        "  < 1e-10        ", "  1e-10 to 1e-8  ", "  1e-8  to 1e-6  ",
        "  1e-6  to 1e-4  ", "  1e-4  to 1e-2  ", "  1e-2  to 1e-1  ",
        "  1e-1  to 1     ", "  1     to 10    ", "  10    to 100   ",
        "  100   to 1000  ", "  ≥ 1000         "
      };
      fprintf(fp, "[gr_euler_recovery] D-magnitude distribution of floored cells:\n");
      for (int b = 0; b < GR_EULER_FLOOR_D_BINS; b++) {
        fprintf(fp, "    %s : %llu (%.2f%%)\n",
          D_labels[b], (unsigned long long)s_floor_D_hist[b],
          100.0 * (double)s_floor_D_hist[b] / (double)total_floor_hits);
      }
      // (E²−M²)/E² histogram — conditioning measure (~1/W² for typical cells).
      const char *s2_labels[GR_EULER_FLOOR_S2_BINS] = {
        "  < 1e-12        ", "  1e-12 to 1e-10 ", "  1e-10 to 1e-8  ",
        "  1e-8  to 1e-6  ", "  1e-6  to 1e-4  ", "  1e-4  to 1e-2  ",
        "  1e-2  to 1     ", "  ≥ 1            "
      };
      fprintf(fp, "[gr_euler_recovery] (E²−M²)/E² distribution of floored cells:\n");
      for (int b = 0; b < GR_EULER_FLOOR_S2_BINS; b++) {
        fprintf(fp, "    %s : %llu (%.2f%%)\n",
          s2_labels[b], (unsigned long long)s_floor_s2_hist[b],
          100.0 * (double)s_floor_s2_hist[b] / (double)total_floor_hits);
      }
      // Ring-buffer sample dump.
      if (s_floor_event_written > 0) {
        fprintf(fp, "[gr_euler_recovery] Floor-event sample (first %d):\n"
          "  idx       D          τ          mom_sq      ρ_raw       p_raw       (E²−M²)/E²\n",
          s_floor_event_written);
        for (int i = 0; i < s_floor_event_written; i++) {
          const struct gkyl_gr_euler_floor_event *e = &s_floor_event_ring[i];
          fprintf(fp, "  %3d  %+10.3e  %+10.3e  %+10.3e  %+10.3e  %+10.3e  %+10.3e\n",
            i, e->D, e->tau, e->mom_sq, e->rho_raw, e->p_raw, e->s_sq_over_E2);
        }
      }
    }
  }

  // RC EOS comparison stats — printed only when RC was actually
  // exercised (i.e., the run used the RYU_CHATTOPADHYAY EOS).
  if (s_rc_attempts > 0) {
    double accept_frac = (double)s_rc_accepted / (double)s_rc_attempts;
    double avg_dW = s_rc_accepted > 0
      ? s_rc_sum_dW_acc / (double)s_rc_accepted : 0.0;
    double escape_frac = s_rc_accepted > 0
      ? (double)s_rc_accepted_escape / (double)s_rc_accepted : 0.0;
    fprintf(fp,
      "[gr_euler_recovery] RC vs TM compare: %llu RC attempts, %llu accepted (%.1f%%)\n"
      "    of accepted, %llu via conditioning escape (%.2f%% of accepted)\n"
      "    accepted |ΔΓ|/Γ_tm:  avg=%.3e  max=%.3e\n"
      "    rejected |ΔΓ|/Γ_tm:  max=%.3e  (rejected = fell back to TM)\n",
      (unsigned long long)s_rc_attempts,
      (unsigned long long)s_rc_accepted,
      100.0 * accept_frac,
      (unsigned long long)s_rc_accepted_escape,
      100.0 * escape_frac,
      avg_dW, s_rc_max_dW_acc,
      s_rc_max_dW_rej);

    // Per-cell snapshot dump of the first N rejected cells. Helps
    // distinguish spurious-root Newton failures (rejection correct)
    // from TM precision-loss at high W (rejection arguably wrong — RC
    // might be the better answer there).
    if (s_rc_reject_written > 0) {
      uint64_t total_rejected = s_rc_attempts - s_rc_accepted;
      fprintf(fp,
        "[gr_euler_recovery] RC rejections (showing first %d of %llu):\n"
        "  idx   |ΔΓ|/Γ_tm      W_tm       W_rc       D          τ          mom_sq     (E²−M²)/E²  θ_tm        p_tm        p_rc\n",
        s_rc_reject_written, (unsigned long long)total_rejected);
      for (int i = 0; i < s_rc_reject_written; i++) {
        const struct gkyl_gr_euler_rc_reject *e = &s_rc_reject_ring[i];
        fprintf(fp,
          "  %3d  %+10.3e  %+9.3e  %+9.3e  %+9.3e  %+9.3e  %+9.3e  %+10.3e  %+10.3e  %+10.3e  %+10.3e\n",
          i, e->dW_rel, e->W_tm, e->W_rc, e->D, e->tau, e->mom_sq,
          e->s_sq_over_E2, e->theta_tm, e->p_tm, e->p_rc);
      }
    }
  }

  // EM-Newton cap-hit diagnostic dump. Each ring-buffer entry is a
  // snapshot of one cell that failed to converge in 100 EM iterations.
  // Useful for understanding which physical regime is pathological —
  // typically cold-flow-near-floor or near-supraluminal.
  if (s_capfail_recorded > 0) {
    fprintf(fp,
      "[gr_euler_recovery] EM-Newton cap-hits: %llu total (showing first %d)\n",
      (unsigned long long)s_capfail_recorded, s_capfail_written);
    fprintf(fp,
      "  idx     γ_eff     D          τ          mom_sq     θ=p/ρ      ρ          p          W       s²/(D+τ)²  γ^xx     γ^yy     γ^zz     γ^xy     γ^xz     γ^yz\n");
    for (int i = 0; i < s_capfail_written; i++) {
      const struct gkyl_gr_euler_newton_capfail *e = &s_capfail_ring[i];
      fprintf(fp,
        "  %3d  %.6f  %+10.3e  %+10.3e  %+10.3e  %+10.3e  %+10.3e  %+10.3e  %7.3f  %+10.3e  %+8.2e %+8.2e %+8.2e %+8.2e %+8.2e %+8.2e\n",
        i, e->gas_gamma, e->D, e->tau, e->mom_sq, e->theta,
        e->rho, e->p, e->W, e->s_sq_over_Dtau_sq,
        e->inv_g_diag[0], e->inv_g_diag[1], e->inv_g_diag[2],
        e->inv_g_off[0],  e->inv_g_off[1],  e->inv_g_off[2]);
    }
  }
}

// Compute the largest α ∈ [0, 1] such that
//   (1-ε)(a + αb)² − (P + 2αQ + α²R) ≥ 0
// where a = D+τ, b = δD+δτ, P = γ^{ij}·S_i·S_j, Q = γ^{ij}·S_i·δS_j,
// R = γ^{ij}·δS_i·δS_j. Equivalently:
//   A α² + 2 B α + C ≥ 0
// with A = (1-ε)b² − R, B = (1-ε)ab − Q, C = (1-ε)a² − P. Assumes
// C ≥ 0 (input state q already satisfies the margin); if not, returns
// 1.0 to let the cascade-repair handle it downstream.
static double
compute_s2_limiter_alpha(const double *prods, double dt,
  const double fluid_old[5], const double S_rate[5])
{
  const double margin = 1.0e-6;

  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) return 1.0;
  double sd = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);
  if (!(sd > 0.0)) return 1.0;
  const double *ig = &prods[GKYL_GR_SP_INV_GIJ];

  // Undensitize q and δq (where δq = dt·S_rate).
  double D    = fluid_old[0] / sd;
  double Sx   = fluid_old[1] / sd;
  double Sy   = fluid_old[2] / sd;
  double Sz   = fluid_old[3] / sd;
  double tau  = fluid_old[4] / sd;
  double dD   = dt * S_rate[0] / sd;
  double dSx  = dt * S_rate[1] / sd;
  double dSy  = dt * S_rate[2] / sd;
  double dSz  = dt * S_rate[3] / sd;
  double dtau = dt * S_rate[4] / sd;

  double a = D + tau;
  double b = dD + dtau;

  // P, Q, R via curved metric contraction (symmetric inv_g).
  double S_v[3]  = { Sx, Sy, Sz };
  double dS_v[3] = { dSx, dSy, dSz };
  double P = 0.0, Q = 0.0, R = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double inv_g_ij = ig[3*i + j];
      P += inv_g_ij * S_v[i]  * S_v[j];
      Q += inv_g_ij * S_v[i]  * dS_v[j];
      R += inv_g_ij * dS_v[i] * dS_v[j];
    }

  double A = (1.0 - margin) * b * b - R;
  double B = (1.0 - margin) * a * b - Q;
  double C = (1.0 - margin) * a * a - P;

  if (C < 0.0) return 1.0;  // q already fails margin; let repair handle

  // Linear-coefficient case.
  if (fabs(A) < 1.0e-30) {
    if (B >= 0.0) return 1.0;  // f non-decreasing, always ≥ 0 from C≥0
    double alpha_root = -0.5 * C / B;
    return fmin(1.0, fmax(0.0, alpha_root));
  }

  double disc = B*B - A*C;
  if (disc < 0.0) {
    // No real roots. f has the sign of A everywhere. If A > 0, f > 0
    // always; if A < 0, f < 0 always — but C ≥ 0 already rules that out.
    return 1.0;
  }
  double sd_disc = sqrt(disc);

  // Safety pullback: scale the computed α by (1 - SAFETY_PULLBACK) so
  // s²_new lands strictly inside the boundary rather than exactly on
  // it. Floating-point error around the boundary would otherwise drop
  // s²_new by ±ε, half the time triggering cascade-repair. With
  // pullback, s²_new lands at (margin + SAFETY_PULLBACK·(D+τ)²)·(D+τ)²
  // worth of slack, comfortably above the cascade threshold.
  // Sweep results (BHL t=15, M=0.3) for source-s² fires (and wp s² in
  // parens), holding the cascade-margin at 1e-6:
  //   1e-4  → 4,897 (2,874)
  //   1e-6  → 4,165 (2,818)
  //   1e-8  → 3,791 (2,472)
  //   1e-10 → 2,587 (1,757)  ← optimum
  //   1e-12 → 4,770 (2,810)
  //   1e-14 → 4,378 (2,312)
  // Below 1e-10 the floating-point landing noise around the cascade
  // boundary dominates the small benefit from a tighter pullback.
  const double SAFETY_PULLBACK = 1.0e-10;

  if (A > 0.0) {
    // Parabola opens up; f < 0 strictly between roots. With C ≥ 0 and
    // A > 0, both roots have the same sign as -B/A's sign. If α_- > 0,
    // both roots > 0 and we cap at α_-. Otherwise both ≤ 0 and α=1 OK.
    double alpha_minus = (-B - sd_disc) / A;
    if (alpha_minus > 0.0)
      return fmin(1.0, alpha_minus * (1.0 - SAFETY_PULLBACK));
    return 1.0;
  } else {
    // A < 0; parabola opens down. f ≥ 0 in [α_lo, α_hi]. With C ≥ 0,
    // 0 ∈ [α_lo, α_hi] so α_lo ≤ 0 ≤ α_hi. Max valid α = min(1, α_hi).
    double r1 = (-B - sd_disc) / A;
    double r2 = (-B + sd_disc) / A;
    double a_hi = fmax(r1, r2);
    if (a_hi <= 0.0) return 1.0;  // shouldn't happen given C ≥ 0
    return fmin(1.0, a_hi * (1.0 - SAFETY_PULLBACK));
  }
}

void
gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
  struct gkyl_gr_euler_eos eos, double t_curr, double dt,
  const double *prods, double *gamma_eff_cell,
  const double fluid_old[5], double fluid_new[5])
{
  (void)t_curr;  // Source math is time-independent.

  // gamma_eff_cell threaded through to compute_source_rate so the Picard
  // initial guess for TM is the previous step's converged value.
  double S_rate[5];
  compute_source_rate(eos, prods, fluid_old, gamma_eff_cell, S_rate);

  // TAU-POSITIVITY LIMITER. δτ = dt·S_τ. Target τ_new = TAU_TARGET
  // (the same floor cascade-repair restores to) so that limited cells
  // land strictly inside A_γ — well-separated from the floating-point
  // ±ε region around 0 that would re-trigger cascade-repair.
  //   α_τ = (τ_old - TAU_TARGET) / -δτ  when δτ < 0 and τ_old > TAU_TARGET
  // If τ_old ≤ TAU_TARGET (already at or below floor), α=0 (no update).
  // If δτ ≥ 0 or τ_new_unlimited ≥ TAU_TARGET, no limit needed.
  const double TAU_TARGET = GR_EULER_TAU_REPAIR_FLOOR;
  double alpha = 1.0;
  double delta_tau = dt * S_rate[4];
  double tau_new_unlimited = fluid_old[4] + delta_tau;
  if (tau_new_unlimited < TAU_TARGET && delta_tau < 0.0) {
    double headroom = fluid_old[4] - TAU_TARGET;
    double alpha_tau = (headroom > 0.0) ? (headroom / -delta_tau) : 0.0;
    if (alpha_tau < 0.0) alpha_tau = 0.0;
    if (alpha_tau > 1.0) alpha_tau = 1.0;
    alpha = fmin(alpha, alpha_tau);
    s_tau_limiter_fires += 1;
  }

  // S²-POSITIVITY LIMITER. Solve the quadratic for max α s.t. the
  // limited source step keeps s²(q_new) ≥ margin·(D+τ)².
  double alpha_s2 = compute_s2_limiter_alpha(prods, dt, fluid_old, S_rate);
  if (alpha_s2 < alpha) {
    alpha = alpha_s2;
    s_s2_limiter_fires += 1;
  }

  double dt_eff = alpha * dt;
  for (int i = 0; i < 5; i++)
    fluid_new[i] = fluid_old[i] + dt_eff * S_rate[i];
}

// ---------------------------------------------------------------------------
// Public API.
// ---------------------------------------------------------------------------

gkyl_moment_spacetime_coupling*
gkyl_moment_spacetime_coupling_new(struct gkyl_moment_spacetime_coupling_inp inp)
{
  // At most one of analytic_spacetime / einstein_eqn must be set.
  assert(!(inp.analytic_spacetime != NULL && inp.einstein_eqn != NULL));

  gkyl_moment_spacetime_coupling *st =
    gkyl_malloc(sizeof(struct gkyl_moment_spacetime_coupling));
  st->grid    = *(inp.grid);
  st->ndim    = st->grid.ndim;
  st->nfluids = inp.nfluids;
  for (int i = 0; i < inp.nfluids; i++)
    st->fluid_param[i] = inp.fluid_param[i];

  // Acquire refs on per-species equation objects so repair_state and
  // check_inv can be driven from the source-update loop.
  for (int i = 0; i < GKYL_MAX_SPECIES; i++) st->eqn[i] = NULL;
  for (int i = 0; i < inp.nfluids; i++) {
    if (inp.eqn[i]) st->eqn[i] = gkyl_wv_eqn_acquire(inp.eqn[i]);
  }

  st->is_static  = inp.is_static;
  st->has_tetrad = inp.has_tetrad;

  st->analytic_spacetime = inp.analytic_spacetime;
  st->spacetime_gauge    = inp.spacetime_gauge;
  st->reinit_freq        = inp.reinit_freq;

  st->einstein_eqn = inp.einstein_eqn;

  return st;
}

void
gkyl_moment_spacetime_coupling_derive_products(
  const gkyl_moment_spacetime_coupling *st,
  double t_curr,
  const struct gkyl_range *update_range,
  const struct gkyl_array *einstein_state,
  struct gkyl_array *prods)
{
  // Phase B (dynamic Bona-Masso) is not yet implemented; this assert will
  // fire if a user tries to drive the path before A4-Phase-B lands.
  assert(st->einstein_eqn == NULL);
  assert(st->analytic_spacetime != NULL);
  (void)einstein_state;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(update_range, iter.idx);
    double *prods_row = gkyl_array_fetch(prods, cidx);

    double xc[GKYL_MAX_DIM] = { 0.0 };
    gkyl_rect_grid_cell_center(&st->grid, iter.idx, xc);

    double x = (st->ndim >= 1) ? xc[0] : 0.0;
    double y = (st->ndim >= 2) ? xc[1] : 0.0;
    double z = (st->ndim >= 3) ? xc[2] : 0.0;

    gkyl_moment_spacetime_coupling_fill_products_analytic(
      st->analytic_spacetime, t_curr, x, y, z, prods_row);
  }
}

void
gkyl_moment_spacetime_coupling_explicit_advance(
  const gkyl_moment_spacetime_coupling *st,
  double t_curr, double dt,
  const struct gkyl_range *update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES],
  const struct gkyl_array *prods,
  struct gkyl_array *gamma_eff_cache)
{
  int nfluids = st->nfluids;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(update_range, iter.idx);
    const double *prods_row = gkyl_array_cfetch(prods, cidx);

    // Per-cell γ_eff warm-start slot. The recovery's Picard iteration
    // reads this as the initial guess and writes back the converged value
    // for the next call. Same array is shared across the three SSP-RK3
    // stages so each stage benefits from the previous stage's refinement.
    double *gamma_eff_cell = (gamma_eff_cache != NULL)
      ? gkyl_array_fetch(gamma_eff_cache, cidx)
      : NULL;

    for (int s = 0; s < nfluids; s++) {
      enum gkyl_eqn_type type = st->fluid_param[s].type;
      // The geometric source terms are identical for the curved and tetrad
      // mod variants — only the wave-structure (flat-flux + GR-correction
      // factorization) differs there, not the ADM source contributions to D,
      // S_i, τ.
      if (type != GKYL_EQN_GR_EULER_MOD && type != GKYL_EQN_GR_EULER_TETRAD_MOD)
        continue;

      struct gkyl_gr_euler_eos eos = st->fluid_param[s].eos;

      double *f = gkyl_array_fetch(fluid[s], cidx);

      // Per-cell setter so check_inv / repair_state can fetch the cell's
      // row in the auxfields products array. Tag the repair context as
      // "source-step" so per-call-site counters in the equation tally
      // here rather than into the wave_prop bucket.
      const struct gkyl_wv_eqn *eqn = st->eqn[s];
      bool can_repair = eqn && eqn->repair_state_func && eqn->check_inv_func;
      if (can_repair) {
        if (eqn->set_cell_idx_func) eqn->set_cell_idx_func(eqn, iter.idx);
        ((struct gkyl_wv_eqn *)eqn)->cur_repair_ctx = 0;
      }

      // Single check + cascade-repair call. The cascade walks all three
      // admissibility constraints in safe order (D, then τ, then S²) in
      // one shot, so a single repair_state call always converges to
      // admissible — no loop needed. Convex-combo updates (f_stage2,
      // final f) do NOT need their own repair pass because the
      // admissibility set is convex (sum of convex constraints, with
      // |S|_γ a norm), so a convex combination of admissible states is
      // automatically admissible. Repair only fires after the three
      // forward-Euler stages.
      #define REPAIR_ONCE(qbuf) do {                                     \
        if (can_repair && !gkyl_wv_eqn_check_inv(eqn, (qbuf)))           \
          gkyl_wv_eqn_repair_state(eqn, (qbuf));                         \
      } while (0)

      // SSP-RK3: three forward Euler stages, repair after each one. The
      // two convex combinations (f_stage2, final f) are guaranteed
      // admissible by convexity of the admissibility set, so no repair
      // pass is needed there.
      double f_old[5], f_new[5], f_stage1[5], f_stage2[5];
      for (int j = 0; j < 5; j++) f_old[j] = f[j];

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        eos, t_curr, dt, prods_row, gamma_eff_cell, f_old, f_new);
      REPAIR_ONCE(f_new);
      for (int j = 0; j < 5; j++) f_stage1[j] = f_new[j];

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        eos, t_curr + dt, dt, prods_row, gamma_eff_cell, f_stage1, f_new);
      REPAIR_ONCE(f_new);
      for (int j = 0; j < 5; j++)
        f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        eos, t_curr + 0.5 * dt, dt, prods_row, gamma_eff_cell, f_stage2, f_new);
      REPAIR_ONCE(f_new);
      for (int j = 0; j < 5; j++)
        f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);

      #undef REPAIR_ONCE
    }
  }
}

void
gkyl_moment_spacetime_coupling_release(gkyl_moment_spacetime_coupling *st)
{
  for (int i = 0; i < GKYL_MAX_SPECIES; i++) {
    if (st->eqn[i]) gkyl_wv_eqn_release(st->eqn[i]);
  }
  gkyl_free(st);
}
