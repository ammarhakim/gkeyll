#pragma once

// Private header for the modular tetrad-basis GR Euler equation object.
// Not for direct inclusion in user-facing code.
//
// Convention-A conserved variables (genuine covariant momentum):
//
//   q[0]   = √γ · ρW                     (D, scalar rest-mass current)
//   q[i+1] = √γ · γ_ij · ρhW² · v^j      (S_i, covariant momentum density)
//   q[4]   = √γ · (ρhW² − p − ρW)        (τ, energy minus rest-mass)
//
// The Banyuls inversion uses |S|² = γ^{ij}·S_i·S_j and recovers the
// contravariant velocity v^i = γ^{ij}·S_j / (ρhW²); both lookups go
// through gkyl_gr_euler_recover_primitives in
// gkyl_wv_gr_euler_prim_priv.h. The packed-tetrad path
// (wv_gr_euler_tetrad.c) and the non-tetrad mod path (wv_gr_euler_mod.c)
// remain on the contravariant-momentum convention; only the modular
// tetrad path defined here uses Convention A.

#include <math.h>
#include <stdint.h>
#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_gr_euler_prim_priv.h>       // struct gkyl_gr_euler_eos
#include <gkyl_wv_gr_euler_tetrad.h>      // enum gkyl_wv_gr_euler_tetrad_rp lives here now

// ---------------------------------------------------------------------------
// Generic tetrad-first scaffolding. Every tetrad-first Riemann solver
// (Lax, HLL, HLLC, Roe) shares the same 6-step Gorard+25 procedure:
// average geometry at the interface, build a Gram-Schmidt-on-γ⁻¹ triad,
// forward-transform L+R states, run a flat-SR Riemann solve in the
// tetrad frame, back-transform waves and speeds. Only step 4 (the SR
// kernel) differs across variants. Each variant's choice of SR kernel +
// num_waves + excision policy is stored on the equation object so the
// shared dispatcher can pick the right pieces at runtime — and so the
// equation can be sent to device for GPU dispatch without per-variant
// branching.
// ---------------------------------------------------------------------------

// SR-kernel signature shared by sr_lax/hll/hllc/roe_minkowski. Takes
// tetrad-frame conservatives, returns waves/speeds in the tetrad frame
// (Cartesian flat, no √γ). Returns the tetrad-frame max-abs speed (CFL
// uses the curved-frame max from the back-transformed speeds).
typedef double (*gkyl_sr_riemann_kernel_t)(
  struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double *waves_tet, double *speeds_tet);

// Excision-boundary policy. Differs across solvers because HLLC's
// star-state construction doesn't degrade gracefully to one-side vacuum.
enum gkyl_tetrad_excision_policy {
  // Lax/HLL/Roe: feed q_tet=0 on the excised side and run the SR kernel
  // using the active cell's geometry. SR kernels' zero-state short-
  // circuits handle the resulting "vacuum" input cleanly.
  GKYL_TETRAD_EXCISION_ZERO_VACUUM,
  // HLLC: any-side excised → return zero waves with sentinel speeds.
  GKYL_TETRAD_EXCISION_SHORT_CIRCUIT,
};

struct wv_gr_euler_tetrad {
  struct gkyl_wv_eqn eqn;                                  // Base equation object.
  struct gkyl_gr_euler_eos eos;                            // Equation of state (IDEAL or APPROXIMATE_SYNGE).
  struct gkyl_range conf_range;                            // For indexing into prods.

  struct gkyl_wv_gr_euler_tetrad_auxfields auxfields;  // Pointer to products + per-callsite status.

  // Riemann-solver configuration. Set once at construction time based on
  // rp_type and read by the generic tetrad-first dispatchers (wave/qfluct).
  // Storing these on the equation object — rather than passing them as
  // arguments to per-variant wrappers — lets the GPU code path send the
  // equation pointer to the device and dispatch the SR kernel from there.
  gkyl_sr_riemann_kernel_t sr_kernel;
  int num_waves;
  enum gkyl_tetrad_excision_policy excision_policy;

  // Cell indices captured by the optional set_interface_idx_func /
  // set_cell_idx_func setters before each callback fires. mutable because
  // wave_prop writes these through a const struct gkyl_wv_eqn*.
  int cur_idxl[GKYL_MAX_DIM];
  int cur_idxr[GKYL_MAX_DIM];
  int cur_cell_idx[GKYL_MAX_DIM];

  // Per-side locally-rotated spacetime scratch buffers, mirroring the regular
  // mod variant. The tetrad equation reads the same spacetime block (lapse,
  // shift, γ_ij, γ^ij, √γ, excision) and does not require stored basis
  // vectors — the "tetrad" name refers to the HIGH_ORDER Riemann-solver
  // strategy (rotate state into a locally-flat orthonormal frame at the
  // interface, solve SR Riemann, back-transform), not to the primitive-
  // variable representation. LOW_ORDER curved Lax and F-wave callbacks
  // use gkyl_gr_euler_banyuls_flux_cell directly.
  double prodl_local[GKYL_GR_SP_NCOMP_BASE];
  double prodr_local[GKYL_GR_SP_NCOMP_BASE];
  int rot_call_parity;
};

// Free function for the reference-count callback.
void gkyl_gr_euler_tetrad_free(const struct gkyl_ref_count *ref);

// ---------------------------------------------------------------------------
// Helpers ported from wv_gr_euler_tetrad.c but driven by the spacetime-products
// layout (see gkyl_moment_spacetime_products.h) instead of the packed 71-comp
// state vector. The math is split into a flat-space SR step and a GR
// correction, matching the packed tetrad factorization. Declared here so
// unit tests can call them directly for equivalence comparison.
// ---------------------------------------------------------------------------

// Recovery helpers. Each accepts an optional prim_status pointer for
// per-callsite instrumentation (NULL skips all accumulation).
GKYL_CU_D
void
gkyl_gr_euler_tetrad_prim_vars(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods,
  struct gkyl_gr_euler_prim_status *stat, double v[5]);

// Direct cell-centered Banyuls flux in the curved coordinate basis. This
// is the canonical reference flux for every flux-jump identity in this
// equation object — both the LOW_ORDER curved-Lax wave_lax_curved and the
// F-wave flux_jump_func use it on each side, and the HIGH_ORDER tetrad
// path's amdq+apdq is checked against it (interface-averaged form) in the
// ctest harness.
//
// Single-pass: recovers (ρ, v^i, p, W, h) with the cell-centered γ⁻¹ and
// √γ, then writes the Banyuls flux directly — α·√γ prefactor, β^x/α
// shift, v_l[i] = γ_ij·v^j for the lowered momentum slot. No SR-flat
// detour, no W_flat / W_curved Valencia ratio.
GKYL_CU_D
void
gkyl_gr_euler_banyuls_flux_cell(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods,
  struct gkyl_gr_euler_prim_status *stat,
  double flux[5]);

GKYL_CU_D
double
gkyl_gr_euler_tetrad_max_abs_speed(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods,
  struct gkyl_gr_euler_prim_status *stat);

// Run a pure SR (Minkowski) Roe wave decomposition in the tetrad frame.
// Inputs: gas_gamma, ql_tet, qr_tet (tetrad-frame conserved variables —
// no √γ, momentum in orthonormal basis so the inner product is Cartesian).
// Outputs: waves_tet[3*5] (three waves × 5 components), speeds[3] (the
// Minkowski-frame wave speeds), and return value is max |s_k|.
//
// The eigenstructure is the same Eulderink-Mellema construction the
// curved-frame Roe uses, but specialized to g^μν = η^μν (g^00 = −1,
// off-diagonals zero, spatial part identity). This is what makes it a
// "true Roe" — the SR Jacobian satisfies A_SR · ∆q_tet = ∆f_SR exactly,
// so ∑ s_k · w_k = ∆f_SR_flat holds at floating-point precision.
//
// IDEAL-GAS ONLY. The c_minus/c_plus/s² coefficients in this eigenstructure
// are derived from the ideal-gas Jacobian (h = 1 + γ/(γ-1)·p/ρ); the
// Mathews-Taub Jacobian has a different algebraic form. The mod-tetrad
// constructor rejects rp_type=ROE when eos.type != IDEAL.
//
// stat (optional, may be NULL): per-callsite prim_status accumulator
// pointer reserved for future instrumentation symmetry with the other
// SR kernels. Roe's inlined primitive recovery does not yet populate it.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_sr_roe_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[3 * 5], double speeds[3]);

// Pure Minkowski SR HLL with Davis/Einfeldt wave-speed bracket.
// Returns the two-wave decomposition with q_HLL as the intermediate
// state. Provably admissibility-preserving in flat space (Mignone-Bodo).
//
// stat (optional, may be NULL): per-callsite prim_status accumulator.
// Both per-side recovery calls write into the same status object.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_sr_hll_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[2 * 5], double speeds[2]);

// Pure Minkowski SR Lax-Friedrichs with symmetric ±amax envelope. More
// diffusive than HLL but admissibility-preserving on admissible inputs.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_sr_lax_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[2 * 5], double speeds[2]);

// Pure Minkowski SR HLLC (Mignone & Bodo 2005, MNRAS 364, 126). Three
// waves (left acoustic, contact, right acoustic) with star-state
// construction. Restores the contact wave that HLL averages over,
// giving sharper resolution of contact discontinuities AND restoring
// τ-positivity from admissible left/right inputs (which HLL does NOT
// guarantee in the τ/D ≪ 1 regime — see SESSION_NOTES_2.md §17).
//
// Implementation notes for the production-grade variant:
//   - Davis (1988) wave-speed bracket per side (eqs 21-23 of MB05).
//   - Stable quadratic for λ* using citardauq form when −b > 0 to
//     avoid catastrophic cancellation near uniform states.
//   - p* from the corrected Rankine-Hugoniot relation
//     p* = (A·λ* − B) / (1 − λ·λ*)  where A = λE − m_x, B = m_x(λ-v_x) − p.
//     This differs from MB05 eq (17) as printed [denominator (1 + λλ*)]
//     in the sign of the λλ* term; the corrected form satisfies
//     p*_L = p*_R = p in the trivial Riemann limit ql = qr.
//   - In degenerate cold-gas regimes (c_s ≈ 0 ⇒ λ_L ≈ λ_R) where the
//     λ* quadratic returns a root outside [λ_L, λ_R] or arbitrarily
//     close to an outer wave (1/(λ_R − λ*) blows up), gracefully
//     degrade to HLL: collapse the contact wave to zero amplitude and
//     let the two outer waves carry the entire HLL-averaged jump.
//
// Caveat: HLL itself is NOT τ-positivity-preserving in our setup
// (verified empirically on the BHL bow-shock — see
// SESSION_NOTES_2.md §17 and ctest test_small_tau_over_D_hll).
// The HLL fallback degrades gracefully but inherits HLL's known
// τ-positivity weakness in the τ/D ≪ 1 regime.
// Diagnostic outparam for sr_hllc_minkowski. Pass NULL to ignore.
// did_fallback: 1 if the cold-gas / degenerate-fan path triggered (HLLC
//   degraded to HLL on this interface), 0 if the full HLLC star-state
//   was used. Useful for diagnosing whether observed τ-violations stem
//   from HLLC's own star-state construction or from the HLL fallback
//   inheriting HLL's known weakness.
// lambda_L, lambda_R: tetrad-frame Davis bracket bounds.
// lambda_star: contact-wave speed (as computed by the quadratic, before
//   the clamp/fallback decision; NaN if the computation never reached
//   that point — e.g. lam_diff < 1e-14 short-circuit).
// fallback_reason: 0 if no fallback, 1 if lam_diff degenerate, 2 if
//   lambda_star not finite, 3 if |lambda_L − lambda_star| < tol
//   (would blow up the U_L* RH inverse), 4 if |lambda_R − lambda_star|
//   < tol (analogous for U_R*). Note 3/4 catch numerical proximity to
//   the bracket edges, NOT lambda_star outside the bracket — the latter
//   is normal for supersonic flow.
// HLLC star-state diagnostics live on the per-callsite prim_status struct
// (stat->hllc). sr_hllc_minkowski writes them when stat != NULL; the other
// SR kernels (Lax/HLL/Roe) leave the hllc sub-struct at zero.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_sr_hllc_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[3 * 5], double speeds[3]);

// Invert a 3×3 symmetric matrix in place. Returns det(g); fills inv_g
// with g^{-1}. Used at the interface to derive a *consistent* (inv_g,
// sqrt(det)) pair from g_iface, avoiding the O((Δγ)²) wave-sum residual
// that arises when these three are independently averaged from L/R cells.
// See TETRAD_REFACTOR_PLAN.md Phase 0(b) Fix 2.
static inline double
gkyl_gr_euler_tetrad_invert_metric_3x3(const double g[3][3],
  double inv_g[3][3])
{
  double det = g[0][0] * (g[1][1]*g[2][2] - g[1][2]*g[1][2])
             - g[0][1] * (g[0][1]*g[2][2] - g[1][2]*g[0][2])
             + g[0][2] * (g[0][1]*g[1][2] - g[1][1]*g[0][2]);
  double inv_det = 1.0 / det;
  inv_g[0][0] =  (g[1][1]*g[2][2] - g[1][2]*g[1][2]) * inv_det;
  inv_g[0][1] = -(g[0][1]*g[2][2] - g[1][2]*g[0][2]) * inv_det;
  inv_g[0][2] =  (g[0][1]*g[1][2] - g[1][1]*g[0][2]) * inv_det;
  inv_g[1][0] = inv_g[0][1];
  inv_g[1][1] =  (g[0][0]*g[2][2] - g[0][2]*g[0][2]) * inv_det;
  inv_g[1][2] = -(g[0][0]*g[1][2] - g[0][1]*g[0][2]) * inv_det;
  inv_g[2][0] = inv_g[0][2];
  inv_g[2][1] = inv_g[1][2];
  inv_g[2][2] =  (g[0][0]*g[1][1] - g[0][1]*g[0][1]) * inv_det;
  return det;
}

// Build a Gram-Schmidt-on-γ⁻¹ triad: e_0 aligned with the contravariant
// x-direction, e_1, e_2 orthogonalized in γ. M[i][a] = e_a^i, M_inv =
// M^T·γ. Eliminates the v_tet^x ↔ v^y, v^z mixing seen with Cholesky
// for non-diagonal γ — see SESSION_NOTES_2.md §12.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_build_triad_contravariant_x(
  const double g_ij[3][3], const double inv_g[3][3],
  double M[3][3], double M_inv[3][3]);

// Forward transform of Convention-A covariant momentum onto the
// contravariant-x triad: S_tet^a = M_inv[a][i]·γ^{ij}·S_j/√γ. For a=0
// this gives (1/√γ^{xx})·S^x, so v_tet^0 = v^x/√γ^{xx} (clean).
GKYL_CU_D
void
gkyl_gr_euler_tetrad_q_to_tetrad_contra(const double q_GR[5],
  double sqrt_det, const double inv_g[3][3], const double M_inv[3][3],
  double q_tet[5]);

// Back-transform of waves to curved-frame Convention A. The coord-x
// momentum slot receives only the a=0 tetrad wave (mirror of the
// forward-clean property in the contravariant-x construction).
GKYL_CU_D
void
gkyl_gr_euler_tetrad_wave_to_curved_contra(const double w_tet[5],
  double sqrt_det, const double M_inv[3][3], double w_GR[5]);

// Speed back-transform: s_coord = α·√γ^{xx}·s_tet − β^x.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_speed_to_curved_contra(double s_tet,
  double lapse, double shift_x, double inv_gxx);

