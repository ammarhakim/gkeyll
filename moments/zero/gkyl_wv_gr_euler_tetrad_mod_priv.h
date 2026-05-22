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
#include <gkyl_wv_gr_euler_tetrad.h>          // shares enum gkyl_wv_gr_euler_tetrad_rp
#include <gkyl_wv_gr_euler_tetrad_mod.h>

struct wv_gr_euler_tetrad_mod {
  struct gkyl_wv_eqn eqn;                                  // Base equation object.
  struct gkyl_gr_euler_eos eos;                            // Equation of state (IDEAL or MATHEWS_TAUB).
  struct gkyl_range conf_range;                            // For indexing into prods.

  struct gkyl_wv_gr_euler_tetrad_mod_auxfields auxfields;  // Pointer to products.

  // Cell indices captured by the optional set_interface_idx_func /
  // set_cell_idx_func setters before each callback fires. mutable because
  // wave_prop writes these through a const struct gkyl_wv_eqn*.
  int cur_idxl[GKYL_MAX_DIM];
  int cur_idxr[GKYL_MAX_DIM];
  int cur_cell_idx[GKYL_MAX_DIM];

  // Per-side locally-rotated spacetime scratch buffers, mirroring the regular
  // mod variant. The tetrad equation reads the same spacetime block (lapse,
  // shift, γ_ij, γ^ij, √γ, excision) and does not require stored basis
  // vectors — the "tetrad" name refers to a flux-factorization strategy
  // (compute flat-space SR flux, then apply a GR correction with √γ·α and
  // shift), not to the primitive-variable representation. See
  // gkyl_gr_euler_tetrad_mod_flux + gkyl_gr_euler_tetrad_mod_flux_correction.
  double prodl_local[GKYL_GR_SP_NCOMP_BASE];
  double prodr_local[GKYL_GR_SP_NCOMP_BASE];
  int rot_call_parity;

  // Diagnostic counters for repair_state — split by call-site context
  // (eqn->cur_repair_ctx: 0 = source-step, 1 = wave_prop). Indexed by
  // gkyl_gr_euler_admissibility_status value (entries 1..3 used; entry 0
  // is OK and never incremented). Mutable because callers reach in
  // through a const struct gkyl_wv_eqn*.
  uint64_t repair_count_source[4];
  uint64_t repair_count_wave_prop[4];
};

// Free function for the reference-count callback.
void gkyl_gr_euler_tetrad_mod_free(const struct gkyl_ref_count *ref);

// ---------------------------------------------------------------------------
// Helpers ported from wv_gr_euler_tetrad.c but driven by the spacetime-products
// layout (see gkyl_moment_spacetime_products.h) instead of the packed 71-comp
// state vector. The math is split into a flat-space SR step and a GR
// correction, matching the packed tetrad factorization. Declared here so
// unit tests can call them directly for equivalence comparison.
// ---------------------------------------------------------------------------

GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_prim_vars(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods, double v[5]);

// Compute the flat-space (special-relativistic) flux from primitives. Uses
// Cartesian dot products on velocity (W_flat = 1/sqrt(1 - vᵢvᵢ)) with no
// metric prefactors. Mirrors gkyl_gr_euler_tetrad_flux in packed.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_flux(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods, double flux_sr[5]);

// Apply the GR coordinate-transformation correction to a flat-space SR flux:
// multiply by α·√γ and replace (vx, W_flat) with (vx - βˣ/α, W_curved).
// Mirrors gkyl_gr_euler_tetrad_flux_correction in packed.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_flux_correction(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods, const double flux_sr[5],
  double flux_gr[5]);

GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_max_abs_speed(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods);

// ---------------------------------------------------------------------------
// Modular tetrad-Roe pipeline. The "tetrad mod" Roe is implemented as a
// composition of four pieces — triad construction, forward state transform,
// pure Minkowski Roe in the tetrad frame, back-transform of waves and
// speeds. Each piece has a clean contract and is unit-tested independently
// (see ctest_wv_gr_euler_tetrad_mod's gr_euler_tetrad_mod_triad_*,
// gr_euler_tetrad_mod_q_to_tetrad_*, gr_euler_tetrad_mod_sr_roe_*, and
// gr_euler_tetrad_mod_wave_to_curved_* tests).
//
// The metric used at the interface is the arithmetic mean of γ_ij on the
// left and right cells (analogously for lapse, shift). This gives a
// symmetric L↔R treatment and matches the standard Roe-at-interface
// convention. The K-weighted average reserved for the Roe primitive
// averages happens inside sr_roe_minkowski; it's a different thing.
// ---------------------------------------------------------------------------

// Build the lower-triangular Cholesky factor L of γ: γ = L L^T, plus L^{-1}.
// This is the "Cholesky-on-γ" decomposition, matching the J.6 construction
// in the GR Kinetics notes (which uses h_ij directly for the momentum
// mapping p_i = √h_rr · p̂_r etc. in eq 915-917). The alternative J.10
// construction (V V^T = γ^{-1}) gives a geometrically equivalent tetrad
// with the same off-diagonal-γ flux-jump residual — that residual is a
// 1D-sweep limitation independent of the triad construction (see the
// wave_to_curved docstring below).
//
// Caller's responsibility to skip excised cells.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_build_triad(const double g_ij[3][3],
  double L[3][3], double L_inv[3][3]);

// Transform a curved-frame conserved state q_GR (Convention A: q[1..3] is
// genuine covariant momentum S_i) to a tetrad-frame conserved state q_tet
// (5-component, no √γ, momentum in the orthonormal-triad basis).
//
// Forward transform uses L_inv = L^{-1} (inverse Cholesky factor):
//   D_tet     = q_GR[0] / √γ
//   τ_tet     = q_GR[4] / √γ
//   S_a_tet   = L_inv[a][i] · (q_GR[i+1] / √γ)   (sum over i)
//
// Derivation: with the orthonormal triad basis vectors ε^i_a = L_inv[a][i]
// satisfying γ_ij·ε^i_a·ε^j_b = δ_ab, the projection of a covariant
// momentum onto the triad is S_a_tet = ε^i_a·S_i = L_inv[a][i]·S_i. Then
// Σ_a (S_a_tet)² = γ^{ij}·S_i·S_j matches the Lorentz scalar in the
// curved frame.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_q_to_tetrad(const double q_GR[5],
  double sqrt_det, const double L_inv[3][3], double q_tet[5]);

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
GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_sr_roe_minkowski(double gas_gamma,
  const double ql_tet[5], const double qr_tet[5],
  double waves_tet[3 * 5], double speeds[3]);

// Pure Minkowski SR HLL with Davis/Einfeldt wave-speed bracket.
// Returns the two-wave decomposition with q_HLL as the intermediate
// state. Provably admissibility-preserving in flat space (Mignone-Bodo).
GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_sr_hll_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  double waves_tet[2 * 5], double speeds[2]);

// Pure Minkowski SR Lax-Friedrichs with symmetric ±amax envelope. More
// diffusive than HLL but admissibility-preserving on admissible inputs.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_sr_lax_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
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
struct gkyl_gr_euler_tetrad_mod_hllc_diag {
  int did_fallback;
  int fallback_reason;
  double lambda_L;
  double lambda_R;
  double lambda_star;
};

GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_sr_hllc_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  double waves_tet[3 * 5], double speeds[3],
  struct gkyl_gr_euler_tetrad_mod_hllc_diag *diag);

// Build a Gram-Schmidt-on-γ⁻¹ triad: e_0 aligned with the contravariant
// x-direction, e_1, e_2 orthogonalized in γ. M[i][a] = e_a^i, M_inv =
// M^T·γ. Eliminates the v_tet^x ↔ v^y, v^z mixing seen with Cholesky
// for non-diagonal γ — see SESSION_NOTES_2.md §12.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_build_triad_contravariant_x(
  const double g_ij[3][3], const double inv_g[3][3],
  double M[3][3], double M_inv[3][3]);

// Forward transform of Convention-A covariant momentum onto the
// contravariant-x triad: S_tet^a = M_inv[a][i]·γ^{ij}·S_j/√γ. For a=0
// this gives (1/√γ^{xx})·S^x, so v_tet^0 = v^x/√γ^{xx} (clean).
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(const double q_GR[5],
  double sqrt_det, const double inv_g[3][3], const double M_inv[3][3],
  double q_tet[5]);

// Back-transform of waves to curved-frame Convention A. The coord-x
// momentum slot receives only the a=0 tetrad wave (mirror of the
// forward-clean property in the contravariant-x construction).
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_wave_to_curved_contra(const double w_tet[5],
  double sqrt_det, const double M_inv[3][3], double w_GR[5]);

// Speed back-transform: s_coord = α·√γ^{xx}·s_tet − β^x.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_speed_to_curved_contra(double s_tet,
  double lapse, double shift_x, double inv_gxx);

// Back-transform: inverse of forward S_a_tet = L_inv[a][i]·S_i. With
// γ = L L^T:
//   w_GR[0]   = √γ · w_tet[0]                  (D-slot)
//   w_GR[i+1] = √γ · L[i][a] · w_tet[a+1]      (sum over a)
//   w_GR[4]   = √γ · w_tet[4]                  (τ-slot)
// With L lower triangular, the coord-x momentum slot only sees tetrad-x;
// coord-y picks up tetrad-x and -y; coord-z all three tetrad axes. This
// cross-coupling carries the off-diagonal-γ contributions back into the
// covariant momentum slots.
//
// The strict flux-jump identity ∑ s·w = Δf_GR does not hold in curved γ for
// this 1D-sweep tetrad-Roe — see Gorard, Hakim, Juno, TenBarge 2025
// (arXiv:2410.02549) Sec. 4 for the structural argument.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_wave_to_curved(const double w_tet[5],
  double sqrt_det, const double L[3][3], double w_GR[5]);

// Transform a tetrad-frame wave speed to a coord-frame wave speed.
// For a wave propagating along the local-x cell normal:
//   s_coord = α · L_inv[0][0] · s_tet − β^1
// where L_inv[0][0] = 1/√γ_xx for diagonal γ is the (0,0) entry of the
// inverse Cholesky factor (= ε^0_0, the basis-vector scaling of e_0
// along ∂_0).
//
// Derivation: the tetrad-frame worldline of a wave at speed s_tet has
// coord motion dx^i = α · ε^i_0 · s_tet · dt − β^i dt. For i=0,
// ε^0_0 = L_inv[0][0] is the only nonzero entry (L_inv lower triangular),
// giving dx^0/dt = α · L_inv[0][0] · s_tet − β^0.
// For Minkowski (α=1, β=0, L_inv=I): s_coord = s_tet ✓.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_speed_to_curved(double s_tet,
  double lapse, double shift_x, const double L_inv[3][3]);
