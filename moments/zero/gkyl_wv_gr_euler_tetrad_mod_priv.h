#pragma once

// Private header for the modular tetrad-basis GR Euler equation object.
// Not for direct inclusion in user-facing code.
//
// Conserved-variable convention (shared with wv_gr_euler.c, wv_gr_euler_mod.c,
// wv_gr_euler_tetrad.c):
//
//   q[0]   = √γ · ρW
//   q[i+1] = √γ · ρhW² · v^i           ← contravariant velocity in slot i
//   q[4]   = √γ · (ρhW² − p − ρW)
//
// The "i" in q[i+1] is a slot label and the stored quantity is contravariant
// momentum density. For static γ_ij this is algebraically equivalent to the
// standard Banyuls covariant-momentum form: the metric is folded into the
// conserved-variable definition, and the flux is written so that v^i appears
// directly in the j-th momentum slot (no metric raise/lower at flux time).
// The recovery polynomial uses the Lorentz scalar |S|² = γ_ij·S^i·S^j; see
// gkyl_gr_euler_tetrad_mod_prim_vars in the .c file.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_gr_euler_tetrad.h>          // shares enum gkyl_wv_gr_euler_tetrad_rp
#include <gkyl_wv_gr_euler_tetrad_mod.h>

struct wv_gr_euler_tetrad_mod {
  struct gkyl_wv_eqn eqn;                                  // Base equation object.
  double gas_gamma;                                        // Adiabatic index.
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
gkyl_gr_euler_tetrad_mod_prim_vars(double gas_gamma, const double q[5],
  const double *prods, double v[5]);

// Compute the flat-space (special-relativistic) flux from primitives. Uses
// Cartesian dot products on velocity (W_flat = 1/sqrt(1 - vᵢvᵢ)) with no
// metric prefactors. Mirrors gkyl_gr_euler_tetrad_flux in packed.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_flux(double gas_gamma, const double q[5],
  const double *prods, double flux_sr[5]);

// Apply the GR coordinate-transformation correction to a flat-space SR flux:
// multiply by α·√γ and replace (vx, W_flat) with (vx - βˣ/α, W_curved).
// Mirrors gkyl_gr_euler_tetrad_flux_correction in packed.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_flux_correction(double gas_gamma, const double q[5],
  const double *prods, const double flux_sr[5], double flux_gr[5]);

GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_max_abs_speed(double gas_gamma, const double q[5],
  const double *prods);

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

// Transform a curved-frame conserved state q_GR to a tetrad-frame conserved
// state q_tet (5-component, no √γ, momentum in the orthonormal-triad basis).
//
// Forward transform uses L^T (transpose of the Cholesky factor of γ):
//   D_tet     = q_GR[0] / √γ
//   τ_tet     = q_GR[4] / √γ
//   S_a_tet   = L[i][a] · (q_GR[i+1] / √γ)   (sum over i)
//
// Derivation: forward of contravariant momentum needs M with M^T M = γ so
// that Σ_a (Mv)_a² = γ_ij v^i v^j is preserved. With γ = L L^T (Cholesky),
// M = L^T. The invariant Σ_a (S_a_tet)² = γ_ij·S^i·S^j (with S^i = ρhW²·v^i)
// then holds in the tetrad frame.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_q_to_tetrad(const double q_GR[5],
  double sqrt_det, const double L[3][3], double q_tet[5]);

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
GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_sr_roe_minkowski(double gas_gamma,
  const double ql_tet[5], const double qr_tet[5],
  double waves_tet[3 * 5], double speeds[3]);

// Back-transform: inverse of forward S_tet = L^T · S^GR. With L_inv = L^{-1}:
//   w_GR[0]   = √γ · w_tet[0]                  (D-slot)
//   w_GR[i+1] = √γ · L_inv[a][i] · w_tet[a+1]  (sum over a)
//   w_GR[4]   = √γ · w_tet[4]                  (τ-slot)
// With L_inv lower triangular, the coord-x momentum slot picks up all three
// tetrad waves through the first column of L_inv; the coord-z slot only sees
// tetrad-z. This cross-coupling carries the off-diagonal-γ contributions
// back into each coord-frame momentum component.
//
// The strict flux-jump identity ∑ s·w = Δf_GR does not hold in curved γ for
// this 1D-sweep tetrad-Roe — see Gorard, Hakim, Juno, TenBarge 2025
// (arXiv:2410.02549) Sec. 4 for the structural argument.
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_wave_to_curved(const double w_tet[5],
  double sqrt_det, const double L_inv[3][3], double w_GR[5]);

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
