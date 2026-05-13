#pragma once

// Private header for the modular tetrad-basis GR Euler equation object.
// Not for direct inclusion in user-facing code.

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

// Build the lower-triangular Cholesky factor L of a symmetric positive-
// definite γ_ij: γ_ij = L_ik L_jk. Also computes its inverse L^{-1}. From
// these the orthonormal triad ε^i_a = (L^{-1})^T satisfies
// γ_ij ε^i_a ε^j_b = δ_ab. We keep both because:
//   forward q_GR → q_tet  uses  S^a_tet = (L^{-1})^a_i (S_i_GR / √γ)
//   back    w_tet → w_GR   uses  w_i_GR_momentum = √γ · L_i^a · w^a_tet
// Caller's responsibility to skip excised cells (where γ_ij may be 0).
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_build_triad(const double g_ij[3][3],
  double L[3][3], double L_inv[3][3]);

// Transform a curved-frame conserved state q_GR (5-component, with √γ
// baked into D, momentum, τ) to a tetrad-frame conserved state q_tet
// (5-component, no √γ, momentum in the orthonormal-triad basis).
//
// The code's conserved momentum is CONTRAVARIANT (q[i+1] = √γ · ρhW² · v^i),
// confirmed by inspection: v_sq = γ_ij · v[i+1] · v[j+1] in prim_vars is the
// curved norm, which requires v[1..3] to be v^i (raised). So the forward
// transform is the inverse-triad applied to the contravariant momentum:
//   D_tet     = q_GR[0] / √γ
//   τ_tet     = q_GR[4] / √γ
//   S^a_tet   = E^a_i · (q_GR[i+1] / √γ)   where  E = L^T (inverse of triad)
//
// In components: q_tet[a+1] = L[i][a] · q_GR[i+1] / √γ.
// L is the Cholesky factor from build_triad (same metric).
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

// Transform a single tetrad-frame wave w_tet (5-component) back to the
// curved frame. With the contravariant-momentum convention (see
// q_to_tetrad), the back-transform is the triad ε = (L^{-1})^T applied
// to the tetrad-frame momentum:
//   w_GR[0]   = √γ · w_tet[0]                 (D-slot)
//   w_GR[i+1] = √γ · L_inv[a][i] · w_tet[a+1] (momentum: tetrad→coord-contra)
//   w_GR[4]   = √γ · w_tet[4]                 (τ-slot)
// L_inv is the inverse Cholesky factor from build_triad (same metric).
GKYL_CU_D
void
gkyl_gr_euler_tetrad_mod_wave_to_curved(const double w_tet[5],
  double sqrt_det, const double L_inv[3][3], double w_GR[5]);

// Transform a tetrad-frame wave speed to a coord-frame wave speed.
// For a wave propagating in the local x-direction (after rotate_to_local):
//   s_coord = α · (L^{-1})[0][0] · s_tet − β^1
// where (L^{-1})[0][0] = ε^1_(1) is the (1,1) entry of the orthonormal
// triad in the rotated frame, α is the lapse, and β^1 is the
// normal-component of the shift in the rotated frame.
//
// Derivation: the tetrad-frame worldline of a wave at speed s_tet is
//   x_tet = s_tet · τ_tet
// Converting to coord-frame coords (Eulerian observer with proper time
// τ_tet = α dt, position dx_coord = ε^1_(1) dx_tet − β^1 dt):
//   dx_coord/dt = ε^1_(1) · (dx_tet/dτ_tet) · (dτ_tet/dt) − β^1
//              = ε^1_(1) · s_tet · α − β^1
// For Minkowski (α=1, β=0, ε=I): s_coord = s_tet ✓.
GKYL_CU_D
double
gkyl_gr_euler_tetrad_mod_speed_to_curved(double s_tet,
  double lapse, double shift_x, const double L_inv[3][3]);
