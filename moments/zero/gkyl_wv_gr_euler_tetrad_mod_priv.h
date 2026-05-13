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
