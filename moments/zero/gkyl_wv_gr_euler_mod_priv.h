#pragma once

// Private header for the modular GR Euler equation object. Not for direct
// inclusion in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wave_spacetime.h>
#include <gkyl_wv_gr_euler.h>          // shares enum gkyl_wv_gr_euler_rp
#include <gkyl_wv_gr_euler_mod.h>

struct wv_gr_euler_mod {
  struct gkyl_wv_eqn eqn;                            // Base equation object.
  double gas_gamma;                                  // Adiabatic index.
  struct gkyl_range conf_range;                      // For indexing into prods.

  struct gkyl_wv_gr_euler_mod_auxfields auxfields;   // Pointer to products.

  // Cell indices captured by the optional set_interface_idx_func /
  // set_cell_idx_func setters before each callback fires. mutable because
  // wave_prop writes these through a const struct gkyl_wv_eqn*.
  int cur_idxl[GKYL_MAX_DIM];
  int cur_idxr[GKYL_MAX_DIM];
  int cur_cell_idx[GKYL_MAX_DIM];

  // Per-side locally-rotated spacetime scratch buffers. Filled by
  // rot_to_local (when the math port is in place) so subsequent waves/qfluct
  // calls see locally-rotated spacetime — matching how the packed equation
  // operates on locally-rotated q[5..66]. rot_call_parity toggles between
  // 0 (next rotate populates left) and 1 (next rotate populates right);
  // reset to 0 by set_interface_idx_func.
  double prodl_local[GKYL_GR_SP_NCOMP_BASE];
  double prodr_local[GKYL_GR_SP_NCOMP_BASE];
  int rot_call_parity;
};

// Free function for the reference-count callback.
void gkyl_gr_euler_mod_free(const struct gkyl_ref_count *ref);

// ---------------------------------------------------------------------------
// Helpers ported from wv_gr_euler.c but driven by the spacetime-products
// layout (see gkyl_moment_spacetime_products.h) instead of the packed 71-comp
// state vector. These take a 5-element hydro state plus a pointer to a row
// of the products array (locally rotated when called from waves/qfluct/flux,
// global when called from per-cell callbacks like max_speed). Declared here
// so unit tests can call them directly for equivalence comparison.
// ---------------------------------------------------------------------------

GKYL_CU_D
void
gkyl_gr_euler_mod_prim_vars(double gas_gamma, const double q[5],
  const double *prods, double v[5]);

GKYL_CU_D
void
gkyl_gr_euler_mod_flux(double gas_gamma, const double q[5],
  const double *prods, double flux[5]);

GKYL_CU_D
double
gkyl_gr_euler_mod_max_abs_speed(double gas_gamma, const double q[5],
  const double *prods);

// ---------------------------------------------------------------------------
// Interface-flux Lax helpers (IFACE_FLUX_PLAN.md §2).
//
// These build the curved Banyuls flux entirely from interface geometry
// (lapse, face-normal shift, sqrt(det γ_iface), γ_ij_iface). The only
// cell-local input is sqrt_det_cell, used to undensitize the conservatives.
// Primitives are recovered with iface γ_ij (Convention B → mom_sq =
// γ_ij·S^i·S^j) and the resulting (ρ, v^i, p, W, h) is iface-frame.
//
// Outputs are face-local-frame fluxes / speeds, ready for direct use in
// wave_prop's Lax-average construction.
// ---------------------------------------------------------------------------

GKYL_CU_D
void
gkyl_gr_euler_mod_banyuls_flux_iface(double gas_gamma,
  const double q[5], double sqrt_det_cell,
  const struct gkyl_wave_spacetime_iface *iface,
  double flux[5]);

GKYL_CU_D
double
gkyl_gr_euler_mod_max_abs_speed_iface(double gas_gamma,
  const double q[5], double sqrt_det_cell,
  const struct gkyl_wave_spacetime_iface *iface);
