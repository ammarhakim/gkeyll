#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_tov_rp {
  WV_GR_TOV_RP_LAX = 0, // Default (Lax fluxes).
  WV_GR_TOV_RP_ROE,
  WV_GR_TOV_RP_HLL,
};

// Input context, packaged as a struct.
struct gkyl_wv_gr_tov_inp {
  double gas_gamma; // Adiabatic index.
  double kappa; // Stress-energy prefactor in the Einstein field equations.

  enum gkyl_wv_gr_tov_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new coupled fluid-Einstein equations object in plane-symmetric spacetimes.
*
* @param gas_gamma Adiabatic index.
* @param kappa Stress-energy prefactor in the Einstein field equations.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the coupled fluid-Einstein equations object in plane-symmetric spacetimes.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_tov_new(double gas_gamma, double kappa, bool use_gpu);

/**
* Create a new coupled fluid-Einstein equations object in plane-symmetric spacetimes, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the coupled fluid-Einstein equations object in plane-symmetric spacetimes.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_tov_inew(const struct gkyl_wv_gr_tov_inp* inp);

/**
* Get adiabatic index.
*
* @param eqn Coupled fluid-Einstein equations object in plane-symmetric spacetimes.
* @return Adiabatic index.
*/
double
gkyl_wv_gr_tov_gas_gamma(const struct gkyl_wv_eqn* eqn);

/**
* Get stress-energy prefactor in the Einstein field equations.
*
* @param eqn Coupled fluid-Einstein equations object in plane-symmetric spacetimes.
* @return Stress-energy prefactor in the Einstein field equations.
*/
double
gkyl_wv_gr_tov_kappa(const struct gkyl_wv_eqn* eqn);