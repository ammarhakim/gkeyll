#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_mhd_rp {
  WV_MHD_RP_LAX = 0, // Default (Lax fluxes).
};

// Input context, packaged as a struct.
struct gkyl_wv_mhd_inp {
  double gas_gamma; // Adiabatic index.

  double light_speed; // Speed of light.
  double b_fact; // Factor of speed of light for magnetic field correction.

  enum gkyl_wv_mhd_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new ideal MHD equations object.
*
* @param gas_gamma Adiabatic index.
* @param light_speed Speed of light.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the ideal MHD equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_mhd_new(double gas_gamma, double light_speed, double b_fact, bool use_gpu);

/**
* Create a new ideal MHD equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the ideal MHD equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_mhd_inew(const struct gkyl_wv_mhd_inp* inp);

/**
* Get adiabatic index.
*
* @param wv Ideal MHD equations object.
* @return Adiabatic index.
*/
double
gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* wv);

/**
* Get speed of light.
*
* @param eqn Ideal MHD equations object.
* @return Speed of light.
*/
double
gkyl_wv_mhd_light_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for magnetic field correction.
*
* @param eqn Ideal MHD equations object.
* @return Factor of speed of light for magnetic field correction.
*/
double
gkyl_wv_mhd_b_fact(const struct gkyl_wv_eqn* eqn);