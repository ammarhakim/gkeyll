#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_mhd_mixture_rp {
  WV_MHD_MIXTURE_RP_LAX, // Default (Lax fluxes).
};

// Input context, packaged as a struct.
struct gkyl_wv_mhd_mixture_inp {
  int num_species; // Number of distinct species in mixture.
  double* gas_gamma_s; // Adiabatic indices for each species in mixture.
  double light_speed; // Speed of light.
  double b_fact; // Factor of speed of light for magnetic field correction.

  enum gkyl_wv_mhd_mixture_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new ideal MHD mixture equations object.
*
* @param num_species Number of distinct species in mixture.
* @param gas_gamma_s Adiabatic indices for each species in mixture.
* @param light_speed Speed of light.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the ideal MHD mixture equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_mhd_mixture_new(int num_species, double* gas_gamma_s, double light_speed, double b_fact, bool use_gpu);

/**
* Create a new ideal MHD mixture equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the ideal MHD mixture equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_mhd_mixture_inew(const struct gkyl_wv_mhd_mixture_inp* inp);

/**
* Get number of distinct species in mixture.
*
* @param wv ideal MHD mixture equations object.
* @return Number of distinct species in mixture.
*/
int
gkyl_wv_mhd_mixture_num_species(const struct gkyl_wv_eqn* wv);

/**
* Get adiabatic indices for each species in mixture.
*
* @param wv ideal MHD mixture equations object.
* @return Adiabatic indices for each species in mixture.
*/
double*
gkyl_wv_mhd_mixture_gas_gamma_s(const struct gkyl_wv_eqn* wv);

/**
* Get speed of light.
*
* @param eqn ideal MHD mixture equations object.
* @return Speed of light.
*/
double
gkyl_wv_mhd_mixture_light_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for magnetic field correction.
*
* @param eqn ideal MHD mixture equations object.
* @return Factor of speed of light for magnetic field correction.
*/
double
gkyl_wv_mhd_mixture_b_fact(const struct gkyl_wv_eqn* eqn);