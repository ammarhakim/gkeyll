#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_mhd_rgfm_rp {
  WV_MHD_RGFM_RP_HLL = 0, // Default (HLL fluxes).
  WV_MHD_RGFM_RP_LAX
};

// Input context, packaged as a struct.
struct gkyl_wv_mhd_rgfm_inp {
  int num_species; // Number of distinct species in the domain.
  double* gas_gamma_s; // Adiabatic indices for each species in the domain.
  int reinit_freq; // Reinitialization frequency for the level set.
  double surface_tension; // Surface tension force for the level set.

  double light_speed; // Speed of light.
  double b_fact; // Factor of speed of light for magnetic field correction.

  enum gkyl_wv_mhd_rgfm_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new ideal MHD Riemann ghost fluid equations object.
*
* @param num_species Number of distinct species in the domain.
* @param gas_gamma_s Adiabatic indices for each species in the domain.
* @param light_speed Speed of light.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param reinit_freq Reinitialization frequency for the level set.
* @param surface_tension Surface tension force for the level set.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the ideal MHD Riemann ghost fluid equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_mhd_rgfm_new(int num_species, double* gas_gamma_s, double light_speed, double b_fact, int reinit_freq, double surface_tension, bool use_gpu);

/**
* Create a new ideal MHD Riemann ghost fluid equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the ideal MHD Riemann ghost fluid equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_mhd_rgfm_inew(const struct gkyl_wv_mhd_rgfm_inp* inp);

/**
* Get number of distinct species in the domain.
*
* @param wv Ideal MHD Riemann ghost fluid equations object.
* @return Number of distinct species in the domain.
*/
int
gkyl_wv_mhd_rgfm_num_species(const struct gkyl_wv_eqn* wv);

/**
* Get adiabatic indices for each species in the domain.
*
* @param wv Ideal MHD Riemann ghost fluid equations object.
* @return Adiabatic indices for each species in the domain.
*/
double*
gkyl_wv_mhd_rgfm_gas_gamma_s(const struct gkyl_wv_eqn* wv);

/**
* Get speed of light.
*
* @param eqn Ideal MHD Riemann ghost fluid equations object.
* @return Speed of light.
*/
double
gkyl_wv_mhd_rgfm_light_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for magnetic field correction.
*
* @param eqn Ideal MHD Riemann ghost fluid equations object.
* @return Factor of speed of light for magnetic field correction.
*/
double
gkyl_wv_mhd_rgfm_b_fact(const struct gkyl_wv_eqn* eqn);

/**
* Get reinitialization frequency for the level set.
*
* @param wv Ideal MHD Riemann ghost fluid equations object.
* @return Reinitialization frequency for the level set.
*/
int
gkyl_wv_mhd_rgfm_reinit_freq(const struct gkyl_wv_eqn* wv);

/**
* Get surface tension force for the level set.
*
* @param wv Ideal MHD Riemann ghost fluid equations object.
* @return Surface tension force for the level set.
*/
double
gkyl_wv_mhd_rgfm_surface_tension(const struct gkyl_wv_eqn* wv);