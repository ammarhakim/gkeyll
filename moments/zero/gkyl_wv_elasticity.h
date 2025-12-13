#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_elasticity_rp {
  WV_ELASTICITY_RP_LAX = 0, // Default (Lax fluxes).
};

// Input context, packaged as a struct.
struct gkyl_wv_elasticity_inp {
  double T_ref; // Reference temperature (unstressed configuration).
  double sound_speed; // Sound speed (speed of pressure waves).
  double shear_speed; // Shear speed (speed of shear waves).
  double heat_capacity; // Specific heat capacity at constant volume.

  double alpha_param; // Alpha parameter (nonlinear elasticity).
  double beta_param; // Beta parameter (nonlinear elasticity).
  double gamma_param; // Gamma parameter (nonlinear elasticity).

  enum gkyl_wv_elasticity_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
*
* @param T_ref Reference temperature (unstressed configuration).
* @param sound_speed Sound speed (speed of pressure waves).
* @param shear_speed Shear speed (speed of shear waves).
* @param heat_capacity Specific heat capacity at constant volume.
* @param alpha_param Alpha parameter (nonlinear elasticity).
* @param beta_param Beta parameter (nonlinear elasticity).
* @param gamma_param Gamma parameter (nonlinear elasticity).
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_elasticity_new(double T_ref, double sound_speed, double shear_speed, double heat_capacity, double alpha_param, double beta_param, double gamma_param,
  bool use_gpu);

/**
* Create a new nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_elasticity_inew(const struct gkyl_wv_elasticity_inp* inp);

/**
* Get reference temperature (unstressed configuration).
*
* @param eqn Nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Reference temperature (unstressed configuration).
*/
double
gkyl_wv_elasticity_T_ref(const struct gkyl_wv_eqn* eqn);

/**
* Get sound speed (speed of pressure waves).
*
* @param eqn Nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Sound speed (speed of pressure waves).
*/
double
gkyl_wv_elasticity_sound_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get shear speed (speed of shear waves).
*
* @param eqn Nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Shear speed (speed of shear waves).
*/
double
gkyl_wv_elasticity_shear_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get specific heat capacity at constant volume.
*
* @param eqn Nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Specific heat capacity at constant volume.
*/
double
gkyl_wv_elasticity_heat_capacity(const struct gkyl_wv_eqn* eqn);

/**
* Get alpha parameter (nonlinear elasticity).
*
* @param eqn Nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Alpha parameter (nonlinear elasticity).
*/
double
gkyl_wv_elasticity_alpha_param(const struct gkyl_wv_eqn* eqn);

/**
* Get beta parameter (nonlinear elasticity).
*
* @param eqn Nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Beta parameter (nonlinear elasticity).
*/
double
gkyl_wv_elasticity_beta_param(const struct gkyl_wv_eqn* eqn);

/**
* Get gamma parameter (nonlinear elasticity).
*
* @param eqn Nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Gamma parameter (nonlinear elasticity).
*/
double
gkyl_wv_elasticity_gamma_param(const struct gkyl_wv_eqn* eqn);