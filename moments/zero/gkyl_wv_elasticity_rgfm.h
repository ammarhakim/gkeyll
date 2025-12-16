#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_elasticity_rgfm_rp {
  WV_ELASTICITY_RGFM_RP_LAX = 0, // Default (Lax fluxes).
};

// Input context, packaged as a struct.
struct gkyl_wv_elasticity_rgfm_inp {
  int num_species; // Number of distinct species in the domain.

  double* rho_ref_s; // Reference densities (unstressed configuration) for each species in the domain.
  double* T_ref_s; // Reference temperatures (unstressed configuration) for each species in the domain.
  double* sound_speed_s; // Sound speeds (speed of pressure waves) for each species in the domain.
  double* shear_speed_s; // Shear speeds (speed of shear waves) for each species in the domain.
  double* heat_capacity_s; // Specific heat capacities at constant volume for each species in the domain.

  double* alpha_param_s; // Alpha parameters (nonlinear elasticity) for each species in the domain.
  double* beta_param_s; // Beta parameters (nonlinear elasticity) for each species in the domain.
  double* gamma_param_s; // Gamma parameters (nonlinear elasticity) for each species in the domain.

  int reinit_freq; // Reinitialization frequency for the level set.
  double surface_tension; // Surface tension force for the level set.

  enum gkyl_wv_elasticity_rgfm_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
*
* @param num_species Number of distinct species in the domain.
* @param rho_ref_s Reference densities (unstressed configuration) for each species in the domain.
* @param T_ref_s Reference temperatures (unstressed configuration) for each species in the domain.
* @param sound_speed_s Sound speeds (speed of pressure waves) for each species in the domain.
* @param shear_speed_s Shear speeds (speed of shear waves) for each species in the domain.
* @param heat_capacity_s Specific heat capacities at constant volume for each species in the domain.
* @param alpha_param_s Alpha parameters (nonlinear elasticity) for each species in the domain.
* @param beta_param_s Beta parameters (nonlinear elasticity) for each species in the domain.
* @param gamma_param_s Gamma parameters (nonlinear elasticity) for each species in the domain.
* @param reinit_freq Reinitialization frequency for the level set.
* @param surface_tension Surface tension force for the level set.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_elasticity_rgfm_new(int num_species, double* rho_ref_s, double* T_ref_s, double* sound_speed_s, double* shear_speed_s, double* heat_capacity_s,
  double* alpha_param_s, double* beta_param_s, double* gamma_param_s, int reinit_freq, double surface_tension, bool use_gpu);

/**
* Create a new nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_elasticity_rgfm_inew(const struct gkyl_wv_elasticity_rgfm_inp* inp);

/**
* Get number of distinct species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Number of distinct species in the domain.
*/
int
gkyl_wv_elasticity_rgfm_num_species(const struct gkyl_wv_eqn* eqn);

/**
* Get reference densities (unstressed configuration) for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Reference densities (unstressed configuration) for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_rho_ref_s(const struct gkyl_wv_eqn* eqn);

/**
* Get reference temperatures (unstressed configuration) for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Reference temperatures (unstressed configuration) for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_T_ref_s(const struct gkyl_wv_eqn* eqn);

/**
* Get sound speeds (speed of pressure waves) for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Sound speeds (speed of pressure waves) for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_sound_speed_s(const struct gkyl_wv_eqn* eqn);

/**
* Get shear speeds (speed of shear waves) for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Shear speeds (speed of shear waves) for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_shear_speed_s(const struct gkyl_wv_eqn* eqn);

/**
* Get specific heat capacities at constant volume for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Specific heat capacities at constant volume for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_heat_capacity_s(const struct gkyl_wv_eqn* eqn);

/**
* Get alpha parameters (nonlinear elasticity) for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Alpha parameters (nonlinear elasticity) for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_alpha_param_s(const struct gkyl_wv_eqn* eqn);

/**
* Get beta parameters (nonlinear elasticity) for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Beta parameters (nonlinear elasticity) for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_beta_param_s(const struct gkyl_wv_eqn* eqn);

/**
* Get gamma parameters (nonlinear elasticity) for each species in the domain.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Gamma parameters (nonlinear elasticity) for each species in the domain.
*/
double*
gkyl_wv_elasticity_rgfm_gamma_param_s(const struct gkyl_wv_eqn* eqn);

/**
* Get reinitialization frequency for the level set.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Reinitialization frequency for the level set.
*/
int
gkyl_wv_elasticity_rgfm_reinit_freq(const struct gkyl_wv_eqn* eqn);

/**
* Get surface tension force for the level set.
*
* @param eqn Nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
* @return Surface tension force for the level set.
*/
double
gkyl_wv_elasticity_rgfm_surface_tension(const struct gkyl_wv_eqn* eqn);