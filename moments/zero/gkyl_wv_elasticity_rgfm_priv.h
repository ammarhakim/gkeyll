#pragma once

// Private header, not for direct use in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_elasticity_rgfm {
  struct gkyl_wv_eqn eqn; // Base equation object.
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
};

/**
* Compute primitive variables given the conserved variables.
*
* @param num_species Number of distinct species in the domain.
* @param T_ref_s Reference temperatures (unstressed configuration) for each species in the domain.
* @param sound_speed_s Sound speeds (speed of pressure waves) for each species in the domain.
* @param shear_speed_s Shear speeds (speed of shear waves) for each species in the domain.
* @param heat_capacity_s Specific heat capacities at constant volume for each species in the domain.
* @param alpha_param_s Alpha parameters (nonlinear elasticity) for each species in the domain.
* @param beta_param_s Beta parameters (nonlinear elasticity) for each species in the domain.
* @param gamma_param_s Gamma parameters (nonlinear elasticity) for each species in the domain.
* @param q Conserved variable vector.
* @param v Primitive variable vector (output).
*/
GKYL_CU_D
void
gkyl_elasticity_rgfm_prim_vars(int num_species, double* T_ref_s, double* sound_speed_s, double* shear_speed_s, double* heat_capacity_s,
  double* alpha_param_s, double* beta_param_s, double* gamma_param_s, const double* q, double* v);

/**
* Compute total inverse deformation gradient tensor given the conserved variables.
*
* @param q Conserved variable vector.
* @param inv_deformation_gradient_total Total inverse deformation gradient tensor (output).
*/
GKYL_CU_D
void
gkyl_elasticity_rgfm_inv_deformation_gradient_total(const double* q, double ***inv_deformation_gradient_total);

/**
* Free nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
*
* @param ref Reference counter for nonlinear elasticity Riemann ghost fluid equations object with Godunov-Romenskii hyperelastic equation of state.
*/
void
gkyl_wv_elasticity_rgfm_free(const struct gkyl_ref_count* ref);