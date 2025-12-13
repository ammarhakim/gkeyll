#pragma once

// Private header, not for direct use in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_elasticity {
  struct gkyl_wv_eqn eqn; // Base equation object.

  double T_ref; // Reference temperature (unstressed configuration).
  double sound_speed; // Sound speed (speed of pressure waves).
  double shear_speed; // Shear speed (speed of shear waves).
  double heat_capacity; // Specific heat capacity at constant volume.

  double alpha_param; // Alpha parameter (nonlinear elasticity).
  double beta_param; // Beta parameter (nonlinear elasticity).
  double gamma_param; // Gamma parameter (nonlinear elasticity).
};

/**
* Compute primitive variables given the conserved variables.
* @param T_ref Reference temperature (unstressed configuration).
* @param sound_speed Sound speed (speed of pressure waves).
* @param shear_speed Shear speed (speed of shear waves).
* @param heat_capacity Specific heat capacity at constant volume.
* @param alpha_param Alpha parameter (nonlinear elasticity).
* @param beta_param Beta parameter (nonlinear elasticity).
* @param gamma_param Gamma parameter (nonlinear elasticity).
* @param q Conserved variable vector.
* @param v Primitive variable vector (output).
*/
GKYL_CU_D
void
gkyl_elasticity_prim_vars(double T_ref, double sound_speed, double shear_speed, double heat_capacity, double alpha_param, double beta_param, double gamma_param,
  const double q[14], double v[14]);

/**
* Compute inverse deformation gradient tensor given the conserved variables.
*
* @param q Conserved variable vector.
* @param inv_deformation_gradient Inverse deformation gradient tensor (output).
*/
GKYL_CU_D
void
gkyl_elasticity_inv_deformation_gradient(const double q[14], double ***inv_spatial_metric);

/**
* Free nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
*
* @param ref Reference counter for nonlinear elasticity equations object with Godunov-Romenskii hyperelastic equation of state.
*/
void
gkyl_wv_elasticity_free(const struct gkyl_ref_count* ref);