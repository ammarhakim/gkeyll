#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_gk_bc_type.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_gk_numerical_diffusion_auxfields { 
  const struct gkyl_array *nu; // Diffusivity.
  const struct gkyl_array *jacobgeo_inv; // Reciprocal of the conf-space Jacobian.
};

/**
 * Create a new gyrokinetic diffusion equation object.
 *
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param conf_range Conf-space range object.
 * @param bc_lower Boundary condition at lower boundaries.
 * @param bc_upper Boundary condition at upper boundaries.
 * @param has_diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param use_gpu Whether to run on host or device.
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_gk_numerical_diffusion_new(const struct gkyl_basis *basis, 
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  const enum gkyl_gyrokinetic_bc_type *bc_lower, const enum gkyl_gyrokinetic_bc_type *bc_upper,
  const bool *has_diff_in_dir, int diff_order, bool use_gpu);

/**
 * Set the auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gk_numerical_diffusion_set_auxfields(const struct gkyl_dg_eqn* eqn,
  struct gkyl_gk_numerical_diffusion_auxfields auxin);

