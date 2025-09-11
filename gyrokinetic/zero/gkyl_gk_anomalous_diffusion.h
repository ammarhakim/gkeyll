#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_gk_anomalous_diffusion_auxfields { 
  const struct gkyl_array *nu; // Diffusivity.
  const struct gkyl_array *jacobgeo_inv; // Reciprocal of the conf-space Jacobian.
};

/**
 * Create a new gyrokinetic anomalous diffusion equation object.
 *
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param conf_range Conf-space range object.
 * @param is_zero_flux_bc Whether each boundary has zero-flux BCs.
 * @param is_absorb_bc Whether each boundary has absorb BCs.
 * @param skip_cell_threshold Threshold which to skip cells
 * @param use_gpu Whether to run on host or device.
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_gk_anomalous_diffusion_new(const struct gkyl_basis *basis, 
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  const bool *is_zero_flux_bc, const bool *is_absorb_bc, 
  double skip_cell_threshold, bool use_gpu);

/**
 * Set the auxiliary fields (e.g. diffusivity).
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gk_anomalous_diffusion_set_auxfields(const struct gkyl_dg_eqn* eqn,
  struct gkyl_gk_anomalous_diffusion_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gk_anomalous_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn *eqn,
  struct gkyl_gk_anomalous_diffusion_auxfields auxin);

#endif
