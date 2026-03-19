#pragma once

#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>


// Input packaged as a struct
struct gkyl_dg_maxwell_inp {
  const struct gkyl_basis *cbasis; // Configuration-space basis functions. 
  const struct gkyl_range *crange; // Configuration-space range for use in indexing forces (EM fields, potentials, etc.). 
  enum gkyl_field_id field_id; // enum to determine what type of maxwell model (e.g., non-relativistic vs. relativistic).
  double lightSpeed; // Speed of light
  double elcErrorSpeedFactor; // Factor multiplying lightSpeed for div E correction
  double mgnErrorSpeedFactor; // Factor multiplying lightSpeed for div B correction
  const struct gkyl_array *conf_flux_surf; // Modal expansion of fluxes at configuration space surfaces. 
  const struct gkyl_array *lapse_vol_nodes; // nodal expansion of lapse at volume nodes
  const struct gkyl_array *shift_vol_nodes; // nodal expansion of shift at volume nodes
  const struct gkyl_array *h_ij_vol_nodes; // nodal expansion of spatial metric at volume nodes
  const struct gkyl_array *det_h_vol_nodes; // nodal expansion of spatial metric determinant at volume nodes
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new Maxwell equation object.
 *
 * @param inp Input parameters defined in gkyl_dg_maxwell_inp struct.
 * @return Pointer to maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_inew(const struct gkyl_dg_maxwell_inp *inp);

/*
 * Create a new Maxwell equation object that lives on NV-GPU.
 *
 * @param inp Input parameters defined in gkyl_dg_maxwell_inp struct.
 * @return Pointer to maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_cu_dev_inew(const struct gkyl_dg_maxwell_inp *inp);
