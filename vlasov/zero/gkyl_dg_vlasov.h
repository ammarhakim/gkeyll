#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Input packaged as a struct
struct gkyl_dg_vlasov_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions. 
  const struct gkyl_basis *phase_basis; // Phase-space basis functions. 
  const struct gkyl_range *hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_range *conf_range; // Configuration-space range for use in indexing forces (EM fields, potentials, etc.). 
  const struct gkyl_range *vel_range; // Velocity-space range for use in indexing velocity-space Jacobian. 
  const struct gkyl_range *phase_range; // Range for indexing velocity-space flux. 
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped. 
  enum gkyl_model_id model_id; // enum to determine what type of Vlasov model (e.g., non-relativistic vs. relativistic).
  bool use_vmap; // bool to determine whether we have a nonuniform velocity map. 
  bool has_E; // bool to determine whether we have electric fields (used for external forces too).
  bool has_phi; // bool to determine whether we have potentials (either electrostatic or gravitational).
  bool has_B; // bool to determine whether we have magnetic fields.
  bool has_rad; // bool to determine whether we have radiation drag forces. 
  const struct gkyl_array *jacob_vel; // Velocity-space Jacobian.  
  const struct gkyl_array *hamil; // Hamiltonian utilized to compute advection in configuration and velocity space. 
  const struct gkyl_array *qmem; // q/m*(E,B) electromagnetic fields (including external electromagnetic fields and forces).
  const struct gkyl_array *pot_tot; // (q/m*(phi + phi_ext) + m*phi_g, q/m*A_ext) total potentials. 
  const struct gkyl_array *rad; // Radiation drag forces, F_rad(v). 
  const struct gkyl_array *vel_flux_surf; // Modal expansion of fluxes at velocity space surfaces. 
  const struct gkyl_array *f_no_J; // Distribution function without velocity-space Jacobian. Used by magnetic field updates.  
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new Vlasov equation object. 
 *
 * @param inp Input parameters defined in gkyl_dg_vlasov_inp struct.
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_inew(const struct gkyl_dg_vlasov_inp *inp);

/**
 * Create a new Vlasov equation object that lives on NV-GPU
 *
 * @param inp Input parameters defined in gkyl_dg_vlasov_inp struct.
 * @return Pointer to Vlasov equation object on device. 
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_cu_dev_inew(const struct gkyl_dg_vlasov_inp *inp);
