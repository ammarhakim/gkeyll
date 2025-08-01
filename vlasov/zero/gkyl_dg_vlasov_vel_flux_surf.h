#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_vlasov_vel_flux_surf gkyl_dg_vlasov_vel_flux_surf;

// Input packaged as a struct
struct gkyl_dg_vlasov_vel_flux_surf_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid. 
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions. 
  const struct gkyl_basis *phase_basis; // Phase-space basis functions. 
  const struct gkyl_range *hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_range *vel_range; // Velocity-space range for use in indexing velocity-space Jacobian. 
  enum gkyl_model_id model_id; // enum to determine what type of Vlasov model (e.g., non-relativistic vs. relativistic).
  bool has_qmem; // bool to determine whether we have electric or magnetic fields (used for external forces too).
  bool has_phi; // bool to determine whether we have potentials (either electrostatic or gravitational).
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new updater to compute the velocity-space fluxes in a modal DG scheme.
 * Updater computes the fluxes at the necessary Gauss-Legendre quadrature points to
 * eliminate aliasing errors by evaluating a total phase space characteristic, e.g., 
 * alpha_v = -grad_x(H) + q/m (E + grad_v(H) x B) - grad_x(q/m*phi + m*phi_g) + q/m*grad_v(H) x curl(A) 
 * at each quadrature point, and then computing F = alpha_v/2*(f^+ + f^-) - |alpha_v|/2*(f^+ - f^-)
 * at quadrature points. Here, f^+/- are the distribution function evaluated just outside/inside 
 * the cell interface and we utilize the full Lax-Friedrichs form to avoid warp divergence on GPUs. 
 * This nodal array is then converted back to a modal expansion on the *lower* surface of each cell
 * for use in the hyper_dg updater for integrating the surface terms in the DG discretization.    
 *
 * @param inp Input parameters defined in gkyl_dg_vlasov_vel_flux_surf_inp struct.
 * @return Pointer to velocity-space surface flux updater. 
 */
struct gkyl_dg_vlasov_vel_flux_surf* 
gkyl_dg_vlasov_vel_flux_surf_inew(const struct gkyl_dg_vlasov_vel_flux_surf_inp *inp);

/**
 * Create new updater to compute the velocity-space fluxes in a modal DG scheme on 
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_vlasov_vel_flux_surf* 
gkyl_dg_vlasov_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_vel_flux_surf_inp *inp);

/**
 * Compute the velocity-space fluxes in a modal DG scheme. 
 *
 * @param up Velocity-space surface flux updater. 
 * @param conf_range Configuration-space range for indexing electromagnetic fields. 
 * @param phase_range Phase-space range for indexing velocity-space flux array. 
 * @param jacob_vel Input velocity-space Jacobian. 
 * @param hamil Input Hamiltonian for computing Hamiltonian forces. 
 * @param qmem Input q/m*(E,B) for computing Lorentz force. 
 * @param pot_tot Input total potentials for computing forces from, e.g., scalar potentials. 
 * @param fin Input distribution function. 
 * @param cflrate Input cflrate. Accumulated to from maximum alpha_v at quadrature points.  
 * @param vel_flux_surf Output modal velocity-space fluxes. 
 */
void 
gkyl_dg_vlasov_vel_flux_surf_advance(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel, const struct gkyl_array *hamil, 
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf);

 /**
 * Host-side wrapper for computing velocity-space fluxes on device. 
 */  
void 
gkyl_dg_vlasov_vel_flux_surf_advance_cu(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel, const struct gkyl_array *hamil, 
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf);