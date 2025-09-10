#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h> 
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_gyrokinetic_auxfields {
  const struct gkyl_array *alpha_surf; // Pointer to surface expansion of phase space flux alpha.
  const struct gkyl_array *sgn_alpha_surf; // Pointer to sign(alpha_surf) at quadrature points.
  const struct gkyl_array *const_sgn_alpha; // Pointer to integer array for if sign(alpha_surf) is single-signed.
  const struct gkyl_array *phi; // Pointer to electrostatic potential.
  const struct gkyl_array *apar; // Pointer to A_\parallel.
  const struct gkyl_array *apardot; // Pointer to d(A_parallel)/dt.
};

/**
 * Create a new Gyrokinetic equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param phase_range Phase space range for use in indexing surface expansion of alpha
 * @param charge Species charge
 * @param mass Species mass
 * @param skip_cell_threshold Threshold for skipping cells in the gyrokinetic equation
 * @param gkmodel_id Model ID for gyrokinetics (e.g., general geometry vs. no toroidal field, see gkyl_eqn_type.h)
 * @param gk_geom Geometry struct
 * @param vel_map Velocity space mapping object.
 * @param use_gpu Boolean to determine if gyrokinetic equation object is on device
 * @return Pointer to Gyrokinetic equation object
 */
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const double charge, const double mass, double skip_cell_threshold, enum gkyl_gkmodel_id gkmodel_id,
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map, bool use_gpu);

/**
 * Set the auxiliary fields (e.g. EM fields) needed in computing
 * gyrokinetic updates.
 *
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gyrokinetic_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin);


