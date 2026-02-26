#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_bc_sheath_gyrokinetic gkyl_bc_sheath_gyrokinetic;

/**
 * Create a new updater to apply conducting sheath BCs in gyrokinetics.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param vel_map Velocity space mapping object.
 * @param cdim Configuration space dimensions.
 * @param q2Dm charge-to-mass ratio times 2.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_sheath_gyrokinetic* gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_basis *basis, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_velocity_map *vel_map, int cdim, double q2Dm, bool use_gpu);

/**
 * Apply the sheath BC with the bc_sheath_gyrokinetic object.
 *
 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param alpha_mu Alpha parameter for sheath BCs.
 * @param distf Distribution function array to apply BC to.
 * @param conf_r Configuration space range (to index phi).
 */
void gkyl_bc_sheath_gyrokinetic_advance(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r);

/**
 * Set the alpha_mu array used in the sheath BC. 
 * This is used to implement a mu-dependent vcut in the sheath BC, 
 * where alpha_mu is the mu dependent multiplying factor in the expression for vcut.
 * 
 * @param up BC updater.
 * @param alpha_mu The alpha_mu array to use in the sheath BC.
 */
void gkyl_bc_sheath_gyrokinetic_set_alpha_mu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *alpha_mu);

/**
 * Get the alpha_mu array used in the sheath BC.
 * 
 * @param up BC updater.
 * @param alpha_mu On output, the alpha_mu array used in the sheath BC.
 */
void gkyl_bc_sheath_gyrokinetic_get_alpha_mu(const struct gkyl_bc_sheath_gyrokinetic *up, struct gkyl_array *alpha_mu);

/**
 * Get the basis of the alpha_mu array used in the sheath BC.
 * 
 * @param up BC updater.
 * @param alpha_mu_basis On output, the basis of the alpha_mu array used in the sheath BC.
 */
void gkyl_bc_sheath_gyrokinetic_get_alpha_mu_basis(const struct gkyl_bc_sheath_gyrokinetic *up, struct gkyl_basis *alpha_mu_basis);

/**
 * Free memory associated with bc_sheath_gyrokinetic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up);
