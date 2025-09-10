#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_gyrokinetic_vars gkyl_dg_calc_gyrokinetic_vars;

/**
 * Create new updater to compute gyrokinetic variables needed in 
 * updates and used for diagnostics. Methods compute:
 * alpha_surf : Surface expansion of phase space flux alpha
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param charge Species charge
 * @param mass Species mass
 * @param gkmodel_id Model ID for gyrokinetics (e.g., general geometry vs. no toroidal field, see gkyl_eqn_type.h)
 * @param gk_geom Geometry struct
 * @param vel_map Velocity space mapping object.
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_gyrokinetic_vars* 
gkyl_dg_calc_gyrokinetic_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass, enum gkyl_gkmodel_id gkmodel_id, 
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map,
  bool use_gpu);

/**
 * Compute surface expansion of phase space flux alpha
 * Computes the Poisson bracket of alpha = dz/dt = {z, H} and then evaluates the resulting polynomial
 * expansion at a surface and projects the evaluated quantity onto the surface basis.
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing gyrokinetic variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param phase_range Phase space range 
 * @param phase_ext_range Extended Phase space range (so we obtain alpha_surf at all the needed surfaces)
 * @param phi Electrostatic potential
 * @param alpha_surf Output surface expansion in a cell on the *lower* edge in each direction 
 * @param sgn_alpha_surf Output sign(alpha) at quadrature points along a surface 
 * @param const_sgn_alpha Output boolean array for if sign(alpha) is a constant on the surface
 *                        If sign(alpha) is a constant, kernels are simpler and we exploit this fact.
 */
void gkyl_dg_calc_gyrokinetic_vars_alpha_surf(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_ext_range, const struct gkyl_array *phi, 
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha);

/**
 * Delete pointer to updater to compute gyrokinetic variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_gyrokinetic_vars_release(struct gkyl_dg_calc_gyrokinetic_vars *up);
