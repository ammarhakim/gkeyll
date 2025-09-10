#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type.
typedef struct gkyl_loss_cone_mask_gyrokinetic gkyl_loss_cone_mask_gyrokinetic;

// Inputs packaged as a struct.
struct gkyl_loss_cone_mask_gyrokinetic_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range.
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations).
  const struct gkyl_range *vel_range; // Velocity space range.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  const struct gkyl_array *bmag; // Magnetic field magnitude.
  const double *bmag_max; // Maximum bmag (on GPU if use_gpu=true).
  double mass; // Species mass.
  double charge; // Species charge.
  int num_quad; // Number of quad points in each direction to use (default: poly_order+1).
  bool use_gpu; // Whether to run on GPU.
};

/**
 * Create new updater that populates an array with the masking function
 *   if (mu > mu_bound)
 *     f = 1
 *   else
 *     f = 0
 * where mu_bound = (0.5*m*vpar^2+q*(phi-phi_m))/(B*(B_max/B-1))
 * is the trapped-passing boundary in vpar-mu space.
 *
 * @param inp Input parameters defined in gkyl_loss_cone_mask_gyrokinetic_inp struct.
 * @return New updater pointer.
 */
struct gkyl_loss_cone_mask_gyrokinetic* 
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp);

/**
 * Compute projection of the loss cone masking function on the phase-space basis.
 *
 * @param up Project on basis updater to run.
 * @param phase_rng Phase-space range.
 * @param conf_rng Configuration-space range.
 * @param phi Electrostatic potential.
 * @param phi_m Electrostatic potential at the mirror throat (on GPU if use_gpu=true).
 * @param mask_out Output masking function.
 */
void gkyl_loss_cone_mask_gyrokinetic_advance(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *phi, const double *phi_m, struct gkyl_array *mask_out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic* up);
