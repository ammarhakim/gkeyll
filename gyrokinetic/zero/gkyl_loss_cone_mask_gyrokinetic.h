#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_velocity_map.h>

// Object type.
typedef struct gkyl_loss_cone_mask_gyrokinetic gkyl_loss_cone_mask_gyrokinetic;

/** Whether a guiding-center trajectory has a loss endpoint in one direction. */
enum gkyl_gk_loss_cone_trajectory_type {
  GKYL_GK_LOSS_CONE_OPEN_TRAJECTORY = 0, // The trajectory reaches an absorbing endpoint.
  GKYL_GK_LOSS_CONE_CLOSED_TRAJECTORY,   // The trajectory has no loss endpoint.
};

// Inputs packaged as a struct.
struct gkyl_loss_cone_mask_gyrokinetic_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_velocity_map *vel_map; // Velocity-space mapping object.
  bool use_gpu; // Flag to indicate if GPU should be used.
  double mass; // Species mass.
  double charge; // Species charge.
  enum gkyl_gk_loss_cone_trajectory_type lower_trajectory; // Lower-z trajectory type.
  enum gkyl_gk_loss_cone_trajectory_type upper_trajectory; // Upper-z trajectory type.
};

/**
 * Create new updater that populates an array with the masking function
 * based on the escape-barrier criterion. Let
 *
 *   H = (1/2) m v_parallel^2 + U(z,mu),  U(z,mu) = mu B(z) + q phi(z).
 *
 * For an open trajectory, the directional escape barrier is the largest U
 * between the particle and its absorbing endpoint. An absorbing plasma edge
 * has no additional endpoint value. If a material-wall potential is supplied,
 * the trajectory extends through the sheath and also includes
 *
 *   U_wall = mu B(z_edge) + q phi_wall.
 *
 * A closed trajectory has an infinite escape barrier. A phase-space corner is
 * trapped only when H is below both directional barriers; a cell is masked as
 * trapped only when all of its corners are trapped.
 *
 * @param inp Input parameters defined in gkyl_loss_cone_mask_gyrokinetic_inp struct.
 * @return New updater pointer.
 */
struct gkyl_loss_cone_mask_gyrokinetic*
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp);

/**
 * Compute the loss-cone mask from phase-space cell corners.
 *
 * The caller supplies the magnetic field magnitude and electrostatic
 * potential arrays. This keeps the updater free of any communication logic;
 * the app is responsible for assembling global data when needed.
 *
 * @param up Loss-cone mask updater to run.
 * @param phase_rng Phase-space range.
 * @param conf_rng Configuration-space range.
 * @param bmag Magnetic field magnitude (cdim DG expansion).
 * @param phi Electrostatic potential throughout the configuration-space range.
 * @param phi_wall_lo Optional lower material-wall potential; NULL for an absorbing plasma edge.
 * @param phi_wall_up Optional upper material-wall potential; NULL for an absorbing plasma edge.
 * @param mask_out Cellwise mask: 1 if every cell corner is trapped, 0 if any
 *                 corner can escape through either boundary.
 */
void gkyl_loss_cone_mask_gyrokinetic_advance(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *bmag, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall_lo, const struct gkyl_array *phi_wall_up,
  struct gkyl_array *mask_out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic *up);
