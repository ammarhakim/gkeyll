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
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range.
  const struct gkyl_velocity_map *vel_map;    // Velocity space mapping object.
  const struct gkyl_array *bmag; // Magnetic field magnitude (cdim DG expansion).
  bool use_gpu; // Flag to indicate if GPU should be used.
  double mass; // Species mass.
  double charge; // Species charge.
};

/**
 * Create new updater that populates an array with the masking function
 * based on the escape-barrier criterion
 *   EB(z,mu) = min(max_{s in [z_L,z]} U(s,mu), max_{s in [z,z_R]} U(s,mu))
 * with U = mu*B + q*phi.
 * A node is trapped if H < EB.
 *
 * @param inp Input parameters defined in gkyl_loss_cone_mask_gyrokinetic_inp struct.
 * @return New updater pointer.
 */
struct gkyl_loss_cone_mask_gyrokinetic*
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp);

/**
 * Compute the loss-cone mask on phase-space cell nodes.
 *
 * @param up Project on basis updater to run.
 * @param phase_rng Phase-space range.
 * @param conf_rng Configuration-space range.
 * @param bmag Magnetic field magnitude (cdim DG expansion).
 * @param phi Electrostatic potential.
 * @param mask_out Output masking function.
 */
void gkyl_loss_cone_mask_gyrokinetic_advance(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *bmag, const struct gkyl_array *phi, struct gkyl_array *mask_out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic *up);
