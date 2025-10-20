#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <kann.h>

struct gkyl_ten_moment_nn_closure_inp
{
  const struct gkyl_rect_grid *grid; // Grid on which to solve equations.
  int poly_order; // Polynomial order of learned DG coefficients.
  double k0; // Damping coefficient.
  kann_t* ann; // Neural network architecture.
};

// Object type.
typedef struct gkyl_ten_moment_nn_closure gkyl_ten_moment_nn_closure;

/**
 * Create a new updater to update the pressure tensor in the ten moment equation system using a machine-learned magnetized closure trained on PKPM simulations,
 * i.e. where the heat flux tensor q_ijk and its derivatives are calculated from predictions for q_par and q_perp on the basis of rho, p_par and p_perp.
 * Returns RHS for accumulation in a forward-Euler method.
 *
 * @param inp Input parameters for updater.
 * @return Pointer to updater.
 */
gkyl_ten_moment_nn_closure*
gkyl_ten_moment_nn_closure_new(struct gkyl_ten_moment_nn_closure_inp inp);

/**
 * Compute the right-hand-side contribution to the ten moment equation system from a machine-learned magnetized closure trained on PKPM simulations.
 * The update_rng MUST be a sub-range of the range on which the array is defined. That is, it must be either the same range as the arary range, or one
 * created using the gkyl_sub_range_init method.
 *
 * @param nnclosure Neural network closure updater object.
 * @param heat_flux_rng Range on which to compute the heat flux tensor (at cell nodes).
 * @param update_rng Range on which to compute update.
 * @param fluid Input array of fluid variables.
 * @param em_tot Total electromagnetic field variables (plasma + external).
 * @param heat_flux Array for storing intermediate computation of heat flux tensor (at cell nodes).
 * @param rhs Right-hand-side output.
 */
void
gkyl_ten_moment_nn_closure_advance(const gkyl_ten_moment_nn_closure *nnclosure, const struct gkyl_range *heat_flux_rng, const struct gkyl_range *update_rng,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot, struct gkyl_array *heat_flux, struct gkyl_array *rhs);

/**
 * Delete updater.
 *
 * @param nnclosure Updater to delete.
 */
 void
 gkyl_ten_moment_nn_closure_release(gkyl_ten_moment_nn_closure *nnclosure);