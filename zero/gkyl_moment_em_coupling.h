#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_moment_em_coupling_data {
  enum gkyl_eqn_type type;
  double charge; // species charge
  double mass; // species mass
  double k0; // closure parameter (default is 0.0, used by 10 moment)
};

struct gkyl_moment_em_coupling_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  int nfluids; // number of fluids
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // species data
  double epsilon0;

  bool has_collision; // has collisions
  // scaling factors for collision frequencies so that nu_sr=nu_base_sr/rho_s
  // nu_rs=nu_base_rs/rho_r, and nu_base_sr=nu_base_rs
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];
  double gas_gamma;
};

// Object type
typedef struct gkyl_moment_em_coupling gkyl_moment_em_coupling;

/**
 * Create new updater to update electromagnetic sources in fluid
 * equations.  Uses implicit time-stepping (Time-centered
 * Crank-Nicholson).
 *
 * @param inp Input parameters to updater
 */
gkyl_moment_em_coupling* gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp);

/**
 * Compute implicit update of the electromagnetic source terms in the
 * multi-fluid system.  The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param mes Moment electromagnetic sources updater object
 * @param update_rng Range on which to solve implicit time-centered update.
 * @param fluid Array of fluid variables (array size: nfluids)
 * @param app_accel Array of applied acceleration terms to fluid equations (for external forces)
 * @param em EM variables
 * @param app_current Applied current array (for external current driving)
 * @param ext_em External EM variables (for EM fields coming from external sources (coils, capacitors, etc.))
 */

void gkyl_moment_em_coupling_advance(const gkyl_moment_em_coupling *mes, double dt,
  const struct gkyl_range *update_rng, 
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], const struct gkyl_array *app_accel[GKYL_MAX_SPECIES],
  struct gkyl_array *em, const struct gkyl_array *app_current, const struct gkyl_array *ext_em);

/**
 * Delete updater.
 *
 * @param mes Updater to delete.
 */
void gkyl_moment_em_coupling_release(gkyl_moment_em_coupling *mes);
