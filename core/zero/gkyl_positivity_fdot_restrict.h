#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

typedef struct gkyl_positivity_fdot_restrict gkyl_positivity_fdot_restrict;

enum gkyl_positivity_fdot_restrict_type {
  GKYL_POSITIVITY_FDOT_RESTRICT_QUAD, // Restrict df/dt at quadrature points.
  GKYL_POSITIVITY_FDOT_RESTRICT_AVG, // Restrict df/dt at cell average (modal coefficient 0).
  GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD, // Diode mode: set df/dt=0 if f<0 and df/dt<0 at quadrature points.
  GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG, // Diode mode: set df/dt=0 if f<0 and df/dt<0 at cell average.
};

struct gkyl_positivity_fdot_restrict_inp {
  struct gkyl_basis *basis; // Basis functions for the discretization.
  enum gkyl_positivity_fdot_restrict_type type; // Mode of restriction (quad or avg).
  double safety_factor; // Safety factor (0 < safety_factor < 1) for target floor in restriction. 0 means no safety factor (i.e. target floor = 0).
  bool use_gpu; // Whether to use GPU kernels if available.
};

/**
 * Create a new df/dt restriction positivity updater.
 *
 * @param inp Input structure with basis, mode, and safety factor.
 * @return New updater struct.
 */
struct gkyl_positivity_fdot_restrict*
gkyl_positivity_fdot_restrict_new(struct gkyl_positivity_fdot_restrict_inp inp);

/**
 * Apply df/dt restriction to maintain positivity.
 *
 * Restricts the rate of change (dfdt) at quadrature points or cell average
 * to ensure the forward-Euler prediction of f remains non-negative.
 *
 * @param up Updater struct.
 * @param range Range of cells to process.
 * @param f Input distribution (current state).
 * @param dfdt Input/output rate of change. Modified in-place.
 * @param dt Timestep to use for prediction.
 */
void
gkyl_positivity_fdot_restrict_advance(gkyl_positivity_fdot_restrict *up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt);

/**
 * Release resources for the df/dt restriction updater.
 *
 * @param up Updater to release.
 */
void
gkyl_positivity_fdot_restrict_release(gkyl_positivity_fdot_restrict *up);
