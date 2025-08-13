#pragma once

// Private header, not for direct use in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_vacuum_einstein {
  struct gkyl_wv_eqn eqn; // Base equation object.

  enum gkyl_spacetime_slicing spacetime_slicing; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution; // Spacetime evolution system.
};

/**
* Compute flux vector. Assumes rotation to local coordinate system.
*
* @param spacetime_slicing Spacetime slicing condition.
* @param spacetime_evolution Spacetime evolution system.
* @param q Conserved variable vector.
* @param flux Flux vector in direction 'dir' (output).
*/
GKYL_CU_D
void
gkyl_vacuum_einstein_flux(enum gkyl_spacetime_slicing spacetime_slicing, enum gkyl_spacetime_evolution spacetime_evolution, const double q[64], double flux[64]);

/**
* Compute inverse spatial metric tensor (in covariant component form) given the conserved variables.
*
* @param q Conserved variable vector.
* @param inv_spatial_metric Inverse spatial metric tensor (output).
*/
GKYL_CU_D
void
gkyl_vacuum_einstein_inv_spatial_metric(const double q[64], double ***inv_spatial_metric);

/**
* Compute maximum absolute wave speed.
*
* @param spacetime_slicing Spacetime slicing condition.
* @param q Conserved variable vector.
* @return Maximum absolute wave speed for a given q.
*/
GKYL_CU_D
static inline double
gkyl_vacuum_einstein_abs_speed(enum gkyl_spacetime_slicing spacetime_slicing, const double q[64]);

/**
* Compute Riemann variables given the conserved variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param qin Conserved variable vector (input).
* @param wout Riemann variable vector (output).
*/
GKYL_CU_D
static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout);

/**
* Compute conserved variables given the Riemann variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param win Riemann variable vector (input).
* @param qout Conserved variable vector (output).
*/
GKYL_CU_D
static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double *qout);

/**
* Free vacuum Einstein equations object using the Bona-Masso formalism.
*
* @param ref Reference counter for vacuum Einstein equations object using the Bona-Masso formalism.
*/
void gkyl_vacuum_einstein_free(const struct gkyl_ref_count* ref);