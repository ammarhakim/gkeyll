#pragma once

#include <gkyl_array.h>
#include <gkyl_gyrokinetic.h>

/**
 * @file gkyl_gk_field_priv.h
 * @brief Private header for gyrokinetic field solver functions.
 *
 * This header contains forward declarations for internal field-related
 * functions used by the gyrokinetic application. These functions handle
 * FEM projections, field solves, boundary conditions, energy diagnostics,
 * FLR corrections, polarization potentials, and biased wall configurations.
 */

// ============================================================================
// FEM Initialization Functions
// ============================================================================

/**
 * Initialize the FEM object for the Boltzmann field solver.
 * Sets up the finite element machinery needed for Boltzmann electron response.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize.
 */
void
gk_field_fem_new_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for 1D (1x) field solves.
 * Used for simulations with one configuration space dimension.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize.
 */
void
gk_field_fem_new_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for 2D or 3D (2x/3x) field solves.
 * Used for simulations with two or three configuration space dimensions.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize.
 */
void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

// ============================================================================
// FEM Projection Functions
// ============================================================================

/**
 * Project a DG field onto the parallel FEM basis to enforce continuity
 * along the parallel (z) direction. This is essential for field solves
 * that require C0 continuity in the parallel direction.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing solver state.
 * @param arr_dg Input DG array to be projected.
 * @param arr_fem Output FEM array (projected result with parallel continuity).
 */
void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);

// ============================================================================
// Charge Density Accumulation Functions
// ============================================================================

/**
 * Accumulate charge density for a Boltzmann electron response.
 * Computes the charge density contribution from a species assuming
 * electrons follow a Boltzmann distribution.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object for storing charge density.
 * @param s Species object providing density moments.
 * @param bflux Boundary flux array (for boundary contributions).
 */
void
gk_field_accumulate_rho_c_boltzmann(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gk_species *s, struct gkyl_array **bflux);

/**
 * Accumulate charge density for an adiabatic species response.
 * Computes the charge density contribution assuming an adiabatic
 * (flux-surface averaged) response for the species.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object for storing charge density.
 * @param s Species object providing density moments.
 * @param bflux Boundary flux array (for boundary contributions).
 */
void
gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gk_species *s, struct gkyl_array **bflux);

/**
 * Accumulate charge density for a full Poisson solve.
 * Computes the charge density from all kinetic species for use
 * in solving the gyrokinetic Poisson equation.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object for storing charge density.
 * @param s Species object providing density moments.
 * @param bflux Boundary flux array (for boundary contributions).
 */
void
gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gk_species *s, struct gkyl_array **bflux);

// ============================================================================
// Field Solve (RHS) Functions
// ============================================================================

/**
 * Compute the right-hand side and solve for the electrostatic potential
 * using a Boltzmann electron response. This assumes electrons are in
 * thermal equilibrium with phi.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the solve state and result.
 */
void
gk_field_boltzmann_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the 1D Poisson equation using parallel FEM projection.
 * Used for 1x simulations where the field equation reduces to a 1D solve.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the solve state and result.
 */
void
gk_field_1x_poisson_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the electrostatic Poisson equation for IWL (Inner Wall Limited)
 * simulations using a deflated approach.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the solve state and result.
 */
void
gk_field_deflate_poisson_es_iwl_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation in 2D or 3D with perpendicular FEM projection.
 * Handles the gyrokinetic Poisson equation for multi-dimensional simulations.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the solve state and result.
 */
void
gk_field_2x3x_poisson_perp_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

// ============================================================================
// FEM Release Functions
// ============================================================================

/**
 * Release resources allocated for the Boltzmann field solver FEM objects.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose FEM resources are to be released.
 */
void
gk_field_fem_release_boltzmann(const gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release resources allocated for the 1x field solver FEM objects.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose FEM resources are to be released.
 */
void
gk_field_fem_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release resources allocated for the 2x/3x field solver FEM objects.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose FEM resources are to be released.
 */
void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f);

// ============================================================================
// Boundary Condition Functions
// ============================================================================

/**
 * Enforce parallel boundary conditions on the field (enabled variant).
 * Applies the appropriate boundary conditions along the parallel direction.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing BC configuration.
 * @param finout Array to apply boundary conditions to (modified in-place).
 */
void
gk_field_enforce_parallel_bc_enabled(const gkyl_gyrokinetic_app *app,
  struct gk_field *field, struct gkyl_array *finout);

/**
 * No-op function for parallel boundary conditions (disabled variant).
 * Used when parallel BCs are not needed or handled elsewhere.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param finout Array (unchanged).
 */
void
gk_field_enforce_parallel_bc_disabled(const gkyl_gyrokinetic_app *app,
  struct gk_field *field, struct gkyl_array *finout);

// ============================================================================
// Energy Diagnostics Functions
// ============================================================================

/**
 * Initialize field energy diagnostic objects.
 * Sets up arrays and reducers needed for computing field energy.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize diagnostics for.
 */
void
gk_field_energy_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release field energy diagnostic resources.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose energy diagnostics are to be released.
 */
void
gk_field_energy_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize time-rate diagnostic objects for the field.
 * Sets up machinery for computing time derivatives of field energy.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize time-rate diagnostics for.
 */
void
gk_field_time_rate_diags_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release time-rate diagnostic resources.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose time-rate diagnostics are to be released.
 */
void
gk_field_time_rate_diags_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Compute field energy and write diagnostics (enabled variant).
 * Calculates the electrostatic field energy and outputs to file.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the potential.
 * @param tm Current simulation time.
 */
void
gk_field_calc_energy_enabled(struct gkyl_gyrokinetic_app *app,
  const struct gk_field *field, double tm);

/**
 * No-op function for field energy calculation (disabled variant).
 * Used when energy diagnostics are turned off.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param tm Current simulation time.
 */
void
gk_field_calc_energy_disabled(struct gkyl_gyrokinetic_app *app,
  const struct gk_field *field, double tm);

/**
 * Compute time derivative of field energy (active variant).
 * Calculates dE/dt for the electrostatic field energy.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the potential.
 * @param dt Time step size.
 * @param energy_reduced Output array for the reduced energy value.
 */
void
gk_field_calc_energy_dt_active(gkyl_gyrokinetic_app *app,
  const struct gk_field *field, double dt, double *energy_reduced);

/**
 * No-op function for time derivative of field energy (inactive variant).
 * Used when time-rate diagnostics are disabled.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param dt Time step size.
 * @param energy_reduced Output array (unchanged).
 */
void
gk_field_calc_energy_dt_none(gkyl_gyrokinetic_app *app,
  const struct gk_field *field, double dt, double *energy_reduced);

// ============================================================================
// Finite Larmor Radius (FLR) Correction Functions
// ============================================================================

/**
 * Initialize FLR correction objects for the field solver.
 * Sets up machinery for gyroaveraging and FLR effects in the potential.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize FLR corrections for.
 */
void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release FLR correction resources.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose FLR resources are to be released.
 */
void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Apply FLR inversion to obtain the gyroaveraged potential.
 * Solves the FLR equation to compute the potential seen by particles.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing solver state.
 * @param phi Electrostatic potential array (input/output).
 */
void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *phi);

/**
 * No-op function for FLR inversion (when FLR corrections are disabled).
 * Used when long-wavelength approximation is valid and FLR effects are negligible.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param phi Electrostatic potential array (unchanged).
 */
void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *phi);

// ============================================================================
// Polarization Potential Functions
// ============================================================================

/**
 * Initialize polarization potential objects.
 * Sets up arrays and solvers for the polarization contribution to the potential.
 *
 * @param f Field object to initialize polarization potential for.
 * @param app Gyrokinetic application object.
 */
void
gk_field_polarization_potential_new(struct gk_field *f, struct gkyl_gyrokinetic_app *app);

/**
 * Release polarization potential resources.
 *
 * @param f Field object whose polarization resources are to be released.
 */
void
gk_field_polarization_potential_release(struct gk_field *f);

// ============================================================================
// Biased Wall Functions
// ============================================================================

/**
 * Initialize biased wall boundary condition objects.
 * Sets up machinery for applying time-dependent wall potentials.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize biased wall for.
 */
void
gk_field_biassed_wall_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release biased wall resources.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose biased wall resources are to be released.
 */
void
gk_field_biassed_wall_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Calculate and apply the wall potential at the current simulation time.
 * Updates phi at the wall boundaries based on the biased wall configuration.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the potential.
 * @param tm Current simulation time.
 */
void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm);
