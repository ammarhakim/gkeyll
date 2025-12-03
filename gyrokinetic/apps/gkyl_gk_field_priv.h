#pragma once

#include <gkyl_array.h>
#include <gkyl_gyrokinetic.h>

// Forward declarations for private field-related functions.

// FEM projection and initialization functions.

/**
 * Project a DG field onto the parallel FEM basis to make it
 * continuous along z (or to solve a Poisson equation in 1x).
 *
 * @param app Application object
 * @param field Field object
 * @param arr_dg Input DG array (to be projected)
 * @param arr_fem Output FEM array (projected result)
 */
void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);

/**
 * Initialize the FEM object for the Boltzmann field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_new_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for the 1x field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_new_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for the 2x or 3x field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Accumulate charge density for a Boltzmann response.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_boltzmann(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux);

/**
 * Accumulate charge density for an adiabatic response.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux);

/**
 * Accumulate charge density for a Poisson solve.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux);

/**
 * Solve for the field using a Boltzmann response.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_boltzmann_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation in 1x with the parallel FEM projection.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_1x_poisson_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation for IWL sims.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_deflate_poisson_es_iwl_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation in 2x and 3x with the parallel FEM projection.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_2x3x_poisson_perp_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Release resources specific to Boltzmann field solver.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_release_boltzmann(const gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release resources specific to 1x field solver.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release resources specific to 2x3x field solver.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f);

// Boundary condition functions.
void
gk_field_enforce_parallel_bc_enabled(const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *finout);

void
gk_field_enforce_parallel_bc_disabled(const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *finout);

// Energy diagnostics functions.
void
gk_field_calc_energy_dt_active(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced);

void
gk_field_calc_energy_dt_none(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced);

void
gk_field_time_rate_diags_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_time_rate_diags_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_calc_energy_enabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm);

void
gk_field_calc_energy_disabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm);

// FLR (Finite Larmor Radius) functions.
void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi);

void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi);

void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

// Polarization potential functions.
void
gk_field_polarization_potential_new(struct gk_field *f, struct gkyl_gyrokinetic_app *app);

void
gk_field_polarization_potential_release(struct gk_field *f);

// Biassed wall functions.
void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm);

void
gk_field_biassed_wall_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_biassed_wall_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

// Energy diagnostics initialization functions.
void
gk_field_energy_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_energy_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

static void
gk_field_calc_ambi_pot_sheath_vals(gkyl_gyrokinetic_app *app, struct gk_field *field);

static void
gk_field_2x3x_add_TSBC_and_SSFG_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f);
