#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fv_proj.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_wv_gr_euler_prim_priv.h>  // struct gkyl_gr_euler_eos

// ---- Source-family parameter bundles -----------------------------------------
// Each operator-split source family is one struct, switched by its 'enabled'
// flag. The same bundles appear in the app input (struct gkyl_moment_species /
// gkyl_moment_field in gkyl_moment.h), the app's internal species/field
// objects, and the coupling input below, so parameters travel by struct
// assignment instead of field-by-field copying.

// Frictional inter-species sources.
struct gkyl_moment_friction {
  bool enabled;
  bool use_explicit; // Integrate frictional sources with an explicit (SSP-RK3) solver.
  double Z; // Ionization number.
  double T_elc; // Electron temperature.
  double Lambda_ee; // Electron-electron collisional term.
};

// Volume-based geometrical sources (expanding/contracting box).
struct gkyl_moment_volume_sources {
  bool enabled;
  double gas_gamma; // Adiabatic index.
  double U0; // Initial comoving plasma velocity.
  double R0; // Initial radial distance from expansion/contraction center.
};

// Reactive sources.
struct gkyl_moment_reactivity {
  bool enabled;
  double gas_gamma; // Adiabatic index.
  double specific_heat_capacity;
  double energy_of_formation;
  double ignition_temperature;
  double reaction_rate;
};

// Coupled fluid-Einstein sources in plane-symmetric spacetimes.
struct gkyl_moment_einstein_medium {
  bool enabled;
  double gas_gamma; // Adiabatic index.
  double kappa; // Stress-energy prefactor.
};

// General relativistic sources for equation systems characterized by a
// single adiabatic index: GR Euler with ultra-relativistic or ideal-gas
// equation of state, and GR MHD.
struct gkyl_moment_gr_sources {
  bool enabled;
  double gas_gamma; // Adiabatic index.
};

// General relativistic two-fluid sources.
struct gkyl_moment_gr_twofluid {
  bool enabled;
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass.
  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.
  double gas_gamma_elc; // Adiabatic index for electrons.
  double gas_gamma_ion; // Adiabatic index for ions.
  double e_fact; // Electric field divergence error propagation speed.
};

// Vacuum Einstein sources in the Bona-Masso formalism (used for both the
// plain and the conformal variants).
struct gkyl_moment_vacuum_einstein {
  bool enabled;
  double excision_threshold; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing slicing; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution evolution; // Spacetime evolution system.
};

struct gkyl_moment_em_coupling_data {
  enum gkyl_eqn_type type; // Equation type.

  double charge; // Species charge.
  double mass; // Species mass.

  double k0; // Closure parameter (for 10-moment equations only; defaults to 0.0).
};

struct gkyl_moment_em_coupling_inp {
  const struct gkyl_rect_grid *grid; // Grid over which the equations are solved.
  int nfluids; // Number of fluid species.
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // Data for each fluid species.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.

  bool static_field; // Is the plasma field static? If true, only J is updated to new time step. 
  double t_ramp_E; // Ramp-up time for the linear ramp function for initializing external electric fields.
  double t_ramp_curr; // Ramp-up time for the linear ramp function for initializing applied currents.

  bool has_collision; // Run with collisions switched on.
  bool use_rel; // Assume special relativistic fluid species.
  
  // Matrix of scaling factors for collision frequencies. Should be symmetric (i.e. nu_base_sr = nu_base_rs).
  // These are defined such that nu_sr = nu_base_sr / rho_s, and nu_rs = nu_base_rs / rho_r.
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];

  bool use_explicit_em_coupling; // Use the explicit source-solver for handling moment-EM coupling (not operational yet).

  bool has_nT_sources; // Run with number density and temperature sources.

  // Source-family parameter bundles (see the definitions at the top of this
  // header). Each family runs when its 'enabled' flag is set.
  struct gkyl_moment_friction friction; // Frictional sources.
  struct gkyl_moment_volume_sources volume_sources; // Volume-based geometrical sources.
  struct gkyl_moment_reactivity reactivity; // Reactive sources.
  struct gkyl_moment_einstein_medium einstein_medium; // Coupled fluid-Einstein sources (plane-symmetric spacetimes).
  struct gkyl_moment_gr_sources gr_ultra_rel; // GR Euler sources, ultra-relativistic equation of state.
  struct gkyl_moment_gr_sources gr_euler; // GR Euler sources, general equation of state.
  struct gkyl_moment_gr_twofluid gr_twofluid; // GR two-fluid sources.
  struct gkyl_moment_vacuum_einstein vacuum_einstein; // Vacuum Einstein sources, Bona-Masso formalism.
  struct gkyl_moment_vacuum_einstein vacuum_einstein_conformal; // Vacuum Einstein sources, conformal Bona-Masso formalism.
  struct gkyl_moment_gr_sources gr_mhd; // GR MHD sources.

  bool has_gr_em_coupling; // Run with explicit (special-)relativistic multi-fluid + Maxwell coupling (SI units, separate fluid/EM inputs).
  struct gkyl_gr_euler_eos gr_em_eos[GKYL_MAX_SPECIES]; // Per-species EOS (IDEAL or APPROXIMATE_SYNGE/RCC) for relativistic-Euler primitive recovery in the GR-EM coupling.
};

// Moment-EM coupling object.
typedef struct gkyl_moment_em_coupling gkyl_moment_em_coupling;

// Status returned by the relativistic explicit GR-EM coupling advance. The
// explicit SSP-RK3 source solve is conditionally stable: it must resolve the
// (relativistic) plasma and cyclotron frequencies. When it cannot at the given
// dt, success is false and dt_suggested gives a stable time-step for the app to
// retry with.
struct gkyl_moment_em_coupling_status {
  bool success; // False if dt under-resolves the plasma/cyclotron frequency.
  double dt_suggested; // Largest stable time-step over the update range.
  double omega_max; // Maximum plasma/cyclotron frequency over the update range.
};

/**
* Create a new moment-EM coupling object, used for integrating the electromagnetic source terms that appear in the multi-fluid equations.
*
* @param inp Input parameters for the moment-EM coupling object.
* @return Moment-EM coupling object.
*/
gkyl_moment_em_coupling*
gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system using an implicit forcing solver (specifically the time-centered
* Crank-Nicolson/implicit Runge-Kutta method). The gkyl_range object to be updated should be a (non-strict) subrange of the range over which the
* fluid variable array is defined (i.e. either the fluid variable array gkyl_range object itself, or a gkyl_range object initialized using the
* gkyl_sub_range_init method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param update_range Range object over which to integrate the electromagnetic sources using an implicit time-centered method.
* @param fluid Array of fluid variables (array size = nfluids).
* @param app_accel Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs Array of RHS/source terms to be applied to the fluid variables  
*              (e.g., Braginskii transport for Euler/Isothermal Euler; gradient-based closure for ten-moment).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources Array of number density and temperature source terms.
*/
void
gkyl_moment_em_coupling_implicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES]);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system using an explicit forcing solver (specifically either the strong
* stability-preserving third-order Runge-Kutta method or the simple first-order forward-Euler method). The gkyl_range object to be updated should be a
* (non-strict) subrange of the range over which the fluid variable array is defined (i.e. either the fluid variable array gkyl_range object itself, or
* a gkyl_range object initialized using the gkyl_sub_range_init method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param update_range Range object over which to integrate the electromagnetic sources using an explicit time integration method.
* @param fluid Array of fluid variables (array size = nfluids).
* @param app_accel Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs Array of RHS/source terms to be applied to the pressure tensor (for the case of 10-moment gradient-based closure only).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param app_current1 Array of stage-1 current terms to be applied to the fluid equations (for stage-1 of external current driving).
* @param app_current2 Array of stage-2 current terms to be applied to the fluid equations (for stage-2 of external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources Array of number density and temperature source terms.
* @param proj_app_curr The finite-volume projection routine for the external current.
* @param nstrang Indicator of which step in the Strang splitting we are currently considering.
*/
void
gkyl_moment_em_coupling_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array *app_current, const struct gkyl_array* app_current1, const struct gkyl_array* app_current2,
  const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES], gkyl_fv_proj* proj_app_curr, int nstrang);

/**
* Integrate the (special-)relativistic multi-fluid + Maxwell coupling source terms using an explicit SSP-RK3 forcing solver, in SI units and with
* the fluid species, geometry, and Maxwell field supplied as SEPARATE inputs (no packed state vector). Each fluid is a modular relativistic-Euler
* state [D, S_x, S_y, S_z, tau]; the field is the standard 8-component Maxwell vector [E_x, E_y, E_z, B_x, B_y, B_z, phi, psi] (E and B are used
* directly, with no vacuum constitutive relation). The geometry is assumed to be Minkowski (flat) in this first version. The solver is conditionally
* stable and reports, via the returned status, whether the supplied dt resolves the relativistic plasma and cyclotron frequencies, along with a
* stable dt_suggested for the app to retry with on failure.
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param update_range Range object over which to integrate the coupling sources.
* @param fluid Array of fluid variables (array size = nfluids).
* @param em Array of electromagnetic variables.
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @return Status indicating whether the explicit step resolved the plasma/cyclotron frequencies, with a suggested stable time-step.
*/
struct gkyl_moment_em_coupling_status
gkyl_moment_em_coupling_gr_em_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], struct gkyl_array* em, const struct gkyl_array* ext_em);

/**
* Delete moment-EM coupling object.
*
* @param mom_em Moment-EM coupling object to delete.
*/
void
gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* mom_em);