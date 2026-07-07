#pragma once

// Private header for the moment species object. Holds struct moment_species --
// one finite-volume fluid (Euler, ten-moment, MHD, GR-Euler, ...) evolved by
// wave-prop or MP/KEP -- and the moment_species_* API the app and time
// steppers consume. Included by gkyl_moment_priv.h so the app struct can hold
// the species array by value; do not include in user-facing headers.

#include <gkyl_array.h>
#include <gkyl_dynvec.h>
#include <gkyl_fv_proj.h>
#include <gkyl_kep_scheme.h>
#include <gkyl_moment.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_mp_scheme.h>
#include <gkyl_ten_moment_nn_closure.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_eqn.h>

// Species data
struct moment_species {
  int ndim;
  char name[128]; // species name
  double charge, mass;

  bool is_static; // is the fluid static?

  // Does this fluid require source terms?
  // Set to true if any term which requires source terms is detected
  // including app->has_field (because the fluids couple to EM fields via sources)
  // and various moment_species inputs such as has_friction, has_volume_sources,
  // has_reactivity, has_app_accel, etc.
  bool update_sources;

  double k0; // Closure parameter (default is 0.0, used by 10 moment).
  bool has_grad_closure; // Has gradient-based closure (only for 10 moment).
  bool has_nn_closure; // Has neural network-based closure (only for 10 moment).
  struct gkyl_kann_net* ann; // Neural network architecture.
  int poly_order; // Polynomial order of learned DG coefficients.
  enum gkyl_braginskii_type type_brag; // which Braginskii equations

  bool has_friction; // Run with frictional sources.
  bool use_explicit_friction; // Use an explicit (SSP-RK3) solver for integrating frictional sources.
  double friction_Z; // Ionization number for frictional sources.
  double friction_T_elc; // Electron temperature for frictional sources.
  double friction_Lambda_ee; // Electron-electron collisional term for frictional sources.

  bool has_volume_sources; // Run with volume-based geometrical sources.
  double volume_gas_gamma; // Adiabatic index for volume-based geometrical sources.
  double volume_U0; // Initial comoving plasma velocity for volume-based geometrical sources.
  double volume_R0; // Initial radial distance from expansion/contraction center for volume-based geometrical sources.

  bool has_reactivity; // Run with reactive sources.
  double reactivity_gas_gamma; // Adiabatic index for reactive sources.
  double reactivity_specific_heat_capacity; // Specific heat capacity for reactive sources.
  double reactivity_energy_of_formation; // Energy of formation for reactive sources.
  double reactivity_ignition_temperature; // Ignition temperature for reactive sources.
  double reactivity_reaction_rate; // Reaction rate for reactive sources.

  bool has_einstein_medium; // Run with coupled fluid-Einstein sources in plane-symmetric spacetimes.
  double medium_gas_gamma; // Adiabatic index for coupled fluid-Einstein sources in plane-symmetric spacetimes.
  double medium_kappa; // Stress-energy prefactor for coupled fluid-Einstein sources in plane-symmetric spacetimes.

  bool has_gr_ultra_rel; // Run with general relativistic source terms (Euler equations, ultra-relativistic equation of state).
  double gr_ultra_rel_gas_gamma; // Adiabatic index for general relativistic Euler equations (ultra-relativistic equation of state).

  bool has_gr_euler; // Run with general relativistic source terms (Euler equations, ideal gas equation of state).
  double gr_euler_gas_gamma; // Adiabatic index for general relativistic Euler equations (ideal gas equation of state).

  bool has_gr_twofluid; // Run with general relativistic two-fluid source terms.
  double gr_twofluid_mass_elc; // Electron mass for general relativistic two-fluid equations.
  double gr_twofluid_mass_ion; // Ion mass for general relativistic two-fluid equations.
  double gr_twofluid_charge_elc; // Electron charge for general relativistic two-fluid equations.
  double gr_twofluid_charge_ion; // Ion charge for general relativistic two-fluid equations.
  double gr_twofluid_gas_gamma_elc; // Adiabatic index for electrons in general relativistic two-fluid equations.
  double gr_twofluid_gas_gamma_ion; // Adiabatic index for ions in general relativistic two-fluid equations.
  double gr_twofluid_e_fact; // Electric field divergence error propagation speed for general relativistic two-fluid equations.

  bool has_vacuum_einstein; // Run with vacuum Einstein sources using the Bona-Masso formalism.
  double vacuum_einstein_excision_threshold; // Excision threshold (lapse) for vacuum Einstein equations using the Bona-Masso formalism.
  enum gkyl_spacetime_slicing vacuum_einstein_spacetime_slicing; // Spacetime slicing condition for vacuum Einstein equations using the Bona-Masso formalism.
  enum gkyl_spacetime_evolution vacuum_einstein_spacetime_evolution; // Spacetime evolution system for vacuum Einstein equations using the Bona-Masso formalism.

  bool has_vacuum_einstein_conformal; // Run with vacuum Einstein sources using the conformal Bona-Masso formalism.
  double vacuum_einstein_conformal_excision_threshold; // Excision threshold (lapse) for vacuum Einstein equations using the conformal Bona-Masso formalism.
  enum gkyl_spacetime_slicing vacuum_einstein_conformal_spacetime_slicing; // Spacetime slicing condition for vacuum Einstein equations using the conformal Bona-Masso formalism.
  enum gkyl_spacetime_evolution vacuum_einstein_conformal_spacetime_evolution; // Spacetime evolution system for vacuum Einstein equations using the conformal Bona-Masso formalism.

  bool has_gr_mhd; // Run with general relativistic source terms (general relativistic magnetohydrodynamics equations).
  double gr_mhd_gas_gamma; // Adiabatic index for general relativistic magnetohydrodynamics equations.

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  bool has_app_accel; // flag to indicate there is applied acceleration
  bool app_accel_evolve; // flag to indicate applied acceleration is time-dependent
  struct gkyl_array *app_accel; // applied acceleration
  gkyl_fv_proj *app_accel_proj; // projector for acceleration

  struct gkyl_array *nT_source; // array for num density and temperature sources
  // projection func for num density and temperature sources
  gkyl_fv_proj *proj_nT_source;
  bool nT_source_set_only_once; // set by user
  bool nT_source_is_set; // to be set at run time

  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  struct gkyl_array *embed_mask;

  enum gkyl_eqn_type eqn_type;  // type ID of equation
  int num_equations;            // number of equations in species
  struct gkyl_wv_eqn *equation; // equation object

  enum gkyl_moment_scheme scheme_type; // scheme to update equations

  // solvers and data to update fluid equations
  union {
    struct {
      gkyl_wave_prop *slvr[3];        // wave-prop solver in each direction
      struct gkyl_array *fdup, *f[4]; // arrays for updates
    };
    struct {
      union {
        gkyl_mp_scheme *mp_slvr;   // monotonicity-preserving scheme
        gkyl_kep_scheme *kep_slvr; // KEP scheme
      };
      struct gkyl_array *f0, *f1, *fnew; // arrays for updates
      struct gkyl_array *cflrate;        // CFL rate in each cell
      struct gkyl_array *alpha;          // for shock detector
    };
  };
  struct gkyl_array *fcurr; // points to current solution (depends on scheme)

  // boundary condition type
  enum gkyl_species_bc_type lower_bct[3], upper_bct[3];
  // boundary condition solvers on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_q; // integrated conserved quantities
  bool is_first_q_write_call; // flag for dynvec written first time

  // Methods installed by moment_species_init from (scheme_type x is_static x
  // equation type); called through the moment_species_* wrappers so callers
  // never branch on scheme or species type. The scheme is a per-species
  // property, so species updated by different schemes can coexist.
  double (*max_dt_func)(const gkyl_moment_app *app,
    const struct moment_species *sp);
  struct gkyl_update_status (*update_func)(gkyl_moment_app *app,
    struct moment_species *sp, double tcurr, double dt); // wave-prop sweep
  double (*rhs_func)(gkyl_moment_app *app, struct moment_species *sp,
    const struct gkyl_array *fin, struct gkyl_array *rhs); // MP/KEP RHS
  // Hook applied after each directional sweep (e.g. the excision scrub for
  // modular GR tetrad species under a spacetime; a no-op otherwise).
  void (*post_update_func)(gkyl_moment_app *app, struct moment_species *sp,
    struct gkyl_array *f);
  void (*copy_func)(const struct moment_species *sp,
    struct gkyl_array *dst, const struct gkyl_array *src);
  void (*release_func)(const struct moment_species *sp);
};

/** moment_species API: implemented in mom_species.c. */

/**
 * Initialize a moment species object from its input struct.
 *
 * @param mom Top-level app input.
 * @param mom_sp Species input.
 * @param app Moment app object.
 * @param sp On output, the initialized species.
 */
void moment_species_init(const struct gkyl_moment *mom,
  const struct gkyl_moment_species *mom_sp,
  struct gkyl_moment_app *app,
  struct moment_species *sp);

/**
 * Apply BCs to species data.
 *
 * @param app Moment app object.
 * @param tcurr Current simulation time.
 * @param sp Species object.
 * @param f Array to apply BCs to.
 */
void moment_species_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_species *sp,
  struct gkyl_array *f);

/**
 * Maximum stable time-step from the species hyperbolic update.
 *
 * @param app Moment app object.
 * @param sp Species object.
 * @return Maximum stable time-step.
 */
double moment_species_max_dt(const gkyl_moment_app *app,
  const struct moment_species *sp);

/**
 * Advance the species by dt (wave-prop dimensional sweep: input in sp->f[0],
 * output in sp->f[ndim]).
 *
 * @param app Moment app object.
 * @param sp Species object.
 * @param tcurr Current simulation time.
 * @param dt Time-step to take.
 * @return Update status (success flag and suggested dt).
 */
struct gkyl_update_status moment_species_update(gkyl_moment_app *app,
  struct moment_species *sp,
  double tcurr, double dt);

/**
 * Compute the RHS of the moment equations (MP/KEP schemes).
 *
 * @param app Moment app object.
 * @param species Species object.
 * @param fin Input conserved variables.
 * @param rhs On output, the RHS from the fluid solver.
 * @return Maximum stable time-step.
 */
double moment_species_rhs(gkyl_moment_app *app, struct moment_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Copy species state arrays. Used for the time-stepper
 * backup/commit/restore.
 *
 * @param sp Species object.
 * @param dst Destination state array.
 * @param src Source state array.
 */
void moment_species_copy(const struct moment_species *sp,
  struct gkyl_array *dst, const struct gkyl_array *src);

/**
 * Release resources allocated by the species.
 *
 * @param sp Species object to release.
 */
void moment_species_release(const struct moment_species *sp);
