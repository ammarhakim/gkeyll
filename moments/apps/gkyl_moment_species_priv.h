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

  // Source-family parameter bundles, copied verbatim from the input (see
  // gkyl_moment_em_coupling.h); moment_coupling_init merges them into the
  // source-solver input.
  struct gkyl_moment_friction friction; // Frictional sources.
  struct gkyl_moment_volume_sources volume_sources; // Volume-based geometrical sources.
  struct gkyl_moment_reactivity reactivity; // Reactive sources.
  struct gkyl_moment_einstein_medium einstein_medium; // Coupled fluid-Einstein sources (plane-symmetric spacetimes).
  struct gkyl_moment_gr_sources gr_ultra_rel; // GR Euler sources, ultra-relativistic equation of state.
  struct gkyl_moment_gr_sources gr_euler; // GR Euler sources, ideal gas equation of state.
  struct gkyl_moment_gr_twofluid gr_twofluid; // GR two-fluid sources.
  struct gkyl_moment_vacuum_einstein vacuum_einstein; // Vacuum Einstein sources, Bona-Masso formalism.
  struct gkyl_moment_vacuum_einstein vacuum_einstein_conformal; // Vacuum Einstein sources, conformal Bona-Masso formalism.
  struct gkyl_moment_gr_sources gr_mhd; // GR MHD sources.

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
