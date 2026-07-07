#pragma once

// Private header for the moment field object. Holds struct moment_field --
// the Maxwell EM field evolved by wave-prop or MP -- and the moment_field_*
// API the app and time steppers consume. moment_field_init() wires the
// *_func dispatch pointers from (scheme x is_static x presence): when no
// field is configured, a no-op method set is installed and nothing is
// allocated, so a fluid-only run pays no Maxwell state. Included by
// gkyl_moment_priv.h so the app struct can hold the field by value; do not
// include in user-facing headers.

#include <gkyl_array.h>
#include <gkyl_dynvec.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment.h>
#include <gkyl_mp_scheme.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_eqn.h>

// Field data
struct moment_field {
  int ndim;
  double epsilon0, mu0;
  bool is_static; // is the field static?

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  struct gkyl_wv_eqn *maxwell; // pointer to Maxwell eqn obj

  bool has_ext_em; // flag to indicate there is external electromagnetic field
  bool ext_em_evolve; // flag to indicate external electromagnetic field is time dependent
  struct gkyl_array *ext_em; // external electromagnetic field
  gkyl_fv_proj *ext_em_proj; // projector for external electromagnetic field
  double t_ramp_E; // linear ramp for turning on external E field

  bool has_app_current; // flag to indicate there is an applied current
  bool app_current_evolve; // flag to indicate applied current is time dependent
  struct gkyl_array *app_current; // applied current
  gkyl_fv_proj *app_current_proj;  // projector for applied current
  double t_ramp_curr; // linear ramp for turning on applied currents

  bool use_explicit_em_coupling; // flag to indicate if em coupling should be explicit, defaults implicit
  bool use_gr_em_coupling; // flag to indicate explicit (special-)relativistic multi-fluid + Maxwell coupling (modular GR-Euler fluids, SI units)
  struct gkyl_array *app_current1; // additional array for applied currents (for use_explicit_em_coupling stages)
  struct gkyl_array *app_current2; // additional array for applied currents (for use_explicit_em_coupling stages)

  // Volume-based geometrical sources (see gkyl_moment_em_coupling.h).
  struct gkyl_moment_volume_sources volume_sources;

  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  struct gkyl_array *embed_mask;

  // scheme to update equations solvers and data to update fluid
  // equations
  enum gkyl_moment_scheme scheme_type;
  union {
    struct {
      gkyl_wave_prop *slvr[3]; // wave-prop solver in each direction
      struct gkyl_array *fdup, *f[4]; // arrays for updates
    };
    struct {
      gkyl_mp_scheme *mp_slvr; // monotonicity-preserving scheme
      struct gkyl_array *f0, *f1, *fnew; // arrays for updates
      struct gkyl_array *cflrate; // CFL rate in each cell
    };
  };
  struct gkyl_array *fcurr; // points to current solution (depends on scheme)

  // boundary condition type
  enum gkyl_field_bc_type lower_bct[3], upper_bct[3];
  // boundary conditions on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_energy; // integrated energy components
  bool is_first_energy_write_call; // flag for dynvec written first time

  // Methods installed by moment_field_init from (scheme_type x is_static x
  // presence); called through the moment_field_* wrappers so callers never
  // branch on scheme or field presence. The no-field wiring is all no-ops.
  double (*max_dt_func)(const gkyl_moment_app *app,
    const struct moment_field *fld);
  struct gkyl_update_status (*update_func)(gkyl_moment_app *app,
    const struct moment_field *fld, double tcurr, double dt); // wave-prop sweep
  double (*rhs_func)(gkyl_moment_app *app, struct moment_field *fld,
    const struct gkyl_array *fin, struct gkyl_array *rhs); // MP RHS
  void (*apply_bc_func)(gkyl_moment_app *app, double tcurr,
    const struct moment_field *field, struct gkyl_array *f);
  void (*copy_func)(const struct moment_field *fld,
    struct gkyl_array *dst, const struct gkyl_array *src);
  void (*release_func)(const struct moment_field *fld);
};

/** moment_field API: implemented in mom_field.c. */

/**
 * Initialize the EM field object from its input struct.
 *
 * @param mom Top-level app input.
 * @param mom_fld Field input.
 * @param app Moment app object.
 * @param fld On output, the initialized field.
 */
void moment_field_init(const struct gkyl_moment *mom,
  const struct gkyl_moment_field *mom_fld,
  struct gkyl_moment_app *app, struct moment_field *fld);

/**
 * Apply BCs to the EM field.
 *
 * @param app Moment app object.
 * @param tcurr Current simulation time.
 * @param field Field object.
 * @param f Array to apply BCs to.
 */
void moment_field_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_field *field,
  struct gkyl_array *f);

/**
 * Maximum stable time-step from the field hyperbolic update.
 *
 * @param app Moment app object.
 * @param fld Field object.
 * @return Maximum stable time-step.
 */
double moment_field_max_dt(const gkyl_moment_app *app,
  const struct moment_field *fld);

/**
 * Advance the EM field by dt (wave-prop dimensional sweep: input in
 * fld->f[0], output in fld->f[ndim]; no-op if the field is static).
 *
 * @param app Moment app object.
 * @param fld Field object.
 * @param tcurr Current simulation time.
 * @param dt Time-step to take.
 * @return Update status (success flag and suggested dt).
 */
struct gkyl_update_status moment_field_update(gkyl_moment_app *app,
  const struct moment_field *fld,
  double tcurr, double dt);

/**
 * Compute the RHS of Maxwell's equations (MP scheme).
 *
 * @param app Moment app object.
 * @param fld Field object.
 * @param fin Input EM field.
 * @param rhs On output, the RHS from the Maxwell solver.
 * @return Maximum stable time-step.
 */
double moment_field_rhs(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Copy field state arrays. Used for the time-stepper backup/commit/restore.
 *
 * @param fld Field object.
 * @param dst Destination state array.
 * @param src Source state array.
 */
void moment_field_copy(const struct moment_field *fld,
  struct gkyl_array *dst, const struct gkyl_array *src);

/**
 * Release resources allocated by the field.
 *
 * @param fld Field object to release.
 */
void moment_field_release(const struct moment_field *fld);
