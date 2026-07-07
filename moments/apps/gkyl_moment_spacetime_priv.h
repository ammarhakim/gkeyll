#pragma once

// Private header for the moment spacetime object. Holds struct
// moment_spacetime -- the GR background that drives the modular GR fluids --
// and the moment_spacetime_* API. The concrete backend (static-analytic
// metric vs dynamic Bona-Masso Einstein evolution) is selected in
// moment_spacetime_init(); each backend installs its own method set on the
// *_func dispatch pointers, and the rest of the app calls the
// moment_spacetime_* wrappers so it never branches on the backend itself.
// Modeled on the Vlasov app's vm_field dispatch. Included by
// gkyl_moment_priv.h so the app struct can hold the spacetime by value; do
// not include in user-facing headers.

#include <gkyl_array.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_moment.h>
#include <gkyl_moment_spacetime_products.h>  // layout of the prods array
#include <gkyl_wave_prop.h>
#include <gkyl_wave_spacetime.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_eqn.h>

// Spacetime data -- parallel to moment_field but for the GR background that
// drives the modular GR fluids. Owns the products array that modular fluid
// equation objects read via auxfields.
struct moment_spacetime {
  int ndim;
  bool is_static;
  bool has_tetrad;
  bool has_einstein_eqn;        // true when einstein_eqn != NULL
  int prods_ncomp;              // GKYL_GR_SP_NCOMP_BASE or *_NCOMP_TETRAD

  // Backend handles (exactly one is non-NULL).
  struct gkyl_wv_eqn *einstein_eqn;
  struct gkyl_gr_spacetime *analytic_spacetime;

  // Static-analytic gauge controls (same controls as the self-contained
  // wv_gr_euler equation).
  enum gkyl_spacetime_gauge spacetime_gauge;
  int reinit_freq;

  // Hyperbolic-update plumbing for the dynamic-Einstein backend; NULL for
  // the static-analytic backend (the wave_prop solvers are constructed only
  // when einstein_eqn is non-NULL).
  gkyl_wave_prop *slvr[3];
  struct gkyl_array *fdup;
  struct gkyl_array *f[4];      // Einstein state arrays (f[0] = current)
  struct gkyl_array *fcurr;     // pointer into f[]
  void *ctx;                    // user IC context
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // The shared spacetime-products array. Modular fluid equation objects
  // acquire pointers into this via gkyl_gr_euler_tetrad_set_auxfields. The
  // coupling object that fills this array and integrates GR sources lives
  // on struct moment_coupling.
  struct gkyl_array *prods;

  // Per-interface tetrad cache. Owned here; created once at app
  // construction by moment_spacetime_create_tetrad_cache (after the first
  // products fill, which the cache is built from) and handed to the tetrad
  // equation objects through moment_coupling_init's auxfield wiring. Its
  // contents are refreshed by every subsequent calc_products. NULL when no
  // tetrad species is present.
  struct gkyl_wave_spacetime *wave_spacetime;

  // BCs (dynamic case).
  enum gkyl_field_bc_type lower_bct[3], upper_bct[3];
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];
  struct gkyl_array *bc_buffer;

  // Embedding/excision mask (1 comp) consumed by the Einstein wave-prop; 1.0
  // everywhere unless the equation carries an embed_geo. NULL in the static
  // case (no hyperbolic update). Mirrors mom_species/mom_field.
  struct gkyl_array *embed_mask;

  // Backend methods installed by moment_spacetime_init (no spacetime vs
  // static-analytic vs dynamic-Einstein); called through the
  // moment_spacetime_* wrappers. The no-spacetime wiring is all no-ops and
  // allocates nothing.
  struct gkyl_update_status (*update_func)(gkyl_moment_app *app,
    struct moment_spacetime *sp, double tcurr, double dt);
  double (*max_dt_func)(const gkyl_moment_app *app,
    const struct moment_spacetime *sp);
  void (*copy_func)(const struct moment_spacetime *sp,
    struct gkyl_array *dst, const struct gkyl_array *src);
  void (*calc_products_func)(gkyl_moment_app *app,
    struct moment_spacetime *sp, double tcurr);
  void (*apply_bc_func)(gkyl_moment_app *app, double tcurr,
    const struct moment_spacetime *sp, struct gkyl_array *f);
  void (*release_func)(const struct moment_spacetime *sp);
};

/** moment_spacetime API: backend-agnostic wrappers that dispatch via the
 *  spacetime's *_func pointers. Implemented in mom_spacetime.c. */

/**
 * Initialize the spacetime object, dispatching on backend: static-analytic
 * metric (analytic_spacetime) vs dynamic Bona-Masso Einstein evolution
 * (einstein_eqn). Exactly one must be set in the input.
 *
 * @param mom Top-level app input.
 * @param mom_st Spacetime input.
 * @param app Moment app object.
 * @param sp On output, the initialized spacetime.
 */
void moment_spacetime_init(const struct gkyl_moment *mom,
  const struct gkyl_moment_spacetime *mom_st,
  struct gkyl_moment_app *app, struct moment_spacetime *sp);

/**
 * Apply BCs to the Einstein state array (no-op for the static-analytic
 * backend: there is no Einstein state to evolve).
 *
 * @param app Moment app object.
 * @param tcurr Current simulation time.
 * @param sp Spacetime object.
 * @param f Einstein state array to apply BCs to.
 */
void moment_spacetime_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_spacetime *sp, struct gkyl_array *f);

/**
 * Maximum stable time-step for the Einstein hyperbolic step (DBL_MAX for the
 * static-analytic backend: no CFL constraint of its own).
 *
 * @param app Moment app object.
 * @param sp Spacetime object.
 * @return Maximum stable time-step.
 */
double moment_spacetime_max_dt(const gkyl_moment_app *app,
  const struct moment_spacetime *sp);

/**
 * Advance the spacetime state by dt: no-op for the static-analytic backend;
 * for the dynamic backend, wave-prop of the Einstein state plus a refresh of
 * the derived geometry (products + tetrad cache) the fluid step consumes.
 *
 * @param app Moment app object.
 * @param sp Spacetime object.
 * @param tcurr Current simulation time.
 * @param dt Time-step to take.
 * @return Update status (success flag and suggested dt).
 */
struct gkyl_update_status moment_spacetime_update(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr, double dt);

/**
 * Copy spacetime state arrays (no-op for the static-analytic backend). Used
 * for the time-stepper backup/commit/restore.
 *
 * @param sp Spacetime object.
 * @param dst Destination state array.
 * @param src Source state array.
 */
void moment_spacetime_copy(const struct moment_spacetime *sp,
  struct gkyl_array *dst, const struct gkyl_array *src);

/**
 * Recompute the derived spacetime quantities the fluid solver consumes
 * (cell-center products + interface tetrad cache) from the current state.
 *
 * @param app Moment app object.
 * @param sp Spacetime object.
 * @param tcurr Current simulation time.
 */
void moment_spacetime_calc_products(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr);

/**
 * Create the per-interface tetrad cache if any species consumes it. Called
 * once from gkyl_moment_app_new after the first products fill; a no-op if
 * the cache exists or no tetrad species is present. The cache pointer is
 * handed to the tetrad equation objects by moment_coupling_init.
 *
 * @param app Moment app object.
 * @param sp Spacetime object.
 * @param tcurr Current simulation time.
 */
void moment_spacetime_create_tetrad_cache(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr);

/**
 * Release resources allocated by the spacetime object.
 *
 * @param sp Spacetime object to release.
 */
void moment_spacetime_release(const struct moment_spacetime *sp);
