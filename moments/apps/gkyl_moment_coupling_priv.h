#pragma once

// Private header for the moment source-coupling objects. Holds struct
// moment_coupling -- the operator-split source solve that ties species,
// field, and spacetime together each Strang half-step (implicit/explicit
// EM coupling, non-ideal closures, Braginskii transport, and the
// spacetime coupling for modular GR fluids) -- and struct mhd_src (MHD
// divergence-constraint sources). Included by gkyl_moment_priv.h; do not
// include in user-facing headers.

#include <gkyl_mhd_src.h>
#include <gkyl_moment.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_moment_spacetime_coupling.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_ten_moment_nn_closure.h>
#include <gkyl_wv_gr_euler_prim_priv.h>  // for gkyl_gr_euler_{prim,repair}_status

// Source data
struct moment_coupling {
  // grid for braginskii variables (braginskii variables located at cell nodes)
  struct gkyl_rect_grid non_ideal_grid;
  // local, local-ext ranges for braginskii variables (loop over nodes)
  struct gkyl_range non_ideal_local, non_ideal_local_ext;

  // Gradient-based closure solver (if present).
  struct gkyl_ten_moment_grad_closure *grad_closure_slvr[GKYL_MAX_SPECIES];
  // Neural network-based closure solver (if present).
  struct gkyl_ten_moment_nn_closure *nn_closure_slvr[GKYL_MAX_SPECIES];
  // Braginskii solver (if present).
  struct gkyl_moment_braginskii *brag_slvr;

  // array for stable time-step from non-ideal terms
  struct gkyl_array *non_ideal_cflrate[GKYL_MAX_SPECIES];
  // array for non-ideal variables
  // Braginskii variables, viscous stress tensor and heat-flux vector for Euler/Isothermal Euler
  // heat-flux tensor for ten-moment
  struct gkyl_array *non_ideal_vars[GKYL_MAX_SPECIES];
  // array for storing RHS of each species from non-ideal term updates
  // Braginskii tranport for Euler/Isothermal Euler
  // Gradient-based closure for ten-moment
  struct gkyl_array *pr_rhs[GKYL_MAX_SPECIES];
  // array for storing RHS of number density and temperature source terms
  struct gkyl_array *nT_sources[GKYL_MAX_SPECIES];

  gkyl_moment_em_coupling *slvr; // source solver function

  // Zero scratch arrays standing in for the field state when the app has no
  // field: several source paths (closures, Braginskii, volume sources,
  // implicit friction) read the EM/external-EM/applied-current arrays
  // unconditionally. NULL when a field is present.
  struct gkyl_array *nofield_em;
  struct gkyl_array *nofield_ext_em;
  struct gkyl_array *nofield_app_current;

  // Spacetime-coupling object for modular GR fluids. Non-NULL only when at
  // least one species is a modular GR variant. Drives both derive_products
  // (which fills app->spacetime.prods) and the GR source-term integration
  // on the modular species. The EM-coupling solver above continues to
  // handle self-contained GR fluids (wv_gr_euler and friends, which carry
  // the spacetime in their state vector); a species is either modular or
  // self-contained, never both.
  gkyl_moment_spacetime_coupling *spacetime_slvr;

  // Per-species GR Euler tetrad instrumentation. Indexed by species index
  // 0..num_species-1; entries for species that are not modular GR are NULL.
  // Wave-prop callsites (handled by the equation object) update the
  // _wave_prop pair; source-step callsites (handled by spacetime_slvr's
  // explicit_advance) update the _source pair.
  struct gkyl_gr_euler_prim_status   *gr_euler_prim_status_wave_prop[GKYL_MAX_SPECIES];
  struct gkyl_gr_euler_repair_status *gr_euler_repair_status_wave_prop[GKYL_MAX_SPECIES];
  struct gkyl_gr_euler_prim_status   *gr_euler_prim_status_source[GKYL_MAX_SPECIES];
  struct gkyl_gr_euler_repair_status *gr_euler_repair_status_source[GKYL_MAX_SPECIES];
};

struct mhd_src {
  gkyl_mhd_src *slvr; // source solver function
};

/** moment_coupling API: implemented in mom_coupling.c. */

/**
 * Initialize the source-coupling solver. Must be called after all species,
 * the field, and the spacetime are initialized (it reads their
 * configuration to assemble the coupled solve).
 *
 * @param app Moment app object.
 * @param src On output, the initialized coupling object.
 */
void moment_coupling_init(const struct gkyl_moment_app *app,
  struct moment_coupling *src);

/**
 * Update the operator-split sources. 'nstrang' is 0 for the first Strang
 * half-step and 1 for the second.
 *
 * @param app Moment app object.
 * @param src Coupling object.
 * @param nstrang Strang half-step index (0 or 1).
 * @param tcurr Current simulation time.
 * @param dt Time-step of the half-step.
 * @return Update status (success flag and suggested dt).
 */
struct gkyl_update_status moment_coupling_update(gkyl_moment_app *app,
  struct moment_coupling *src, int nstrang, double tcurr, double dt);

/**
 * Release resources allocated by the coupling object.
 *
 * @param app Moment app object.
 * @param src Coupling object to release.
 */
void moment_coupling_release(const struct gkyl_moment_app *app,
  const struct moment_coupling *src);

/** mhd_src API: implemented in mom_species.c. */

/**
 * Initialize the MHD divergence-constraint source solver.
 *
 * @param app Moment app object.
 * @param sp MHD species input.
 * @param src On output, the initialized MHD source object.
 */
void mhd_src_init(const struct gkyl_moment_app *app,
  const struct gkyl_moment_species *sp, struct mhd_src *src);

/**
 * Update the MHD sources. 'nstrang' is 0 for the first Strang half-step and
 * 1 for the second.
 *
 * @param app Moment app object.
 * @param src MHD source object.
 * @param nstrang Strang half-step index (0 or 1).
 * @param tcurr Current simulation time.
 * @param dt Time-step of the half-step.
 */
void mhd_src_update(gkyl_moment_app *app, struct mhd_src *src, int nstrang,
  double tcurr, double dt);

/**
 * Release resources allocated by the MHD source object.
 *
 * @param src MHD source object to release.
 */
void mhd_src_release(const struct mhd_src *src);
