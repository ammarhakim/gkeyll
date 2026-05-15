#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>

// Per-fluid configuration for the spacetime-coupling updater. Only the
// modular GR fluid types (GKYL_EQN_GR_EULER_MOD, GKYL_EQN_GR_EULER_TETRAD_MOD)
// are valid here; other equation types are owned by moment_em_coupling.
struct gkyl_moment_spacetime_coupling_data {
  enum gkyl_eqn_type type;  // Equation type (mod variant).
  double gas_gamma;         // Adiabatic index.
};

// Input for moment_spacetime_coupling. Mirrors gkyl_moment_em_coupling_inp's
// shape but holds only the fields that the spacetime path actually needs.
struct gkyl_moment_spacetime_coupling_inp {
  const struct gkyl_rect_grid *grid;  // Grid (used to compute cell centers).
  int nfluids;                         // Number of fluid species.
  struct gkyl_moment_spacetime_coupling_data fluid_param[GKYL_MAX_SPECIES];

  // Per-species equation objects. The coupling acquires a ref on each
  // (via gkyl_wv_eqn_acquire) at construction so it can call equation
  // function pointers (set_cell_idx_func, check_inv_func, repair_state_func)
  // from the source-update loop. Entries for non-mod species may be NULL.
  struct gkyl_wv_eqn *eqn[GKYL_MAX_SPECIES];

  bool is_static;       // True when the spacetime never evolves; the app
                        // calls derive_products once at IC and never again.
  bool has_tetrad;      // True when at least one species needs the tetrad
                        // add-on block (e.g. wv_gr_euler_tetrad_mod). The
                        // size of the products array determines this.

  // Phase A — static analytic spacetime: derive_products evaluates the
  // gr_spacetime callbacks at each cell center. Exactly one of
  // analytic_spacetime and einstein_eqn must be non-NULL.
  struct gkyl_gr_spacetime *analytic_spacetime;
  enum gkyl_spacetime_gauge spacetime_gauge;
  int reinit_freq;

  // Phase B — dynamic Bona-Masso evolver: derive_products reads from the
  // evolved Einstein-state array passed at derive time. NULL in Phase A.
  struct gkyl_wv_eqn *einstein_eqn;
};

// Opaque handle.
typedef struct gkyl_moment_spacetime_coupling gkyl_moment_spacetime_coupling;

/**
 * Construct a spacetime-coupling updater. The returned object owns no fluid
 * arrays; callers pass the fluid arrays and products array at advance time.
 */
gkyl_moment_spacetime_coupling*
gkyl_moment_spacetime_coupling_new(struct gkyl_moment_spacetime_coupling_inp inp);

/**
 * Fill the spacetime-products array from the current spacetime state. For
 * an analytic background this evaluates gr_spacetime callbacks at each
 * cell center. For a dynamic Bona-Masso evolver (Phase B) the einstein_state
 * array supplies the lapse/shift/metric/K block, which is unpacked into the
 * products layout (derived inverse-metric and determinant are recomputed
 * each time so the products array stays self-consistent).
 *
 * The is_static gating lives at the app layer: for static spacetimes the
 * caller invokes this once at IC and never again; this updater always
 * produces fresh products when invoked.
 *
 * @param st            Coupling object.
 * @param t_curr        Current simulation time (passed through to callbacks).
 * @param update_range  Range over which to fill products.
 * @param einstein_state NULL for analytic; for dynamic spacetimes the
 *                       evolved Einstein state array.
 * @param prods          Spacetime-products array to fill. Must have ncomp
 *                       matching the configured products layout.
 */
void
gkyl_moment_spacetime_coupling_derive_products(
  const gkyl_moment_spacetime_coupling *st,
  double t_curr,
  const struct gkyl_range *update_range,
  const struct gkyl_array *einstein_state,
  struct gkyl_array *prods);

/**
 * Integrate the geometric source terms on the per-species hydro fluid
 * arrays using an SSP-RK3 forward-Euler scheme that mirrors
 * explicit_gr_euler_source_update in sources_explicit.c, but reads
 * spacetime from the products array instead of from packed q[5..66].
 *
 * The mod fluid arrays carry only the 5 hydro components; the products
 * array supplies all spacetime needed by the source-term math.
 *
 * @param st            Coupling object.
 * @param t_curr        Current simulation time.
 * @param dt            Time step.
 * @param update_range  Range over which to integrate sources.
 * @param fluid         Per-species writable hydro arrays.
 * @param prods         Current spacetime-products array.
 */
void
gkyl_moment_spacetime_coupling_explicit_advance(
  const gkyl_moment_spacetime_coupling *st,
  double t_curr, double dt,
  const struct gkyl_range *update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES],
  const struct gkyl_array *prods);

/**
 * Release the coupling object.
 */
void
gkyl_moment_spacetime_coupling_release(gkyl_moment_spacetime_coupling *st);
