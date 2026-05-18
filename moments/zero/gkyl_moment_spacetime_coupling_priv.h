#pragma once

// Private header for moment_spacetime_coupling. Not for direct inclusion in
// user-facing code.

#include <stdint.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_moment_spacetime_coupling.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>

struct gkyl_moment_spacetime_coupling {
  struct gkyl_rect_grid grid;
  int ndim;
  int nfluids;
  struct gkyl_moment_spacetime_coupling_data fluid_param[GKYL_MAX_SPECIES];

  // Acquired refs on per-species equation objects. Used by the source-
  // update loop to invoke set_cell_idx_func, check_inv_func, and
  // repair_state_func. NULL for non-mod species.
  struct gkyl_wv_eqn *eqn[GKYL_MAX_SPECIES];

  bool is_static;
  bool has_tetrad;

  // Analytic-spacetime backend.
  struct gkyl_gr_spacetime *analytic_spacetime;
  enum gkyl_spacetime_gauge spacetime_gauge;
  int reinit_freq;

  // Dynamic Bona-Masso backend (Phase B).
  struct gkyl_wv_eqn *einstein_eqn;
};

/**
 * Per-cell helper: fill one row of the products array by evaluating the
 * gr_spacetime callbacks at (t, x, y, z). Exposed for unit testing the
 * derive-products path without needing a gkyl_array.
 */
void
gkyl_moment_spacetime_coupling_fill_products_analytic(
  struct gkyl_gr_spacetime *spacetime,
  double t_curr, double x, double y, double z, double *prods);

/**
 * Per-cell forward-Euler source-term integrator for the modular GR Euler
 * equation. Mirrors explicit_gr_euler_source_update_euler in
 * sources_explicit.c but reads spacetime from the per-cell products row
 * instead of from packed q[5..66]. Operates on the 5-element hydro state.
 */
void
gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
  double gas_gamma, double t_curr, double dt,
  const double *prods,
  const double fluid_old[5], double fluid_new[5]);

// Diagnostic: count of times the τ-positivity limiter has scaled down
// a source step to keep τ ≥ 0. Reset is not exposed (process-lifetime
// monotonic counter).
uint64_t
gkyl_moment_spacetime_coupling_tau_limiter_fires(void);

// Diagnostic: count of times the s²-positivity limiter has scaled down
// a source step (further than the τ-limiter would have) to keep
// s²(q_new) ≥ margin·(D+τ)² in the curved metric.
uint64_t
gkyl_moment_spacetime_coupling_s2_limiter_fires(void);
