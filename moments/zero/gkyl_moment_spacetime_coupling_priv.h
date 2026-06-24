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
#include <gkyl_wv_gr_euler_prim_priv.h>

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

