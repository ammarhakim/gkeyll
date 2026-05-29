#pragma once

// Private header for wave_spacetime — not for direct use in user code.

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_util.h>
#include <gkyl_wave_spacetime.h>

// Free function bound into ref_count_init.
void
gkyl_wave_spacetime_free(const struct gkyl_ref_count *ref);

// CPU builder. Walks @a range and fills the per-cell iface[GKYL_MAX_CDIM]
// entries against the current prods+spacetime state. Each cell owns the
// face on its lower side in each direction (wave_geom convention), so the
// iface[d] entry for cell idx is built from cell at idx-δ_d (L) and cell at
// idx (R).
void
gkyl_wave_spacetime_build_cpu(struct gkyl_wave_spacetime *ws,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_wave_geom *wg,
  const struct gkyl_gr_spacetime *spacetime,
  const struct gkyl_array *prods,
  double t_curr);
