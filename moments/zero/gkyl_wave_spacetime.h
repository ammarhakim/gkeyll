#pragma once

#include <gkyl_array.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>

// Excision tag attached to each interface struct. Drives the tetrad-first
// Riemann scaffolding's branch logic without re-checking flags every call.
enum gkyl_ws_iface_kind {
  GKYL_WS_IFACE_INTERIOR = 0,    // L and R cells both active; arithmetic mean.
  GKYL_WS_IFACE_EXCISION = 1,    // exactly one side excised; active-cell values.
  GKYL_WS_IFACE_BOTH_EXCISED = 2 // both sides excised; Riemann solver skips.
};

// Per-interface tetrad cache: everything wave_tetrad_high_order needs to do
// the Gorard+25 transform without re-running the per-call averaging+triad
// build. All quantities are in the face-local frame (face normal = x), i.e.
// the same frame the wave-propagation scheme rotates conservatives into.
//
// `M_inv` is the Gram-Schmidt-on-γ⁻¹ back-transform; `inv_g_iface` is the
// interface-averaged inverse spatial metric (used by q_to_tetrad_contra);
// `sqrt_det_iface` is sqrt of det(g_iface) computed consistently with
// `inv_g_iface` (Phase 0 Fix 2); `alpha` and `shift_n` are the interface
// lapse and face-normal shift used in the speed back-transform.
struct gkyl_wave_spacetime_iface {
  double M_inv[3][3];
  double g_iface[3][3];      // forward spatial metric γ_ij at the face
                             // (face-local frame; face normal = x).
  double inv_g_iface[3][3];  // inverse spatial metric γ^{ij} at the face.
  double sqrt_det_iface;
  double alpha;
  double shift_n;
  enum gkyl_ws_iface_kind kind;
};

// Per-cell entry. Each cell owns the lower face in each direction (matching
// wave_geom's convention). Allocate up to GKYL_MAX_CDIM directions even if
// the grid is lower-dimensional; the unused entries cost ~200 bytes/cell
// and keep the indexing uniform.
struct gkyl_wave_spacetime_cell {
  struct gkyl_wave_spacetime_iface iface[GKYL_MAX_CDIM];
};

// Object header — mirrors gkyl_wave_geom (range + data array + flags +
// refcount + on_dev). `is_static` is informational; Phase 2 always rebuilds
// on refresh, Phase 3 will gate this.
struct gkyl_wave_spacetime {
  int ndim;
  struct gkyl_range range;
  struct gkyl_array *data;   // gkyl_wave_spacetime_cell per cell over `range`.

  bool is_static;
  double last_build_time;

  uint32_t flags;
  struct gkyl_ref_count ref_count;
  struct gkyl_wave_spacetime *on_dev;
};

/**
 * Construct a wave_spacetime cache. The constructor immediately runs the
 * builder once at time @a t_curr; callers can subsequently call
 * gkyl_wave_spacetime_refresh to rebuild.
 *
 * @param grid      Grid the cache is defined on.
 * @param range     Range over which to construct entries (typically the
 *                  configuration-space range with ghost layers).
 * @param wg        Wave geometry providing the per-face tau1/tau2/norm
 *                  rotation used to push prods from the global frame into
 *                  each face-local frame before averaging.
 * @param spacetime Analytic spacetime callbacks (for face-position eval at
 *                  excision-adjacent faces). May be NULL in pure-interior
 *                  setups (no excision).
 * @param prods     Cell-centered spacetime products array; the cache reads
 *                  lapse/shift/γ/inv-γ/excision flags from here. Layout per
 *                  gkyl_moment_spacetime_products.h.
 * @param t_curr    Time at construction (used only as a label).
 * @param is_static Hint for future Phase 3 short-circuiting; ignored in
 *                  Phase 2 (refresh always rebuilds).
 * @param use_gpu   GPU mirror flag.
 */
struct gkyl_wave_spacetime *
gkyl_wave_spacetime_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_wave_geom *wg,
  const struct gkyl_gr_spacetime *spacetime,
  const struct gkyl_array *prods,
  double t_curr,
  bool is_static,
  bool use_gpu);

/**
 * Rebuild the cache against the current spacetime/prods state. No-op only
 * when @a is_static was set true at construction AND the requested time
 * matches `last_build_time` exactly (Phase 3 short-circuit; in Phase 2
 * always rebuilds and updates `last_build_time`).
 */
void
gkyl_wave_spacetime_refresh(struct gkyl_wave_spacetime *ws,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_wave_geom *wg,
  const struct gkyl_gr_spacetime *spacetime,
  const struct gkyl_array *prods,
  double t_curr);

/**
 * Inline accessor returning the per-cell entry for cell index @a idx.
 * Mirrors gkyl_wave_geom_get; runs on host and device.
 */
GKYL_CU_DH
static inline const struct gkyl_wave_spacetime_cell *
gkyl_wave_spacetime_get(const struct gkyl_wave_spacetime *ws, const int *idx)
{
  return (const struct gkyl_wave_spacetime_cell *)
    gkyl_array_cfetch(ws->data, gkyl_range_idx(&ws->range, idx));
}

/**
 * Acquire a borrowed reference (refcount++). Release with
 * gkyl_wave_spacetime_release.
 */
struct gkyl_wave_spacetime *
gkyl_wave_spacetime_acquire(const struct gkyl_wave_spacetime *ws);

/**
 * Release a reference. Frees the object when the refcount reaches zero.
 */
void
gkyl_wave_spacetime_release(const struct gkyl_wave_spacetime *ws);
