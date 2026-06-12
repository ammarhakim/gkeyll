#pragma once

#include <gkyl_array.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_moment_spacetime_products.h>
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
  double inv_g_iface[3][3];  // inverse spatial metric γ^{ij} at the face
                             // (face-local frame; face normal = x).
  double sqrt_det_iface;
  double alpha;
  double shift_n;
  enum gkyl_ws_iface_kind kind;
};

// Per-cell entry. Each cell owns the lower face in each direction (matching
// wave_geom's convention). Allocate up to GKYL_MAX_CDIM directions even if
// the grid is lower-dimensional; the unused entries cost ~200 bytes/cell
// and keep the indexing uniform.
//
// `cell_prods_local[d]` is the cell's OWN spacetime-products row rotated
// into the face-local frame of direction d (face normal = x), in the
// standard GKYL_GR_SP_* row layout so consumers keep their
// `const double *prods` signatures. This is the per-(cell, dir) geometry
// that used to be smuggled into the equation object through the
// rotation-parity contract (prodl_local/prodr_local) — consumers now
// fetch it by (idx, dir) like every other cached quantity
// (WAVE_SPACETIME_PARITY_PLAN.md). Unused dirs (d ≥ ndim) are
// zero-filled.
//
// Storage contract for the lower-face convention at the edge of the build
// range: cells at the BOTTOM edge of the range have no lower neighbor, so
// their iface[d] cannot be built from an L/R pair — the only choice
// available on the local extended range is a one-sided fill from the
// cell's own prods (tagged INTERIOR). wave_prop never consumes that
// entry: its edge loop spans [lower-1, upper+2] over the update range,
// so with the standard two-ghost-layer extended range the deepest face
// it requests belongs to the first ghost cell, whose lower neighbor (the
// second ghost layer) is inside the build range. The one-sided fill
// guarantees that every active-cell entry in the cache is nonetheless
// valid (sqrt_det_iface > 0) regardless of consumer loop extents; only
// faces whose owning cell is itself excised carry the zero-filled
// BOTH_EXCISED entry, and those are short-circuited by every consumer
// via the excision flags before the cache is read.
struct gkyl_wave_spacetime_cell {
  struct gkyl_wave_spacetime_iface iface[GKYL_MAX_CDIM];
  double cell_prods_local[GKYL_MAX_CDIM][GKYL_GR_SP_NCOMP_BASE];
};

// Object header — mirrors gkyl_wave_geom (range + data array + flags +
// refcount + on_dev). Dynamism is the APP's concern: the app calls
// refresh only for Einstein-evolved (non-static) spacetimes, and
// refresh always rebuilds — the cache carries no static/dynamic state.
struct gkyl_wave_spacetime {
  int ndim;
  struct gkyl_range range;
  struct gkyl_array *data;   // gkyl_wave_spacetime_cell per cell over `range`.

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
 * @param spacetime Analytic spacetime callbacks. When non-NULL,
 *                  excision-adjacent faces are built by evaluating the
 *                  metric AT THE FACE CENTROID (the excision geometry is
 *                  known exactly). When NULL (Einstein-evolved runs),
 *                  those faces fall back to the active cell's center
 *                  values — the slot an apparent-horizon finder fills.
 *                  Interior faces always use L/R averaging (for a
 *                  dynamic spacetime there is no other option).
 * @param prods     Cell-centered spacetime products array; the cache reads
 *                  lapse/shift/γ/inv-γ/excision flags from here. Layout per
 *                  gkyl_moment_spacetime_products.h.
 * @param t_curr    Time at construction (used only as a label).
 * @param use_gpu   GPU mirror flag.
 */
struct gkyl_wave_spacetime *
gkyl_wave_spacetime_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_wave_geom *wg,
  const struct gkyl_gr_spacetime *spacetime,
  const struct gkyl_array *prods,
  double t_curr,
  bool use_gpu);

/**
 * Rebuild the cache against the current spacetime/prods state. Always
 * rebuilds and updates `last_build_time`; the app invokes this only for
 * Einstein-evolved (dynamic) spacetimes — static runs build once at
 * construction.
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

// ---------------------------------------------------------------------------
// Per-interface tetrad geometry helpers. These are the operations the
// cache builder performs per face — and that the tetrad equation's
// cache-less fallback path performs per call to reproduce the cache
// bit-identically — so they are owned here, not by any equation object.
// ---------------------------------------------------------------------------

/**
 * Rotate a rank-1 spatial object from the global frame into the
 * face-local frame (rows: norm, tau1, tau2).
 */
static inline void
gkyl_wave_spacetime_rotate_rank1(const double *tau1, const double *tau2,
  const double *norm, const double *in, double *out)
{
  out[0] = in[0]*norm[0] + in[1]*norm[1] + in[2]*norm[2];
  out[1] = in[0]*tau1[0] + in[1]*tau1[1] + in[2]*tau1[2];
  out[2] = in[0]*tau2[0] + in[1]*tau2[1] + in[2]*tau2[2];
}

/**
 * Rotate a rank-2 spatial object (row-major 9 components) into the
 * face-local frame: out = R·in·Rᵀ with R rows (norm, tau1, tau2).
 */
static inline void
gkyl_wave_spacetime_rotate_rank2(const double *tau1, const double *tau2,
  const double *norm, const double *in, double *out)
{
  const double *R[3] = { norm, tau1, tau2 };
  double tmp[3][3];
  for (int a = 0; a < 3; a++)
    for (int i = 0; i < 3; i++)
      tmp[a][i] = in[3*i+0]*R[a][0] + in[3*i+1]*R[a][1] + in[3*i+2]*R[a][2];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      out[3*a + b] = tmp[a][0]*R[b][0] + tmp[a][1]*R[b][1] + tmp[a][2]*R[b][2];
}

/**
 * Rotate a rank-3 spatial object (27 components, row-major over
 * (k, i, j)) into the face-local frame: out = R·R·R·in with R rows
 * (norm, tau1, tau2). Used for the metric-derivative block ∂_k γ_ij.
 */
static inline void
gkyl_wave_spacetime_rotate_rank3(const double *tau1, const double *tau2,
  const double *norm, const double *in, double *out)
{
  const double *R[3] = { norm, tau1, tau2 };
  double tmp1[3][3][3];
  for (int a = 0; a < 3; a++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        tmp1[a][j][k] =
          in[0*9 + 3*j + k]*R[a][0] +
          in[1*9 + 3*j + k]*R[a][1] +
          in[2*9 + 3*j + k]*R[a][2];
  double tmp2[3][3][3];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      for (int k = 0; k < 3; k++)
        tmp2[a][b][k] =
          tmp1[a][0][k]*R[b][0] +
          tmp1[a][1][k]*R[b][1] +
          tmp1[a][2][k]*R[b][2];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      for (int c = 0; c < 3; c++)
        out[9*a + 3*b + c] =
          tmp2[a][b][0]*R[c][0] +
          tmp2[a][b][1]*R[c][1] +
          tmp2[a][b][2]*R[c][2];
}

/**
 * Rotate a full spacetime-products row (GKYL_GR_SP_* layout) from the
 * global frame into a face-local frame: scalars copy through; shift,
 * ∂α, and cell-center rotate rank-1; γ, K, ∂β, γ⁻¹ rotate rank-2;
 * ∂γ rotates rank-3. Used by the cache builder to fill
 * cell_prods_local[d] and defines the face-local row semantics every
 * consumer sees.
 */
static inline void
gkyl_wave_spacetime_rotate_prods_row(const double *tau1, const double *tau2,
  const double *norm, const double *in_p, double *out_p)
{
  out_p[GKYL_GR_SP_LAPSE]        = in_p[GKYL_GR_SP_LAPSE];
  out_p[GKYL_GR_SP_EXCISION]     = in_p[GKYL_GR_SP_EXCISION];
  out_p[GKYL_GR_SP_SPATIAL_DET]  = in_p[GKYL_GR_SP_SPATIAL_DET];

  gkyl_wave_spacetime_rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_SHIFT],   &out_p[GKYL_GR_SP_SHIFT]);
  gkyl_wave_spacetime_rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_DALPHA],  &out_p[GKYL_GR_SP_DALPHA]);
  gkyl_wave_spacetime_rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_CELLCTR], &out_p[GKYL_GR_SP_CELLCTR]);

  gkyl_wave_spacetime_rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_GIJ],     &out_p[GKYL_GR_SP_GIJ]);
  gkyl_wave_spacetime_rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_KIJ],     &out_p[GKYL_GR_SP_KIJ]);
  gkyl_wave_spacetime_rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_DBETA],   &out_p[GKYL_GR_SP_DBETA]);
  gkyl_wave_spacetime_rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_INV_GIJ], &out_p[GKYL_GR_SP_INV_GIJ]);

  gkyl_wave_spacetime_rotate_rank3(tau1, tau2, norm, &in_p[GKYL_GR_SP_DGIJ],    &out_p[GKYL_GR_SP_DGIJ]);
}

/**
 * Invert a 3×3 symmetric matrix. Returns det(g); fills inv_g with g^{-1}.
 * Used at the interface to derive a *consistent* (inv_g, sqrt(det)) pair
 * from g_iface, avoiding the O((Δγ)²) wave-sum residual that arises when
 * these are independently averaged from L/R cells
 * (TETRAD_REFACTOR_PLAN.md Phase 0 Fix 2).
 */
static inline double
gkyl_wave_spacetime_invert_metric_3x3(const double g[3][3],
  double inv_g[3][3])
{
  double det = g[0][0] * (g[1][1]*g[2][2] - g[1][2]*g[1][2])
             - g[0][1] * (g[0][1]*g[2][2] - g[1][2]*g[0][2])
             + g[0][2] * (g[0][1]*g[1][2] - g[1][1]*g[0][2]);
  double inv_det = 1.0 / det;
  inv_g[0][0] =  (g[1][1]*g[2][2] - g[1][2]*g[1][2]) * inv_det;
  inv_g[0][1] = -(g[0][1]*g[2][2] - g[1][2]*g[0][2]) * inv_det;
  inv_g[0][2] =  (g[0][1]*g[1][2] - g[1][1]*g[0][2]) * inv_det;
  inv_g[1][0] = inv_g[0][1];
  inv_g[1][1] =  (g[0][0]*g[2][2] - g[0][2]*g[0][2]) * inv_det;
  inv_g[1][2] = -(g[0][0]*g[1][2] - g[0][1]*g[0][2]) * inv_det;
  inv_g[2][0] = inv_g[0][2];
  inv_g[2][1] = inv_g[1][2];
  inv_g[2][2] =  (g[0][0]*g[1][1] - g[0][1]*g[0][1]) * inv_det;
  return det;
}

/**
 * Build a Gram-Schmidt-on-γ⁻¹ triad: e_0 aligned with the contravariant
 * x-direction (∂^x = γ^{xj}·∂_j), e_1, e_2 orthogonalized in γ.
 * M[i][a] = e_a^i; M_inv = Mᵀ·γ. Eliminates the v_tet^x ↔ v^y, v^z
 * mixing seen with Cholesky-on-γ for non-diagonal γ
 * (SESSION_NOTES_2.md §12).
 *
 * EXTERN by design (defined in wave_spacetime.c): under -O3 -ffast-math,
 * call-boundary placement of this function affects FP reordering and
 * therefore the bit-stability of cached M_inv values and production
 * baselines. Do not convert to inline without a baseline-refresh
 * checkpoint.
 */
GKYL_CU_D
void
gkyl_wave_spacetime_build_triad_contravariant_x(
  const double g_ij[3][3], const double inv_g[3][3],
  double M[3][3], double M_inv[3][3]);

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
