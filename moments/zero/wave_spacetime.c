#include <assert.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_spacetime.h>
#include <gkyl_wave_spacetime_priv.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>  // invert_metric_3x3, build_triad_contravariant_x

// ---------------------------------------------------------------------------
// Local rotation helpers — mirror the static-inline copies in
// wv_gr_euler_tetrad.c so the cache builder can produce face-local-frame
// quantities without taking a dependency on internal symbols of the
// equation object. Behaviour is identical (the equation's rot_call_parity
// path uses the same algebra).
// ---------------------------------------------------------------------------

static inline void
ws_rotate_rank1(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  out[0] = in[0]*norm[0] + in[1]*norm[1] + in[2]*norm[2];
  out[1] = in[0]*tau1[0] + in[1]*tau1[1] + in[2]*tau1[2];
  out[2] = in[0]*tau2[0] + in[1]*tau2[1] + in[2]*tau2[2];
}

static inline void
ws_rotate_rank2(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
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

// Fill iface in BOTH_EXCISED state — sentinel values consistent with the
// Riemann solver's short-circuit branch (it never reads the matrices when
// it sees this kind).
static inline void
ws_fill_both_excised(struct gkyl_wave_spacetime_iface *iface)
{
  memset(iface, 0, sizeof(*iface));
  iface->kind = GKYL_WS_IFACE_BOTH_EXCISED;
}

// Build the iface struct for a single face given the active prods data
// already in face-local frame (rotated to put the face normal on x).
//
// If @a use_avg is true, gL_local/gR_local are the rotated metrics from
// the L and R cells respectively and we form g_iface = ½(gL + gR); the
// active-side lapse/shift averaging is done outside via @a alpha_in /
// @a shift_n_in already containing the average.
//
// If @a use_avg is false, gL_local is the rotated metric from the single
// active cell and gR_local is ignored.
static void
ws_build_iface_from_rotated(struct gkyl_wave_spacetime_iface *iface,
  enum gkyl_ws_iface_kind kind,
  const double gL_local[3][3], const double gR_local[3][3],
  bool use_avg, double alpha_in, double shift_n_in)
{
  double g_iface[3][3];
  if (use_avg) {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        g_iface[i][j] = 0.5 * (gL_local[i][j] + gR_local[i][j]);
  } else {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        g_iface[i][j] = gL_local[i][j];
  }

  // Phase 0 Fix 2: derive (inv_g_iface, det) consistently from a single
  // matrix inversion rather than from independent averages.
  double inv_g[3][3];
  double det = gkyl_gr_euler_tetrad_invert_metric_3x3(g_iface, inv_g);

  double M[3][3], M_inv[3][3];
  gkyl_gr_euler_tetrad_build_triad_contravariant_x(g_iface, inv_g, M, M_inv);

  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++) {
      iface->M_inv[a][b]       = M_inv[a][b];
      iface->g_iface[a][b]     = g_iface[a][b];
      iface->inv_g_iface[a][b] = inv_g[a][b];
    }
  iface->sqrt_det_iface = sqrt(det);
  iface->alpha          = alpha_in;
  iface->shift_n        = shift_n_in;
  iface->kind           = kind;
}

void
gkyl_wave_spacetime_build_cpu(struct gkyl_wave_spacetime *ws,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_wave_geom *wg,
  const struct gkyl_gr_spacetime *spacetime,
  const struct gkyl_array *prods,
  double t_curr)
{
  (void)grid;
  (void)spacetime;  // Phase 2: face-position eval at excision boundary is
                    // a Phase-2+ extension; identity-mapped pass uses the
                    // active-cell values that the current Riemann path does.

  const int ndim = ws->ndim;
  const struct gkyl_range *range = &ws->range;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long roff = gkyl_range_idx(range, iter.idx);
    struct gkyl_wave_spacetime_cell *cell = gkyl_array_fetch(ws->data, roff);

    // Pre-fetch wave_geom for this cell (holds the per-face rotation).
    const struct gkyl_wave_cell_geom *wcg = gkyl_wave_geom_get(wg, iter.idx);

    const double *pR = gkyl_array_cfetch(prods, roff);
    double exc_R = pR[GKYL_GR_SP_EXCISION];

    for (int d = 0; d < ndim; d++) {
      struct gkyl_wave_spacetime_iface *iface = &cell->iface[d];

      // L = neighbor at idx - δ_d. Skip if out of range (lower-boundary
      // ghost layer); the wave-prop inner loop never reads those ifaces.
      int idxL[GKYL_MAX_CDIM];
      for (int k = 0; k < GKYL_MAX_CDIM; k++) idxL[k] = iter.idx[k];
      idxL[d] -= 1;
      if (idxL[d] < range->lower[d]) {
        ws_fill_both_excised(iface);
        continue;
      }

      const double *pL = gkyl_array_cfetch(prods, gkyl_range_idx(range, idxL));
      double exc_L = pL[GKYL_GR_SP_EXCISION];

      bool excised_L = exc_L < 1.0e-8;
      bool excised_R = exc_R < 1.0e-8;

      if (excised_L && excised_R) {
        ws_fill_both_excised(iface);
        continue;
      }

      const double *tau1 = wcg->tau1[d];
      const double *tau2 = wcg->tau2[d];
      const double *norm = wcg->norm[d];

      double gL_loc[3][3] = {{0}}, gR_loc[3][3] = {{0}};
      double shiftL_loc[3] = {0}, shiftR_loc[3] = {0};

      if (!excised_L) {
        ws_rotate_rank2(tau1, tau2, norm,
          &pL[GKYL_GR_SP_GIJ], (double *)gL_loc);
        ws_rotate_rank1(tau1, tau2, norm,
          &pL[GKYL_GR_SP_SHIFT], shiftL_loc);
      }
      if (!excised_R) {
        ws_rotate_rank2(tau1, tau2, norm,
          &pR[GKYL_GR_SP_GIJ], (double *)gR_loc);
        ws_rotate_rank1(tau1, tau2, norm,
          &pR[GKYL_GR_SP_SHIFT], shiftR_loc);
      }

      if (excised_L) {
        // Active = R; mirror current Riemann path which uses only the
        // active cell's values at excision-adjacent faces.
        ws_build_iface_from_rotated(iface, GKYL_WS_IFACE_EXCISION,
          gR_loc, gR_loc, /*use_avg=*/false,
          pR[GKYL_GR_SP_LAPSE], shiftR_loc[0]);
      } else if (excised_R) {
        ws_build_iface_from_rotated(iface, GKYL_WS_IFACE_EXCISION,
          gL_loc, gL_loc, /*use_avg=*/false,
          pL[GKYL_GR_SP_LAPSE], shiftL_loc[0]);
      } else {
        double alpha_iface   = 0.5 * (pL[GKYL_GR_SP_LAPSE] + pR[GKYL_GR_SP_LAPSE]);
        double shift_n_iface = 0.5 * (shiftL_loc[0] + shiftR_loc[0]);
        ws_build_iface_from_rotated(iface, GKYL_WS_IFACE_INTERIOR,
          gL_loc, gR_loc, /*use_avg=*/true,
          alpha_iface, shift_n_iface);
      }
    }

    // Zero-fill unused higher-dim ifaces so they're not read uninitialized.
    for (int d = ndim; d < GKYL_MAX_CDIM; d++)
      ws_fill_both_excised(&cell->iface[d]);
  }

  ws->last_build_time = t_curr;
}

static bool
wave_spacetime_is_cu_dev(const struct gkyl_wave_spacetime *ws)
{
  return GKYL_IS_CU_ALLOC(ws->flags);
}

void
gkyl_wave_spacetime_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wave_spacetime *ws =
    container_of(ref, struct gkyl_wave_spacetime, ref_count);
  gkyl_array_release(ws->data);
  if (wave_spacetime_is_cu_dev(ws))
    gkyl_cu_free(ws->on_dev);
  gkyl_free(ws);
}

struct gkyl_wave_spacetime *
gkyl_wave_spacetime_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_wave_geom *wg,
  const struct gkyl_gr_spacetime *spacetime,
  const struct gkyl_array *prods,
  double t_curr,
  bool is_static,
  bool use_gpu)
{
  // GPU mirror deferred to a later phase; CPU-only for now matches the
  // current state of the mod-tetrad equation object.
  assert(!use_gpu);

  struct gkyl_wave_spacetime *ws = gkyl_malloc(sizeof(*ws));
  ws->ndim  = grid->ndim;
  ws->range = *range;
  ws->data  = gkyl_array_new(GKYL_USER,
    sizeof(struct gkyl_wave_spacetime_cell), range->volume);
  // GKYL_USER arrays don't support gkyl_array_clear; the builder pass that
  // runs below populates every iface entry.

  ws->is_static       = is_static;
  ws->last_build_time = -1.0;

  ws->flags = 0;
  GKYL_CLEAR_CU_ALLOC(ws->flags);
  ws->ref_count = gkyl_ref_count_init(gkyl_wave_spacetime_free);
  ws->on_dev    = ws;

  gkyl_wave_spacetime_build_cpu(ws, grid, wg, spacetime, prods, t_curr);
  return ws;
}

void
gkyl_wave_spacetime_refresh(struct gkyl_wave_spacetime *ws,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_wave_geom *wg,
  const struct gkyl_gr_spacetime *spacetime,
  const struct gkyl_array *prods,
  double t_curr)
{
  // Phase 2: always rebuild. Phase 3 will gate this on is_static +
  // last_build_time matching t_curr.
  gkyl_wave_spacetime_build_cpu(ws, grid, wg, spacetime, prods, t_curr);
}

struct gkyl_wave_spacetime *
gkyl_wave_spacetime_acquire(const struct gkyl_wave_spacetime *ws)
{
  gkyl_ref_count_inc(&ws->ref_count);
  return (struct gkyl_wave_spacetime *)ws;
}

void
gkyl_wave_spacetime_release(const struct gkyl_wave_spacetime *ws)
{
  gkyl_ref_count_dec(&ws->ref_count);
}
