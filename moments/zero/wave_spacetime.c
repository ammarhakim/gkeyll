#include <assert.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_spacetime.h>
#include <gkyl_wave_spacetime_priv.h>

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
  double det = gkyl_wave_spacetime_invert_metric_3x3(g_iface, inv_g);

  double M[3][3], M_inv[3][3];
  gkyl_wave_spacetime_build_triad_contravariant_x(g_iface, inv_g, M, M_inv);

  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++) {
      iface->M_inv[a][b]       = M_inv[a][b];
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

    // Cell's own products rotated into each face-local frame (the
    // per-(cell, dir) geometry consumers fetch by index — see the
    // cell_prods_local storage contract in the header). Filled for
    // every cell including excised ones: the excision flag travels
    // with the row and every consumer short-circuits on it.
    for (int d = 0; d < ndim; d++)
      gkyl_wave_spacetime_rotate_prods_row(wcg->tau1[d], wcg->tau2[d],
        wcg->norm[d], pR, cell->cell_prods_local[d]);
    for (int d = ndim; d < GKYL_MAX_CDIM; d++)
      for (int k = 0; k < GKYL_GR_SP_NCOMP_BASE; k++)
        cell->cell_prods_local[d][k] = 0.0;

    for (int d = 0; d < ndim; d++) {
      struct gkyl_wave_spacetime_iface *iface = &cell->iface[d];

      // L = neighbor at idx - δ_d. At the bottom edge of the build range
      // no lower neighbor exists, so the only choice available on the
      // local extended range is to fill this face one-sided from the
      // cell's OWN prods. wave_prop never consumes this entry — its edge
      // loop spans [lower-1, upper+2] over the update range, so with the
      // standard two-ghost-layer extended range the deepest face it
      // requests is owned by the FIRST ghost cell, whose lower neighbor
      // (the second ghost layer) is inside this range. The one-sided
      // fill exists so that every entry in the cache is nonetheless
      // VALID: a zero fill here would place sqrt_det_iface = 0 within
      // reach of a division should any consumer's loop extents ever
      // change. (If the edge cell itself is excised, fall back to the
      // BOTH_EXCISED zero fill — its metric may be unusable and every
      // consumer already short-circuits on the excision flags.)
      int idxL[GKYL_MAX_CDIM];
      for (int k = 0; k < GKYL_MAX_CDIM; k++) idxL[k] = iter.idx[k];
      idxL[d] -= 1;
      if (idxL[d] < range->lower[d]) {
        if (exc_R < 1.0e-8) {
          ws_fill_both_excised(iface);
          continue;
        }
        double gR_loc[3][3] = {{0}};
        double shiftR_loc[3] = {0};
        gkyl_wave_spacetime_rotate_rank2(wcg->tau1[d], wcg->tau2[d],
          wcg->norm[d], &pR[GKYL_GR_SP_GIJ], (double *)gR_loc);
        gkyl_wave_spacetime_rotate_rank1(wcg->tau1[d], wcg->tau2[d],
          wcg->norm[d], &pR[GKYL_GR_SP_SHIFT], shiftR_loc);
        ws_build_iface_from_rotated(iface, GKYL_WS_IFACE_INTERIOR,
          gR_loc, gR_loc, /*use_avg=*/false,
          pR[GKYL_GR_SP_LAPSE], shiftR_loc[0]);
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
        gkyl_wave_spacetime_rotate_rank2(tau1, tau2, norm,
          &pL[GKYL_GR_SP_GIJ], (double *)gL_loc);
        gkyl_wave_spacetime_rotate_rank1(tau1, tau2, norm,
          &pL[GKYL_GR_SP_SHIFT], shiftL_loc);
      }
      if (!excised_R) {
        gkyl_wave_spacetime_rotate_rank2(tau1, tau2, norm,
          &pR[GKYL_GR_SP_GIJ], (double *)gR_loc);
        gkyl_wave_spacetime_rotate_rank1(tau1, tau2, norm,
          &pR[GKYL_GR_SP_SHIFT], shiftR_loc);
      }

      if (excised_L || excised_R) {
        if (spacetime != NULL) {
          // Excision-face geometry is KNOWN — evaluate the analytic
          // metric AT THE FACE CENTROID rather than substituting the
          // active cell's center values (a half-cell offset toward
          // weaker curvature, largest exactly where gradients peak).
          // Interior faces stay on the averaging policy: for a dynamic
          // spacetime averaging is the only option, and the analytic
          // shortcut must not make the two regimes diverge.
          double xc[GKYL_MAX_CDIM] = { 0.0, 0.0, 0.0 };
          gkyl_rect_grid_cell_center(grid, iter.idx, xc);
          xc[d] -= 0.5 * grid->dx[d];   // this cell's lower-face centroid
          double pF[GKYL_GR_SP_NCOMP_BASE];
          gkyl_moment_spacetime_coupling_fill_products_analytic(
            (struct gkyl_gr_spacetime *)spacetime, t_curr,
            xc[0], xc[1], xc[2], pF);
          double gF_loc[3][3] = {{0}};
          double shiftF_loc[3] = {0};
          gkyl_wave_spacetime_rotate_rank2(tau1, tau2, norm,
            &pF[GKYL_GR_SP_GIJ], (double *)gF_loc);
          gkyl_wave_spacetime_rotate_rank1(tau1, tau2, norm,
            &pF[GKYL_GR_SP_SHIFT], shiftF_loc);
          ws_build_iface_from_rotated(iface, GKYL_WS_IFACE_EXCISION,
            gF_loc, gF_loc, /*use_avg=*/false,
            pF[GKYL_GR_SP_LAPSE], shiftF_loc[0]);
        }
        else if (excised_L) {
          // No analytic spacetime (Einstein-evolved run): active-cell
          // values. This branch is the slot an apparent-horizon finder
          // fills with its declared excision-boundary geometry.
          ws_build_iface_from_rotated(iface, GKYL_WS_IFACE_EXCISION,
            gR_loc, gR_loc, /*use_avg=*/false,
            pR[GKYL_GR_SP_LAPSE], shiftR_loc[0]);
        }
        else {
          ws_build_iface_from_rotated(iface, GKYL_WS_IFACE_EXCISION,
            gL_loc, gL_loc, /*use_avg=*/false,
            pL[GKYL_GR_SP_LAPSE], shiftL_loc[0]);
        }
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
  // Always rebuild; the app gates calls on spacetime dynamism.
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

// Build an orthonormal triad whose FIRST basis vector is aligned with the
// CONTRAVARIANT x-direction (∂^x = γ^{xj}·∂_j) rather than the coordinate
// x-direction. The other two basis vectors come from Gram-Schmidt-in-γ
// starting from (0,1,0) and (0,0,1). This eliminates the v_tet^x ↔ v^y,v^z
// mixing seen with Cholesky-on-γ for non-diagonal γ.
//
// Outputs:
//   M[i][a]    = e_a^i  (matrix of basis-vector components; columns = e_a)
//   M_inv[a][i] = (e^a)_i = γ_ij · M[j][a]  (dual co-vector components)
//
// Properties:
//   γ_ij·M[i][a]·M[j][b] = δ_ab          (orthonormality)
//   M_inv·M = M·M_inv = I (matrix product)
//   M[i][0] = γ^{xi}/√γ^{xx}              (contravariant x alignment)
//   M_inv[0][i] = δ^x_i/√γ^{xx} = (1/√γ^{xx}, 0, 0)
//
// For diagonal γ this reduces to the Cholesky construction.
void
gkyl_wave_spacetime_build_triad_contravariant_x(
  const double g_ij[3][3], const double inv_g[3][3],
  double M[3][3], double M_inv[3][3])
{
  // First basis vector: contravariant x-direction, normalized in γ.
  //   e_0^i = γ^{xi} / √γ^{xx}
  double sqrt_inv_gxx = sqrt(inv_g[0][0]);
  for (int i = 0; i < 3; i++) M[i][0] = inv_g[0][i] / sqrt_inv_gxx;

  // Helper: γ-inner-product of two contravariant vectors.
  // <u, v>_γ = γ_ij · u^i · v^j
  #define GAMMA_DOT(u, v)                                              \
    (g_ij[0][0]*(u)[0]*(v)[0] + g_ij[1][1]*(u)[1]*(v)[1]               \
   + g_ij[2][2]*(u)[2]*(v)[2]                                          \
   + g_ij[0][1]*((u)[0]*(v)[1] + (u)[1]*(v)[0])                        \
   + g_ij[0][2]*((u)[0]*(v)[2] + (u)[2]*(v)[0])                        \
   + g_ij[1][2]*((u)[1]*(v)[2] + (u)[2]*(v)[1]))

  // Second basis vector: Gram-Schmidt-in-γ of (0,1,0) against e_0.
  double e0[3] = { M[0][0], M[1][0], M[2][0] };
  double f1[3] = { 0.0, 1.0, 0.0 };
  double c01 = GAMMA_DOT(f1, e0);  // projection coefficient
  double e1[3] = { f1[0] - c01*e0[0], f1[1] - c01*e0[1], f1[2] - c01*e0[2] };
  double norm1 = sqrt(GAMMA_DOT(e1, e1));
  for (int i = 0; i < 3; i++) M[i][1] = e1[i] / norm1;

  // Third basis vector: Gram-Schmidt of (0,0,1) against e_0 and e_1.
  double e1_norm[3] = { M[0][1], M[1][1], M[2][1] };
  double f2[3] = { 0.0, 0.0, 1.0 };
  double c02 = GAMMA_DOT(f2, e0);
  double c12 = GAMMA_DOT(f2, e1_norm);
  double e2[3] = {
    f2[0] - c02*e0[0] - c12*e1_norm[0],
    f2[1] - c02*e0[1] - c12*e1_norm[1],
    f2[2] - c02*e0[2] - c12*e1_norm[2],
  };
  double norm2 = sqrt(GAMMA_DOT(e2, e2));
  for (int i = 0; i < 3; i++) M[i][2] = e2[i] / norm2;

  #undef GAMMA_DOT

  // M_inv = M^T · γ (since γ-orthonormality says M^T·γ·M = I, so
  // M^{-1} = M^T·γ).
  for (int a = 0; a < 3; a++) {
    for (int i = 0; i < 3; i++) {
      M_inv[a][i] = M[0][a]*g_ij[0][i]
                  + M[1][a]*g_ij[1][i]
                  + M[2][a]*g_ij[2][i];
    }
  }
}
