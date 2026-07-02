#include <assert.h>
#include <math.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_spacetime_coupling.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_gr_euler_prim_priv.h>

// Single-call cascade: each repair_state invocation fixes all three
// admissibility constraints (D > 0, τ ≥ 0, s² > 0) in safe order, so the
// post-call state is always admissible. No iteration needed.

// ---------------------------------------------------------------------------
// Internal helpers.
// ---------------------------------------------------------------------------

// Compute the inverse of a 3x3 symmetric matrix using the same closed-form
// expression used by gkyl_gr_euler_inv_spatial_metric in wv_gr_euler.c. This
// is duplicated here (rather than calling the packed helper) so the
// spacetime-coupling object stays decoupled from the packed equation; once
// the packed wv_gr_euler is removed this remains the canonical path.
static inline void
compute_inv_spatial_metric(const double g[3][3], double inv_g[3][3])
{
  double spatial_det =
    g[0][0] * (g[1][1] * g[2][2] - g[2][1] * g[1][2])
  - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
  + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);

  double trace = g[0][0] + g[1][1] + g[2][2];

  double g_sq[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_sq[i][j] = 0.0;
      for (int k = 0; k < 3; k++) g_sq[i][j] += g[i][k] * g[k][j];
    }

  double sq_trace = g_sq[0][0] + g_sq[1][1] + g_sq[2][2];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double euclid = (i == j) ? 1.0 : 0.0;
      inv_g[i][j] = (1.0 / spatial_det) *
        ((0.5 * (trace * trace - sq_trace) * euclid)
         - (trace * g[i][j])
         + g_sq[i][j]);
    }
  }
}

void
gkyl_moment_spacetime_coupling_fill_products_analytic(
  struct gkyl_gr_spacetime *spacetime,
  double t_curr, double x, double y, double z, double *prods)
{
  double spatial_det, lapse;
  bool in_excision_region;

  // Allocate the heap arrays the gr_spacetime callbacks fill. Packed code
  // does the same allocation pattern; matching it keeps memory ownership
  // stories identical, and the cost is per-cell-per-step (cheap for typical
  // grid sizes).
  double *shift = gkyl_malloc(sizeof(double[3]));
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i]    = gkyl_malloc(sizeof(double[3]));
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }
  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));
  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++)
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  // Match the finite-difference epsilon used by gr_euler_impose_gauge in
  // level_set.c so derive_products numerically agrees with what the packed
  // path reinitialises into q[5..66] each sweep.
  const double eps = pow(10.0, -8.0);

  spacetime->lapse_function_func(spacetime, t_curr, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, t_curr, x, y, z, &shift);
  spacetime->spatial_metric_det_func(spacetime, t_curr, x, y, z, &spatial_det);
  spacetime->excision_region_func(spacetime, t_curr, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, t_curr, x, y, z, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, t_curr, x, y, z, eps, eps, eps, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, t_curr, x, y, z, eps, eps, eps, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, t_curr, x, y, z, eps, eps, eps, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, t_curr, x, y, z, eps, eps, eps, &spatial_metric_der);

  prods[GKYL_GR_SP_LAPSE]    = lapse;
  prods[GKYL_GR_SP_SHIFT + 0] = shift[0];
  prods[GKYL_GR_SP_SHIFT + 1] = shift[1];
  prods[GKYL_GR_SP_SHIFT + 2] = shift[2];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_GIJ + 3*i + j] = spatial_metric[i][j];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_KIJ + 3*i + j] = extrinsic_curvature[i][j];

  prods[GKYL_GR_SP_EXCISION] = in_excision_region ? -1.0 : 1.0;

  for (int i = 0; i < 3; i++)
    prods[GKYL_GR_SP_DALPHA + i] = lapse_der[i];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_DBETA + 3*i + j] = shift_der[i][j];

  for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        prods[GKYL_GR_SP_DGIJ + 9*k + 3*i + j] = spatial_metric_der[k][i][j];

  // Derived blocks: inverse spatial metric and determinant. Computing them
  // once here saves recomputation in every source-update / wave call.
  double g[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g[i][j] = spatial_metric[i][j];

  double inv_g[3][3];
  compute_inv_spatial_metric(g, inv_g);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      prods[GKYL_GR_SP_INV_GIJ + 3*i + j] = inv_g[i][j];

  prods[GKYL_GR_SP_SPATIAL_DET] = spatial_det;

  prods[GKYL_GR_SP_CELLCTR + 0] = x;
  prods[GKYL_GR_SP_CELLCTR + 1] = y;
  prods[GKYL_GR_SP_CELLCTR + 2] = z;

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]);
    gkyl_free(shift_der[i]);
    for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]);
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(extrinsic_curvature);
  gkyl_free(shift_der);
  gkyl_free(spatial_metric_der);
  gkyl_free(shift);
  gkyl_free(lapse_der);
}

// ---------------------------------------------------------------------------
// Newton-iteration histogram bin-edges. Referenced by the static-inline
// gkyl_gr_euler_recover_primitives via the extern declaration in
// gkyl_wv_gr_euler_prim_priv.h. The accumulators themselves live on
// per-callsite prim_status struct buckets owned by the app — no
// process-global storage.
// ---------------------------------------------------------------------------

const int gkyl_gr_euler_newton_bin_edges[GR_EULER_NEWTON_HIST_NBINS] = {
  4, 8, 16, 32, 64, 99, 100
};



// ---------------------------------------------------------------------------
// Generic per-cell explicit source step over a wv_eqn. The source physics
// lives on each equation object (eqn->source_func for the RHS, the optional
// eqn->source_limiter_func for a positivity limiter, eqn->repair_state_func for
// admissibility), so this integrator is equation-agnostic and is used
// uniformly for the GR-Euler fluid species and the Einstein spacetime state.
// ---------------------------------------------------------------------------

#define GKYL_MOMENT_SOURCE_MAX_NEQN 128  // max state size sourced here (conformal Einstein = 77)

// One explicit Euler step: out = q_in + alpha*dt*S(q_in), where S is the
// equation's source RHS and alpha its optional positivity limiter (1 if none);
// then repair if the equation flags q_out inadmissible (anchored on q_in).
static void
source_forward_euler(const struct gkyl_wv_eqn *eqn, double dt,
  const double *q_in, double *q_out)
{
  int n = eqn->num_equations;
  double S[GKYL_MOMENT_SOURCE_MAX_NEQN];
  eqn->source_func(eqn, q_in, S);
  double alpha = eqn->source_limiter_func
    ? eqn->source_limiter_func(eqn, q_in, S, dt) : 1.0;
  for (int j = 0; j < n; j++) q_out[j] = q_in[j] + (alpha * dt * S[j]);

  if (eqn->repair_state_func && eqn->check_inv_func
      && !gkyl_wv_eqn_check_inv(eqn, q_out))
    gkyl_wv_eqn_repair_state(eqn, q_in, q_out);
}

// RK convex combination: out = c0*a + c1*b over the equation's components.
static void
source_combine(const struct gkyl_wv_eqn *eqn, double c0, double c1,
  const double *a, const double *b, double *out)
{
  int n = eqn->num_equations;
  for (int j = 0; j < n; j++) out[j] = (c0 * a[j]) + (c1 * b[j]);
}

// ---------------------------------------------------------------------------
// Public API.
// ---------------------------------------------------------------------------

gkyl_moment_spacetime_coupling*
gkyl_moment_spacetime_coupling_new(struct gkyl_moment_spacetime_coupling_inp inp)
{
  // At most one of analytic_spacetime / einstein_eqn must be set.
  assert(!(inp.analytic_spacetime != NULL && inp.einstein_eqn != NULL));

  gkyl_moment_spacetime_coupling *st =
    gkyl_malloc(sizeof(struct gkyl_moment_spacetime_coupling));
  st->grid    = *(inp.grid);
  st->ndim    = st->grid.ndim;
  st->nfluids = inp.nfluids;
  for (int i = 0; i < inp.nfluids; i++)
    st->fluid_param[i] = inp.fluid_param[i];

  // Acquire refs on per-species equation objects so repair_state and
  // check_inv can be driven from the source-update loop.
  for (int i = 0; i < GKYL_MAX_SPECIES; i++) st->eqn[i] = NULL;
  for (int i = 0; i < inp.nfluids; i++) {
    if (inp.eqn[i]) st->eqn[i] = gkyl_wv_eqn_acquire(inp.eqn[i]);
  }

  st->is_static  = inp.is_static;
  st->has_tetrad = inp.has_tetrad;

  st->analytic_spacetime = inp.analytic_spacetime;
  st->spacetime_gauge    = inp.spacetime_gauge;
  st->reinit_freq        = inp.reinit_freq;

  st->einstein_eqn = inp.einstein_eqn;

  return st;
}

void
gkyl_moment_spacetime_coupling_derive_products(
  const gkyl_moment_spacetime_coupling *st,
  double t_curr,
  const struct gkyl_range *update_range,
  const struct gkyl_array *einstein_state,
  struct gkyl_array *prods)
{
  // Phase B (dynamic Bona-Masso) is not yet implemented; this assert will
  // fire if a user tries to drive the path before A4-Phase-B lands.
  assert(st->einstein_eqn == NULL);
  assert(st->analytic_spacetime != NULL);
  (void)einstein_state;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(update_range, iter.idx);
    double *prods_row = gkyl_array_fetch(prods, cidx);

    double xc[GKYL_MAX_DIM] = { 0.0 };
    gkyl_rect_grid_cell_center(&st->grid, iter.idx, xc);

    double x = (st->ndim >= 1) ? xc[0] : 0.0;
    double y = (st->ndim >= 2) ? xc[1] : 0.0;
    double z = (st->ndim >= 3) ? xc[2] : 0.0;

    gkyl_moment_spacetime_coupling_fill_products_analytic(
      st->analytic_spacetime, t_curr, x, y, z, prods_row);
  }
}

void
gkyl_moment_spacetime_coupling_explicit_advance(
  const gkyl_moment_spacetime_coupling *st,
  double t_curr, double dt,
  const struct gkyl_range *update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES],
  struct gkyl_array *spacetime)
{
  int nfluids = st->nfluids;
  const struct gkyl_wv_eqn *einstein_eqn = st->einstein_eqn;

  // SSP-RK3 stage coefficients: each stage advances every state by
  //   u <- c0*u^n + c1*(u_in + alpha*dt L(u_in)).
  static const double c0[3] = { 0.0, 0.75, 1.0 / 3.0 };
  static const double c1[3] = { 1.0, 0.25, 2.0 / 3.0 };

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(update_range, iter.idx);

    // Per-cell setup: copy each sourced state's start-of-step value (u^n) into a
    // scratch buffer for the SSP-RK3 combines, and arm each equation's cell
    // index (so source_func/check_inv/repair can fetch the cell's products row)
    // and source-step repair context. A species with no source equation
    // (st->eqn[s] == NULL) is skipped throughout.
    double f_old[GKYL_MAX_SPECIES][5], f_cur[GKYL_MAX_SPECIES][5];
    for (int s = 0; s < nfluids; s++) {
      const struct gkyl_wv_eqn *eqn = st->eqn[s];
      if (!eqn) continue;
      if (eqn->set_cell_idx_func) eqn->set_cell_idx_func(eqn, iter.idx);
      if (eqn->repair_state_func && eqn->check_inv_func)
        ((struct gkyl_wv_eqn *)eqn)->cur_repair_ctx = 0;  // route to source-side repair_status
      const double *f = gkyl_array_cfetch(fluid[s], cidx);
      for (int j = 0; j < 5; j++) { f_old[s][j] = f[j]; f_cur[s][j] = f[j]; }
    }

    double st_old[GKYL_MOMENT_SOURCE_MAX_NEQN];  // u^n (like f_old)
    double st_cur[GKYL_MOMENT_SOURCE_MAX_NEQN];  // working state (like f_cur)
    if (einstein_eqn) {
      const double *q = gkyl_array_cfetch(spacetime, cidx);
      for (int j = 0; j < einstein_eqn->num_equations; j++) { st_old[j] = q[j]; st_cur[j] = q[j]; }
    }

    // SSP-RK3: three stages stepping the fluid(s) AND the spacetime state on the
    // same stages, each via the same generic forward-Euler + combine. The source
    // physics (rate, positivity limiter, admissibility repair) lives on each
    // equation object, so this loop is equation-agnostic.
    for (int stage = 0; stage < 3; stage++) {
      for (int s = 0; s < nfluids; s++) {
        const struct gkyl_wv_eqn *eqn = st->eqn[s];
        if (!eqn) continue;
        double f_new[5];
        source_forward_euler(eqn, dt, f_cur[s], f_new);
        source_combine(eqn, c0[stage], c1[stage], f_old[s], f_new, f_cur[s]);
      }
      if (einstein_eqn) {
        double st_new[GKYL_MOMENT_SOURCE_MAX_NEQN];
        source_forward_euler(einstein_eqn, dt, st_cur, st_new);
        source_combine(einstein_eqn, c0[stage], c1[stage], st_old, st_new, st_cur);
      }
    }

    // Commit results.
    for (int s = 0; s < nfluids; s++) {
      if (!st->eqn[s]) continue;
      double *f = gkyl_array_fetch(fluid[s], cidx);
      for (int j = 0; j < 5; j++) f[j] = f_cur[s][j];
    }
    if (einstein_eqn) {
      double *q = gkyl_array_fetch(spacetime, cidx);
      for (int j = 0; j < einstein_eqn->num_equations; j++) q[j] = st_cur[j];
    }
  }
}

void
gkyl_moment_spacetime_coupling_release(gkyl_moment_spacetime_coupling *st)
{
  for (int i = 0; i < GKYL_MAX_SPECIES; i++) {
    if (st->eqn[i]) gkyl_wv_eqn_release(st->eqn[i]);
  }
  gkyl_free(st);
}
