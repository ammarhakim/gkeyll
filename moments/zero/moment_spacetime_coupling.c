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

// Compute the geometric source rate vector S(q) for the modular GR Euler
// equation. dq/dt = S(q) for the Banyuls source terms; integrators wrap
// this to take a time step. Excised cells return zero.
//
// Returns the rate of change of the densitized conservative state
// d(√γ·U)/dt = α√γ·S(w), so callers get the same scaling the original
// inlined source step used.
static void
compute_source_rate(double gas_gamma, const double *prods,
  const double q[5], double S_rate[5])
{
  bool in_excision_region = prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (in_excision_region) {
    for (int i = 0; i < 5; i++) S_rate[i] = 0.0;
    return;
  }

  double lapse   = prods[GKYL_GR_SP_LAPSE];
  double shift_x = prods[GKYL_GR_SP_SHIFT + 0];
  double shift_y = prods[GKYL_GR_SP_SHIFT + 1];
  double shift_z = prods[GKYL_GR_SP_SHIFT + 2];
  double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];

  // Convention A primitive recovery via the shared helper.
  double D    = q[0] / sqrt(spatial_det);
  double momx = q[1] / sqrt(spatial_det);
  double momy = q[2] / sqrt(spatial_det);
  double momz = q[3] / sqrt(spatial_det);
  double Etot = q[4] / sqrt(spatial_det);

  const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
  double inv_g[3][3] = {
    { ig[0], ig[1], ig[2] },
    { ig[3], ig[4], ig[5] },
    { ig[6], ig[7], ig[8] },
  };

  struct gkyl_gr_euler_prim prim;
  gkyl_gr_euler_recover_primitives(gas_gamma,
    D, momx, momy, momz, Etot, inv_g, &prim);

  double rho = prim.rho;
  double vx  = prim.v[0];
  double vy  = prim.v[1];
  double vz  = prim.v[2];
  double p   = prim.p;
  double W   = prim.W;
  double h   = prim.h;

  // Spacetime 4-velocity and contravariant 4-metric.
  double u4[4];
  u4[0] = W / lapse;
  u4[1] = (W * vx) - (shift_x * (W / lapse));
  u4[2] = (W * vy) - (shift_y * (W / lapse));
  u4[3] = (W * vz) - (shift_z * (W / lapse));

  double shift[3] = { shift_x, shift_y, shift_z };

  double inv_g4[4][4];
  inv_g4[0][0] = - (1.0 / (lapse * lapse));
  for (int i = 0; i < 3; i++) {
    inv_g4[0][i + 1] = (1.0 / (lapse * lapse)) * shift[i];
    inv_g4[i + 1][0] = (1.0 / (lapse * lapse)) * shift[i];
  }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      inv_g4[i + 1][j + 1] = prods[GKYL_GR_SP_INV_GIJ + 3*i + j]
        - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);

  // Perfect-fluid stress-energy tensor.
  double T[4][4];
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      T[i][j] = (rho * h * u4[i] * u4[j]) + (p * inv_g4[i][j]);

  // Spacetime-derivative lookups from products.
  double lapse_der[3] = {
    prods[GKYL_GR_SP_DALPHA + 0],
    prods[GKYL_GR_SP_DALPHA + 1],
    prods[GKYL_GR_SP_DALPHA + 2],
  };

  double shift_der[3][3];
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
      shift_der[j][i] = prods[GKYL_GR_SP_DBETA + 3*j + i];

  double spatial_metric_der[3][3][3];
  for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        spatial_metric_der[k][i][j] = prods[GKYL_GR_SP_DGIJ + 9*k + 3*i + j];

  double extrinsic_curvature[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      extrinsic_curvature[i][j] = prods[GKYL_GR_SP_KIJ + 3*i + j];

  const double *g_ij = &prods[GKYL_GR_SP_GIJ];

  double shift_lower[3] = { 0.0, 0.0, 0.0 };
  for (int m = 0; m < 3; m++)
    for (int k = 0; k < 3; k++)
      shift_lower[m] += g_ij[3*m + k] * shift[k];

  for (int i = 0; i < 5; i++) S_rate[i] = 0.0;

  double prefac = sqrt(spatial_det) * lapse;  // = α√γ = √(-g)

  // Energy density source — see source_euler for derivation.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      S_rate[4] += prefac * (T[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
      S_rate[4] += prefac * (2.0 * T[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
      S_rate[4] += prefac * (T[i + 1][j + 1] * extrinsic_curvature[i][j]);
    }
    S_rate[4] -= prefac * (T[0][0] * shift[i] * lapse_der[i]);
    S_rate[4] -= prefac * (T[0][i + 1] * lapse_der[i]);
  }

  // Momentum density sources — full Banyuls (1/2)·T^{μν}·∂_j g_{μν}
  // expansion. See source_euler for the per-term key.
  for (int j = 0; j < 3; j++) {
    S_rate[1 + j] -= prefac * (T[0][0] * lapse * lapse_der[j]);

    for (int k = 0; k < 3; k++) {
      for (int l = 0; l < 3; l++) {
        S_rate[1 + j] += prefac * (0.5 * T[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
        S_rate[1 + j] += prefac * (0.5 * T[k + 1][l + 1] * spatial_metric_der[j][k][l]);
      }
      S_rate[1 + j] += prefac * (T[0][0] * shift_lower[k] * shift_der[j][k]);

      for (int i = 0; i < 3; i++) {
        S_rate[1 + j] += prefac * (T[0][i + 1] * g_ij[3*i + k] * shift_der[j][k]);
        S_rate[1 + j] += prefac * (T[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
      }
    }
  }
}

// Per-cell forward-Euler source step for the modular GR Euler equation.
// Thin wrapper over compute_source_rate: q_new = q_old + dt·S(q_old).
void
gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
  double gas_gamma, double t_curr, double dt,
  const double *prods,
  const double fluid_old[5], double fluid_new[5])
{
  (void)t_curr;  // Source math is time-independent.

  double S_rate[5];
  compute_source_rate(gas_gamma, prods, fluid_old, S_rate);

  for (int i = 0; i < 5; i++)
    fluid_new[i] = fluid_old[i] + dt * S_rate[i];
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
  const struct gkyl_array *prods)
{
  int nfluids = st->nfluids;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(update_range, iter.idx);
    const double *prods_row = gkyl_array_cfetch(prods, cidx);

    for (int s = 0; s < nfluids; s++) {
      enum gkyl_eqn_type type = st->fluid_param[s].type;
      // The geometric source terms are identical for the curved and tetrad
      // mod variants — only the wave-structure (flat-flux + GR-correction
      // factorization) differs there, not the ADM source contributions to D,
      // S_i, τ.
      if (type != GKYL_EQN_GR_EULER_MOD && type != GKYL_EQN_GR_EULER_TETRAD_MOD)
        continue;

      double gas_gamma = st->fluid_param[s].gas_gamma;

      double *f = gkyl_array_fetch(fluid[s], cidx);

      // Per-cell setter so check_inv / repair_state can fetch the cell's
      // row in the auxfields products array. Tag the repair context as
      // "source-step" so per-call-site counters in the equation tally
      // here rather than into the wave_prop bucket.
      const struct gkyl_wv_eqn *eqn = st->eqn[s];
      bool can_repair = eqn && eqn->repair_state_func && eqn->check_inv_func;
      if (can_repair) {
        if (eqn->set_cell_idx_func) eqn->set_cell_idx_func(eqn, iter.idx);
        ((struct gkyl_wv_eqn *)eqn)->cur_repair_ctx = 0;
      }

      // Single check + cascade-repair call. The cascade walks all three
      // admissibility constraints in safe order (D, then τ, then S²) in
      // one shot, so a single repair_state call always converges to
      // admissible — no loop needed. Convex-combo updates (f_stage2,
      // final f) do NOT need their own repair pass because the
      // admissibility set is convex (sum of convex constraints, with
      // |S|_γ a norm), so a convex combination of admissible states is
      // automatically admissible. Repair only fires after the three
      // forward-Euler stages.
      #define REPAIR_ONCE(qbuf) do {                                     \
        if (can_repair && !gkyl_wv_eqn_check_inv(eqn, (qbuf)))           \
          gkyl_wv_eqn_repair_state(eqn, (qbuf));                         \
      } while (0)

      // SSP-RK3: three forward Euler stages, repair after each one. The
      // two convex combinations (f_stage2, final f) are guaranteed
      // admissible by convexity of the admissibility set, so no repair
      // pass is needed there.
      double f_old[5], f_new[5], f_stage1[5], f_stage2[5];
      for (int j = 0; j < 5; j++) f_old[j] = f[j];

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        gas_gamma, t_curr, dt, prods_row, f_old, f_new);
      REPAIR_ONCE(f_new);
      for (int j = 0; j < 5; j++) f_stage1[j] = f_new[j];

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        gas_gamma, t_curr + dt, dt, prods_row, f_stage1, f_new);
      REPAIR_ONCE(f_new);
      for (int j = 0; j < 5; j++)
        f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        gas_gamma, t_curr + 0.5 * dt, dt, prods_row, f_stage2, f_new);
      REPAIR_ONCE(f_new);
      for (int j = 0; j < 5; j++)
        f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);

      #undef REPAIR_ONCE
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
