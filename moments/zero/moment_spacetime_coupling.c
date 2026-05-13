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

// Per-cell forward-Euler source step for the modular GR Euler equation.
// Mirrors explicit_gr_euler_source_update_euler in sources_explicit.c
// (lines 767-988), substituting reads of fluid_old[5..66] with reads of the
// per-cell products row. Operates entirely on 5-element hydro state.
void
gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
  double gas_gamma, double t_curr, double dt,
  const double *prods,
  const double fluid_old[5], double fluid_new[5])
{
  (void)t_curr;  // Source math is time-independent; t_curr matches packed signature.

  double lapse   = prods[GKYL_GR_SP_LAPSE];
  double shift_x = prods[GKYL_GR_SP_SHIFT + 0];
  double shift_y = prods[GKYL_GR_SP_SHIFT + 1];
  double shift_z = prods[GKYL_GR_SP_SHIFT + 2];

  bool in_excision_region = prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);

  // Inside the excision region the source contribution is zero — just pass
  // the hydro state through.
  if (in_excision_region) {
    for (int i = 0; i < 5; i++) fluid_new[i] = fluid_old[i];
    return;
  }

  double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];

  // Primitive recovery (port of the Newton iteration in packed lines
  // 840-891). The same iteration is used by gkyl_gr_euler_mod_prim_vars,
  // but inlined here so the source step stays a single self-contained
  // function — matching packed structure for an apples-to-apples diff.
  double D    = fluid_old[0] / sqrt(spatial_det);
  double momx = fluid_old[1] / sqrt(spatial_det);
  double momy = fluid_old[2] / sqrt(spatial_det);
  double momz = fluid_old[3] / sqrt(spatial_det);
  double Etot = fluid_old[4] / sqrt(spatial_det);

  double s_sq = ((Etot + D) * (Etot + D)) - ((momx*momx) + (momy*momy) + (momz*momz));
  double C, C0;
  if (s_sq < pow(10.0, -8.0)) {
    C  = D / sqrt(pow(10.0, -8.0));
    C0 = (D + Etot) / sqrt(pow(10.0, -8.0));
  } else {
    C  = D / sqrt(s_sq);
    C0 = (D + Etot) / sqrt(s_sq);
  }

  double alpha0 = -1.0 / (gas_gamma * gas_gamma);
  double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
  double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0 -
    (C*C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
  double alpha4 = (C0*C0) - 1.0;
  double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);

  double guess = 1.0;
  int iter = 0;
  while (iter < 100) {
    double poly = (alpha4 * (guess*guess*guess) * (guess - eta))
                + (alpha2 * (guess*guess)) + (alpha1 * guess) + alpha0;
    double poly_der = alpha1 + (2.0 * alpha2 * guess)
                    + (4.0 * alpha4 * (guess*guess*guess))
                    - (3.0 * eta * alpha4 * (guess*guess));
    double guess_new = guess - (poly / poly_der);
    if (fabs(guess - guess_new) < pow(10.0, -8.0)) {
      iter = 100;
    } else {
      iter += 1;
      guess = guess_new;
    }
  }

  double W = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0) / gas_gamma) *
    ((1.0 - (C * guess)) / ((C0*C0) * (guess*guess))))));
  double h = 1.0 / (C * guess);

  double rho = D / W;
  double vx = momx / (rho * h * (W*W));
  double vy = momy / (rho * h * (W*W));
  double vz = momz / (rho * h * (W*W));
  double p = (rho * h * (W*W)) - D - Etot;

  if (rho < pow(10.0, -8.0)) rho = pow(10.0, -8.0);
  if (p   < pow(10.0, -8.0)) p   = pow(10.0, -8.0);

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

  // For non-symmetric tensors the products layout follows the packed
  // convention: ∂_j β^i is stored row-major as q[31 + 3*j + i] in packed,
  // i.e. prods[GKYL_GR_SP_DBETA + 3*j + i]. shift_der[j][i] mirrors that.
  double shift_der[3][3];
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
      shift_der[j][i] = prods[GKYL_GR_SP_DBETA + 3*j + i];

  // ∂_k γ_ij is stored row-major as q[40 + 9*k + 3*i + j] in packed.
  double spatial_metric_der[3][3][3];
  for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        spatial_metric_der[k][i][j] = prods[GKYL_GR_SP_DGIJ + 9*k + 3*i + j];

  double extrinsic_curvature[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      extrinsic_curvature[i][j] = prods[GKYL_GR_SP_KIJ + 3*i + j];

  // Three-momentum, matching packed line 945-947. Packed uses (rho + p) here
  // rather than the canonical rho*h; preserving the expression keeps mod
  // bit-equivalent to the existing packed path. If/when packed is removed
  // this should be revisited.
  double mom[3];
  mom[0] = (rho + p) * (W * W) * vx;
  mom[1] = (rho + p) * (W * W) * vy;
  mom[2] = (rho + p) * (W * W) * vz;

  // Initial copy: density and spacetime are unchanged; only momentum and
  // energy receive source contributions.
  for (int i = 0; i < 5; i++) fluid_new[i] = fluid_old[i];

  // Energy density source — port of packed lines 953-963.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fluid_new[4] += dt * (T[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
      fluid_new[4] += dt * (2.0 * T[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
      fluid_new[4] += dt * (T[i + 1][j + 1] * extrinsic_curvature[i][j]);
    }
    fluid_new[4] -= dt * (T[0][0] * shift[i] * lapse_der[i]);
    fluid_new[4] -= dt * (T[0][i + 1] * lapse_der[i]);
  }

  // Momentum density sources — port of packed lines 965-981.
  for (int j = 0; j < 3; j++) {
    fluid_new[1 + j] -= dt * (T[0][0] * lapse * lapse_der[j]);

    for (int k = 0; k < 3; k++) {
      for (int l = 0; l < 3; l++) {
        fluid_new[1 + j] += dt * (0.5 * T[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
        fluid_new[1 + j] += dt * (0.5 * T[k + 1][l + 1] * spatial_metric_der[j][k][l]);
      }
      fluid_new[1 + j] += dt * ((mom[k] / lapse) * shift_der[j][k]);

      for (int i = 0; i < 3; i++) {
        fluid_new[1 + j] += dt * (T[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
      }
    }
  }
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
      // Phase A covers only GR Euler mod; tetrad mod uses the same source
      // math (only the wave structure differs) and will be added later.
      if (type != GKYL_EQN_GR_EULER_MOD && type != GKYL_EQN_GR_EULER_TETRAD_MOD)
        continue;

      double gas_gamma = st->fluid_param[s].gas_gamma;

      double *f = gkyl_array_fetch(fluid[s], cidx);

      // SSP-RK3 mirroring explicit_gr_euler_source_update (sources_explicit.c
      // line 990). Spacetime products do not change across stages (they are
      // either static or refreshed externally between source halves).
      double f_old[5], f_new[5], f_stage1[5], f_stage2[5];
      for (int j = 0; j < 5; j++) f_old[j] = f[j];

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        gas_gamma, t_curr, dt, prods_row, f_old, f_new);
      for (int j = 0; j < 5; j++) f_stage1[j] = f_new[j];

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        gas_gamma, t_curr + dt, dt, prods_row, f_stage1, f_new);
      for (int j = 0; j < 5; j++)
        f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);

      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        gas_gamma, t_curr + 0.5 * dt, dt, prods_row, f_stage2, f_new);
      for (int j = 0; j < 5; j++)
        f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
gkyl_moment_spacetime_coupling_release(gkyl_moment_spacetime_coupling *st)
{
  gkyl_free(st);
}
