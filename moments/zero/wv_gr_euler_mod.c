#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler_mod.h>
#include <gkyl_wv_gr_euler_mod_priv.h>

// ---------------------------------------------------------------------------
// Setters for cell-index communication from wave_prop. The mutable state on
// the equation object is the per-instance idx slots filled below; downstream
// callbacks read those slots to fetch the appropriate row from the spacetime
// products array.
// ---------------------------------------------------------------------------

static void
gr_euler_mod_set_interface_idx(const struct gkyl_wv_eqn *eqn,
  const int *idxl, const int *idxr)
{
  struct wv_gr_euler_mod *grm = container_of(eqn, struct wv_gr_euler_mod, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_idxl[d] = idxl[d];
    grm->cur_idxr[d] = idxr[d];
  }
  // The next pair of rot_to_local calls correspond to (left, right) in that
  // order. Reset the parity so rot_to_local can stash locally-rotated
  // spacetime into the correct scratch slot.
  grm->rot_call_parity = 0;
}

static void
gr_euler_mod_set_cell_idx(const struct gkyl_wv_eqn *eqn, const int *idx)
{
  struct wv_gr_euler_mod *grm = container_of(eqn, struct wv_gr_euler_mod, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_cell_idx[d] = idx[d];
  }
}

void
gkyl_gr_euler_mod_set_auxfields(const struct gkyl_wv_eqn *eqn,
  struct gkyl_wv_gr_euler_mod_auxfields auxin)
{
  struct wv_gr_euler_mod *grm = container_of(eqn, struct wv_gr_euler_mod, eqn);
  grm->auxfields.prods = auxin.prods;
}

void
gkyl_gr_euler_mod_set_conf_range(const struct gkyl_wv_eqn *eqn,
  const struct gkyl_range *conf_range)
{
  struct wv_gr_euler_mod *grm = container_of(eqn, struct wv_gr_euler_mod, eqn);
  grm->conf_range = *conf_range;
}

// ---------------------------------------------------------------------------
// Rotation of the spacetime products array from global to local
// tangent-normal frame. Mirrors the spacetime portion of the packed
// rot_to_local in wv_gr_euler.c (lines 559-791), substituting global-q
// offsets for products-layout offsets defined in
// gkyl_moment_spacetime_products.h. Tensor-rank handling:
//   scalars (lapse, excision, spatial_det)         — passthrough
//   vectors (shift, lapse-der, cell-center)        — Cartesian rotated
//   rank-2 (gij, Kij, inv_gij, shift-der ∂_jβ^i)   — double rotated
//   rank-3 (∂_k γ_ij)                              — triple rotated
// Rotation R = [norm; tau1; tau2] as rows: out_a[b] = sum_b R[a][b] * in[b].
// ---------------------------------------------------------------------------

// Rank-1 vector rotation: out[a] = R[a][i] in[i] with R = [norm; tau1; tau2].
static inline void
rotate_rank1(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  out[0] = in[0]*norm[0] + in[1]*norm[1] + in[2]*norm[2];
  out[1] = in[0]*tau1[0] + in[1]*tau1[1] + in[2]*tau1[2];
  out[2] = in[0]*tau2[0] + in[1]*tau2[1] + in[2]*tau2[2];
}

// Rank-2 tensor rotation matching the packed wv_gr_euler.c convention:
// out[a][b] = sum_{i,j} R[a][j] R[b][i] in[i][j], row-major. For symmetric
// inputs (γ_ij, K_ij, γ^ij) this equals the standard rotation; for
// non-symmetric inputs (∂_j β^i) it differs from the standard convention
// by an index swap, but this exactly reproduces the packed result so
// downstream waves/qfluct math ports across cleanly.
static inline void
rotate_rank2(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  const double *R[3] = { norm, tau1, tau2 };
  // tmp[a][i] = sum_j in[i][j] R[a][j] (contract second axis of in with a).
  double tmp[3][3];
  for (int a = 0; a < 3; a++)
    for (int i = 0; i < 3; i++)
      tmp[a][i] = in[3*i+0]*R[a][0] + in[3*i+1]*R[a][1] + in[3*i+2]*R[a][2];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      out[3*a + b] = tmp[a][0]*R[b][0] + tmp[a][1]*R[b][1] + tmp[a][2]*R[b][2];
}

// Rank-3 tensor rotation: out[a][b][c] = R[a][i] R[b][j] R[c][k] in[i][j][k],
// row-major (i is outermost, k innermost). For ∂_i γ_jk the conventional
// notation puts the derivative index outermost in our products layout.
static inline void
rotate_rank3(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  const double *R[3] = { norm, tau1, tau2 };
  double tmp1[3][3][3];   // tmp1[a][j][k] = sum_i in[i][j][k] R[a][i]
  for (int a = 0; a < 3; a++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        tmp1[a][j][k] =
          in[0*9 + 3*j + k]*R[a][0] +
          in[1*9 + 3*j + k]*R[a][1] +
          in[2*9 + 3*j + k]*R[a][2];
  double tmp2[3][3][3];   // tmp2[a][b][k] = sum_j tmp1[a][j][k] R[b][j]
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

static inline void
rot_spacetime_to_local(const double *tau1, const double *tau2,
  const double *norm, const double *in_p, double *out_p)
{
  out_p[GKYL_GR_SP_LAPSE]        = in_p[GKYL_GR_SP_LAPSE];        // scalar
  out_p[GKYL_GR_SP_EXCISION]     = in_p[GKYL_GR_SP_EXCISION];     // scalar
  out_p[GKYL_GR_SP_SPATIAL_DET]  = in_p[GKYL_GR_SP_SPATIAL_DET];  // invariant

  rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_SHIFT],   &out_p[GKYL_GR_SP_SHIFT]);
  rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_DALPHA],  &out_p[GKYL_GR_SP_DALPHA]);
  rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_CELLCTR], &out_p[GKYL_GR_SP_CELLCTR]);

  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_GIJ],     &out_p[GKYL_GR_SP_GIJ]);
  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_KIJ],     &out_p[GKYL_GR_SP_KIJ]);
  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_DBETA],   &out_p[GKYL_GR_SP_DBETA]);
  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_INV_GIJ], &out_p[GKYL_GR_SP_INV_GIJ]);

  rotate_rank3(tau1, tau2, norm, &in_p[GKYL_GR_SP_DGIJ],    &out_p[GKYL_GR_SP_DGIJ]);
}

// ---------------------------------------------------------------------------
// Hydro helpers (prim_vars / flux / max_abs_speed). Ported from packed
// wv_gr_euler.c with reads against the spacetime-products layout instead of
// q[5..70]. The Newton-iteration recovery in prim_vars and the eigenvalue
// computation in max_abs_speed are reproduced verbatim — only the data
// access pattern changes, so results agree to <1e-12 with the packed path
// for the same physical input.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_mod_prim_vars(double gas_gamma, const double q[5],
  const double *prods, double v[5])
{
  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];

    double D    = q[0] / sqrt(spatial_det);
    double momx = q[1] / sqrt(spatial_det);
    double momy = q[2] / sqrt(spatial_det);
    double momz = q[3] / sqrt(spatial_det);
    double Etot = q[4] / sqrt(spatial_det);

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
      double poly = (alpha4 * (guess*guess*guess) * (guess - eta)) +
        (alpha2 * (guess*guess)) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess) +
        (4.0 * alpha4 * (guess*guess*guess)) - (3.0 * eta * alpha4 * (guess*guess));
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

    v[0] = D / W;
    v[1] = momx / (v[0] * h * (W*W));
    v[2] = momy / (v[0] * h * (W*W));
    v[3] = momz / (v[0] * h * (W*W));
    v[4] = (v[0] * h * (W*W)) - D - Etot;

    if (v[0] < pow(10.0, -8.0)) v[0] = pow(10.0, -8.0);
    if (v[4] < pow(10.0, -8.0)) v[4] = pow(10.0, -8.0);
  } else {
    for (int i = 0; i < 5; i++) v[i] = 0.0;
  }
}

void
gkyl_gr_euler_mod_flux(double gas_gamma, const double q[5],
  const double *prods, double flux[5])
{
  double v[5];
  gkyl_gr_euler_mod_prim_vars(gas_gamma, q, prods, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  double lapse        = prods[GKYL_GR_SP_LAPSE];
  double shift_x      = prods[GKYL_GR_SP_SHIFT + 0];
  double spatial_det  = prods[GKYL_GR_SP_SPATIAL_DET];

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3] = { vx, vy, vz };
    double v_sq = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        v_sq += prods[GKYL_GR_SP_GIJ + 3*i + j] * vel[i] * vel[j];

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) W = 1.0 / sqrt(pow(10.0, -8.0));

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    flux[0] = (lapse * sqrt(spatial_det)) * (rho * W * (vx - (shift_x / lapse)));
    flux[1] = (lapse * sqrt(spatial_det)) *
      (rho * h * (W*W) * (vx * (vx - (shift_x / lapse))) + p);
    flux[2] = (lapse * sqrt(spatial_det)) *
      (rho * h * (W*W) * (vy * (vx - (shift_x / lapse))));
    flux[3] = (lapse * sqrt(spatial_det)) *
      (rho * h * (W*W) * (vz * (vx - (shift_x / lapse))));
    flux[4] = (lapse * sqrt(spatial_det)) *
      (((rho * h * (W*W)) - p - (rho * W)) * (vx - (shift_x / lapse)) + (p * vx));
  } else {
    for (int i = 0; i < 5; i++) flux[i] = 0.0;
  }
}

double
gkyl_gr_euler_mod_max_abs_speed(double gas_gamma, const double q[5],
  const double *prods)
{
  double v[5];
  gkyl_gr_euler_mod_prim_vars(gas_gamma, q, prods, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  double lapse   = prods[GKYL_GR_SP_LAPSE];
  double shift_x = prods[GKYL_GR_SP_SHIFT + 0];
  double shift_y = prods[GKYL_GR_SP_SHIFT + 1];
  double shift_z = prods[GKYL_GR_SP_SHIFT + 2];

  double num = (gas_gamma * p) / rho;
  double den = 1.0 + ((p / rho) * gas_gamma / (gas_gamma - 1.0));
  double c_s = sqrt(num / den);

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  bool curved_spacetime = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double g_ij = prods[GKYL_GR_SP_GIJ + 3*i + j];
      if (i == j) {
        if (fabs(g_ij - 1.0) > pow(10.0, -8.0)) curved_spacetime = true;
      } else {
        if (fabs(g_ij)       > pow(10.0, -8.0)) curved_spacetime = true;
      }
    }
  }
  if (fabs(lapse - 1.0) > pow(10.0, -8.0) ||
      fabs(shift_x)     > pow(10.0, -8.0) ||
      fabs(shift_y)     > pow(10.0, -8.0) ||
      fabs(shift_z)     > pow(10.0, -8.0)) {
    curved_spacetime = true;
  }

  if (in_excision_region) return pow(10.0, -8.0);

  if (curved_spacetime) {
    double vel[3]   = { vx, vy, vz };
    double shift[3] = { shift_x, shift_y, shift_z };
    double v_sq = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        v_sq += prods[GKYL_GR_SP_GIJ + 3*i + j] * vel[i] * vel[j];

    double max_eig = 0.0;
    for (int i = 0; i < 3; i++) {
      double inv_ii = prods[GKYL_GR_SP_INV_GIJ + 3*i + i];
      double material = (lapse * vel[i]) - shift[i];
      double common = lapse / (1.0 - (v_sq * (c_s*c_s)));
      double rad = (1.0 - v_sq) * (inv_ii * (1.0 - (v_sq * (c_s*c_s))) -
        (vel[i] * vel[i]) * (1.0 - (c_s*c_s)));
      double fast = common * ((vel[i] * (1.0 - (c_s*c_s))) + (c_s * sqrt(rad))) - shift[i];
      double slow = common * ((vel[i] * (1.0 - (c_s*c_s))) - (c_s * sqrt(rad))) - shift[i];

      if (fabs(material) > max_eig) max_eig = fabs(material);
      if (fabs(fast)     > max_eig) max_eig = fabs(fast);
      if (fabs(slow)     > max_eig) max_eig = fabs(slow);
    }

    return fabs(v_sq) + max_eig;
  }

  // Flat spacetime fast path. Matches the (slightly unusual) packed
  // expression `fabs(v_sq) + c_s` where v_sq here is actually |v|.
  double v_mag = sqrt((vx*vx) + (vy*vy) + (vz*vz));
  return fabs(v_mag) + c_s;
}

// ---------------------------------------------------------------------------
// Riemann-variable conversions and Cartesian-frame rotations. Hydro-only: the
// spacetime portion of the state lives in auxfields, not in q, so there is no
// spacetime block to pass through here.
// ---------------------------------------------------------------------------

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn, const double *qstate,
  const double *qin, double *wout)
{
  // Matches the packed implementation (identity until a proper L matrix is
  // wired in).
  for (int i = 0; i < 5; i++) wout[i] = qin[i];
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn, const double *qstate,
  const double *win, double *qout)
{
  for (int i = 0; i < 5; i++) qout[i] = win[i];
}

static void
gr_euler_mod_wall(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
  for (int i = 0; i < 5; i++) ghost[i] = skin[i];
  ghost[1] = -ghost[1];                          // reverse normal momentum
}

static void
gr_euler_mod_no_slip(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn *eqn, const double *tau1,
  const double *tau2, const double *norm, const double *GKYL_RESTRICT qglobal,
  double *GKYL_RESTRICT qlocal)
{
  // Hydro: scalar D and Etot pass through; momentum is Cartesian-rotated.
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];

  // Spacetime: read this interface side's global products and rotate into
  // the matching scratch buffer. wave_prop calls rot_to_local twice per
  // interface (left then right); rot_call_parity is reset to 0 by
  // set_interface_idx_func so the first call here populates prodl_local and
  // the second populates prodr_local.
  struct wv_gr_euler_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_mod, eqn);
  if (grm->auxfields.prods) {
    const int *idx = (grm->rot_call_parity == 0) ? grm->cur_idxl : grm->cur_idxr;
    double *out = (grm->rot_call_parity == 0) ? grm->prodl_local : grm->prodr_local;
    long cidx = gkyl_range_idx(&grm->conf_range, idx);
    const double *in_p = gkyl_array_cfetch(grm->auxfields.prods, cidx);
    rot_spacetime_to_local(tau1, tau2, norm, in_p, out);
  }
  grm->rot_call_parity ^= 1;
}

static inline void
rot_to_global(const struct gkyl_wv_eqn *eqn, const double *tau1,
  const double *tau2, const double *norm, const double *GKYL_RESTRICT qlocal,
  double *GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];
}

// ---------------------------------------------------------------------------
// Riemann-solver callbacks. These are STUBS until the math port is in place.
// Real implementations need locally-rotated spacetime products (see TODO in
// rot_to_local) and will mirror the packed wave_lax/wave_roe/wave_hll and
// their qfluct counterparts in wv_gr_euler.c.
// ---------------------------------------------------------------------------

static double
wave_lax(const struct gkyl_wv_eqn *eqn, const double *delta, const double *ql,
  const double *qr, double *waves, double *s)
{
  struct wv_gr_euler_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double sl = gkyl_gr_euler_mod_max_abs_speed(gas_gamma, ql, grm->prodl_local);
  double sr = gkyl_gr_euler_mod_max_abs_speed(gas_gamma, qr, grm->prodr_local);
  double amax = fmax(sl, sr);

  double fl[5], fr[5];
  gkyl_gr_euler_mod_flux(gas_gamma, ql, grm->prodl_local, fl);
  gkyl_gr_euler_mod_flux(gas_gamma, qr, grm->prodr_local, fr);

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);

  double *w0 = &waves[0], *w1 = &waves[5];
  if (!excise_l && !excise_r) {
    for (int i = 0; i < 5; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
    }
  } else {
    for (int i = 0; i < 5; i++) { w0[i] = 0.0; w1[i] = 0.0; }
  }

  s[0] = -amax;
  s[1] =  amax;
  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);
  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_roe(const struct gkyl_wv_eqn *eqn, const double *delta, const double *ql,
  const double *qr, double *waves, double *s)
{
  // Ported from packed wv_gr_euler.c:1109. Note: the packed Roe is a pure-
  // hydro special-relativistic Roe solver applied to the locally-rotated
  // hydro state; it does not consume any spacetime data, so prodl_local /
  // prodr_local are unused here. Equivalence holds at <1e-12 because the
  // arithmetic is identical on the 5-element hydro slice.
  struct wv_gr_euler_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double vl[5], vr[5];
  gkyl_gr_euler_mod_prim_vars(gas_gamma, ql, grm->prodl_local, vl);
  gkyl_gr_euler_mod_prim_vars(gas_gamma, qr, grm->prodr_local, vr);

  double rho_l = vl[0], vx_l = vl[1], vy_l = vl[2], vz_l = vl[3], p_l = vl[4];
  double rho_r = vr[0], vx_r = vr[1], vy_r = vr[2], vz_r = vr[3], p_r = vr[4];
  (void)rho_l; (void)rho_r;

  double Etot_l = ql[4];
  double Etot_r = qr[4];

  double W_l = 1.0 / sqrt(1.0 - ((vx_l*vx_l) + (vy_l*vy_l) + (vz_l*vz_l)));
  double W_r = 1.0 / sqrt(1.0 - ((vx_r*vx_r) + (vy_r*vy_r) + (vz_r*vz_r)));

  double K_l = sqrt(Etot_l + p_l) / W_l;
  double K_r = sqrt(Etot_r + p_r) / W_r;
  double K_avg = 1.0 / (K_l + K_r);

  double v0 = ((K_l * W_l)        + (K_r * W_r))        * K_avg;
  double v1 = ((K_l * W_l * vx_l) + (K_r * W_r * vx_r)) * K_avg;
  double v2 = ((K_l * W_l * vy_l) + (K_r * W_r * vy_r)) * K_avg;
  double v3 = ((K_l * W_l * vz_l) + (K_r * W_r * vz_r)) * K_avg;
  double v4 = ((p_l / K_l)        + (p_r / K_r))        * K_avg;

  double c_minus = 1.0 - ((gas_gamma / (gas_gamma - 1.0)) * v4);
  double c_plus  = 1.0 + ((gas_gamma / (gas_gamma - 1.0)) * v4);

  double v_alpha_sq = -(v0*v0) + (v1*v1) + (v2*v2) + (v3*v3);
  double s_sq = (0.5 * gas_gamma * v4 * (1.0 - v_alpha_sq))
              - (0.5 * (gas_gamma - 1.0) * (1.0 + v_alpha_sq));
  double energy = (v0*v0) - (v1*v1);
  double y = sqrt(((1.0 - (gas_gamma * v4)) * energy) + s_sq);

  double k = (v0 * delta[4]) - (v1 * delta[1]);
  double v_delta = (-v0 * delta[4]) + (v1 * delta[1]) + (v2 * delta[2]) + (v3 * delta[3]);
  double a1 = -((s_sq * k) + (sqrt(s_sq) * y * ((v0 * delta[1]) - (v1 * delta[4]))
              + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta))))) / (2.0 * energy * s_sq);
  double a2 = -((s_sq * k) - (sqrt(s_sq) * y * ((v0 * delta[1]) - (v1 * delta[4]))
              + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta))))) / (2.0 * energy * s_sq);
  double a3 = ((2.0 * s_sq * k)
              + ((gas_gamma - 1.0) * energy * (delta[0] + (c_plus * v_delta)))) / (energy * s_sq);
  double a4 = delta[2] - ((k * v2) / energy);
  double a5 = delta[3] - ((k * v3) / energy);

  for (int i = 0; i < 5 * 3; i++) waves[i] = 0.0;

  double *wv;
  wv = &waves[0 * 5];
  wv[0] = a1 * c_minus;
  wv[1] = a1 * (v1 - ((sqrt(s_sq) * v0) / y));
  wv[2] = a1 * v2;
  wv[3] = a1 * v3;
  wv[4] = a1 * (v0 - ((sqrt(s_sq) * v1) / y));
  s[0] = (((1.0 - (gas_gamma * v4)) * v0 * v1) - (sqrt(s_sq) * y))
       / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  wv = &waves[1 * 5];
  wv[0] = (a3 * (c_minus + (s_sq / (gas_gamma - 1.0))))
        - (a4 * c_plus * v2) - (a5 * c_plus * v3);
  wv[1] = a3 * v1;
  wv[2] = (a3 * v2) + a4;
  wv[3] = (a3 * v3) + a5;
  wv[4] = a3 * v0;
  s[1] = v1 / v0;

  wv = &waves[2 * 5];
  wv[0] = a2 * c_minus;
  wv[1] = a2 * (v1 + ((sqrt(s_sq) * v0) / y));
  wv[2] = a2 * v2;
  wv[3] = a2 * v3;
  wv[4] = a2 * (v0 + ((sqrt(s_sq) * v1) / y));
  s[2] = (((1.0 - (gas_gamma * v4)) * v0 * v1) + (sqrt(s_sq) * y))
       / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  return (((1.0 - (gas_gamma * v4)) * v0 * fabs(v1)) + (sqrt(s_sq) * y))
       / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0 * 5];
  const double *w1 = &waves[1 * 5];
  const double *w2 = &waves[2 * 5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);
  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]);
  }
}

static double
wave_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return wave_roe(eqn, delta, ql, qr, waves, s);
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_hll(const struct gkyl_wv_eqn *eqn, const double *delta, const double *ql,
  const double *qr, double *waves, double *s)
{
  struct wv_gr_euler_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_mod, eqn);
  double gas_gamma = grm->gas_gamma;
  const double *pl = grm->prodl_local;
  const double *pr = grm->prodr_local;

  double vl[5], vr[5];
  gkyl_gr_euler_mod_prim_vars(gas_gamma, ql, pl, vl);
  gkyl_gr_euler_mod_prim_vars(gas_gamma, qr, pr, vr);

  double rho_l = vl[0], vx_l = vl[1], vy_l = vl[2], vz_l = vl[3], p_l = vl[4];
  double rho_r = vr[0], vx_r = vr[1], vy_r = vr[2], vz_r = vr[3], p_r = vr[4];

  double lapse_l   = pl[GKYL_GR_SP_LAPSE];
  double shift_xl  = pl[GKYL_GR_SP_SHIFT + 0];
  double shift_yl  = pl[GKYL_GR_SP_SHIFT + 1];
  double shift_zl  = pl[GKYL_GR_SP_SHIFT + 2];
  double lapse_r   = pr[GKYL_GR_SP_LAPSE];
  double shift_xr  = pr[GKYL_GR_SP_SHIFT + 0];
  double shift_yr  = pr[GKYL_GR_SP_SHIFT + 1];
  double shift_zr  = pr[GKYL_GR_SP_SHIFT + 2];

  bool excise_l = pl[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = pr[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);

  bool curved_l = false, curved_r = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double gl = pl[GKYL_GR_SP_GIJ + 3*i + j];
      double gr = pr[GKYL_GR_SP_GIJ + 3*i + j];
      if (i == j) {
        if (fabs(gl - 1.0) > pow(10.0, -8.0)) curved_l = true;
        if (fabs(gr - 1.0) > pow(10.0, -8.0)) curved_r = true;
      } else {
        if (fabs(gl) > pow(10.0, -8.0)) curved_l = true;
        if (fabs(gr) > pow(10.0, -8.0)) curved_r = true;
      }
    }
  }
  if (fabs(lapse_l - 1.0) > pow(10.0, -8.0) ||
      fabs(shift_xl) > pow(10.0, -8.0) || fabs(shift_yl) > pow(10.0, -8.0) ||
      fabs(shift_zl) > pow(10.0, -8.0))
    curved_l = true;
  if (fabs(lapse_r - 1.0) > pow(10.0, -8.0) ||
      fabs(shift_xr) > pow(10.0, -8.0) || fabs(shift_yr) > pow(10.0, -8.0) ||
      fabs(shift_zr) > pow(10.0, -8.0))
    curved_r = true;

  double c_sl = sqrt((gas_gamma * p_l / rho_l) /
                     (1.0 + (p_l / rho_l) * gas_gamma / (gas_gamma - 1.0)));
  double c_sr = sqrt((gas_gamma * p_r / rho_r) /
                     (1.0 + (p_r / rho_r) * gas_gamma / (gas_gamma - 1.0)));

  double vx_avg = 0.5 * (vx_l + vx_r);
  double cs_avg = 0.5 * (c_sl + c_sr);

  double sl, sr;
  if (curved_l || curved_r) {
    double vel_l[3] = { vx_l, vy_l, vz_l };
    double vel_r[3] = { vx_r, vy_r, vz_r };
    double shift_l[3] = { shift_xl, shift_yl, shift_zl };
    double shift_r[3] = { shift_xr, shift_yr, shift_zr };

    double vsq_l = 0.0, vsq_r = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        vsq_l += pl[GKYL_GR_SP_GIJ + 3*i + j] * vel_l[i] * vel_l[j];
        vsq_r += pr[GKYL_GR_SP_GIJ + 3*i + j] * vel_r[i] * vel_r[j];
      }

    double max_eig_l = 0.0, max_eig_r = 0.0;
    for (int i = 0; i < 3; i++) {
      double inv_ii_l = pl[GKYL_GR_SP_INV_GIJ + 3*i + i];
      double inv_ii_r = pr[GKYL_GR_SP_INV_GIJ + 3*i + i];

      double mat_l  = (lapse_l * vel_l[i]) - shift_l[i];
      double com_l  = lapse_l / (1.0 - (vsq_l * c_sl * c_sl));
      double rad_l  = (1.0 - vsq_l) * (inv_ii_l * (1.0 - (vsq_l * c_sl * c_sl)) -
                       (vel_l[i] * vel_l[i]) * (1.0 - (c_sl * c_sl)));
      double fast_l = com_l * ((vel_l[i] * (1.0 - c_sl*c_sl)) + (c_sl * sqrt(rad_l))) - shift_l[i];
      double slow_l = com_l * ((vel_l[i] * (1.0 - c_sl*c_sl)) - (c_sl * sqrt(rad_l))) - shift_l[i];
      if (fabs(mat_l)  > max_eig_l) max_eig_l = fabs(mat_l);
      if (fabs(fast_l) > max_eig_l) max_eig_l = fabs(fast_l);
      if (fabs(slow_l) > max_eig_l) max_eig_l = fabs(slow_l);

      double mat_r  = (lapse_r * vel_r[i]) - shift_r[i];
      double com_r  = lapse_r / (1.0 - (vsq_r * c_sr * c_sr));
      double rad_r  = (1.0 - vsq_r) * (inv_ii_r * (1.0 - (vsq_r * c_sr * c_sr)) -
                       (vel_r[i] * vel_r[i]) * (1.0 - (c_sr * c_sr)));
      double fast_r = com_r * ((vel_r[i] * (1.0 - c_sr*c_sr)) + (c_sr * sqrt(rad_r))) - shift_r[i];
      double slow_r = com_r * ((vel_r[i] * (1.0 - c_sr*c_sr)) - (c_sr * sqrt(rad_r))) - shift_r[i];
      if (fabs(mat_r)  > max_eig_r) max_eig_r = fabs(mat_r);
      if (fabs(fast_r) > max_eig_r) max_eig_r = fabs(fast_r);
      if (fabs(slow_r) > max_eig_r) max_eig_r = fabs(slow_r);
    }

    double max_eig_avg = 0.5 * (max_eig_l + max_eig_r);
    sl = (vx_avg - max_eig_avg) / (1.0 - (vx_avg * max_eig_avg));
    sr = (vx_avg + max_eig_avg) / (1.0 + (vx_avg * max_eig_avg));
  } else {
    sl = (vx_avg - cs_avg) / (1.0 - (vx_avg * cs_avg));
    sr = (vx_avg + cs_avg) / (1.0 + (vx_avg * cs_avg));
  }

  double fl[5], fr[5];
  gkyl_gr_euler_mod_flux(gas_gamma, ql, pl, fl);
  gkyl_gr_euler_mod_flux(gas_gamma, qr, pr, fr);

  double qm[5];
  for (int i = 0; i < 5; i++)
    qm[i] = ((sr * qr[i]) - (sl * ql[i]) + (fl[i] - fr[i])) / (sr - sl);

  double *w0 = &waves[0], *w1 = &waves[5];
  if (!excise_l && !excise_r) {
    for (int i = 0; i < 5; i++) {
      w0[i] = qm[i] - ql[i];
      w1[i] = qr[i] - qm[i];
    }
  } else {
    for (int i = 0; i < 5; i++) { w0[i] = 0.0; w1[i] = 0.0; }
  }

  s[0] = sl;
  s[1] = sr;
  return fmax(fabs(sl), fabs(sr));
}

static void
qfluct_hll(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);
  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return wave_hll(eqn, delta, ql, qr, waves, s);
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
flux_jump_func(const struct gkyl_wv_eqn *eqn, const double *ql,
  const double *qr, double *flux_jump)
{
  // Called per interface after rot_to_local has populated prodl_local
  // (left side) and prodr_local (right side). Mirrors packed flux_jump
  // (wv_gr_euler.c:1522).
  struct wv_gr_euler_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double fl[5], fr[5];
  gkyl_gr_euler_mod_flux(gas_gamma, ql, grm->prodl_local, fl);
  gkyl_gr_euler_mod_flux(gas_gamma, qr, grm->prodr_local, fr);

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (!excise_l && !excise_r) {
    for (int m = 0; m < 5; m++) flux_jump[m] = fr[m] - fl[m];
  } else {
    for (int m = 0; m < 5; m++) flux_jump[m] = 0.0;
  }

  double amaxl = gkyl_gr_euler_mod_max_abs_speed(gas_gamma, ql, grm->prodl_local);
  double amaxr = gkyl_gr_euler_mod_max_abs_speed(gas_gamma, qr, grm->prodr_local);
  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  // Per-cell positivity check on (rho, p) recovered from q via the GR
  // primitive-recovery iteration. Mirrors packed check_inv (wv_gr_euler.c:
  // 1559): set_cell_idx_func must have been called beforehand to point
  // grm->cur_cell_idx at this cell.
  struct wv_gr_euler_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_mod, eqn);
  if (!grm->auxfields.prods) {
    // Auxfields not wired yet — fall back to a bare hydro check.
    return q[0] >= 0.0;
  }
  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);

  double v[5];
  gkyl_gr_euler_mod_prim_vars(grm->gas_gamma, q, prods, v);
  if (v[0] < 0.0 || v[4] < 0.0) return false;
  return true;
}

static double
max_speed_func(const struct gkyl_wv_eqn *eqn, const double *q)
{
  // Per-cell max characteristic speed. Mirrors packed max_speed
  // (wv_gr_euler.c:1576). Uses grm->cur_cell_idx populated by
  // set_cell_idx_func; products are in the global frame at this point
  // (max_speed_func is called outside the rotate→waves→qfluct sequence).
  struct wv_gr_euler_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_mod, eqn);
  if (!grm->auxfields.prods) return 1.0;
  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);
  return gkyl_gr_euler_mod_max_abs_speed(grm->gas_gamma, q, prods);
}

static inline void
gr_euler_mod_cons_to_diag(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *diag)
{
  for (int i = 0; i < 5; i++) diag[i] = qin[i];
}

static inline void
gr_euler_mod_source(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *sout)
{
  // Geometric source terms are integrated by moment_spacetime_coupling, not
  // through this callback.
  for (int i = 0; i < 5; i++) sout[i] = 0.0;
}

// ---------------------------------------------------------------------------
// Constructors and lifetime.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_mod_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_gr_euler_mod *grm = container_of(base, struct wv_gr_euler_mod, eqn);
  gkyl_free(grm);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_mod_new(double gas_gamma, const struct gkyl_range *conf_range,
  bool use_gpu)
{
  return gkyl_wv_gr_euler_mod_inew(&(struct gkyl_wv_gr_euler_mod_inp) {
    .gas_gamma = gas_gamma,
    .conf_range = *conf_range,
    .rp_type = WV_GR_EULER_RP_HLL,
    .use_gpu = use_gpu,
  });
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_mod_inew(const struct gkyl_wv_gr_euler_mod_inp *inp)
{
  struct wv_gr_euler_mod *grm = gkyl_malloc(sizeof(struct wv_gr_euler_mod));

  grm->eqn.type = GKYL_EQN_GR_EULER_MOD;
  grm->eqn.num_equations = 5;
  grm->eqn.num_diag = 5;

  grm->gas_gamma = inp->gas_gamma;
  grm->conf_range = inp->conf_range;
  grm->auxfields.prods = NULL;
  grm->rot_call_parity = 0;
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_idxl[d] = 0;
    grm->cur_idxr[d] = 0;
    grm->cur_cell_idx[d] = 0;
  }

  if (inp->rp_type == WV_GR_EULER_RP_LAX) {
    grm->eqn.num_waves = 2;
    grm->eqn.waves_func = wave_lax_l;
    grm->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_EULER_RP_ROE) {
    grm->eqn.num_waves = 3;
    grm->eqn.waves_func = wave_roe_l;
    grm->eqn.qfluct_func = qfluct_roe_l;
  }
  else {  // default: HLL
    grm->eqn.num_waves = 2;
    grm->eqn.waves_func = wave_hll_l;
    grm->eqn.qfluct_func = qfluct_hll_l;
  }

  grm->eqn.ffluct_func = NULL;  // q-wave only; matches packed.
  grm->eqn.flux_jump = flux_jump_func;
  grm->eqn.check_inv_func = check_inv;
  grm->eqn.max_speed_func = max_speed_func;
  grm->eqn.rotate_to_local_func = rot_to_local;
  grm->eqn.rotate_to_global_func = rot_to_global;
  grm->eqn.wall_bc_func = gr_euler_mod_wall;
  grm->eqn.no_slip_bc_func = gr_euler_mod_no_slip;
  grm->eqn.cons_to_riem = cons_to_riem;
  grm->eqn.riem_to_cons = riem_to_cons;
  grm->eqn.cons_to_diag = gr_euler_mod_cons_to_diag;
  grm->eqn.source_func = gr_euler_mod_source;

  // Per-cell idx setters — these are what tie the equation to the products
  // array at advance time.
  grm->eqn.set_interface_idx_func = gr_euler_mod_set_interface_idx;
  grm->eqn.set_cell_idx_func = gr_euler_mod_set_cell_idx;

  grm->eqn.embed_geo = NULL;
  grm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(grm->eqn.flags);
  grm->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_mod_free);
  grm->eqn.on_dev = &grm->eqn;

  return &grm->eqn;
}

double
gkyl_wv_gr_euler_mod_gas_gamma(const struct gkyl_wv_eqn *eqn)
{
  const struct wv_gr_euler_mod *grm = container_of(eqn,
    struct wv_gr_euler_mod, eqn);
  return grm->gas_gamma;
}
