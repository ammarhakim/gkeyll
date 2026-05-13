#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler_tetrad_mod.h>
#include <gkyl_wv_gr_euler_tetrad_mod_priv.h>

// ---------------------------------------------------------------------------
// Setters for cell-index communication from wave_prop. Identical to the
// regular mod variant; the per-instance idx slots populated below are read
// by downstream callbacks to fetch the appropriate row from the spacetime
// products array.
// ---------------------------------------------------------------------------

static void
gr_euler_tetrad_mod_set_interface_idx(const struct gkyl_wv_eqn *eqn,
  const int *idxl, const int *idxr)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_idxl[d] = idxl[d];
    grm->cur_idxr[d] = idxr[d];
  }
  grm->rot_call_parity = 0;
}

static void
gr_euler_tetrad_mod_set_cell_idx(const struct gkyl_wv_eqn *eqn, const int *idx)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_cell_idx[d] = idx[d];
  }
}

void
gkyl_gr_euler_tetrad_mod_set_auxfields(const struct gkyl_wv_eqn *eqn,
  struct gkyl_wv_gr_euler_tetrad_mod_auxfields auxin)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  grm->auxfields.prods = auxin.prods;
}

void
gkyl_gr_euler_tetrad_mod_set_conf_range(const struct gkyl_wv_eqn *eqn,
  const struct gkyl_range *conf_range)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  grm->conf_range = *conf_range;
}

// ---------------------------------------------------------------------------
// Rotation of the spacetime products array from global to local
// tangent-normal frame. Shared verbatim with the regular mod variant — the
// tetrad equation reads the same spacetime block; rotation conventions are
// identical.
// ---------------------------------------------------------------------------

static inline void
rotate_rank1(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  out[0] = in[0]*norm[0] + in[1]*norm[1] + in[2]*norm[2];
  out[1] = in[0]*tau1[0] + in[1]*tau1[1] + in[2]*tau1[2];
  out[2] = in[0]*tau2[0] + in[1]*tau2[1] + in[2]*tau2[2];
}

static inline void
rotate_rank2(const double *tau1, const double *tau2, const double *norm,
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

static inline void
rotate_rank3(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
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

static inline void
rot_spacetime_to_local(const double *tau1, const double *tau2,
  const double *norm, const double *in_p, double *out_p)
{
  out_p[GKYL_GR_SP_LAPSE]        = in_p[GKYL_GR_SP_LAPSE];
  out_p[GKYL_GR_SP_EXCISION]     = in_p[GKYL_GR_SP_EXCISION];
  out_p[GKYL_GR_SP_SPATIAL_DET]  = in_p[GKYL_GR_SP_SPATIAL_DET];

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
// Hydro helpers (prim_vars / flux / flux_correction / max_abs_speed).
// prim_vars and max_abs_speed match the regular mod variant — they don't
// see the tetrad split. The flux is what differs:
//   gkyl_gr_euler_tetrad_mod_flux:            flat-space SR flux (Cartesian W,
//                                             no lapse/shift/√γ)
//   gkyl_gr_euler_tetrad_mod_flux_correction: applies α·√γ and replaces
//                                             (vx, W_flat) with (vx - βˣ/α,
//                                             W_curved) to recover the
//                                             curved-space flux.
// Mirrors the packed factorization in wv_gr_euler_tetrad.c at lines 10-124.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_tetrad_mod_prim_vars(double gas_gamma, const double q[5],
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
gkyl_gr_euler_tetrad_mod_flux(double gas_gamma, const double q[5],
  const double *prods, double flux_sr[5])
{
  double v[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, q, prods, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double v_dot = (vx*vx) + (vy*vy) + (vz*vz);
    double W = 1.0 / sqrt(1.0 - v_dot);
    if (v_dot > 1.0 - pow(10.0, -8.0)) W = 1.0 / sqrt(pow(10.0, -8.0));

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    flux_sr[0] = rho * W * vx;
    flux_sr[1] = (rho * h * (W*W) * (vx * vx)) + p;
    flux_sr[2] = rho * h * (W*W) * (vy * vx);
    flux_sr[3] = rho * h * (W*W) * (vz * vx);
    flux_sr[4] = ((rho * h * (W*W)) - (rho * W)) * vx;
  } else {
    for (int i = 0; i < 5; i++) flux_sr[i] = 0.0;
  }
}

void
gkyl_gr_euler_tetrad_mod_flux_correction(double gas_gamma, const double q[5],
  const double *prods, const double flux_sr[5], double flux_gr[5])
{
  double v[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, q, prods, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];
  (void)vy; (void)vz;  // used via vel[] below

  double lapse        = prods[GKYL_GR_SP_LAPSE];
  double shift_x      = prods[GKYL_GR_SP_SHIFT + 0];
  double spatial_det  = prods[GKYL_GR_SP_SPATIAL_DET];

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel[3] = { v[1], v[2], v[3] };
    double v_sq = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        v_sq += prods[GKYL_GR_SP_GIJ + 3*i + j] * vel[i] * vel[j];

    double v_dot = (v[1]*v[1]) + (v[2]*v[2]) + (v[3]*v[3]);
    double W_flat = 1.0 / sqrt(1.0 - v_dot);
    if (v_dot > 1.0 - pow(10.0, -8.0)) W_flat = 1.0 / sqrt(pow(10.0, -8.0));

    double W_curved = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) W_curved = 1.0 / sqrt(pow(10.0, -8.0));

    // Numerical-safety clamp on vx for the divisions below. Mirrors packed.
    if (fabs(vx) < pow(10.0, -8.0)) {
      vx = (vx > 0.0) ? pow(10.0, -8.0) : -pow(10.0, -8.0);
    }

    double prefac = lapse * sqrt(spatial_det);
    double vmsh   = vx - (shift_x / lapse);

    flux_gr[0] = prefac * ((flux_sr[0] * vmsh * W_curved) / (vx * W_flat));
    flux_gr[1] = prefac * ((((flux_sr[1] - p) * vmsh * (W_curved*W_curved))
                         / (vx * (W_flat*W_flat))) + p);
    flux_gr[2] = prefac * ((flux_sr[2] * vmsh * (W_curved*W_curved))
                         / (vx * (W_flat*W_flat)));
    flux_gr[3] = prefac * ((flux_sr[3] * vmsh * (W_curved*W_curved))
                         / (vx * (W_flat*W_flat)));
    flux_gr[4] = prefac * (((((flux_sr[4] + (rho * vx * W_flat)) * (W_curved*W_curved))
                          / (vx * (W_flat*W_flat))) - p - (rho * W_curved)) * vmsh + (p * vx));
  } else {
    for (int i = 0; i < 5; i++) flux_gr[i] = 0.0;
  }
}

double
gkyl_gr_euler_tetrad_mod_max_abs_speed(double gas_gamma, const double q[5],
  const double *prods)
{
  double v[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, q, prods, v);
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

  double v_mag = sqrt((vx*vx) + (vy*vy) + (vz*vz));
  return fabs(v_mag) + c_s;
}

// ---------------------------------------------------------------------------
// Riemann-variable conversions and Cartesian-frame rotations.
// ---------------------------------------------------------------------------

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn, const double *qstate,
  const double *qin, double *wout)
{
  for (int i = 0; i < 5; i++) wout[i] = qin[i];
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn, const double *qstate,
  const double *win, double *qout)
{
  for (int i = 0; i < 5; i++) qout[i] = win[i];
}

static void
gr_euler_tetrad_mod_wall(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
  for (int i = 0; i < 5; i++) ghost[i] = skin[i];
  ghost[1] = -ghost[1];
}

static void
gr_euler_tetrad_mod_no_slip(const struct gkyl_wv_eqn *eqn, double t, int nc,
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
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];

  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
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
// Riemann-solver callbacks. Wave-speed estimates and Roe eigenstructure are
// shared with the regular mod variant; the tetrad-specific piece is the
// flat-flux + correction step used inside wave_lax / wave_hll / flux_jump.
// ---------------------------------------------------------------------------

static double
wave_lax(const struct gkyl_wv_eqn *eqn, const double *delta, const double *ql,
  const double *qr, double *waves, double *s)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double sl = gkyl_gr_euler_tetrad_mod_max_abs_speed(gas_gamma, ql, grm->prodl_local);
  double sr = gkyl_gr_euler_tetrad_mod_max_abs_speed(gas_gamma, qr, grm->prodr_local);
  double amax = fmax(sl, sr);

  double fl_sr[5], fr_sr[5];
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, ql, grm->prodl_local, fl_sr);
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qr, grm->prodr_local, fr_sr);

  double fl_gr[5], fr_gr[5];
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, ql, grm->prodl_local, fl_sr, fl_gr);
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, qr, grm->prodr_local, fr_sr, fr_gr);

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);

  double *w0 = &waves[0], *w1 = &waves[5];
  if (!excise_l && !excise_r) {
    for (int i = 0; i < 5; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr_gr[i] - fl_gr[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr_gr[i] - fl_gr[i]) / amax);
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
  // Identical to the regular mod variant's wave_roe: the special-relativistic
  // Roe decomposition operates only on the hydro slice and never sees the
  // spacetime block, so the tetrad split has no effect here.
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double vl[5], vr[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, ql, grm->prodl_local, vl);
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, qr, grm->prodr_local, vr);

  double rho_l = vl[0], vx_l = vl[1], vy_l = vl[2], vz_l = vl[3], p_l = vl[4];
  double rho_r = vr[0], vx_r = vr[1], vy_r = vr[2], vz_r = vr[3], p_r = vr[4];

  double Etot_l = ql[4];
  double Etot_r = qr[4];
  (void)Etot_l; (void)Etot_r;

  double W_l = 1.0 / sqrt(1.0 - ((vx_l*vx_l) + (vy_l*vy_l) + (vz_l*vz_l)));
  double W_r = 1.0 / sqrt(1.0 - ((vx_r*vx_r) + (vy_r*vy_r) + (vz_r*vz_r)));

  // Eulderink–Mellema (1995) Roe averaging weights. K_i² = ρ_i h_i W_i² is
  // built from SR-conserved quantities; the conserved state ql/qr carries an
  // explicit sqrt(γ) prefactor for any spacetime, so we strip it off here
  // before forming K. This keeps the τ-convention (q[4] = ρhW² − p − ρW with
  // rest mass subtracted, so the non-rel limit τ → p/(γ−1) + ½ρ|v|² is
  // smooth) consistent with how K reconstructs ρhW²: explicitly add D and p
  // back rather than carrying √γ-coupled terms around the algebra.
  //
  // With this K, the Roe-averaged primitives
  //   ṽ_α = (K_L z_α,L + K_R z_α,R) / (K_L + K_R)
  // for z = (W, Wu, Wv, Ww, ε) with ε = p/(ρh) yield an A_Roe(ṽ) satisfying
  // A_Roe·∆q = ∆f exactly (the defining Roe property). The previous code
  // used K = sqrt(τ+p)/W and v4 = avg(p/K), neither matching EM, which broke
  // both the wave-sum and flux-jump identities (verified empirically; see
  // the Roe-properties test before/after this fix).
  double h_l = 1.0 + ((p_l / rho_l) * gas_gamma / (gas_gamma - 1.0));
  double h_r = 1.0 + ((p_r / rho_r) * gas_gamma / (gas_gamma - 1.0));
  double eps_l = p_l / (rho_l * h_l);
  double eps_r = p_r / (rho_r * h_r);

  double sqrt_det_l = sqrt(grm->prodl_local[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_r = sqrt(grm->prodr_local[GKYL_GR_SP_SPATIAL_DET]);
  double D_sr_l   = ql[0] / sqrt_det_l;             // = ρW (SR rest-mass)
  double D_sr_r   = qr[0] / sqrt_det_r;
  double tau_sr_l = ql[4] / sqrt_det_l;             // = ρhW² − p − ρW (SR τ)
  double tau_sr_r = qr[4] / sqrt_det_r;

  // Paper Eq. (10.3): K² = √(−g) ρh. For Minkowski √(−g)=1, so K = √(ρh).
  // Since D + τ + p = ρhW² (rest-mass-subtracted τ + D + p collapses to
  // the full energy density), we have K = √((D+τ+p)/W²) = √(D+τ+p) / W.
  double K_l = sqrt(D_sr_l + tau_sr_l + p_l) / W_l;
  double K_r = sqrt(D_sr_r + tau_sr_r + p_r) / W_r;
  double K_avg = 1.0 / (K_l + K_r);

  double v0 = ((K_l * W_l)        + (K_r * W_r))        * K_avg;
  double v1 = ((K_l * W_l * vx_l) + (K_r * W_r * vx_r)) * K_avg;
  double v2 = ((K_l * W_l * vy_l) + (K_r * W_r * vy_r)) * K_avg;
  double v3 = ((K_l * W_l * vz_l) + (K_r * W_r * vz_r)) * K_avg;
  double v4 = ((K_l * eps_l)      + (K_r * eps_r))      * K_avg;

  double c_minus = 1.0 - ((gas_gamma / (gas_gamma - 1.0)) * v4);
  double c_plus  = 1.0 + ((gas_gamma / (gas_gamma - 1.0)) * v4);

  double v_alpha_sq = -(v0*v0) + (v1*v1) + (v2*v2) + (v3*v3);
  double s_sq = (0.5 * gas_gamma * v4 * (1.0 - v_alpha_sq))
              - (0.5 * (gas_gamma - 1.0) * (1.0 + v_alpha_sq));
  double energy = (v0*v0) - (v1*v1);
  double y = sqrt(((1.0 - (gas_gamma * v4)) * energy) + s_sq);

  // Wave amplitudes derived from inverting the 5×5 eigenvector basis with
  // the corrected eigenvector forms (see derivation in
  // wv_gr_euler_roe_derivation.md, sections "Translated to code ordering"
  // and "Corrected wave amplitudes"). The eigenvectors below carry the
  // rest-mass-subtraction terms in the τ slot that the original code
  // omitted; the amplitudes here are the unique solutions of
  // Σ a_k r_k = ∆q for that basis.
  //
  // The combined Row-0 + Row-4 of the inverse system collapses (the τ
  // slot cancels the D slot's c_-, a3·s²/(Γ−1) and c_+·(a4 v_2 + a5 v_3)
  // contributions) into the paper's E-slot equation:
  //   v^0 A + (s/y) v^1 B = δ[0] + δ[4]
  // paired with Row 1 (Sx):
  //   v^1 A + (s/y) v^0 B = δ[1]
  // gives A = a1+a2+a3, B = a2−a1.
  double sum04 = delta[0] + delta[4];
  double A_sum  = ((v0 * sum04) - (v1 * delta[1])) / energy;
  double B_diff = y * ((v0 * delta[1]) - (v1 * sum04)) / (sqrt(s_sq) * energy);

  double a4 = delta[2] - (v2 * A_sum);
  double a5 = delta[3] - (v3 * A_sum);

  double a3 = ((gas_gamma - 1.0) / s_sq) *
              (delta[0] - (A_sum * c_minus) + (c_plus * ((v2 * a4) + (v3 * a5))));

  double a1 = 0.5 * (A_sum - a3 - B_diff);
  double a2 = 0.5 * (A_sum - a3 + B_diff);

  for (int i = 0; i < 5 * 3; i++) waves[i] = 0.0;

  // Eigenvectors with the rest-mass-subtraction τ-slot terms (EM 1995
  // Eq. 10.15 translated from (D, E, Sx, Sy, Sz) into the code ordering
  // (D, Sx, Sy, Sz, τ) via τ = E − D). See wv_gr_euler_roe_derivation.md.
  // The previous implementation kept the bare paper components in the τ
  // slot, which broke the flux-jump identity.
  double *wv;
  wv = &waves[0 * 5];
  wv[0] = a1 * c_minus;
  wv[1] = a1 * (v1 - ((sqrt(s_sq) * v0) / y));
  wv[2] = a1 * v2;
  wv[3] = a1 * v3;
  wv[4] = a1 * (v0 - ((sqrt(s_sq) * v1) / y) - c_minus);
  s[0] = (((1.0 - (gas_gamma * v4)) * v0 * v1) - (sqrt(s_sq) * y))
       / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  wv = &waves[1 * 5];
  wv[0] = (a3 * (c_minus + (s_sq / (gas_gamma - 1.0))))
        - (a4 * c_plus * v2) - (a5 * c_plus * v3);
  wv[1] = a3 * v1;
  wv[2] = (a3 * v2) + a4;
  wv[3] = (a3 * v3) + a5;
  wv[4] = (a3 * (v0 - c_minus - (s_sq / (gas_gamma - 1.0))))
        + (a4 * c_plus * v2) + (a5 * c_plus * v3);
  s[1] = v1 / v0;

  wv = &waves[2 * 5];
  wv[0] = a2 * c_minus;
  wv[1] = a2 * (v1 + ((sqrt(s_sq) * v0) / y));
  wv[2] = a2 * v2;
  wv[3] = a2 * v3;
  wv[4] = a2 * (v0 + ((sqrt(s_sq) * v1) / y) - c_minus);
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
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;
  const double *pl = grm->prodl_local;
  const double *pr = grm->prodr_local;

  double vl[5], vr[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, ql, pl, vl);
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, qr, pr, vr);

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

  // Tetrad split: compute flat-space SR fluxes, then apply GR correction.
  double fl_sr[5], fr_sr[5];
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, ql, pl, fl_sr);
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qr, pr, fr_sr);

  double fl[5], fr[5];
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, ql, pl, fl_sr, fl);
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, qr, pr, fr_sr, fr);

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
  // Mirrors the regular mod variant's flux_jump_func but uses the corrected
  // (flat + GR-correction) flux to match packed wv_gr_euler_tetrad.c:1581.
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double fl_sr[5], fr_sr[5];
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, ql, grm->prodl_local, fl_sr);
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qr, grm->prodr_local, fr_sr);

  double fl[5], fr[5];
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, ql, grm->prodl_local, fl_sr, fl);
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, qr, grm->prodr_local, fr_sr, fr);

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (!excise_l && !excise_r) {
    for (int m = 0; m < 5; m++) flux_jump[m] = fr[m] - fl[m];
  } else {
    for (int m = 0; m < 5; m++) flux_jump[m] = 0.0;
  }

  double amaxl = gkyl_gr_euler_tetrad_mod_max_abs_speed(gas_gamma, ql, grm->prodl_local);
  double amaxr = gkyl_gr_euler_tetrad_mod_max_abs_speed(gas_gamma, qr, grm->prodr_local);
  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  if (!grm->auxfields.prods) {
    return q[0] >= 0.0;
  }
  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);

  double v[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(grm->gas_gamma, q, prods, v);
  if (v[0] < 0.0 || v[4] < 0.0) return false;
  return true;
}

static double
max_speed_func(const struct gkyl_wv_eqn *eqn, const double *q)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  if (!grm->auxfields.prods) return 1.0;
  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);
  return gkyl_gr_euler_tetrad_mod_max_abs_speed(grm->gas_gamma, q, prods);
}

static inline void
gr_euler_tetrad_mod_cons_to_diag(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *diag)
{
  for (int i = 0; i < 5; i++) diag[i] = qin[i];
}

static inline void
gr_euler_tetrad_mod_source(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *sout)
{
  // Integrated by moment_spacetime_coupling.
  for (int i = 0; i < 5; i++) sout[i] = 0.0;
}

// ---------------------------------------------------------------------------
// Constructors and lifetime.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_tetrad_mod_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_gr_euler_tetrad_mod *grm = container_of(base,
    struct wv_gr_euler_tetrad_mod, eqn);
  gkyl_free(grm);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_mod_new(double gas_gamma,
  const struct gkyl_range *conf_range, bool use_gpu)
{
  return gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp) {
      .gas_gamma = gas_gamma,
      .conf_range = *conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_HLL,
      .use_gpu = use_gpu,
    });
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_mod_inew(
  const struct gkyl_wv_gr_euler_tetrad_mod_inp *inp)
{
  struct wv_gr_euler_tetrad_mod *grm =
    gkyl_malloc(sizeof(struct wv_gr_euler_tetrad_mod));

  grm->eqn.type = GKYL_EQN_GR_EULER_TETRAD_MOD;
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

  if (inp->rp_type == WV_GR_EULER_TETRAD_RP_LAX) {
    grm->eqn.num_waves = 2;
    grm->eqn.waves_func = wave_lax_l;
    grm->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_EULER_TETRAD_RP_ROE) {
    grm->eqn.num_waves = 3;
    grm->eqn.waves_func = wave_roe_l;
    grm->eqn.qfluct_func = qfluct_roe_l;
  }
  else {  // default: HLL
    grm->eqn.num_waves = 2;
    grm->eqn.waves_func = wave_hll_l;
    grm->eqn.qfluct_func = qfluct_hll_l;
  }

  grm->eqn.ffluct_func = NULL;
  grm->eqn.flux_jump = flux_jump_func;
  grm->eqn.check_inv_func = check_inv;
  grm->eqn.max_speed_func = max_speed_func;
  grm->eqn.rotate_to_local_func = rot_to_local;
  grm->eqn.rotate_to_global_func = rot_to_global;
  grm->eqn.wall_bc_func = gr_euler_tetrad_mod_wall;
  grm->eqn.no_slip_bc_func = gr_euler_tetrad_mod_no_slip;
  grm->eqn.cons_to_riem = cons_to_riem;
  grm->eqn.riem_to_cons = riem_to_cons;
  grm->eqn.cons_to_diag = gr_euler_tetrad_mod_cons_to_diag;
  grm->eqn.source_func = gr_euler_tetrad_mod_source;

  grm->eqn.set_interface_idx_func = gr_euler_tetrad_mod_set_interface_idx;
  grm->eqn.set_cell_idx_func = gr_euler_tetrad_mod_set_cell_idx;

  grm->eqn.embed_geo = NULL;
  grm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(grm->eqn.flags);
  grm->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_tetrad_mod_free);
  grm->eqn.on_dev = &grm->eqn;

  return &grm->eqn;
}

double
gkyl_wv_gr_euler_tetrad_mod_gas_gamma(const struct gkyl_wv_eqn *eqn)
{
  const struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  return grm->gas_gamma;
}
