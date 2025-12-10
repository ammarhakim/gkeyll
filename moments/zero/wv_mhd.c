#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_mhd_priv.h>

void
gkyl_mhd_prim_vars(double gas_gamma, const double q[9], double v[9])
{
  double rho = q[0];
  double momx = q[1];
  double momy = q[2];
  double momz = q[3];
  double Etot = q[4];
  double Bx = q[5];
  double By = q[6];
  double Bz = q[7];
  double psi = q[8];

  double vx = momx / rho;
  double vy = momy / rho;
  double vz = momz / rho;

  double p = (gas_gamma - 1.0) * (Etot - (0.5 * rho * ((vx * vx) + (vy * vy) + (vz * vz))) - (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))));

  v[0] = rho;
  v[1] = vx;
  v[2] = vy;
  v[3] = vz;
  v[4] = p;
  v[5] = Bx;
  v[6] = By;
  v[7] = Bz;
  v[8] = psi;
}

static inline double
gkyl_mhd_max_abs_speed(double gas_gamma, const double q[9])
{
  double v[9] = { 0.0 };
  gkyl_mhd_prim_vars(gas_gamma, q, v);

  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];
  double Bx = v[5];
  double By = v[6];
  double Bz = v[7];

  double v_mag = sqrt((vx * vx) + (vy * vy) + (vz * vz));
  double B_mag = sqrt((Bx * Bx) + (By * By) + (Bz * Bz));

  double max_abs_speed = 0.0;
  double alfven_eig = Bx / sqrt(rho);
  double slow_magnetosonic_eig = sqrt(((gas_gamma * p) + (B_mag * B_mag) - sqrt(((gas_gamma * p) + (B_mag * B_mag)) * ((gas_gamma * p) + (B_mag * B_mag))
    - (4.0 * gas_gamma * p * (Bx * Bx)))) / (2.0 * rho));
  double fast_magnetosonic_eig = sqrt(((gas_gamma * p) + (B_mag * B_mag) + sqrt(((gas_gamma * p) + (B_mag * B_mag)) * ((gas_gamma * p) + (B_mag * B_mag))
    - (4.0 * gas_gamma * p * (Bx * Bx)))) / (2.0 * rho));

  if (fabs(v_mag) + fabs(alfven_eig) > max_abs_speed) {
    max_abs_speed = fabs(v_mag) + fabs(alfven_eig);
  }
  if (fabs(v_mag) + fabs(slow_magnetosonic_eig) > max_abs_speed) {
    max_abs_speed = fabs(v_mag) + fabs(slow_magnetosonic_eig);
  }
  if (fabs(v_mag) + fabs(fast_magnetosonic_eig) > max_abs_speed) {
    max_abs_speed = fabs(v_mag) + fabs(fast_magnetosonic_eig);
  }

  return max_abs_speed;
}

void
gkyl_mhd_flux(double gas_gamma, double light_speed, double b_fact, const double q[9], double flux[9])
{
  double v[9] = { 0.0 };
  gkyl_mhd_prim_vars(gas_gamma, q, v);

  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];
  double E = q[4];
  double Bx = v[5];
  double By = v[6];
  double Bz = v[7];
  double psi = v[8];

  flux[0] = rho * vx;
  flux[1] = (rho * (vx * vx)) + (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz)))) - (Bx * Bx);
  flux[2] = (rho * (vx * vy)) - (Bx * By);
  flux[3] = (rho * (vx * vz)) - (Bx * Bz);
  flux[4] = (E * vx) + (vx * (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))))) - (Bx * ((vx * Bx) + (vy * By) + (vz * Bz)));

  flux[5] = b_fact * psi;
  flux[6] = (By * vx) - (Bx * vy);
  flux[7] = (Bz * vx) - (Bx * vz);
  flux[8] = b_fact * (light_speed * light_speed) * Bx;
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 9; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 9; i++) {
    qout[i] = win[i];
  }
}

static void
mhd_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 9; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
  ghost[5] = -ghost[5];
}

static void
mhd_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 9; i++) {
    if ((i > 0 && i < 4) || (i > 4 && i < 8)) {
      ghost[i] = -skin[i];
    }
    else {
      ghost[i] = skin[i];
    }
  }
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];

  qlocal[5] = (qglobal[5] * norm[0]) + (qglobal[6] * norm[1]) + (qglobal[7] * norm[2]);
  qlocal[6] = (qglobal[5] * tau1[0]) + (qglobal[6] * tau1[1]) + (qglobal[7] * tau1[2]);
  qlocal[7] = (qglobal[5] * tau2[0]) + (qglobal[6] * tau2[1]) + (qglobal[7] * tau2[2]);
  qlocal[8] = qglobal[8];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];

  qglobal[5] = (qlocal[5] * norm[0]) + (qlocal[6] * tau1[0]) + (qlocal[7] * tau2[0]);
  qglobal[6] = (qlocal[5] * norm[1]) + (qlocal[6] * tau1[1]) + (qlocal[7] * tau2[1]);
  qglobal[7] = (qlocal[5] * norm[2]) + (qlocal[6] * tau1[2]) + (qlocal[7] * tau2[2]);
  qglobal[8] = qlocal[8];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);

  double gas_gamma = mhd->gas_gamma;
  double light_speed = mhd->light_speed;
  double b_fact = mhd->b_fact;

  double sl = gkyl_mhd_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_mhd_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[9], fr[9];
  gkyl_mhd_flux(gas_gamma, light_speed, b_fact, ql, fl);
  gkyl_mhd_flux(gas_gamma, light_speed, b_fact, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[9];
  for (int i = 0; i < 9; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[9];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 9; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_hll(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);

  double gas_gamma = mhd->gas_gamma;
  double light_speed = mhd->light_speed;
  double b_fact = mhd->b_fact;

  double vl[9] = { 0.0 };
  double vr[9] = { 0.0 };
  gkyl_mhd_prim_vars(gas_gamma, ql, vl);
  gkyl_mhd_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[0];
  double vx_l = vl[1];
  double vy_l = vl[2];
  double vz_l = vl[3];
  double p_l = vl[4];
  double Bx_l = vl[5];
  double By_l = vl[6];
  double Bz_l = vl[7];

  double B_mag_l = sqrt((Bx_l * Bx_l) + (By_l * By_l) + (Bz_l * Bz_l));

  double max_eig_l = 0.0;
  double alfven_eig_l = Bx_l / sqrt(rho_l);
  double slow_magnetosonic_eig_l = sqrt(((gas_gamma * p_l) + (B_mag_l * B_mag_l) - sqrt(((gas_gamma * p_l) + (B_mag_l * B_mag_l)) * ((gas_gamma * p_l) + (B_mag_l * B_mag_l))
    - (4.0 * gas_gamma * p_l * (Bx_l * Bx_l)))) / (2.0 * rho_l));
  double fast_magnetosonic_eig_l = sqrt(((gas_gamma * p_l) + (B_mag_l * B_mag_l) + sqrt(((gas_gamma * p_l) + (B_mag_l * B_mag_l)) * ((gas_gamma * p_l) + (B_mag_l * B_mag_l))
    - (4.0 * gas_gamma * p_l * (Bx_l * Bx_l)))) / (2.0 * rho_l));

  if (fabs(alfven_eig_l) > max_eig_l) {
    max_eig_l = fabs(alfven_eig_l);
  }
  if (fabs(slow_magnetosonic_eig_l) > max_eig_l) {
    max_eig_l = fabs(slow_magnetosonic_eig_l);
  }
  if (fabs(fast_magnetosonic_eig_l) > max_eig_l) {
    max_eig_l = fabs(fast_magnetosonic_eig_l);
  }

  double rho_r = vr[0];
  double vx_r = vr[1];
  double vy_r = vr[2];
  double vz_r = vr[3];
  double p_r = vr[4];
  double Bx_r = vr[5];
  double By_r = vr[6];
  double Bz_r = vr[7];

  double B_mag_r = sqrt((Bx_r * Bx_r) + (By_r * By_r) + (Bz_r * Bz_r));

  double max_eig_r = 0.0;
  double alfven_eig_r = Bx_r / sqrt(rho_r);
  double slow_magnetosonic_eig_r = sqrt(((gas_gamma * p_r) + (B_mag_r * B_mag_r) - sqrt(((gas_gamma * p_r) + (B_mag_r * B_mag_r)) * ((gas_gamma * p_r) + (B_mag_r * B_mag_r))
    - (4.0 * gas_gamma * p_r * (Bx_r * Bx_r)))) / (2.0 * rho_r));
  double fast_magnetosonic_eig_r = sqrt(((gas_gamma * p_r) + (B_mag_r * B_mag_r) + sqrt(((gas_gamma * p_r) + (B_mag_r * B_mag_r)) * ((gas_gamma * p_r) + (B_mag_r * B_mag_r))
    - (4.0 * gas_gamma * p_r * (Bx_r * Bx_r)))) / (2.0 * rho_r));

  if (fabs(alfven_eig_r) > max_eig_r) {
    max_eig_r = fabs(alfven_eig_r);
  }
  if (fabs(slow_magnetosonic_eig_r) > max_eig_r) {
    max_eig_r = fabs(slow_magnetosonic_eig_r);
  }
  if (fabs(fast_magnetosonic_eig_r) > max_eig_r) {
    max_eig_r = fabs(fast_magnetosonic_eig_r);
  }

  double vx_avg = 0.5 * (vx_l + vx_r);
  double max_eig_avg = 0.5 * (max_eig_l + max_eig_r);

  double sl = fmin(vx_l - max_eig_l, vx_r - max_eig_r);
  double sr = fmax(vx_l + max_eig_l, vx_r + max_eig_r);

  double fl[9] = { 0.0 };
  double fr[9] = { 0.0 };
  gkyl_mhd_flux(gas_gamma, light_speed, b_fact, ql, fl);
  gkyl_mhd_flux(gas_gamma, light_speed, b_fact, qr, fr);

  double qm[9] = { 0.0 };
  double *w0 = &waves[0], *w1 = &waves[9];
  for (int i = 0; i < 9; i++) {
    qm[i] = ((sr * qr[i]) - (sl * ql[i]) + (fl[i] - fr[i])) / (sr - sl);

    w0[i] = qm[i] - ql[i];
    w1[i] = qr[i] - qm[i];
  }

  s[0] = sl;
  s[1] = sr;

  return fmax(fabs(sl), fabs(sr));
}

static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[9];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 9; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_hll(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);

  double gas_gamma = mhd->gas_gamma;
  double light_speed = mhd->light_speed;
  double b_fact = mhd->b_fact;

  double fl[9], fr[9];
  gkyl_mhd_flux(gas_gamma, light_speed, b_fact, ql, fl);
  gkyl_mhd_flux(gas_gamma, light_speed, b_fact, qr, fr);

  for (int m = 0; m < 9; m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_mhd_max_abs_speed(gas_gamma, ql);
  double amaxr = gkyl_mhd_max_abs_speed(gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double gas_gamma = mhd->gas_gamma;

  double v[9] = { 0.0 };
  gkyl_mhd_prim_vars(gas_gamma, q, v);

  if (v[0] < 0.0 || v[4] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double gas_gamma = mhd->gas_gamma;

  return gkyl_mhd_max_abs_speed(gas_gamma, q);
}

static inline void
mhd_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 9; i++) {
    diag[i] = qin[i];
  }
}

static inline void
mhd_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 9; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_mhd_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_mhd *mhd = container_of(base->on_dev, struct wv_mhd, eqn);
    gkyl_cu_free(mhd);
  }

  struct wv_mhd *mhd = container_of(base, struct wv_mhd, eqn);
  gkyl_free(mhd);
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_new(double gas_gamma, double light_speed, double b_fact, bool use_gpu)
{
  return gkyl_wv_mhd_inew(&(struct gkyl_wv_mhd_inp) {
      .gas_gamma = gas_gamma,
      .light_speed = light_speed,
      .b_fact = b_fact,
      .rp_type = WV_MHD_RP_HLL,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_inew(const struct gkyl_wv_mhd_inp* inp)
{
  struct wv_mhd *mhd = gkyl_malloc(sizeof(struct wv_mhd));

  mhd->eqn.type = GKYL_EQN_MHD;
  mhd->eqn.num_equations = 9;
  mhd->eqn.num_diag = 9;

  mhd->gas_gamma = inp->gas_gamma;
  mhd->light_speed = inp->light_speed;
  mhd->b_fact = inp->b_fact;

  if (inp->rp_type == WV_MHD_RP_LAX) {
    mhd->eqn.num_waves = 2;
    mhd->eqn.waves_func = wave_lax_l;
    mhd->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_MHD_RP_HLL) {
    mhd->eqn.num_waves = 2;
    mhd->eqn.waves_func = wave_hll_l;
    mhd->eqn.qfluct_func = qfluct_hll_l;
  }

  mhd->eqn.flux_jump = flux_jump;
  mhd->eqn.check_inv_func = check_inv;
  mhd->eqn.max_speed_func = max_speed;
  mhd->eqn.rotate_to_local_func = rot_to_local;
  mhd->eqn.rotate_to_global_func = rot_to_global;
  
  mhd->eqn.wall_bc_func = mhd_wall;
  mhd->eqn.no_slip_bc_func = mhd_no_slip;

  mhd->eqn.cons_to_riem = cons_to_riem;
  mhd->eqn.riem_to_cons = riem_to_cons;

  mhd->eqn.cons_to_diag = mhd_cons_to_diag;

  mhd->eqn.source_func = mhd_source;

  mhd->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mhd->eqn.flags);
  mhd->eqn.ref_count = gkyl_ref_count_init(gkyl_mhd_free);
  mhd->eqn.on_dev = &mhd->eqn; // On the CPU, the equation object points ot itself.

  mhd->eqn.embed_geo = NULL;

  return &mhd->eqn;
}

double
gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double gas_gamma = mhd->gas_gamma;

  return gas_gamma;
}

double
gkyl_wv_mhd_light_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double light_speed = mhd->light_speed;

  return light_speed;
}

double
gkyl_wv_mhd_b_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double b_fact = mhd->b_fact;

  return b_fact;
}