#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_mhd_rgfm.h>
#include <gkyl_wv_mhd_rgfm_priv.h>

void
gkyl_mhd_rgfm_prim_vars(int num_species, double* gas_gamma_s, const double* q, double* v)
{
  double rho_total = q[0];
  double momx_total = q[1];
  double momy_total = q[2];
  double momz_total = q[3];
  double E_total = q[4];
  double Bx_total = q[5];
  double By_total = q[6];
  double Bz_total = q[7];
  double reinit_param = q[7 + (2 * num_species)];

  double *level_set_cons_s = gkyl_malloc(sizeof(double[num_species - 1]));
  for (int i = 0; i < num_species - 1; i++) {
    level_set_cons_s[i] = q[8 + i];
  }

  double *rho_cons_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_cons_s[i] = q[7 + num_species + i];
  }

  double vx_total = momx_total / rho_total;
  double vy_total = momy_total / rho_total;
  double vz_total = momz_total / rho_total;

  double *level_set_s = gkyl_malloc(sizeof(double[num_species]));
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_s[i] = level_set_cons_s[i] / rho_total;
    level_set_total += level_set_s[i];
  }
  level_set_s[num_species - 1] = 1.0 - level_set_total;

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = rho_cons_s[i] / level_set_s[i];
  }

  double *p_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    p_s[i] = (gas_gamma_s[i] - 1.0) * (E_total - (0.5 * rho_total * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))) -
      (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))));
  }

  double p_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    p_total += level_set_s[i] * p_s[i];
  }

  v[0] = rho_total;
  v[1] = vx_total;
  v[2] = vy_total;
  v[3] = vz_total;
  v[4] = p_total;
  v[5] = Bx_total;
  v[6] = By_total;
  v[7] = Bz_total;
  for (int i = 0; i < num_species - 1; i++) {
    v[8 + i] = level_set_s[i];
  }
  for (int i = 0; i < num_species; i++) {
    v[7 + num_species + i] = rho_s[i];
  }
  v[7 + (2 * num_species)] = reinit_param;

  gkyl_free(level_set_cons_s);
  gkyl_free(rho_cons_s);
  gkyl_free(level_set_s);
  gkyl_free(rho_s);
  gkyl_free(p_s);
}

static inline double
gkyl_mhd_rgfm_max_abs_speed(int num_species, double* gas_gamma_s, const double* q)
{
  double *v = gkyl_malloc(sizeof(double[8 + (2 * num_species)]));
  gkyl_mhd_rgfm_prim_vars(num_species, gas_gamma_s, q, v);

  double vx_total = v[1];
  double vy_total = v[2];
  double vz_total = v[3];
  double p_total = v[4];
  double Bx_total = v[5];
  double By_total = v[6];
  double Bz_total = v[7];

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = v[7 + num_species + i];
  }

  double v_mag = sqrt((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total));
  double B_mag = sqrt((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total));

  double max_abs_speed = 0.0;
  for (int i = 0; i < num_species; i++) {
    double alfven_eig = B_mag / sqrt(rho_s[i]);
    double slow_magnetosonic_eig = 0.5 * ((((gas_gamma_s[i] * p_total) + (B_mag * B_mag)) / rho_s[i]) -
      sqrt(((((gas_gamma_s[i] * p_total) + (B_mag * B_mag)) / rho_s[i]) * (((gas_gamma_s[i] * p_total) + (B_mag * B_mag)) / rho_s[i]))) -
      (4.0 * ((gas_gamma_s[i] * p_total) / rho_s[i]) * ((B_mag * B_mag) / rho_s[i])));
    double fast_magnetosonic_eig = 0.5 * ((((gas_gamma_s[i] * p_total) + (B_mag * B_mag)) / rho_s[i]) +
      sqrt(((((gas_gamma_s[i] * p_total) + (B_mag * B_mag)) / rho_s[i]) * (((gas_gamma_s[i] * p_total) + (B_mag * B_mag)) / rho_s[i]))) -
      (4.0 * ((gas_gamma_s[i] * p_total) / rho_s[i]) * ((B_mag * B_mag) / rho_s[i])));

    if (fabs(v_mag) + fabs(alfven_eig) > max_abs_speed) {
      max_abs_speed = fabs(v_mag) + fabs(alfven_eig);
    }
    if (fabs(v_mag) + fabs(slow_magnetosonic_eig) > max_abs_speed) {
      max_abs_speed = fabs(v_mag) + fabs(slow_magnetosonic_eig);
    }
    if (fabs(v_mag) + fabs(fast_magnetosonic_eig) > max_abs_speed) {
      max_abs_speed = fabs(v_mag) + fabs(fast_magnetosonic_eig);
    }
  }

  gkyl_free(v);
  gkyl_free(rho_s);

  return max_abs_speed;
}

void
gkyl_mhd_rgfm_flux(int num_species, double* gas_gamma_s, const double* q, double* flux)
{
  double *v = gkyl_malloc(sizeof(double[8 + (2 * num_species)]));
  gkyl_mhd_rgfm_prim_vars(num_species, gas_gamma_s, q, v);

  double rho_total = v[0];
  double vx_total = v[1];
  double vy_total = v[2];
  double vz_total = v[3];
  double p_total = v[4];
  double E_total = q[4];
  double Bx_total = v[5];
  double By_total = v[6];
  double Bz_total = v[7];
  double reinit_param = v[7 + (2 * num_species)];

  double *level_set_s = gkyl_malloc(sizeof(double[num_species]));
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_s[i] = v[8 + i];
    level_set_total += level_set_s[i];
  }
  level_set_s[num_species - 1] = 1.0 - level_set_total;

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = v[7 + num_species + i];
  }

  flux[0] = rho_total * vx_total;
  flux[1] = (rho_total * (vx_total * vx_total)) + (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)))) + (Bx_total * Bx_total);
  flux[2] = (rho_total * (vx_total * vy_total)) - (Bx_total * By_total);
  flux[3] = (rho_total * (vx_total * vz_total)) - (Bx_total * Bz_total);
  flux[4] = (E_total * vx_total) + (vx_total * (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))))) -
    (Bx_total * ((vx_total * Bx_total) + (vy_total * By_total) + (Bz_total * Bz_total)));

  flux[5] = 0.0;
  flux[6] = (By_total * vx_total) - (Bx_total * vy_total);
  flux[7] = (Bz_total * vx_total) - (Bx_total * vz_total);
  
  for (int i = 0; i < num_species - 1; i++) {
    flux[8 + i] = rho_total * (vx_total * level_set_s[i]);
  }
  for (int i = 0; i < num_species; i++) {
    flux[7 + num_species + i] = level_set_s[i] * (vx_total * rho_s[i]);
  }

  flux[7 + (2 * num_species)] = 0.0;

  gkyl_free(v);
  gkyl_free(level_set_s);
  gkyl_free(rho_s);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 8 + (2 * num_species); i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 8 + (2 * num_species); i++) {
    qout[i] = win[i];
  }
}

static void
mhd_rgfm_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  for (int i = 0; i < 8 + (2 * num_species); i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
  ghost[5] = -ghost[5];
}

static void
mhd_rgfm_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  for (int i = 0; i < 8 + (2 * num_species); i++) {
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
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  for (int i = 0; i < 8 + (2 * num_species); i++) {
    qlocal[i] = qglobal[i];
  }

  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);

  qlocal[5] = (qglobal[5] * norm[0]) + (qglobal[6] * norm[1]) + (qglobal[7] * norm[2]);
  qlocal[6] = (qglobal[5] * tau1[0]) + (qglobal[6] * tau1[1]) + (qglobal[7] * tau1[2]);
  qlocal[7] = (qglobal[5] * tau2[0]) + (qglobal[6] * tau2[1]) + (qglobal[7] * tau2[2]);
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  for (int i = 0; i < 8 + (2 * num_species); i++) {
    qglobal[i] = qlocal[i];
  }

  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);

  qglobal[5] = (qlocal[5] * norm[0]) + (qlocal[6] * tau1[0]) + (qlocal[7] * tau2[0]);
  qglobal[6] = (qlocal[5] * norm[1]) + (qlocal[6] * tau1[1]) + (qlocal[7] * tau2[1]);
  qglobal[7] = (qlocal[5] * norm[2]) + (qlocal[6] * tau1[2]) + (qlocal[7] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;
  double* gas_gamma_s = mhd_rgfm->gas_gamma_s;
  int reinit_freq = mhd_rgfm->reinit_freq;

  double sl = gkyl_mhd_rgfm_max_abs_speed(num_species, gas_gamma_s, ql);
  double sr = gkyl_mhd_rgfm_max_abs_speed(num_species, gas_gamma_s, qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double[8 + (2 * num_species)]));
  double *fr = gkyl_malloc(sizeof(double[8 + (2 * num_species)]));
  gkyl_mhd_rgfm_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_mhd_rgfm_flux(num_species, gas_gamma_s, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[8 + (2 * num_species)];
  for (int i = 0; i < 8 + (2 * num_species); i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  gkyl_free(fl);
  gkyl_free(fr);

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  const double *w0 = &waves[0], *w1 = &waves[8 + (2 * num_species)];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 8 + (2 * num_species); i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr,
  const double phil, const double phir, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr,
  const double phil, const double phir, const double* waves, const double* s, double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;
  double *gas_gamma_s = mhd_rgfm->gas_gamma_s;
  int reinit_freq = mhd_rgfm->reinit_freq;

  double *fr = gkyl_malloc(sizeof(double[8 + (2 * num_species)]));
  double *fl = gkyl_malloc(sizeof(double[8 + (2 * num_species)]));
  gkyl_mhd_rgfm_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_mhd_rgfm_flux(num_species, gas_gamma_s, qr, fr);

  for (int m = 0; m < 8 + (2 * num_species); m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_mhd_rgfm_max_abs_speed(num_species, gas_gamma_s, ql);
  double amaxr = gkyl_mhd_rgfm_max_abs_speed(num_species, gas_gamma_s, qr);
  
  gkyl_free(fr);
  gkyl_free(fl);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;
  double *gas_gamma_s = mhd_rgfm->gas_gamma_s;

  double *v = gkyl_malloc(sizeof(double[8 + (2 * num_species)]));
  gkyl_mhd_rgfm_prim_vars(num_species, gas_gamma_s, q, v);

  for (int i = 0; i < num_species; i++) {
    if (v[8 + num_species + i] < 0.0) {
      gkyl_free(v);
      return false;
    }
  }
  
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_total += v[8 + i];
  }
  if (level_set_total > 1.0) {
    gkyl_free(v);
    return false;
  }

  if (v[0] < 0.0 || v[4] < 0.0) {
    gkyl_free(v);
    return false;
  }
  else {
    gkyl_free(v);
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;
  double* gas_gamma_s = mhd_rgfm->gas_gamma_s;

  return gkyl_mhd_rgfm_max_abs_speed(num_species, gas_gamma_s, q);
}

static inline void
mhd_rgfm_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  for (int i = 0; i < 8 + (2 * num_species); i++) {
    diag[i] = qin[i];
  }
}

static inline void
mhd_rgfm_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  for (int i = 0; i < 8 + (2 * num_species); i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_mhd_rgfm_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_mhd_rgfm *mhd_rgfm = container_of(base->on_dev, struct wv_mhd_rgfm, eqn);
    gkyl_cu_free(mhd_rgfm);
  }

  struct wv_mhd_rgfm *mhd_rgfm = container_of(base, struct wv_mhd_rgfm, eqn);
  gkyl_free(mhd_rgfm);
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_rgfm_new(int num_species, double* gas_gamma_s, int reinit_freq, bool use_gpu)
{
  return gkyl_wv_mhd_rgfm_inew(&(struct gkyl_wv_mhd_rgfm_inp) {
      .num_species = num_species,
      .gas_gamma_s = gas_gamma_s,
      .reinit_freq = reinit_freq,
      .rp_type = WV_MHD_RGFM_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_rgfm_inew(const struct gkyl_wv_mhd_rgfm_inp* inp)
{
  struct wv_mhd_rgfm *mhd_rgfm = gkyl_malloc(sizeof(struct wv_mhd_rgfm));

  mhd_rgfm->eqn.type = GKYL_EQN_MHD_RGFM;
  mhd_rgfm->eqn.num_equations = 8 + (2 * inp->num_species);
  mhd_rgfm->eqn.num_diag = 8 + (2 * inp->num_species);

  mhd_rgfm->num_species = inp->num_species;
  mhd_rgfm->gas_gamma_s = inp->gas_gamma_s;
  mhd_rgfm->reinit_freq = inp->reinit_freq;

  if (inp->rp_type == WV_MHD_RGFM_RP_LAX) {
    mhd_rgfm->eqn.num_waves = 2;
    mhd_rgfm->eqn.waves_func = wave_lax_l;
    mhd_rgfm->eqn.qfluct_func = qfluct_lax_l;
  }

  mhd_rgfm->eqn.flux_jump = flux_jump;
  mhd_rgfm->eqn.check_inv_func = check_inv;
  mhd_rgfm->eqn.max_speed_func = max_speed;
  mhd_rgfm->eqn.rotate_to_local_func = rot_to_local;
  mhd_rgfm->eqn.rotate_to_global_func = rot_to_global;
  
  mhd_rgfm->eqn.wall_bc_func = mhd_rgfm_wall;
  mhd_rgfm->eqn.no_slip_bc_func = mhd_rgfm_no_slip;

  mhd_rgfm->eqn.cons_to_riem = cons_to_riem;
  mhd_rgfm->eqn.riem_to_cons = riem_to_cons;

  mhd_rgfm->eqn.cons_to_diag = mhd_rgfm_cons_to_diag;

  mhd_rgfm->eqn.source_func = mhd_rgfm_source;

  mhd_rgfm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mhd_rgfm->eqn.flags);
  mhd_rgfm->eqn.ref_count = gkyl_ref_count_init(gkyl_mhd_rgfm_free);
  mhd_rgfm->eqn.on_dev = &mhd_rgfm->eqn; // On the CPU, the equation object points ot itself.

  mhd_rgfm->eqn.embed_geo = NULL;

  return &mhd_rgfm->eqn;
}

int
gkyl_wv_mhd_rgfm_num_species(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int num_species = mhd_rgfm->num_species;

  return num_species;
}

double*
gkyl_wv_mhd_rgfm_gas_gamma_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  double *gas_gamma_s = mhd_rgfm->gas_gamma_s;

  return gas_gamma_s;
}

int
gkyl_wv_mhd_rgfm_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd_rgfm *mhd_rgfm = container_of(eqn, struct wv_mhd_rgfm, eqn);
  int reinit_freq = mhd_rgfm->reinit_freq;

  return reinit_freq;
}
