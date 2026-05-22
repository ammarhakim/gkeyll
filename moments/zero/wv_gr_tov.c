#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_tov.h>
#include <gkyl_wv_gr_tov_priv.h>


void
gkyl_gr_tov_prim_vars(double gas_gamma, const double q[8], double v[8])
{
  for (int i = 0; i < 8; i++) {
    v[i] = 0.0;
  }


  double Etot = q[1];
  double mom_r = q[2];
  double Phi = q[3]; // potential (in the metric term dt)
  double m = q[4];
  double r = q[5];

  double beta = 0.25 * (2 - gas_gamma);
  double p_term = (4.0 * beta * beta * Etot * Etot) + ((gas_gamma - 1.0) * ((Etot * Etot) - (mom_r * mom_r)));
  double p = -(2.0 * beta * Etot) + sqrt(fmax(0.0, p_term));

  double rho = p / (gas_gamma - 1.0);

  double vel = 0.0; // v^r
  double denom = Etot + p;
  if (denom > 0.0) {
    vel = mom_r / denom;
  }
  else {
    vel = 0.0;
  }

  v[0] = 0.0;

  v[1] = Etot; 
  v[2] = mom_r;
  v[3] = Phi; 
  v[4] = m;
  v[5] = r; 
  
  v[6] = rho;
  v[7] = vel; 
}

static inline double
gkyl_gr_tov_max_abs_speed(double gas_gamma, const double q[8])
{
  double v[8] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);

  double vel = v[7];
  double c_s = sqrt(gas_gamma - 1.0);

  //return fabs(vel) + c_s;
  return 1.0; // Return speed of light.
}

void
gkyl_gr_tov_flux(double gas_gamma, double kappa, const double q[8], double flux1[8], double flux2[8])
{
  for (int i = 0; i < 8; i++) {
    flux1[i] = 0.0;
  }

  for (int i = 0; i < 8; i++) {
    flux2[i] = 0.0;
  }

  double v[8] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);
  
  double mom_r = v[2];
  double Phi = v[3];
  double m = v[4];
  double r = v[5]; 

  double rho = v[6];
  double vel = v[7];

  double p = (gas_gamma - 1.0) * rho;

  double lapse = exp(Phi);
  double a = 1.0 / sqrt(1.0 - (2.0 * m / r));

  flux1[0] = 0.0;
  flux1[1] = r * r * (lapse / a ) * mom_r;
  flux1[2] = r * r * (lapse / a ) * mom_r * vel;
  flux1[3] = 0.0; // is there flux on Phi (lapse)?

  flux2[1] = 0.0;
  flux2[2] = (lapse / a ) * p;
}

static void
gkyl_gr_tov_flux_total(double gas_gamma, double kappa, const double q[8], double flux[8])
{
  double f1[8], f2[8];
  gkyl_gr_tov_flux(gas_gamma, kappa, q, f1, f2);

  double r = q[5];
  double r2 = r * r;

  for (int i = 0; i < 8; i++) {
    flux[i] = f2[i];
  }

  flux[1] += f1[1] / r2;
  flux[2] += f1[2] / r2;
}

static void
gkyl_gr_tov_flux_at_radius(double gas_gamma, double kappa, const double q[8], double r_face, double flux[8])
{
  for (int i = 0; i < 8; i++) {
    flux[i] = 0.0;
  }

  double v[8] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);

  double mom_r = v[2];
  double Phi = v[3];
  double m = v[4];
  double r_cell = v[5];
  double rho = v[6];
  double vel = v[7];

  double p = (gas_gamma - 1.0) * rho;
  double lapse = exp(Phi);
  double a = 1.0 / sqrt(1.0 - (2.0 * m / r_cell));
  double X = lapse / a;
  double r_face2 = r_face * r_face;

  double f1_E = r_face2 * X * mom_r;
  double f1_S = r_face2 * X * mom_r * vel;
  double f2_E = 0.0;
  double f2_S = X * p;

  flux[1] = f2_E + f1_E / r_face2;
  flux[2] = f2_S + f1_S / r_face2;
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 8; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 8; i++) {
    qout[i] = win[i];
  }
}

static void
gr_tov_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 8; i++) {
    ghost[i] = skin[i];
  }
  
  ghost[2] = -skin[2]; // radial momentum is odd across the spherical center.
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  for (int i = 0; i < 8; i++) {
    qlocal[i] = qglobal[i];
  }
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  for (int i = 0; i < 8; i++) {
    qglobal[i] = qlocal[i];
  }
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;
  double kappa = gr_tov->kappa;

  double sl = gkyl_gr_tov_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_tov_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  // The only change was to the physical flux used to construct the waves: 
  // instead of using a single flux f(q), I used the Appendix-B split spherical flux evaluated at the interface
  // => the waves are now built from the split flux at r_face instead of the old total flux using each state’s own q[5]
  double r_face = 0.5 * (ql[5] + qr[5]);
  double fl[8], fr[8];
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, fr);

  double *w0 = &waves[0], *w1 = &waves[8];
  for (int i = 0; i < 8; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }
  for (int i = 0; i < 8; i++) {
    if (i != 1 && i != 2) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[8];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 8; i++) {
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
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;
  double kappa = gr_tov->kappa;

  double vl[8], vr[8];
  gkyl_gr_tov_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_tov_prim_vars(gas_gamma, qr, vr);

  double c_sl = sqrt(gas_gamma - 1.0);
  double c_sr = sqrt(gas_gamma - 1.0);
  double vel_l = vl[7];
  double vel_r = vr[7];

  double lapse_l = exp(ql[3]);
  double lapse_r = exp(qr[3]);
  double a_l = 1.0 / sqrt(1.0 - (2.0 * ql[4] / ql[5]));
  double a_r = 1.0 / sqrt(1.0 - (2.0 * qr[4] / qr[5]));
  double metric_speed_l = lapse_l / a_l;
  double metric_speed_r = lapse_r / a_r;

  // vel is the locally measured radial velocity
  double slow_acoustic_eig_l = metric_speed_l * (vel_l - c_sl) / (1.0 - (vel_l * c_sl));
  double fast_acoustic_eig_l = metric_speed_l * (vel_l + c_sl) / (1.0 + (vel_l * c_sl));
  double slow_acoustic_eig_r = metric_speed_r * (vel_r - c_sr) / (1.0 - (vel_r * c_sr));
  double fast_acoustic_eig_r = metric_speed_r * (vel_r + c_sr) / (1.0 + (vel_r * c_sr));

  double sl = fmin(slow_acoustic_eig_l, slow_acoustic_eig_r); //most negative wave speed (after the acoustic eigenvalues estiamtion)
  double sr = fmax(fast_acoustic_eig_l, fast_acoustic_eig_r); //most postivie wave speed
  if (sr <= sl) {
    double amax = fmax(gkyl_gr_tov_max_abs_speed(gas_gamma, ql), gkyl_gr_tov_max_abs_speed(gas_gamma, qr));
    sl = -amax;
    sr = amax;
  }

  double fl[8], fr[8];
  double r_face = 0.5 * (ql[5] + qr[5]);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, fr);

  double qm[8];
  for (int i = 0; i < 8; i++) {
    qm[i] = ((sr * qr[i]) - (sl * ql[i]) + (fl[i] - fr[i])) / (sr - sl);
  }

  double *w0 = &waves[0], *w1 = &waves[8];
  for (int i = 0; i < 8; i++) {
    w0[i] = qm[i] - ql[i];
    w1[i] = qr[i] - qm[i];
  }
  for (int i = 0; i < 8; i++) {
    if (i != 1 && i != 2) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = sl;
  s[1] = sr;

  return fmax(fabs(sl), fabs(sr));
}

static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[8];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 8; i++) {
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

  return 0.0;
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
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;
  double kappa = gr_tov->kappa;

  double fr[8], fl[8];
  //gkyl_gr_tov_flux_total(gas_gamma, kappa, ql, fl);
  //gkyl_gr_tov_flux_total(gas_gamma, kappa, qr, fr);
  double r_face = 0.5 * (ql[5] + qr[5]);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, fr);

  for (int m = 0; m < 8; m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_gr_tov_max_abs_speed(gas_gamma, ql);
  double amaxr = gkyl_gr_tov_max_abs_speed(gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;

  double v[8] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);

  if (v[6] < 0.0 || fabs(v[7]) > 1.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;

  return gkyl_gr_tov_max_abs_speed(gas_gamma, q);
}

static inline void
gr_tov_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 8; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_tov_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
}

void
gkyl_gr_tov_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_tov *gr_tov = container_of(base->on_dev, struct wv_gr_tov, eqn);
    gkyl_cu_free(gr_tov);
  }

  struct wv_gr_tov *gr_tov = container_of(base, struct wv_gr_tov, eqn);
  gkyl_free(gr_tov);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_tov_new(double gas_gamma, double kappa, bool use_gpu)
{
  return gkyl_wv_gr_tov_inew(&(struct gkyl_wv_gr_tov_inp) {
      .gas_gamma = gas_gamma,
      .kappa = kappa,
      .rp_type = WV_GR_TOV_RP_HLL,//WV_GR_TOV_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_tov_inew(const struct gkyl_wv_gr_tov_inp* inp)
{
  struct wv_gr_tov *gr_tov = gkyl_malloc(sizeof(struct wv_gr_tov));

  gr_tov->eqn.type = GKYL_EQN_GR_TOV;
  gr_tov->eqn.num_equations = 8;
  gr_tov->eqn.num_diag = 8;

  gr_tov->gas_gamma = inp->gas_gamma;
  gr_tov->kappa = inp->kappa;

  if (inp->rp_type == WV_GR_TOV_RP_LAX) {
    gr_tov->eqn.num_waves = 2;
    gr_tov->eqn.waves_func = wave_lax_l;
    gr_tov->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_TOV_RP_HLL) {
    gr_tov->eqn.num_waves = 2;
    gr_tov->eqn.waves_func = wave_hll_l;
    gr_tov->eqn.qfluct_func = qfluct_hll_l;
  }

  gr_tov->eqn.flux_jump = flux_jump;
  gr_tov->eqn.check_inv_func = check_inv;
  gr_tov->eqn.max_speed_func = max_speed;
  gr_tov->eqn.rotate_to_local_func = rot_to_local;
  gr_tov->eqn.rotate_to_global_func = rot_to_global;

  gr_tov->eqn.wall_bc_func = gr_tov_wall;
  
  gr_tov->eqn.cons_to_riem = cons_to_riem;
  gr_tov->eqn.riem_to_cons = riem_to_cons;

  gr_tov->eqn.cons_to_diag = gr_tov_cons_to_diag;

  gr_tov->eqn.source_func = gr_tov_source;

  gr_tov->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_tov->eqn.flags);
  gr_tov->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_tov_free);
  gr_tov->eqn.on_dev = &gr_tov->eqn; // On the CPU, the equation object points to itself.

  gr_tov->eqn.embed_geo = NULL;

  return &gr_tov->eqn;
}

double
gkyl_wv_gr_tov_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;

  return gas_gamma;
}

double
gkyl_wv_gr_tov_kappa(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double kappa = gr_tov->kappa;

  return kappa;
}
