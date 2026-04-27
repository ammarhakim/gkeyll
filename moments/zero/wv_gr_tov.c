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


  double PI = q[1]; // tau + mom_r
  double PHI = q[2]; // tau - mom_r
  double Phi = q[3]; // potential (in the metric term dt)
  double m = q[4];
  double r = q[5];

  double beta = 0.25 * (2 - gas_gamma);
  double p_term = (beta * beta * (PI + PHI) * (PI + PHI) + (gas_gamma - 1.0) * PI * PHI);
  double p = - beta * (PI + PHI) + sqrt(fmax(0.0, p_term));

  double rho = p / (gas_gamma - 1.0);

  double vel = 0.0; // v^r
  double denom = PI + PHI + 2.0 * p;
  if (denom > 0.0) {
    vel = (PI - PHI) / denom;
  }
  else {
    vel = 0.0;
  }

  v[0] = 0.0;

  v[1] = PI; 
  v[2] = PHI;
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
gkyl_gr_tov_flux(double gas_gamma, double kappa, const double q[8], double flux[8])
{
  for (int i = 0; i < 8; i++) {
    flux[i] = 0.0;
  }

  double v[8] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);
  
  double PI = v[1];
  double PHI = v[2];
  double Phi = v[3];

  double rho = v[6];
  double vel = v[7];

  double p = (gas_gamma - 1.0) * rho;

  flux[0] = 0.0;

  flux[1] = 0.5 * (PI - PHI) * (1 + vel) + p; 
  flux[2] = 0.5 * (PI - PHI) * (1 - vel) - p;
  flux[3] = 0.0; // is there flux on Phi (lapse)?
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
  for (int i = 0; i < 7; i++) {
    ghost[i] = skin[i];
  }
  
  ghost[7] = -skin[7];
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

  double fl[8], fr[8];
  gkyl_gr_tov_flux(gas_gamma, kappa, ql, fl);
  gkyl_gr_tov_flux(gas_gamma, kappa, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[8];
  for (int i = 0; i < 8; i++) {
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
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;
  double kappa = gr_tov->kappa;

  double fr[8], fl[8];
  gkyl_gr_tov_flux(gas_gamma, kappa, ql, fl);
  gkyl_gr_tov_flux(gas_gamma, kappa, qr, fr);

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
      .rp_type = WV_GR_TOV_RP_LAX,
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
