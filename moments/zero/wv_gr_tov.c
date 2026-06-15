// Primitive variable recovery for the ideal-gas (Gamma-law) EOS:
//   p = (Gamma - 1) * rho0 * eps
// Entropy is now dynamical, so eps is an independent thermodynamic degree of freedom (it must be recovered from the energy equation, not from a barotrope)
//
// Conventions :
// q[0] = r^2 * D with D = rho0 * W
// q[1] = r^2 * tau with tau  = Etot - D = rho0 * h * W^2 - p - rho0 * W
// q[2] = r^2 * S with S = rho0 * h * W^2 * v
// q[3] = Phi, q[4] = m, q[5] = r (geometry, passed through)
//
// Pressure is a genuine independent primitive for a Gamma-law gas -- (rho0, vel) no longer determine it

#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_tov.h>
#include <gkyl_wv_gr_tov_priv.h>
#include "tov_solver.h"

static inline double
safe_r2(double r) // Numerical guard so we don't end up dividing by r = 0
{
  double r2 = r * r;
  if (r2 > 1.0e-300) {
    return r2;
  }
  else {
    return 1.0e-300;
  }
}

void
gkyl_gr_tov_prim_vars(double gas_gamma, const double q[8], double v[9])
{
  // v[0] = D       conserved rest-mass density  (D = rho0 * W)   <- not a primitive
  // v[1] = Etot    conserved energy             (tau + D)        <- not a primitive
  // v[2] = mom_r   conserved momentum S         (rho0 h W^2 v)   <- not a primitive
  // v[3] = Phi     geometry (passed through)
  // v[4] = m       geometry
  // v[5] = r       geometry
  // v[6] = rho     <- the actual primitive rest-mass density
  // v[7] = vel     <- primitive
  // v[8] = p       <- primitive
  for (int i = 0; i < 9; i++) {
    v[i] = 0.0;
  }

  double Phi = q[3]; // potential (in the metric term dt)
  double m = q[4];
  double r = q[5];
  double r2 = safe_r2(r);

  double D = q[0] / r2;
  double tau   = q[1] / r2;   // tau = Etot - D (Valencia energy, rest mass removed)
  double mom_r = q[2] / r2;

  double rho_floor = 1.0e-15;
  double p_floor = 1.0e-300;
  double D_safe = fmax(D, rho_floor);

  // |v| < 1 requires Q = Etot + p = tau + D + p > |S|. => 
  // set a hard lower bound on p so it can never leave the subluminal region.
  double p_min = fmax(p_floor, fabs(mom_r) - tau - D_safe + 1.0e-16);

  // Initial guess: (assume W~1), so rho0*eps ~ tau and p ~ (Gamma-1)*tau <- exact for a static star
  double p = fmax((gas_gamma - 1.0) * tau, p_min);

  double vel = 0.0, W = 1.0, rho = D_safe, eps = 0.0;

  for (int iter = 0; iter < 100; iter++) {
    double Q = tau + D_safe + p; // Q = rho0 h W^2, the relativistic inertia density
    vel = mom_r / Q;     // Q > |S|because p >= p_min so |vel| < 1
    if (vel * vel > 1.0 - 1.0e-12) { // velocity limiter, clamp |v| just below lightspeed sqrt(1.0 - 1.0e-12), keeping the original direction with copysign
      vel = copysign(sqrt(1.0 - 1.0e-12), vel);
    }
    double W2 = 1.0 / (1.0 - vel * vel);
    W = sqrt(W2);
    rho = D_safe / W;

    // Internal energy from the energy balance: rho0 eps W = tau - D(W-1) - p(W^2-1)
    eps = (tau - D_safe * (W - 1.0) - p * (W2 - 1.0)) / (D_safe * W);
    if (eps < 0.0) {
      eps = 0.0;
    }

    // residual (EOS pressure implied by the current guess, minus the guess)
    double p_eos = (gas_gamma - 1.0) * rho * eps;
    double f_val = p_eos - p;

    // Newton derivative df/dp = v^2 * cs^2 - 1
    double h = 1.0 + eps + p / rho;
    double cs2 = gas_gamma * p / (rho * h);
    if (cs2 < 0.0) {
      cs2 = 0.0;
    }
    if (cs2 > 1.0 - 1.0e-12) {
      cs2 = 1.0 - 1.0e-12;
    }
    double f_prime = vel * vel * cs2 - 1.0;

    double dp = -f_val / f_prime;
    p += dp;
    if (p < p_min) {
      p = p_min;
    }

    if (fabs(dp) < 1.0e-12 * fmax(fabs(p), p_floor)) {
      break;
    }
  }

  // Recompute the final state
  double Q = tau + D_safe + p;
  vel = mom_r / Q;
  if (vel * vel > 1.0 - 1.0e-12) {
    vel = copysign(sqrt(1.0 - 1.0e-12), vel);
  }
  double W2 = 1.0 / (1.0 - vel * vel);
  W = sqrt(W2);
  rho = D_safe / W;
  eps = (tau - D_safe * (W - 1.0) - p * (W2 - 1.0)) / (D_safe * W);
  if (eps < 0.0) {
    eps = 0.0;
    p = fmax((gas_gamma - 1.0) * rho * eps, p_min); // = p_min; cold fallback
  }

  double Etot = tau + D; 
  v[0] = D;
  v[1] = Etot;
  v[2] = mom_r;
  v[3] = Phi;
  v[4] = m;
  v[5] = r;
  v[6] = rho;
  v[7] = vel;
  v[8] = p;
}

static inline double
gkyl_gr_tov_max_abs_speed(double gas_gamma, const double q[8])
{
  double v[9] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);

  double Phi = v[3];
  double m = v[4];
  double r = v[5]; 
  double rho = v[6];
  double vel = v[7];
  double p = v[8];
  double eps = p / ((gas_gamma - 1.0) * rho);
  double h = 1.0 + eps + p / rho;
  double cs2 = gas_gamma * p / (rho * h);
  if (cs2 < 0.0) {
    cs2 = 0.0;
  }
  if (cs2 > 1.0 - 1.0e-12) {
    cs2 = 1.0 - 1.0e-12;
  }
  double c_s = sqrt(cs2);
  double lapse = exp(Phi);
  double a = 1.0 / sqrt(1.0 - (2.0 * m / r));

  double lambda_minus = (lapse / a) * (vel - c_s) / (1.0 - vel * c_s);
  double lambda_plus  = (lapse / a) * (vel + c_s) / (1.0 + vel * c_s);

  return fmax(fabs(lambda_minus), fabs(lambda_plus));
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

  double v[9] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);
  
  double D = v[0];
  double mom_r = v[2];
  double Phi = v[3];
  double m = v[4];
  double r = v[5]; 

  double vel = v[7];
  double p = v[8];

  double lapse = exp(Phi);
  double a = 1.0 / sqrt(1.0 - (2.0 * m / r));
  double X = lapse / a;

  flux1[0] = r * r * X * D * vel;
  flux1[1] = r * r * X * (mom_r - D * vel); // tau-flux = energy-flux minus mass-flux
  flux1[2] = r * r * X * mom_r * vel;
  flux1[3] = 0.0;

  flux2[0] = 0.0;
  flux2[1] = 0.0;
  flux2[2] = X * p;
}

static void
gkyl_gr_tov_flux_total(double gas_gamma, double kappa, const double q[8], double flux[8])
{
  double f1[8], f2[8];
  gkyl_gr_tov_flux(gas_gamma, kappa, q, f1, f2);

  double r = q[5];
  double r2 = safe_r2(r);

  for (int i = 0; i < 8; i++) {
    flux[i] = f1[i] + (r2 * f2[i]);
  }
}

// Flux of the r^2-weighted conserved variables evaluated AT an interface
static void
gkyl_gr_tov_flux_at_radius(double gas_gamma, double kappa, const double q[8], double r_face,
  double Phi_face, double m_face, double flux[8])
{
  for (int i = 0; i < 8; i++) {
    flux[i] = 0.0;
  }

  double v[9] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);

  double D = v[0];
  double mom_r = v[2];
  double vel = v[7];
  double p = v[8];
  double lapse = exp(Phi_face);
  // Regularity at the center: m ~ r^3 so 2m/r -> 0 as r -> 0 (a -> 1). At the center face the wall BC gives r_face = m_face = 0, so guard the 0/0 ratio.
  double compactness = 0.0;
  if (fabs(r_face) > 1.0e-300) {
    compactness = 2.0 * m_face / r_face;
  }
  double a = 1.0 / sqrt(1.0 - compactness);
  double X = lapse / a;
  double r_face2 = safe_r2(r_face);

  flux[0] = r_face2 * X * D * vel;
  flux[1] = r_face2 * X * (mom_r - D * vel); // tau-flux = energy-flux minus mass-flux
  flux[2] = r_face2 * X * (mom_r * vel + p);
}

// Well-balancing for the static TOV star. Returns deq[8], the jump in the r^2-weighted conserved variables of the "frozen initial TOV equilibrium" between the two cell-center radii r_l = ql[5], r_r = qr[5].
// The HLL dissipation acts on Delta_q = qr - ql. Subtracting deq makes it act on delta_q = Delta_q - deq. For the static star the discrete state equals the reference, so delta_q -> 0 and the spurious dissipative mass/energy fluxes
// vanish at every interface. 
// Conservation is conserved: the wave construction keeps the physical face-flux jump fr - fl, and deq only redistributes dissipation between amdq/apdq (their sum stays fr - fl) for any deq.
//
// The reference is read from the frozen conserved slots q[6]=r^2*p_eq, q[7]=r^2*D_eq (snapshotted at t=0 by the IC projection), so deq matches Delta_q to machine precision at t=0.
// Only the mass (0) and energy (1) channels get an equilibrium jump; at rest S = 0 so deq[2] = 0 and momentum stays fully HLL.

static void
gr_tov_wb_equilibrium_jump(const struct wv_gr_tov *gr_tov, const double* ql, const double* qr,
  double r_face, double Phi_face, double m_face, double deq[8], double *dfeq2)
{
  for (int i = 0; i < 8; i++) {
    deq[i] = 0.0;
  }
  *dfeq2 = 0.0;

  double r_l = ql[5], r_r = qr[5];
  if (r_l < 0.0 || r_r < 0.0) {
    return; // center ghost (reflected, r < 0): no equilibrium lookup, let the wall BC act
  }

  // We read the equilibrium straight out of the frozen conserved slots q[6] = r^2*p_eq and q[7] = r^2*D_eq, snapshotted at t=0 by the IC projection.
  // These undergo the IDENTICAL r^2-averaging as q[0] and q[1], so the cancellation is exact (machine precision) at t=0. 
  // q[6] and q[7] are never modified by the flux waves, the source update, or the limiter, and never read by prims, so they remain a good frozen reference.
  double gas_gamma = gr_tov->gas_gamma;
  double r2p_l_eq = ql[6], r2p_r_eq = qr[6]; // r^2 * p_eq
  double r2D_l_eq = ql[7], r2D_r_eq = qr[7]; // r^2 * D_eq

  deq[0] = r2D_r_eq - r2D_l_eq; // jump in r^2*D_eq
  deq[1] = (r2p_r_eq - r2p_l_eq) / (gas_gamma - 1.0); // tau_eq = p_eq/(g-1), so r^2*tau_eq = q[6]/(g-1)
  deq[2] = 0.0; // S = 0 at rest; momentum dissipation stays HLL.

  // Equilibrium momentum flux jump, using the SAME face geometry as the actual flux (gkyl_gr_tov_flux_at_radius) so the cancellation is exact at t = 0.
  // The face equilibrium pressure jump is reconstructed from the frozen slots, converting r^2*p_eq back to a face-centered r_face^2 * X * p_eq.
  double lapse = exp(Phi_face);
  double compactness = 0.0;
  if (fabs(r_face) > 1.0e-300) {
    compactness = 2.0 * m_face / r_face;
  }
  double a = 1.0 / sqrt(1.0 - compactness);
  double X = lapse / a;
  double r_l2 = safe_r2(r_l), r_r2 = safe_r2(r_r);
  double p_l_eq = r2p_l_eq / r_l2;
  double p_r_eq = r2p_r_eq / r_r2;
  double r_face2 = safe_r2(r_face);
  *dfeq2 = r_face2 * X * (p_r_eq - p_l_eq);
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
  ghost[4] = -skin[4]; // enclosed mass is odd under the regular r -> -r extension.
  ghost[5] = -skin[5]; // this places the lower interface at r = 0 for shell fluxes.
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
  double Phi_face = 0.5 * (ql[3] + qr[3]);
  double m_face = 0.5 * (ql[4] + qr[4]);
  double fl[8], fr[8];
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, Phi_face, m_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, Phi_face, m_face, fr);

  double *w0 = &waves[0], *w1 = &waves[8];
  for (int i = 0; i < 8; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }
  for (int i = 0; i < 8; i++) {
    if (i != 0 && i != 1 && i != 2) {
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

  double vl[9], vr[9];
  gkyl_gr_tov_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_tov_prim_vars(gas_gamma, qr, vr);


  double rho_l = vl[6];
  double rho_r = vr[6];
  double p_l = vl[8];
  double p_r = vr[8];
  double eps_l = p_l / ((gas_gamma - 1.0) * rho_l);
  double eps_r = p_r / ((gas_gamma - 1.0) * rho_r);
  double h_l = 1.0 + eps_l + p_l / rho_l;
  double h_r = 1.0 + eps_r + p_r / rho_r;
  double cs2_l = gas_gamma * p_l / (rho_l * h_l);
  double cs2_r = gas_gamma * p_r / (rho_r * h_r);
  if (cs2_l < 0.0) {
    cs2_l = 0.0;
  }
  if (cs2_l > 1.0 - 1.0e-12) {
    cs2_l = 1.0 - 1.0e-12;
  }
  if (cs2_r < 0.0) {
    cs2_r = 0.0;
  }
  if (cs2_r > 1.0 - 1.0e-12) {
    cs2_r = 1.0 - 1.0e-12;
  }
  double c_sl = sqrt(cs2_l);
  double c_sr = sqrt(cs2_r);
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
  double Phi_face = 0.5 * (ql[3] + qr[3]);
  double m_face = 0.5 * (ql[4] + qr[4]);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, Phi_face, m_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, Phi_face, m_face, fr);

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
    if (i != 0 && i != 1 && i != 2) {
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

// Well-balanced HLL (Kappeli & Mishra 2014). Identical signal speeds and physical face fluxes as wave_hll; the ONLY change is that the dissipative jump fed to the wave construction is the equilibrium "deviation" delta_q = (qr - ql) - deq rather
// than the raw jump qr - ql. The waves are set directly so that  sl*w0 + sr*w1 = fr - fl   (the physical face-flux jump is unchanged), so amdq + apdq = fr - fl is preserved exactly for any deq - I wanted to make it strictly conservative. 
// deq only redistributes dissipation between amdq/apdq, and vanishes at a static equilibrium so the spurious mass/energy dissipation disappears there.
static double
wave_hll_wb(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;
  double kappa = gr_tov->kappa;

  double vl[9], vr[9];
  gkyl_gr_tov_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_tov_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[6], rho_r = vr[6];
  double p_l = vl[8], p_r = vr[8];
  double eps_l = p_l / ((gas_gamma - 1.0) * rho_l);
  double eps_r = p_r / ((gas_gamma - 1.0) * rho_r);
  double h_l = 1.0 + eps_l + p_l / rho_l;
  double h_r = 1.0 + eps_r + p_r / rho_r;
  double cs2_l = gas_gamma * p_l / (rho_l * h_l);
  double cs2_r = gas_gamma * p_r / (rho_r * h_r);
  if (cs2_l < 0.0) { cs2_l = 0.0; }
  if (cs2_l > 1.0 - 1.0e-12) { cs2_l = 1.0 - 1.0e-12; }
  if (cs2_r < 0.0) { cs2_r = 0.0; }
  if (cs2_r > 1.0 - 1.0e-12) { cs2_r = 1.0 - 1.0e-12; }
  double c_sl = sqrt(cs2_l);
  double c_sr = sqrt(cs2_r);
  double vel_l = vl[7];
  double vel_r = vr[7];

  double lapse_l = exp(ql[3]);
  double lapse_r = exp(qr[3]);
  double a_l = 1.0 / sqrt(1.0 - (2.0 * ql[4] / ql[5]));
  double a_r = 1.0 / sqrt(1.0 - (2.0 * qr[4] / qr[5]));
  double metric_speed_l = lapse_l / a_l;
  double metric_speed_r = lapse_r / a_r;

  double slow_acoustic_eig_l = metric_speed_l * (vel_l - c_sl) / (1.0 - (vel_l * c_sl));
  double fast_acoustic_eig_l = metric_speed_l * (vel_l + c_sl) / (1.0 + (vel_l * c_sl));
  double slow_acoustic_eig_r = metric_speed_r * (vel_r - c_sr) / (1.0 - (vel_r * c_sr));
  double fast_acoustic_eig_r = metric_speed_r * (vel_r + c_sr) / (1.0 + (vel_r * c_sr));

  double sl = fmin(slow_acoustic_eig_l, slow_acoustic_eig_r);
  double sr = fmax(fast_acoustic_eig_l, fast_acoustic_eig_r);
  if (sr <= sl) {
    double amax = fmax(gkyl_gr_tov_max_abs_speed(gas_gamma, ql), gkyl_gr_tov_max_abs_speed(gas_gamma, qr));
    sl = -amax;
    sr = amax;
  }

  double fl[8], fr[8];
  double r_face = 0.5 * (ql[5] + qr[5]);
  double Phi_face = 0.5 * (ql[3] + qr[3]);
  double m_face = 0.5 * (ql[4] + qr[4]);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, Phi_face, m_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, Phi_face, m_face, fr);

  // Equilibrium-subtracted jump: delta_q = (qr - ql) - deq; and the frozen equilibrium momentum flux jump dfeq2 (subtracted from df[2] only).
  double deq[8], dfeq2 = 0.0;
  gr_tov_wb_equilibrium_jump(gr_tov, ql, qr, r_face, Phi_face, m_face, deq, &dfeq2);

  double *w0 = &waves[0], *w1 = &waves[8];
  double inv = 1.0 / (sr - sl);
  for (int i = 0; i < 8; i++) {
    double dq = (qr[i] - ql[i]) - deq[i];
    double df = fr[i] - fl[i];
    if (i == 2) {
      df -= dfeq2; // well-balance the momentum central flux: df[2] - df_eq[2] -> 0 at rest
    }
    w0[i] = (sr * dq - df) * inv; // left-going wave
    w1[i] = (df - sl * dq) * inv; // right-going wave; sl*w0 + sr*w1 = df (- dfeq2 for mom)
  }
  for (int i = 0; i < 8; i++) {
    if (i != 0 && i != 1 && i != 2) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = sl;
  s[1] = sr;

  return fmax(fabs(sl), fabs(sr));
}

static double
wave_hll_wb_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_hll_wb(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }
}

// HLLC Riemann solver (Mignone & Bodo 2005, MNRAS 364:126 for the relativistic
// contact closure; Toro/Spruce/Speares 1994 for the base 3-wave structure).
static double
wave_hllc(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_tov *gr_tov = container_of(eqn, struct wv_gr_tov, eqn);
  double gas_gamma = gr_tov->gas_gamma;
  double kappa = gr_tov->kappa;

  double *w0 = &waves[0], *w1 = &waves[8], *w2 = &waves[16];
  for (int i = 0; i < 24; i++) {
    waves[i] = 0.0;
  }

  // Common interface geometry (shared by both states for telescoping).
  double r_face = 0.5 * (ql[5] + qr[5]);
  double Phi_face = 0.5 * (ql[3] + qr[3]);
  double m_face = 0.5 * (ql[4] + qr[4]);
  double lapse_face = exp(Phi_face);
  double compactness = 0.0;
  if (fabs(r_face) > 1.0e-300) {
    compactness = 2.0 * m_face / r_face;
  }
  double a_face = 1.0 / sqrt(1.0 - compactness);
  double X_face = lapse_face / a_face;

  // Actual coordinate fluxes at the interface (carry X_face * r_face^2).
  double fl[8], fr[8];
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, Phi_face, m_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, Phi_face, m_face, fr);

  // Primitive recovery for SR wave speeds and the contact closure.
  double vl[9], vr[9];
  gkyl_gr_tov_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_tov_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[6], rho_r = vr[6];
  double vel_l = vl[7], vel_r = vr[7];
  double p_l = vl[8], p_r = vr[8];

  double eps_l = p_l / ((gas_gamma - 1.0) * rho_l);
  double eps_r = p_r / ((gas_gamma - 1.0) * rho_r);
  double h_l = 1.0 + eps_l + p_l / rho_l;
  double h_r = 1.0 + eps_r + p_r / rho_r;
  double cs2_l = gas_gamma * p_l / (rho_l * h_l);
  double cs2_r = gas_gamma * p_r / (rho_r * h_r);
  if (cs2_l < 0.0) cs2_l = 0.0;
  if (cs2_l > 1.0 - 1.0e-12) cs2_l = 1.0 - 1.0e-12;
  if (cs2_r < 0.0) cs2_r = 0.0;
  if (cs2_r > 1.0 - 1.0e-12) cs2_r = 1.0 - 1.0e-12;
  double cs_l = sqrt(cs2_l);
  double cs_r = sqrt(cs2_r);

  // SR characteristic speeds (no metric factor).
  double slow_l = (vel_l - cs_l) / (1.0 - vel_l * cs_l);
  double fast_l = (vel_l + cs_l) / (1.0 + vel_l * cs_l);
  double slow_r = (vel_r - cs_r) / (1.0 - vel_r * cs_r);
  double fast_r = (vel_r + cs_r) / (1.0 + vel_r * cs_r);
  double sl_sr = fmin(slow_l, slow_r);
  double sr_sr = fmax(fast_l, fast_r);

  // Degenerate fan: fall back to plain HLL stored in the 3-wave layout.
  bool use_hll = (sr_sr <= sl_sr);

  // HLL averages of total energy E and momentum m, in q-space BUT with only the
  // common frame factor X_face stripped (not r^2). Keeping the shared r^2 factor
  // is essential: at the center face r_face -> 0, so an explicit 1/r^2 strip would
  // overflow. The contact quadratic A*lambda^2 - B*lambda + C = 0 has A, B, C all
  // carrying the same r^2, so lambda* (a velocity) comes out r^2-invariant, while
  // p_star_r2 is an r^2-weighted pressure that cancels cleanly against r2*p in the
  // star states below. Energy flux = momentum (f_E = m); momentum flux = m*v + p.
  double r2 = safe_r2(r_face);
  double Eh_l = ql[0] + ql[1]; // r^2 * (D + tau) = r^2 * E_l
  double Eh_r = qr[0] + qr[1];
  double mh_l = ql[2]; // r^2 * S_l
  double mh_r = qr[2];
  double fEh_l = (fl[0] + fl[1]) / X_face; // r^2 * S_l (energy flux / X)
  double fEh_r = (fr[0] + fr[1]) / X_face;
  double fmh_l = fl[2] / X_face; // r^2 * (S*v + p)_l
  double fmh_r = fr[2] / X_face;

  double lambda_star = 0.0, p_star_r2 = 0.0;
  if (!use_hll) {
    double dsr = sr_sr - sl_sr;
    double Eh_hll  = (sr_sr * Eh_r - sl_sr * Eh_l + fEh_l - fEh_r) / dsr;
    double mh_hll  = (sr_sr * mh_r - sl_sr * mh_l + fmh_l - fmh_r) / dsr;
    double fEh_hll = (sr_sr * fEh_l - sl_sr * fEh_r + sr_sr * sl_sr * (Eh_r - Eh_l)) / dsr;
    double fmh_hll = (sr_sr * fmh_l - sl_sr * fmh_r + sr_sr * sl_sr * (mh_r - mh_l)) / dsr;

    // Contact speed solves fEh_hll*lambda^2 - (Eh_hll + fmh_hll)*lambda + mh_hll = 0.
    double A = fEh_hll;
    double B = Eh_hll + fmh_hll;
    double C = mh_hll;
    if (fabs(A) < 1.0e-300) {
      // Linear degeneracy: B*lambda = C.
      if (fabs(B) < 1.0e-300) {
        use_hll = true;
      }
      else {
        lambda_star = C / B;
        p_star_r2 = fmh_hll - fEh_hll * lambda_star;
      }
    }
    else {
      double disc = B * B - 4.0 * A * C;
      if (disc < 0.0) {
        use_hll = true;
      }
      else {
        lambda_star = (B - sqrt(disc)) / (2.0 * A); // physical (minus) root
        p_star_r2 = fmh_hll - fEh_hll * lambda_star;
      }
    }

    // Consistency / robustness: the contact must lie inside the outer fan,
    // sl <= lambda* <= sr. In the near-uniform atmosphere the fan collapses to a
    // razor-thin sliver (sl_sr ~ sr_sr ~ tiny) and the quadratic above becomes
    // ill-conditioned, returning a contact speed far outside the fan. That single
    // spurious speed would set the global CFL and crater dt. When the 3-wave model
    // is inconsistent like this, fall back to HLL for the face (no contact to
    // resolve there anyway). The negated compound also catches a NaN root.
    if (!(lambda_star >= sl_sr && lambda_star <= sr_sr)) {
      use_hll = true;
    }
  }

  if (!use_hll) {
    double denomL = sl_sr - lambda_star;
    double denomR = sr_sr - lambda_star;
    if (fabs(denomL) < 1.0e-12 || fabs(denomR) < 1.0e-12) {
      use_hll = true;
    }
    else {
      // Left star state (q-space, r^2-weighted), from SR RH across sl.
      // p_star_r2 and r2*p_l both carry the same r2 = safe_r2(r_face).
      double qL0 = ql[0] * (sl_sr - vel_l) / denomL;
      double qL2 = (ql[2] * (sl_sr - vel_l) + (p_star_r2 - r2 * p_l)) / denomL;
      double qLE = (sl_sr * (ql[0] + ql[1]) - ql[2] + lambda_star * p_star_r2) / denomL;
      double qL1 = qLE - qL0;

      // Right star state, from SR RH across sr.
      double qR0 = qr[0] * (sr_sr - vel_r) / denomR;
      double qR2 = (qr[2] * (sr_sr - vel_r) + (p_star_r2 - r2 * p_r)) / denomR;
      double qRE = (sr_sr * (qr[0] + qr[1]) - qr[2] + lambda_star * p_star_r2) / denomR;
      double qR1 = qRE - qR0;

      w0[0] = qL0 - ql[0]; w0[1] = qL1 - ql[1]; w0[2] = qL2 - ql[2];
      w1[0] = qR0 - qL0;   w1[1] = qR1 - qL1;   w1[2] = qR2 - qL2;
      w2[0] = qr[0] - qR0; w2[1] = qr[1] - qR1; w2[2] = qr[2] - qR2;

      // Coordinate wave speeds.
      s[0] = X_face * sl_sr;
      s[1] = X_face * lambda_star;
      s[2] = X_face * sr_sr;

      return fmax(fabs(s[0]), fabs(s[2]));
    }
  }

  // HLL fallback in 3-wave layout: middle wave carries no jump.
  double sl_c, sr_c;
  if (sr_sr <= sl_sr) {
    double amax = fmax(gkyl_gr_tov_max_abs_speed(gas_gamma, ql), gkyl_gr_tov_max_abs_speed(gas_gamma, qr));
    sl_c = -amax;
    sr_c = amax;
  }
  else {
    sl_c = X_face * sl_sr;
    sr_c = X_face * sr_sr;
  }
  double qm[8];
  for (int i = 0; i < 8; i++) {
    qm[i] = ((sr_c * qr[i]) - (sl_c * ql[i]) + (fl[i] - fr[i])) / (sr_c - sl_c);
  }
  for (int i = 0; i < 3; i++) {
    w0[i] = qm[i] - ql[i];
    w1[i] = 0.0;
    w2[i] = qr[i] - qm[i];
  }
  s[0] = sl_c;
  s[1] = 0.5 * (sl_c + sr_c);
  s[2] = sr_c;

  return fmax(fabs(sl_c), fabs(sr_c));
}

static void
qfluct_hllc(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[8], *w2 = &waves[16];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i = 0; i < 8; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]);
  }
}

static double
wave_hllc_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_hllc(eqn, delta, ql, qr, waves, s);
  }
  else {
    // Lax fallback writes only 2 waves; zero the third HLLC slot so the limiter
    // never reads uninitialized data.
    double smax = wave_lax(eqn, delta, ql, qr, waves, s);
    double *w2 = &waves[16];
    for (int i = 0; i < 8; i++) {
      w2[i] = 0.0;
    }
    s[2] = s[1];
    return smax;
  }
}

static void
qfluct_hllc_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_hllc(eqn, ql, qr, waves, s, amdq, apdq);
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
  double r_face = 0.5 * (ql[5] + qr[5]);
  double Phi_face = 0.5 * (ql[3] + qr[3]);
  double m_face = 0.5 * (ql[4] + qr[4]);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, ql, r_face, Phi_face, m_face, fl);
  gkyl_gr_tov_flux_at_radius(gas_gamma, kappa, qr, r_face, Phi_face, m_face, fr);

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

  double v[9] = { 0.0 };
  gkyl_gr_tov_prim_vars(gas_gamma, q, v);

  if (v[6] < 0.0 || v[8] < 0.0 || fabs(v[7]) > 1.0) {
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
  else if (inp->rp_type == WV_GR_TOV_RP_HLLC) {
    gr_tov->eqn.num_waves = 3;
    gr_tov->eqn.waves_func = wave_hllc_l;
    gr_tov->eqn.qfluct_func = qfluct_hllc_l;
  }
  else if (inp->rp_type == WV_GR_TOV_RP_HLL_WB) {
    gr_tov->eqn.num_waves = 2;
    gr_tov->eqn.waves_func = wave_hll_wb_l;
    gr_tov->eqn.qfluct_func = qfluct_hll_l; // waves set directly; same fluctuation formula
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
