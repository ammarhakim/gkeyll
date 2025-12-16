#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <float.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_embed_geo.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#define RHOU 1
#define RHOV 2
#define RHOW 3

struct wv_coldfluid {
  struct gkyl_wv_eqn eqn; // base object
};


/**
 * Free Coldfluid eqn object.
 *
 * @param ref Reference counter for Coldfluid eqn
 */
void gkyl_coldfluid_free(const struct gkyl_ref_count *ref);

/**
 * Compute flux. Assumes rotation to local coordinate system.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
GKYL_CU_DH
static void
gkyl_coldfluid_flux(const double q[4], double flux[4])
{
  double u = q[1]/q[0];
  flux[0] = q[1]; // rho*u
  flux[1] = q[1]*u; // rho*u*u 
  flux[2] = q[2]*u; // rho*v*u
  flux[3] = q[3]*u; // rho*w*u
}

GKYL_CU_DH
static inline void
coldfluid_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // density and moment as copied as-is
  for (int i=0; i<4; ++i) diag[i] = qin[i];
  double ke = 0.5*(qin[1]*qin[1] + qin[2]*qin[2] + qin[3]*qin[3])/qin[0];
  diag[4] = ke;
}

GKYL_CU_DH
static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    wout[i] = qin[i];
}

GKYL_CU_DH
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    qout[i] = win[i];
}

// Cold fluid perfectly reflecting wall
GKYL_CU_DH
static void
coldfluid_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density 
  ghost[0] = skin[0];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = skin[2];
  ghost[3] = skin[3];
}

// Cold fluid no-slip wall
GKYL_CU_DH
static void
coldfluid_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density 
  ghost[0] = skin[0];

  // zero momentum
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];
  ghost[3] = -skin[3];
}

GKYL_CU_DH
static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
}

GKYL_CU_DH
static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
}

GKYL_CU_DH
static void
wave_embed_absorb(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = DBL_EPSILON;
  qphi[1] = 0.0;
  qphi[2] = 0.0;
  qphi[3] = 0.0;
}

GKYL_CU_DH
static void
wave_embed_reflect(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = q[0];
  qphi[1] = -q[1];
  qphi[2] = q[2];
  qphi[3] = q[3];
}

GKYL_CU_DH
static double
wave_embedded(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, const double phil,
  const double phir, double *waves, double *s)
{
  const struct wv_coldfluid *coldfluid = container_of(eqn, struct wv_coldfluid, eqn);
  double qphi[4] = {0.0};
  double deltaphi[4] = {0.0};
  double amax = 0.0;

  double fl[4], fr[4] = {0.0};
  if ((phil < 0.0) && (phir > 0.0)) {
    double rhor = qr[0];
    double ur = qr[1]/qr[0];
    double sr = fabs(ur);
    amax = sr;
   
    eqn->embed_geo->embed_func(qr, qphi, deltaphi, eqn->embed_geo->ctx);

    gkyl_coldfluid_flux(qphi, fl);
    gkyl_coldfluid_flux(qr, fr);

    double *w0 = &waves[0], *w1 = &waves[4];
    if (amax > 0.0) {
      for (int i=0; i<4; ++i) {
        w1[i] = 0.5*((qr[i]-qphi[i]) + (fr[i]-fl[i])/amax);
      }
    }
    else {
      for (int i=0; i<4; ++i) {
        w1[i] = 0.0;
      }      
    }
  }
  
  if ((phir < 0.0) && (phil > 0.0)) {
    double rhol = ql[0];
    double ul = ql[1]/ql[0];
    double sl = fabs(ul);
    amax = sl;
    
    eqn->embed_geo->embed_func(ql, qphi, deltaphi, eqn->embed_geo->ctx);

    gkyl_coldfluid_flux(ql, fl);
    gkyl_coldfluid_flux(qphi, fr);
    
    double *w0 = &waves[0], *w1 = &waves[4];
    if (amax > 0.0) {
      for (int i=0; i<4; ++i) {
        w0[i] = 0.5*((qphi[i]-ql[i]) - (fr[i]-fl[i])/amax);
      }
    }
    else {
      for (int i=0; i<4; ++i) {
        w0[i] = 0.0;
      }      
    }
  }

  if (amax > 0.0) {
    s[0] = -amax;
    s[1] = amax;
  }
  else {
    s[0] = 0.0;
    s[1] = 0.0;
  }

  return s[1];
}

// Waves and speeds using Lax fluxes
GKYL_CU_DH
static double
wave_lax(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_coldfluid *coldfluid = container_of(eqn, struct wv_coldfluid, eqn);

  double rhol = ql[0], rhor = qr[0];
  double ul = ql[1]/ql[0], ur = qr[1]/qr[0];
  double sl = fabs(ul), sr = fabs(ur);
  double amax = fmax(sl, sr);

  double fl[4], fr[4];
  gkyl_coldfluid_flux(ql, fl);
  gkyl_coldfluid_flux(qr, fr);

  double *w0 = &waves[0], *w1 = &waves[4];
  if (amax > 0.0) {
    for (int i=0; i<4; ++i) {
      w0[i] = 0.5*((qr[i]-ql[i]) - (fr[i]-fl[i])/amax);
      w1[i] = 0.5*((qr[i]-ql[i]) + (fr[i]-fl[i])/amax);
    }
    s[0] = -amax;
    s[1] = amax;
  }
  else {
    for (int i=0; i<4; ++i) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
    s[0] = 0.0;
    s[1] = 0.0;
  }
  return s[1];
}

GKYL_CU_DH
static void
qfluct_lax(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[4];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<4; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

GKYL_CU_DH
static double
wave_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double phil,
  double phir, double *waves, double *s)
{
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);
}

GKYL_CU_DH
static void
qfluct_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, double phil,
  double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// Waves and speeds using Roe averaging
GKYL_CU_DH
static double
wave_roe(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  double f[4];
  double ur = qr[RHOU]/qr[0], ul = ql[RHOU]/ql[0];

  double *wv = 0;

  if ((ul < 0) && (0 < ur)) { // vacuum intermediate state will be formed
    gkyl_coldfluid_flux(ql, f);
    wv = &waves[0];
    for(int m=0; m<4; ++m) wv[m] = -f[m];
    s[0] = ul;

    gkyl_coldfluid_flux(qr, f);
    wv = &waves[4];
    for(int m=0; m<4; ++m) wv[m] = f[m];
    s[1] = ur;
  }
  else {
    // no vacuum state
    double rl = ql[0];
    double rr = qr[0];
    // compute Roe averaged speed
    double uav = (sqrt(rl)*ul + sqrt(rr)*ur)/(sqrt(rl)+sqrt(rr));
            
    if(uav<0) {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = 0.0;
    }
    else {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = 0;

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];
    }
    s[0] = uav;
    s[1] = uav;
  }

  return fmax(fabs(s[0]), fabs(s[1]));
}

GKYL_CU_DH
static void
qfluct_roe(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[4];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<4; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

GKYL_CU_DH
static double
wave_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double phil,
  double phir, double *waves, double *s)
{
  // clear waves and wave speeds
  int mwaves = (type == GKYL_WV_HIGH_ORDER_FLUX) ? eqn->num_waves : 2;
  int meqn = eqn->num_equations;
  for (int i=0; i<mwaves; ++i) {
    double *w = &waves[i*meqn];
    for (int j=0; j<meqn; ++j)
      w[j] = 0.0;
    s[i] = 0.0;
  }
  // Check if inside an embedded surface
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else {
    if (type == GKYL_WV_HIGH_ORDER_FLUX)
      return wave_roe(eqn, delta, ql, qr, waves, s);
    else
      return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // can't happen
}

GKYL_CU_DH
static void
qfluct_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, double phil, double phir, const double *waves,
  const double *s, double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX && (phil > 0.0) && (phir > 0.0))
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

GKYL_CU_DH
static void
ffluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = 4, mwaves = 2;
  
  for (int m=0; m<meqn; ++m) {
    amdq[m] = 0.0; apdq[m] = 0.0;

    for (int mw=0; mw<mwaves; ++mw) {
      const double *wv = &waves[mw*meqn];
      
      if (s[mw] < 0.0) {
        amdq[m] += wv[m];
      }
      else if (s[mw] > 0.0) {
        apdq[m] += wv[m];
      }
      else {
        amdq[m] += 0.5*wv[m];
        apdq[m] += 0.5*wv[m];
      }
    }
  }
}

GKYL_CU_DH
static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  double fr[4], fl[4];
  gkyl_coldfluid_flux(ql, fl);
  gkyl_coldfluid_flux(qr, fr);

  for (int m=0; m<4; ++m) flux_jump[m] = fr[m]-fl[m];

  double amaxl = ql[RHOU]/ql[0];
  double amaxr = qr[RHOU]/qr[0];

  return fmax(amaxl, amaxr);
}

GKYL_CU_DH
static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return q[0] > 0.0;
}

GKYL_CU_DH
static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_coldfluid *coldfluid = container_of(eqn, struct wv_coldfluid, eqn);
  return fabs(q[RHOU]/q[0]);
}

GKYL_CU_DH
static inline void
coldfluid_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 4; i++) {
    sout[i] = 0.0;
  }
}
