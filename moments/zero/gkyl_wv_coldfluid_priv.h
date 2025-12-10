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


struct wv_coldfluid {
  struct gkyl_wv_eqn eqn; // base object
};


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
    for (int i=0; i<4; ++i) {
      w1[i] = 0.5*((qr[i]-qphi[i]) + (fr[i]-fl[i])/amax);
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
    for (int i=0; i<4; ++i) {
      w0[i] = 0.5*((qphi[i]-ql[i]) - (fr[i]-fl[i])/amax);
    }
  }

  s[0] = -amax;
  s[1] = amax;
  
  return s[1];
}
