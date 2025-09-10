#include <gkyl_canonical_pb_kernels.h> 
double canonical_pb_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxdvInv0 = 1.0/(dxv[0]*dxv[1]); 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 2.0*(6.708203932499369*f[3]*hamil[5]+6.708203932499369*f[2]*hamil[4]+3.0*f[1]*hamil[3]+3.0*f[0]*hamil[2])*dxdvInv0; 
  out[2] += -0.4*(15.0*f[4]*hamil[5]+15.0*f[2]*hamil[3]+15.0*f[0]*hamil[1])*dxdvInv0; 
  out[3] += 2.0*(3.0*f[5]*hamil[5]+6.708203932499369*f[1]*hamil[5]+6.0*f[4]*hamil[4]+6.708203932499369*f[0]*hamil[4]+3.0*f[2]*hamil[2]-3.0*f[1]*hamil[1])*dxdvInv0; 
  out[4] += -0.4*(30.0*f[2]*hamil[5]+30.0*hamil[3]*f[4]+33.54101966249685*f[0]*hamil[3]+33.54101966249685*hamil[1]*f[2])*dxdvInv0; 
  out[5] += -0.4*(15.0*hamil[3]*f[5]-30.0*f[2]*hamil[4]-15.0*hamil[2]*f[4]+33.54101966249684*f[1]*hamil[3]+33.54101966249684*hamil[1]*f[3])*dxdvInv0; 
  return 0.; 
} 
