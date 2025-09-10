#include <gkyl_canonical_pb_kernels.h> 
double canonical_pb_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxdvInv0 = 1.0/(dxv[0]*dxv[1]); 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  double dvInv1 = 1.0/dxv[2]; 
  out[1] += 4.242640687119286*(f[5]*hamil[7]+f[3]*hamil[6]+f[1]*hamil[4]+f[0]*hamil[2])*dxdvInv0; 
  out[2] += -4.242640687119286*(f[6]*hamil[7]+f[3]*hamil[5]+f[2]*hamil[4]+f[0]*hamil[1])*dxdvInv0; 
  out[4] += 4.242640687119286*(f[6]*hamil[6]-1.0*f[5]*hamil[5]+f[2]*hamil[2]-1.0*f[1]*hamil[1])*dxdvInv0; 
  out[5] += 4.242640687119286*(f[1]*hamil[7]+f[0]*hamil[6]+hamil[4]*f[5]+hamil[2]*f[3])*dxdvInv0; 
  out[6] += -4.242640687119286*(f[2]*hamil[7]+hamil[4]*f[6]+f[0]*hamil[5]+hamil[1]*f[3])*dxdvInv0; 
  out[7] += 4.242640687119286*(f[2]*hamil[6]+hamil[2]*f[6]-1.0*f[1]*hamil[5]-1.0*hamil[1]*f[5])*dxdvInv0; 
  return 0.; 
} 
