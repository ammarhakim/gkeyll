#include <gkyl_canonical_pb_kernels.h> 
double canonical_pb_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxdvInv0 = 1.0/(dxv[0]*dxv[1]); 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 6.0*(f[1]*hamil[3]+f[0]*hamil[2])*dxdvInv0; 
  out[2] += -6.0*(f[2]*hamil[3]+f[0]*hamil[1])*dxdvInv0; 
  out[3] += 6.0*(f[2]*hamil[2]-1.0*f[1]*hamil[1])*dxdvInv0; 
  return 0.; 
} 
