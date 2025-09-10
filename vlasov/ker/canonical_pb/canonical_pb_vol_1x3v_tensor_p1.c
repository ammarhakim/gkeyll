#include <gkyl_canonical_pb_kernels.h> 
double canonical_pb_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxdvInv0 = 1.0/(dxv[0]*dxv[1]); 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  double dvInv1 = 1.0/dxv[2]; 
  double dvInv2 = 1.0/dxv[3]; 
  out[1] += 3.0*(f[13]*hamil[15]+f[10]*hamil[14]+f[8]*hamil[12]+f[6]*hamil[11]+f[4]*hamil[9]+f[3]*hamil[7]+f[1]*hamil[5]+f[0]*hamil[2])*dxdvInv0; 
  out[2] += -3.0*(f[14]*hamil[15]+f[10]*hamil[13]+f[9]*hamil[12]+f[7]*hamil[11]+f[4]*hamil[8]+f[3]*hamil[6]+f[2]*hamil[5]+f[0]*hamil[1])*dxdvInv0; 
  out[5] += 3.0*(f[14]*hamil[14]-1.0*f[13]*hamil[13]+f[9]*hamil[9]-1.0*f[8]*hamil[8]+f[7]*hamil[7]-1.0*f[6]*hamil[6]+f[2]*hamil[2]-1.0*f[1]*hamil[1])*dxdvInv0; 
  out[6] += 3.0*(f[8]*hamil[15]+f[4]*hamil[14]+hamil[12]*f[13]+f[1]*hamil[11]+hamil[9]*f[10]+f[0]*hamil[7]+hamil[5]*f[6]+hamil[2]*f[3])*dxdvInv0; 
  out[7] += -3.0*(f[9]*hamil[15]+hamil[12]*f[14]+f[4]*hamil[13]+f[2]*hamil[11]+hamil[8]*f[10]+hamil[5]*f[7]+f[0]*hamil[6]+hamil[1]*f[3])*dxdvInv0; 
  out[8] += 3.0*(f[6]*hamil[15]+f[3]*hamil[14]+hamil[11]*f[13]+f[1]*hamil[12]+hamil[7]*f[10]+f[0]*hamil[9]+hamil[5]*f[8]+hamil[2]*f[4])*dxdvInv0; 
  out[9] += -3.0*(f[7]*hamil[15]+hamil[11]*f[14]+f[3]*hamil[13]+f[2]*hamil[12]+hamil[6]*f[10]+hamil[5]*f[9]+f[0]*hamil[8]+hamil[1]*f[4])*dxdvInv0; 
  out[11] += 3.0*(f[9]*hamil[14]+hamil[9]*f[14]-1.0*f[8]*hamil[13]-1.0*hamil[8]*f[13]+f[2]*hamil[7]+hamil[2]*f[7]-1.0*f[1]*hamil[6]-1.0*hamil[1]*f[6])*dxdvInv0; 
  out[12] += 3.0*(f[7]*hamil[14]+hamil[7]*f[14]-1.0*f[6]*hamil[13]-1.0*hamil[6]*f[13]+f[2]*hamil[9]+hamil[2]*f[9]-1.0*f[1]*hamil[8]-1.0*hamil[1]*f[8])*dxdvInv0; 
  out[13] += 3.0*(f[1]*hamil[15]+f[0]*hamil[14]+hamil[5]*f[13]+f[6]*hamil[12]+f[8]*hamil[11]+hamil[2]*f[10]+f[3]*hamil[9]+f[4]*hamil[7])*dxdvInv0; 
  out[14] += -3.0*(f[2]*hamil[15]+hamil[5]*f[14]+f[0]*hamil[13]+f[7]*hamil[12]+f[9]*hamil[11]+hamil[1]*f[10]+f[3]*hamil[8]+f[4]*hamil[6])*dxdvInv0; 
  out[15] += 3.0*(f[2]*hamil[14]+hamil[2]*f[14]-1.0*f[1]*hamil[13]-1.0*hamil[1]*f[13]+f[7]*hamil[9]+hamil[7]*f[9]-1.0*f[6]*hamil[8]-1.0*hamil[6]*f[8])*dxdvInv0; 
  return 0.; 
} 
