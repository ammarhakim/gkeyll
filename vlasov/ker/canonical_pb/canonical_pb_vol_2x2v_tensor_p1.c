#include <gkyl_canonical_pb_kernels.h> 
double canonical_pb_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *hamil, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxdvInv0 = 1.0/(dxv[0]*dxv[2]); 
  double dxdvInv1 = 1.0/(dxv[1]*dxv[3]); 
  double dxInv0 = 1.0/dxv[0]; 
  double dxInv1 = 1.0/dxv[1]; 
  double dvInv0 = 1.0/dxv[2]; 
  double dvInv1 = 1.0/dxv[3]; 
  out[1] += 3.0*(f[12]*hamil[15]+f[9]*hamil[14]+f[8]*hamil[13]+f[5]*hamil[11]+f[4]*hamil[10]+f[2]*hamil[7]+f[1]*hamil[6]+f[0]*hamil[3])*dxdvInv0; 
  out[2] += 3.0*(f[11]*hamil[15]+f[7]*hamil[14]+f[6]*hamil[13]+f[5]*hamil[12]+f[3]*hamil[10]+f[2]*hamil[9]+f[1]*hamil[8]+f[0]*hamil[4])*dxdvInv1; 
  out[3] += -3.0*(f[14]*hamil[15]+f[10]*hamil[13]+f[9]*hamil[12]+f[7]*hamil[11]+f[4]*hamil[8]+f[3]*hamil[6]+f[2]*hamil[5]+f[0]*hamil[1])*dxdvInv0; 
  out[4] += -3.0*(f[13]*hamil[15]+f[10]*hamil[14]+f[8]*hamil[12]+f[6]*hamil[11]+f[4]*hamil[9]+f[3]*hamil[7]+f[1]*hamil[5]+f[0]*hamil[2])*dxdvInv1; 
  out[5] += 3.0*(f[7]*hamil[15]*dxdvInv1+f[11]*hamil[14]*dxdvInv1+f[3]*hamil[13]*dxdvInv1+f[2]*hamil[12]*dxdvInv1+f[6]*hamil[10]*dxdvInv1+f[5]*hamil[9]*dxdvInv1+f[0]*hamil[8]*dxdvInv1+f[1]*hamil[4]*dxdvInv1+f[8]*hamil[15]*dxdvInv0+f[4]*hamil[14]*dxdvInv0+f[12]*hamil[13]*dxdvInv0+f[1]*hamil[11]*dxdvInv0+f[9]*hamil[10]*dxdvInv0+f[0]*hamil[7]*dxdvInv0+f[5]*hamil[6]*dxdvInv0+f[2]*hamil[3]*dxdvInv0); 
  out[6] += 3.0*(f[14]*hamil[14]-1.0*f[12]*hamil[12]+f[10]*hamil[10]-1.0*f[8]*hamil[8]+f[7]*hamil[7]-1.0*f[5]*hamil[5]+f[3]*hamil[3]-1.0*f[1]*hamil[1])*dxdvInv0; 
  out[7] += 3.0*(f[5]*hamil[15]*dxdvInv1+f[2]*hamil[14]*dxdvInv1+f[1]*hamil[13]*dxdvInv1+f[11]*hamil[12]*dxdvInv1+f[0]*hamil[10]*dxdvInv1+f[7]*hamil[9]*dxdvInv1+f[6]*hamil[8]*dxdvInv1+f[3]*hamil[4]*dxdvInv1-1.0*f[10]*hamil[15]*dxdvInv0-1.0*hamil[13]*f[14]*dxdvInv0-1.0*f[4]*hamil[12]*dxdvInv0-1.0*f[3]*hamil[11]*dxdvInv0-1.0*hamil[8]*f[9]*dxdvInv0-1.0*hamil[6]*f[7]*dxdvInv0-1.0*f[0]*hamil[5]*dxdvInv0-1.0*hamil[1]*f[2]*dxdvInv0); 
  out[8] += -3.0*(f[10]*hamil[15]*dxdvInv1+f[13]*hamil[14]*dxdvInv1+f[4]*hamil[12]*dxdvInv1+f[3]*hamil[11]*dxdvInv1+f[8]*hamil[9]*dxdvInv1+f[6]*hamil[7]*dxdvInv1+f[0]*hamil[5]*dxdvInv1+f[1]*hamil[2]*dxdvInv1-1.0*f[5]*hamil[15]*dxdvInv0-1.0*f[2]*hamil[14]*dxdvInv0-1.0*f[1]*hamil[13]*dxdvInv0-1.0*hamil[11]*f[12]*dxdvInv0-1.0*f[0]*hamil[10]*dxdvInv0-1.0*hamil[7]*f[9]*dxdvInv0-1.0*hamil[6]*f[8]*dxdvInv0-1.0*hamil[3]*f[4]*dxdvInv0); 
  out[9] += 3.0*(f[13]*hamil[13]-1.0*f[11]*hamil[11]+f[10]*hamil[10]+f[8]*hamil[8]-1.0*f[7]*hamil[7]-1.0*f[5]*hamil[5]+f[4]*hamil[4]-1.0*f[2]*hamil[2])*dxdvInv1; 
  out[10] += -3.0*(f[8]*hamil[15]*dxdvInv1+f[4]*hamil[14]*dxdvInv1+hamil[12]*f[13]*dxdvInv1+f[1]*hamil[11]*dxdvInv1+hamil[9]*f[10]*dxdvInv1+f[0]*hamil[7]*dxdvInv1+hamil[5]*f[6]*dxdvInv1+hamil[2]*f[3]*dxdvInv1+f[7]*hamil[15]*dxdvInv0+hamil[11]*f[14]*dxdvInv0+f[3]*hamil[13]*dxdvInv0+f[2]*hamil[12]*dxdvInv0+hamil[6]*f[10]*dxdvInv0+hamil[5]*f[9]*dxdvInv0+f[0]*hamil[8]*dxdvInv0+hamil[1]*f[4]*dxdvInv0); 
  out[11] += 3.0*(f[2]*hamil[15]*dxdvInv1+f[5]*hamil[14]*dxdvInv1+f[0]*hamil[13]*dxdvInv1+f[7]*hamil[12]*dxdvInv1+hamil[9]*f[11]*dxdvInv1+f[1]*hamil[10]*dxdvInv1+f[3]*hamil[8]*dxdvInv1+hamil[4]*f[6]*dxdvInv1+f[10]*hamil[14]*dxdvInv0+hamil[10]*f[14]*dxdvInv0-1.0*f[8]*hamil[12]*dxdvInv0-1.0*hamil[8]*f[12]*dxdvInv0+f[3]*hamil[7]*dxdvInv0+hamil[3]*f[7]*dxdvInv0-1.0*f[1]*hamil[5]*dxdvInv0-1.0*hamil[1]*f[5]*dxdvInv0); 
  out[12] += 3.0*(f[10]*hamil[13]*dxdvInv1+hamil[10]*f[13]*dxdvInv1-1.0*f[7]*hamil[11]*dxdvInv1-1.0*hamil[7]*f[11]*dxdvInv1+f[4]*hamil[8]*dxdvInv1+hamil[4]*f[8]*dxdvInv1-1.0*f[2]*hamil[5]*dxdvInv1-1.0*hamil[2]*f[5]*dxdvInv1+f[1]*hamil[15]*dxdvInv0+f[0]*hamil[14]*dxdvInv0+f[5]*hamil[13]*dxdvInv0+hamil[6]*f[12]*dxdvInv0+f[8]*hamil[11]*dxdvInv0+f[2]*hamil[10]*dxdvInv0+hamil[3]*f[9]*dxdvInv0+f[4]*hamil[7]*dxdvInv0); 
  out[13] += -3.0*(f[4]*hamil[15]*dxdvInv1+f[8]*hamil[14]*dxdvInv1+hamil[9]*f[13]*dxdvInv1+f[10]*hamil[12]*dxdvInv1+f[0]*hamil[11]*dxdvInv1+f[1]*hamil[7]*dxdvInv1+hamil[2]*f[6]*dxdvInv1+f[3]*hamil[5]*dxdvInv1-1.0*f[7]*hamil[14]*dxdvInv0-1.0*hamil[7]*f[14]*dxdvInv0+f[5]*hamil[12]*dxdvInv0+hamil[5]*f[12]*dxdvInv0-1.0*f[3]*hamil[10]*dxdvInv0-1.0*hamil[3]*f[10]*dxdvInv0+f[1]*hamil[8]*dxdvInv0+hamil[1]*f[8]*dxdvInv0); 
  out[14] += 3.0*(f[8]*hamil[13]*dxdvInv1+hamil[8]*f[13]*dxdvInv1-1.0*f[5]*hamil[11]*dxdvInv1-1.0*hamil[5]*f[11]*dxdvInv1+f[4]*hamil[10]*dxdvInv1+hamil[4]*f[10]*dxdvInv1-1.0*f[2]*hamil[7]*dxdvInv1-1.0*hamil[2]*f[7]*dxdvInv1-1.0*f[3]*hamil[15]*dxdvInv0-1.0*hamil[6]*f[14]*dxdvInv0-1.0*f[7]*hamil[13]*dxdvInv0-1.0*f[0]*hamil[12]*dxdvInv0-1.0*f[10]*hamil[11]*dxdvInv0-1.0*hamil[1]*f[9]*dxdvInv0-1.0*f[2]*hamil[8]*dxdvInv0-1.0*f[4]*hamil[5]*dxdvInv0); 
  out[15] += 3.0*(f[4]*hamil[13]*dxdvInv1+hamil[4]*f[13]*dxdvInv1-1.0*f[2]*hamil[11]*dxdvInv1-1.0*hamil[2]*f[11]*dxdvInv1+f[8]*hamil[10]*dxdvInv1+hamil[8]*f[10]*dxdvInv1-1.0*f[5]*hamil[7]*dxdvInv1-1.0*hamil[5]*f[7]*dxdvInv1+f[3]*hamil[14]*dxdvInv0+hamil[3]*f[14]*dxdvInv0-1.0*f[1]*hamil[12]*dxdvInv0-1.0*hamil[1]*f[12]*dxdvInv0+f[7]*hamil[10]*dxdvInv0+hamil[7]*f[10]*dxdvInv0-1.0*f[5]*hamil[8]*dxdvInv0-1.0*hamil[5]*f[8]*dxdvInv0); 
  return 0.; 
} 
