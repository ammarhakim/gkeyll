#include <gkyl_advection_kernels.h> 
GKYL_CU_DH double advection_vol_3x_ser_p1(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u[NDIM]:   Advection velocity.
  // f:         Input function.
  // out:       Incremented output.
  const double rdx2 = 2.0/dxv[0]; 
  const double rdy2 = 2.0/dxv[1]; 
  const double rdz2 = 2.0/dxv[2]; 
  double cflFreq_mid = 0.0; 
  cflFreq_mid += fabs(0.5303300858899105*u[0]*rdx2); 
  cflFreq_mid += fabs(0.5303300858899105*u[8]*rdy2); 
  cflFreq_mid += fabs(0.5303300858899105*u[16]*rdz2); 

  out[1] += 0.6123724356957944*(f[7]*u[7]+f[6]*u[6]+f[5]*u[5]+f[4]*u[4]+f[3]*u[3]+f[2]*u[2]+f[1]*u[1]+f[0]*u[0])*rdx2; 
  out[2] += 0.6123724356957944*(f[7]*u[15]+f[6]*u[14]+f[5]*u[13]+f[4]*u[12]+f[3]*u[11]+f[2]*u[10]+f[1]*u[9]+f[0]*u[8])*rdy2; 
  out[3] += 0.6123724356957944*(f[7]*u[23]+f[6]*u[22]+f[5]*u[21]+f[4]*u[20]+f[3]*u[19]+f[2]*u[18]+f[1]*u[17]+f[0]*u[16])*rdz2; 
  out[4] += 0.6123724356957944*((f[6]*u[15]+f[7]*u[14]+f[3]*u[13]+f[2]*u[12]+f[5]*u[11]+f[4]*u[10]+f[0]*u[9]+f[1]*u[8])*rdy2+(f[5]*u[7]+u[5]*f[7]+f[3]*u[6]+u[3]*f[6]+f[1]*u[4]+u[1]*f[4]+f[0]*u[2]+u[0]*f[2])*rdx2); 
  out[5] += 0.6123724356957944*((f[6]*u[23]+f[7]*u[22]+f[3]*u[21]+f[2]*u[20]+f[5]*u[19]+f[4]*u[18]+f[0]*u[17]+f[1]*u[16])*rdz2+(f[4]*u[7]+u[4]*f[7]+f[2]*u[6]+u[2]*f[6]+f[1]*u[5]+u[1]*f[5]+f[0]*u[3]+u[0]*f[3])*rdx2); 
  out[6] += 0.6123724356957944*((f[5]*u[23]+f[3]*u[22]+f[7]*u[21]+f[1]*u[20]+f[6]*u[19]+f[0]*u[18]+f[4]*u[17]+f[2]*u[16])*rdz2+(f[4]*u[15]+f[2]*u[14]+f[1]*u[13]+f[7]*u[12]+f[0]*u[11]+f[6]*u[10]+f[5]*u[9]+f[3]*u[8])*rdy2); 
  out[7] += 0.6123724356957944*((f[3]*u[23]+f[5]*u[22]+f[6]*u[21]+f[0]*u[20]+f[7]*u[19]+f[1]*u[18]+f[2]*u[17]+f[4]*u[16])*rdz2+(f[2]*u[15]+f[4]*u[14]+f[0]*u[13]+f[6]*u[12]+f[1]*u[11]+f[7]*u[10]+f[3]*u[9]+f[5]*u[8])*rdy2+(f[1]*u[7]+u[1]*f[7]+f[0]*u[6]+u[0]*f[6]+f[4]*u[5]+u[4]*f[5]+f[2]*u[3]+u[2]*f[3])*rdx2); 

  return cflFreq_mid; 
} 
