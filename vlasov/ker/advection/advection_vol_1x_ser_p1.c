#include <gkyl_advection_kernels.h> 
GKYL_CU_DH double advection_vol_1x_ser_p1(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u[NDIM]:   Advection velocity.
  // f:         Input function.
  // out:       Incremented output.
  const double rdx2 = 2.0/dxv[0]; 
  double cflFreq_mid = 0.0; 
  cflFreq_mid += fabs(1.060660171779821*u[0]*rdx2); 

  out[1] += 1.224744871391589*(f[1]*u[1]+f[0]*u[0])*rdx2; 

  return cflFreq_mid; 
} 
