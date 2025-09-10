#include <gkyl_advection_kernels.h> 
GKYL_CU_DH double advection_vol_1x_ser_p2(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u[NDIM]:   Advection velocity.
  // f:         Input function.
  // out:       Incremented output.
  const double rdx2 = 2.0/dxv[0]; 
  double cflFreq_mid = 0.0; 
  cflFreq_mid += fabs((1.767766952966369*u[0]-1.976423537605237*u[2])*rdx2); 

  out[1] += 1.224744871391589*(f[2]*u[2]+f[1]*u[1]+f[0]*u[0])*rdx2; 
  out[2] += (2.449489742783178*(f[1]*u[2]+u[1]*f[2])+2.738612787525831*(f[0]*u[1]+u[0]*f[1]))*rdx2; 

  return cflFreq_mid; 
} 
