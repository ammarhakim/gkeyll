#include <gkyl_conf_diffusion_kernels.h>

GKYL_CU_DH double conf_diffusion_vol_2x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates.
  // dx[NDIM]: Cell spacing.
  // K: Full tensor K^{ij}=J D g^{ij}, stored by tensor component.
  // jacobgeo_inv: 1/J, used in the CFL estimate.
  // f: Physical scalar field N/J.
  // out: Incremented output.

  const double rdx20 = 2.0/dx[0];
  const double rdx2sq0 = pow(2.0/dx[0],2.0);
  const double rdx21 = 2.0/dx[1];
  const double rdx2sq1 = pow(2.0/dx[1],2.0);

  out[1] += -(1.5*K[2]*f[3]*rdx2sq0)-1.5*K[0]*f[1]*rdx2sq0-1.5*f[3]*K[5]*rdx20*rdx21-1.5*f[2]*K[4]*rdx20*rdx21; 
  out[2] += -(1.5*f[3]*K[13]*rdx2sq1)-1.5*f[2]*K[12]*rdx2sq1-1.5*f[3]*K[10]*rdx20*rdx21-1.5*f[1]*K[8]*rdx20*rdx21; 
  out[3] += -(1.5*f[2]*K[13]*rdx2sq1)-1.5*f[3]*K[12]*rdx2sq1-1.5*K[0]*f[3]*rdx2sq0-1.5*f[1]*K[2]*rdx2sq0-1.5*f[3]*K[11]*rdx20*rdx21-1.5*f[1]*K[9]*rdx20*rdx21-1.5*f[3]*K[7]*rdx20*rdx21-1.5*f[2]*K[6]*rdx20*rdx21; 

  return jacobgeo_inv[0]*K[12]*rdx2sq1+K[0]*jacobgeo_inv[0]*rdx2sq0;
}
