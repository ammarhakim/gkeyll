#include <gkyl_conf_diffusion_kernels.h>

GKYL_CU_DH double conf_diffusion_vol_1x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates.
  // dx[NDIM]: Cell spacing.
  // K: Full tensor K^{ij}=J D g^{ij}, stored by tensor component.
  // jacobgeo_inv: 1/J, used in the CFL estimate.
  // f: Physical scalar field N/J.
  // out: Incremented output.

  const double rdx20 = 2.0/dx[0];
  const double rdx2sq0 = pow(2.0/dx[0],2.0);

  out[1] += -(2.1213203435596424*K[0]*f[1]*rdx2sq0); 

  return 2.0*K[0]*jacobgeo_inv[0]*rdx2sq0;
}
