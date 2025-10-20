#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_vol_1x_ser_p2_constcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 81.0*coeff[0]*pow(rdx2, 4.0); 
}

GKYL_CU_DH double dg_diffusion_vlasov_order4_vol_1x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order4_vol_1x_ser_p2_constcoeff_diffx(w, dx, coeff, q, out);

  return cflFreq;
}

