#include <gkyl_gk_numerical_diffusion_kernels.h>

static GKYL_CU_DH double gk_numerical_diffusion_order4_vol_1x1v_ser_p1_varnu_diffx(const double *w, const double *dx, const double *nu, const double *fin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity.
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // fin: Input field times conf-space Jacobian
  // out: Incremented output

  const double rdx2fac = pow(2.0/dx[0],4); 

  return 2.8284271247461907*nu[0]*rdx2fac; 
}

GKYL_CU_DH double gk_numerical_diffusion_order4_vol_1x1v_ser_p1_varnu_diffdirsx(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[6];
  fin[0] = 0.7071067811865476*(Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.7071067811865476*(Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.7071067811865476*(jacobgeo_inv[1]*Jfin[3]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.7071067811865476*(jacobgeo_inv[0]*Jfin[3]+jacobgeo_inv[1]*Jfin[2]); 
  fin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*Jfin[5]+21.213203435596427*jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_inv[0]*Jfin[5]+21.21320343559643*jacobgeo_inv[1]*Jfin[4]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order4_vol_1x1v_ser_p1_varnu_diffx(w, dx, nu, fin, out);

  return cflFreq;
}

