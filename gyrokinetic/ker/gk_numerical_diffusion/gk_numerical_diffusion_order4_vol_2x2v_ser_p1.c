#include <gkyl_gk_numerical_diffusion_kernels.h>

static GKYL_CU_DH double gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffx(const double *w, const double *dx, const double *nu, const double *fin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity.
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // fin: Input field times conf-space Jacobian
  // out: Incremented output

  const double rdx2fac = pow(2.0/dx[0],4); 

  return 2.0*nu[0]*rdx2fac; 
}

static GKYL_CU_DH double gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffz(const double *w, const double *dx, const double *nu, const double *fin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity.
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // fin: Input field times conf-space Jacobian
  // out: Incremented output

  const double rdx2fac = pow(2.0/dx[1],4); 

  return 2.0*nu[4]*rdx2fac; 
}

GKYL_CU_DH double gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsx(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[24];
  fin[0] = 0.5*(jacobgeo_inv[3]*Jfin[5]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.5*(jacobgeo_inv[2]*Jfin[5]+Jfin[2]*jacobgeo_inv[3]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.5*(jacobgeo_inv[1]*Jfin[5]+Jfin[1]*jacobgeo_inv[3]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.5*(jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[7]+jacobgeo_inv[1]*Jfin[6]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.5*(jacobgeo_inv[3]*Jfin[12]+jacobgeo_inv[2]*Jfin[9]+jacobgeo_inv[1]*Jfin[8]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.5*(jacobgeo_inv[0]*Jfin[5]+Jfin[0]*jacobgeo_inv[3]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[6] = 0.5*(jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+jacobgeo_inv[1]*Jfin[3]); 
  fin[7] = 0.5*(jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[0]*Jfin[7]+jacobgeo_inv[3]*Jfin[6]+jacobgeo_inv[2]*Jfin[3]); 
  fin[8] = 0.5*(jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[3]*Jfin[9]+jacobgeo_inv[0]*Jfin[8]+jacobgeo_inv[1]*Jfin[4]); 
  fin[9] = 0.5*(jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[3]*Jfin[8]+jacobgeo_inv[2]*Jfin[4]); 
  fin[10] = 0.5*(jacobgeo_inv[3]*Jfin[15]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[0]*Jfin[10]); 
  fin[11] = 0.5*(jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[1]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]); 
  fin[12] = 0.5*(jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[2]*Jfin[8]+jacobgeo_inv[3]*Jfin[4]); 
  fin[13] = 0.5*(jacobgeo_inv[2]*Jfin[15]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[1]*Jfin[10]); 
  fin[14] = 0.5*(jacobgeo_inv[1]*Jfin[15]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[2]*Jfin[10]); 
  fin[15] = 0.5*(jacobgeo_inv[0]*Jfin[15]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[3]*Jfin[10]); 
  fin[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*Jfin[20]+15.000000000000002*(jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[1]*Jfin[17])+15.0*jacobgeo_inv[0]*Jfin[16]); 
  fin[17] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*Jfin[20]+15.0*(jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[0]*Jfin[17])+15.000000000000002*jacobgeo_inv[1]*Jfin[16]); 
  fin[18] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*Jfin[20]+15.0*(jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[3]*Jfin[17])+15.000000000000002*jacobgeo_inv[2]*Jfin[16]); 
  fin[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*Jfin[23]+15.000000000000002*(jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21])+15.0*jacobgeo_inv[0]*Jfin[19]); 
  fin[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*Jfin[20]+15.000000000000002*(jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[2]*Jfin[17])+15.0*jacobgeo_inv[3]*Jfin[16]); 
  fin[21] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*Jfin[23]+15.0*(jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[0]*Jfin[21])+15.000000000000002*jacobgeo_inv[1]*Jfin[19]); 
  fin[22] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*Jfin[23]+15.0*(jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[3]*Jfin[21])+15.000000000000002*jacobgeo_inv[2]*Jfin[19]); 
  fin[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*Jfin[23]+15.000000000000002*(jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21])+15.0*jacobgeo_inv[3]*Jfin[19]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffx(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsxz(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[24];
  fin[0] = 0.5*(jacobgeo_inv[3]*Jfin[5]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.5*(jacobgeo_inv[2]*Jfin[5]+Jfin[2]*jacobgeo_inv[3]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.5*(jacobgeo_inv[1]*Jfin[5]+Jfin[1]*jacobgeo_inv[3]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.5*(jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[7]+jacobgeo_inv[1]*Jfin[6]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.5*(jacobgeo_inv[3]*Jfin[12]+jacobgeo_inv[2]*Jfin[9]+jacobgeo_inv[1]*Jfin[8]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.5*(jacobgeo_inv[0]*Jfin[5]+Jfin[0]*jacobgeo_inv[3]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[6] = 0.5*(jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+jacobgeo_inv[1]*Jfin[3]); 
  fin[7] = 0.5*(jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[0]*Jfin[7]+jacobgeo_inv[3]*Jfin[6]+jacobgeo_inv[2]*Jfin[3]); 
  fin[8] = 0.5*(jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[3]*Jfin[9]+jacobgeo_inv[0]*Jfin[8]+jacobgeo_inv[1]*Jfin[4]); 
  fin[9] = 0.5*(jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[3]*Jfin[8]+jacobgeo_inv[2]*Jfin[4]); 
  fin[10] = 0.5*(jacobgeo_inv[3]*Jfin[15]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[0]*Jfin[10]); 
  fin[11] = 0.5*(jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[1]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]); 
  fin[12] = 0.5*(jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[2]*Jfin[8]+jacobgeo_inv[3]*Jfin[4]); 
  fin[13] = 0.5*(jacobgeo_inv[2]*Jfin[15]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[1]*Jfin[10]); 
  fin[14] = 0.5*(jacobgeo_inv[1]*Jfin[15]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[2]*Jfin[10]); 
  fin[15] = 0.5*(jacobgeo_inv[0]*Jfin[15]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[3]*Jfin[10]); 
  fin[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*Jfin[20]+15.000000000000002*(jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[1]*Jfin[17])+15.0*jacobgeo_inv[0]*Jfin[16]); 
  fin[17] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*Jfin[20]+15.0*(jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[0]*Jfin[17])+15.000000000000002*jacobgeo_inv[1]*Jfin[16]); 
  fin[18] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*Jfin[20]+15.0*(jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[3]*Jfin[17])+15.000000000000002*jacobgeo_inv[2]*Jfin[16]); 
  fin[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*Jfin[23]+15.000000000000002*(jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21])+15.0*jacobgeo_inv[0]*Jfin[19]); 
  fin[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*Jfin[20]+15.000000000000002*(jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[2]*Jfin[17])+15.0*jacobgeo_inv[3]*Jfin[16]); 
  fin[21] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*Jfin[23]+15.0*(jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[0]*Jfin[21])+15.000000000000002*jacobgeo_inv[1]*Jfin[19]); 
  fin[22] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*Jfin[23]+15.0*(jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[3]*Jfin[21])+15.000000000000002*jacobgeo_inv[2]*Jfin[19]); 
  fin[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*Jfin[23]+15.000000000000002*(jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21])+15.0*jacobgeo_inv[3]*Jfin[19]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffx(w, dx, nu, fin, out);
  cflFreq += gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffz(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsz(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[24];
  fin[0] = 0.5*(jacobgeo_inv[3]*Jfin[5]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.5*(jacobgeo_inv[2]*Jfin[5]+Jfin[2]*jacobgeo_inv[3]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.5*(jacobgeo_inv[1]*Jfin[5]+Jfin[1]*jacobgeo_inv[3]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.5*(jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[7]+jacobgeo_inv[1]*Jfin[6]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.5*(jacobgeo_inv[3]*Jfin[12]+jacobgeo_inv[2]*Jfin[9]+jacobgeo_inv[1]*Jfin[8]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.5*(jacobgeo_inv[0]*Jfin[5]+Jfin[0]*jacobgeo_inv[3]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[6] = 0.5*(jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+jacobgeo_inv[1]*Jfin[3]); 
  fin[7] = 0.5*(jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[0]*Jfin[7]+jacobgeo_inv[3]*Jfin[6]+jacobgeo_inv[2]*Jfin[3]); 
  fin[8] = 0.5*(jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[3]*Jfin[9]+jacobgeo_inv[0]*Jfin[8]+jacobgeo_inv[1]*Jfin[4]); 
  fin[9] = 0.5*(jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[3]*Jfin[8]+jacobgeo_inv[2]*Jfin[4]); 
  fin[10] = 0.5*(jacobgeo_inv[3]*Jfin[15]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[0]*Jfin[10]); 
  fin[11] = 0.5*(jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[1]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]); 
  fin[12] = 0.5*(jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[2]*Jfin[8]+jacobgeo_inv[3]*Jfin[4]); 
  fin[13] = 0.5*(jacobgeo_inv[2]*Jfin[15]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[1]*Jfin[10]); 
  fin[14] = 0.5*(jacobgeo_inv[1]*Jfin[15]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[2]*Jfin[10]); 
  fin[15] = 0.5*(jacobgeo_inv[0]*Jfin[15]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[3]*Jfin[10]); 
  fin[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*Jfin[20]+15.000000000000002*(jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[1]*Jfin[17])+15.0*jacobgeo_inv[0]*Jfin[16]); 
  fin[17] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*Jfin[20]+15.0*(jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[0]*Jfin[17])+15.000000000000002*jacobgeo_inv[1]*Jfin[16]); 
  fin[18] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*Jfin[20]+15.0*(jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[3]*Jfin[17])+15.000000000000002*jacobgeo_inv[2]*Jfin[16]); 
  fin[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*Jfin[23]+15.000000000000002*(jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21])+15.0*jacobgeo_inv[0]*Jfin[19]); 
  fin[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*Jfin[20]+15.000000000000002*(jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[2]*Jfin[17])+15.0*jacobgeo_inv[3]*Jfin[16]); 
  fin[21] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*Jfin[23]+15.0*(jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[0]*Jfin[21])+15.000000000000002*jacobgeo_inv[1]*Jfin[19]); 
  fin[22] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*Jfin[23]+15.0*(jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[3]*Jfin[21])+15.000000000000002*jacobgeo_inv[2]*Jfin[19]); 
  fin[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*Jfin[23]+15.000000000000002*(jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21])+15.0*jacobgeo_inv[3]*Jfin[19]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffz(w, dx, nu, fin, out);

  return cflFreq;
}

