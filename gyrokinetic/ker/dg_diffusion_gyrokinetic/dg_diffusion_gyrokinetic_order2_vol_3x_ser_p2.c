#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 9.0*coeff[0]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 9.0*coeff[1]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[2]; 
  return 9.0*coeff[2]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 1.125*(coeff[19]*jacobgeo_inv[19]+coeff[18]*jacobgeo_inv[18]+coeff[17]*jacobgeo_inv[17]+coeff[16]*jacobgeo_inv[16]+coeff[15]*jacobgeo_inv[15]+coeff[14]*jacobgeo_inv[14]+coeff[13]*jacobgeo_inv[13]+coeff[12]*jacobgeo_inv[12]+coeff[11]*jacobgeo_inv[11]+coeff[10]*jacobgeo_inv[10]+coeff[9]*jacobgeo_inv[9]+coeff[8]*jacobgeo_inv[8]+coeff[7]*jacobgeo_inv[7]+coeff[6]*jacobgeo_inv[6]+coeff[5]*jacobgeo_inv[5]+coeff[4]*jacobgeo_inv[4]+coeff[3]*jacobgeo_inv[3]+coeff[2]*jacobgeo_inv[2]+coeff[1]*jacobgeo_inv[1]+coeff[0]*jacobgeo_inv[0])*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 1.125*(jacobgeo_inv[19]*coeff[39]+jacobgeo_inv[18]*coeff[38]+jacobgeo_inv[17]*coeff[37]+jacobgeo_inv[16]*coeff[36]+jacobgeo_inv[15]*coeff[35]+jacobgeo_inv[14]*coeff[34]+jacobgeo_inv[13]*coeff[33]+jacobgeo_inv[12]*coeff[32]+jacobgeo_inv[11]*coeff[31]+jacobgeo_inv[10]*coeff[30]+jacobgeo_inv[9]*coeff[29]+jacobgeo_inv[8]*coeff[28]+jacobgeo_inv[7]*coeff[27]+jacobgeo_inv[6]*coeff[26]+jacobgeo_inv[5]*coeff[25]+jacobgeo_inv[4]*coeff[24]+jacobgeo_inv[3]*coeff[23]+jacobgeo_inv[2]*coeff[22]+jacobgeo_inv[1]*coeff[21]+jacobgeo_inv[0]*coeff[20])*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[2]; 
  return 1.125*(jacobgeo_inv[19]*coeff[59]+jacobgeo_inv[18]*coeff[58]+jacobgeo_inv[17]*coeff[57]+jacobgeo_inv[16]*coeff[56]+jacobgeo_inv[15]*coeff[55]+jacobgeo_inv[14]*coeff[54]+jacobgeo_inv[13]*coeff[53]+jacobgeo_inv[12]*coeff[52]+jacobgeo_inv[11]*coeff[51]+jacobgeo_inv[10]*coeff[50]+jacobgeo_inv[9]*coeff[49]+jacobgeo_inv[8]*coeff[48]+jacobgeo_inv[7]*coeff[47]+jacobgeo_inv[6]*coeff[46]+jacobgeo_inv[5]*coeff[45]+jacobgeo_inv[4]*coeff[44]+jacobgeo_inv[3]*coeff[43]+jacobgeo_inv[2]*coeff[42]+jacobgeo_inv[1]*coeff[41]+jacobgeo_inv[0]*coeff[40])*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

