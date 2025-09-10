#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 4.0*coeff[0]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 4.0*coeff[1]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[2]; 
  return 4.0*coeff[2]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 0.5*(coeff[7]*jacobgeo_inv[7]+coeff[6]*jacobgeo_inv[6]+coeff[5]*jacobgeo_inv[5]+coeff[4]*jacobgeo_inv[4]+coeff[3]*jacobgeo_inv[3]+coeff[2]*jacobgeo_inv[2]+coeff[1]*jacobgeo_inv[1]+coeff[0]*jacobgeo_inv[0])*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 0.5*(jacobgeo_inv[7]*coeff[15]+jacobgeo_inv[6]*coeff[14]+jacobgeo_inv[5]*coeff[13]+jacobgeo_inv[4]*coeff[12]+jacobgeo_inv[3]*coeff[11]+jacobgeo_inv[2]*coeff[10]+jacobgeo_inv[1]*coeff[9]+jacobgeo_inv[0]*coeff[8])*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[2]; 
  return 0.5*(jacobgeo_inv[7]*coeff[23]+jacobgeo_inv[6]*coeff[22]+jacobgeo_inv[5]*coeff[21]+jacobgeo_inv[4]*coeff[20]+jacobgeo_inv[3]*coeff[19]+jacobgeo_inv[2]*coeff[18]+jacobgeo_inv[1]*coeff[17]+jacobgeo_inv[0]*coeff[16])*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffx(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffy(w, dx, coeff, jacobgeo_inv, q, out);
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffz(w, dx, coeff, jacobgeo_inv, q, out);

  return cflFreq;
}

