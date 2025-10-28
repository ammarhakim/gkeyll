#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 9.0*coeff[0]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 9.0*coeff[1]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[2]; 
  return 9.0*coeff[2]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffy(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_constcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[0]; 
  return 3.181980515339463*coeff[0]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[1]; 
  return 3.181980515339463*coeff[20]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // coeff: Diffusion coefficient vector
  // q: Input field
  // out: Incremented output

  const double rdx2 = 2.0/dx[2]; 
  return 3.181980515339463*coeff[40]*pow(rdx2, 2.0); 
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffx(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffdirsy(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffy(w, dx, coeff, q, out);
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *q, double* GKYL_RESTRICT out) 
{ 
  double cflFreq = 0.;
  
  cflFreq += dg_diffusion_vlasov_order2_vol_3x_ser_p2_varcoeff_diffz(w, dx, coeff, q, out);

  return cflFreq;
}

