#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  out[0] += -0.0625*(8.660254037844386*coeff[1]*qr[2]-8.660254037844386*coeff[1]*ql[2]+((-9.0*qr[0])-9.0*ql[0]+18.0*qc[0])*coeff[1])*Jfac; 
  out[1] += -0.0625*(8.660254037844386*coeff[1]*qr[3]-8.660254037844386*coeff[1]*ql[3]-9.0*coeff[1]*qr[1]-9.0*coeff[1]*ql[1]+18.0*coeff[1]*qc[1])*Jfac; 
  out[2] += -0.0625*(7.0*coeff[1]*qr[2]+7.0*coeff[1]*ql[2]+46.0*coeff[1]*qc[2]+(8.660254037844386*ql[0]-8.660254037844386*qr[0])*coeff[1])*Jfac; 
  out[3] += -0.0625*(7.0*coeff[1]*qr[3]+7.0*coeff[1]*ql[3]+46.0*coeff[1]*qc[3]-8.660254037844386*coeff[1]*qr[1]+8.660254037844386*coeff[1]*ql[1])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  out[0] += -0.03125*((15.0*qr[3]+15.0*ql[3]+30.0*qc[3]-15.58845726811989*qr[1]+15.58845726811989*ql[1])*coeff[7]+(15.0*qr[2]+15.0*ql[2]+30.0*qc[2]-15.58845726811989*qr[0]+15.58845726811989*ql[0])*coeff[6]+(8.660254037844386*qr[3]-8.660254037844386*ql[3]-9.0*qr[1]-9.0*ql[1]+18.0*qc[1])*coeff[5]+(8.660254037844386*qr[2]-8.660254037844386*ql[2]-9.0*qr[0]-9.0*ql[0]+18.0*qc[0])*coeff[4])*Jfac; 
  out[1] += -0.03125*((15.0*qr[2]+15.0*ql[2]+30.0*qc[2]-15.58845726811989*qr[0]+15.58845726811989*ql[0])*coeff[7]+(15.0*qr[3]+15.0*ql[3]+30.0*qc[3]-15.58845726811989*qr[1]+15.58845726811989*ql[1])*coeff[6]+(8.660254037844386*qr[2]-8.660254037844386*ql[2]-9.0*qr[0]-9.0*ql[0]+18.0*qc[0])*coeff[5]+(8.660254037844386*qr[3]-8.660254037844386*ql[3]-9.0*qr[1]-9.0*ql[1]+18.0*qc[1])*coeff[4])*Jfac; 
  out[2] += -0.03125*((12.12435565298214*qr[3]-12.12435565298214*ql[3]-15.0*qr[1]-15.0*ql[1]+30.0*qc[1])*coeff[7]+(12.12435565298214*qr[2]-12.12435565298214*ql[2]-15.0*qr[0]-15.0*ql[0]+30.0*qc[0])*coeff[6]+(7.0*qr[3]+7.0*ql[3]+46.0*qc[3]-8.660254037844386*qr[1]+8.660254037844386*ql[1])*coeff[5]+(7.0*qr[2]+7.0*ql[2]+46.0*qc[2]-8.660254037844386*qr[0]+8.660254037844386*ql[0])*coeff[4])*Jfac; 
  out[3] += -0.03125*((12.12435565298214*qr[2]-12.12435565298214*ql[2]-15.0*qr[0]-15.0*ql[0]+30.0*qc[0])*coeff[7]+(12.12435565298214*qr[3]-12.12435565298214*ql[3]-15.0*qr[1]-15.0*ql[1]+30.0*qc[1])*coeff[6]+(7.0*qr[2]+7.0*ql[2]+46.0*qc[2]-8.660254037844386*qr[0]+8.660254037844386*ql[0])*coeff[5]+(7.0*qr[3]+7.0*ql[3]+46.0*qc[3]-8.660254037844386*qr[1]+8.660254037844386*ql[1])*coeff[4])*Jfac; 

  return 0.;

}

