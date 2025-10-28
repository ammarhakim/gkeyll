#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0625*(563.489130329947*coeff[0]*qr[2]+563.489130329947*coeff[0]*ql[2]-1126.978260659894*coeff[0]*qc[2]-545.5960043841961*coeff[0]*qr[1]+545.5960043841961*coeff[0]*ql[1]+315.0*coeff[0]*qr[0]+315.0*coeff[0]*ql[0]-630.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*qr[2]-6587.944671898812*coeff[0]*ql[2]-7245.0*coeff[0]*qr[1]-7245.0*coeff[0]*ql[1]-15750.0*coeff[0]*qc[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.0078125*(405.0*coeff[0]*qr[2]+405.0*coeff[0]*ql[2]+18090.0*coeff[0]*qc[2]+1568.558255214003*coeff[0]*qr[1]-1568.558255214003*coeff[0]*ql[1]-1609.968943799849*coeff[0]*qr[0]-1609.968943799849*coeff[0]*ql[0]+3219.937887599698*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

