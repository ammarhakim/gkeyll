#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],6.);

  out[0] += 0.0625*(563.489130329947*coeff[1]*qr[5]+563.489130329947*coeff[1]*ql[5]-1126.978260659894*coeff[1]*qc[5]-545.5960043841961*coeff[1]*qr[2]+545.5960043841961*coeff[1]*ql[2]+(315.0*qr[0]+315.0*ql[0]-630.0*qc[0])*coeff[1])*Jfac; 
  out[1] += 0.0625*(563.4891303299469*coeff[1]*qr[7]+563.4891303299469*coeff[1]*ql[7]-1126.978260659894*coeff[1]*qc[7]-545.5960043841961*coeff[1]*qr[3]+545.5960043841961*coeff[1]*ql[3]+315.0*coeff[1]*qr[1]+315.0*coeff[1]*ql[1]-630.0*coeff[1]*qc[1])*Jfac; 
  out[2] += 0.0078125*(6587.944671898812*coeff[1]*qr[5]-6587.944671898812*coeff[1]*ql[5]-7245.0*coeff[1]*qr[2]-7245.0*coeff[1]*ql[2]-15750.0*coeff[1]*qc[2]+(4364.768035073569*qr[0]-4364.768035073569*ql[0])*coeff[1])*Jfac; 
  out[3] += 0.0078125*(6587.944671898817*coeff[1]*qr[7]-6587.944671898817*coeff[1]*ql[7]-7245.0*coeff[1]*qr[3]-7245.0*coeff[1]*ql[3]-15750.0*coeff[1]*qc[3]+4364.768035073569*coeff[1]*qr[1]-4364.768035073569*coeff[1]*ql[1])*Jfac; 
  out[4] += -0.0625*(545.5960043841964*coeff[1]*qr[6]-545.5960043841964*coeff[1]*ql[6]-315.0*coeff[1]*qr[4]-315.0*coeff[1]*ql[4]+630.0*coeff[1]*qc[4])*Jfac; 
  out[5] += -0.0078125*(405.0*coeff[1]*qr[5]+405.0*coeff[1]*ql[5]+18090.0*coeff[1]*qc[5]+1568.558255214003*coeff[1]*qr[2]-1568.558255214003*coeff[1]*ql[2]+((-1609.968943799849*qr[0])-1609.968943799849*ql[0]+3219.937887599698*qc[0])*coeff[1])*Jfac; 
  out[6] += -0.0078125*(7245.0*coeff[1]*qr[6]+7245.0*coeff[1]*ql[6]+15750.0*coeff[1]*qc[6]-4364.768035073571*coeff[1]*qr[4]+4364.768035073571*coeff[1]*ql[4])*Jfac; 
  out[7] += -0.0078125*(405.0*coeff[1]*qr[7]+405.0*coeff[1]*ql[7]+18090.0*coeff[1]*qc[7]+1568.558255214004*coeff[1]*qr[3]-1568.558255214004*coeff[1]*ql[3]-1609.968943799848*coeff[1]*qr[1]-1609.968943799848*coeff[1]*ql[1]+3219.937887599697*coeff[1]*qc[1])*Jfac; 

  return 0.;

}

