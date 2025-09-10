#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],6.);

  out[0] += 0.0625*(563.489130329947*coeff[0]*qr[4]+563.489130329947*coeff[0]*ql[4]-1126.978260659894*coeff[0]*qc[4]-545.5960043841961*coeff[0]*qr[1]+545.5960043841961*coeff[0]*ql[1]+315.0*coeff[0]*qr[0]+315.0*coeff[0]*ql[0]-630.0*coeff[0]*qc[0])*rdx2Sq; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*qr[4]-6587.944671898812*coeff[0]*ql[4]-7245.0*coeff[0]*qr[1]-7245.0*coeff[0]*ql[1]-15750.0*coeff[0]*qc[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*rdx2Sq; 
  out[2] += 0.0625*(563.4891303299469*coeff[0]*qr[6]+563.4891303299469*coeff[0]*ql[6]-1126.978260659894*coeff[0]*qc[6]-545.5960043841961*coeff[0]*qr[3]+545.5960043841961*coeff[0]*ql[3]+315.0*coeff[0]*qr[2]+315.0*coeff[0]*ql[2]-630.0*coeff[0]*qc[2])*rdx2Sq; 
  out[3] += 0.0078125*(6587.944671898817*coeff[0]*qr[6]-6587.944671898817*coeff[0]*ql[6]-7245.0*coeff[0]*qr[3]-7245.0*coeff[0]*ql[3]-15750.0*coeff[0]*qc[3]+4364.768035073569*coeff[0]*qr[2]-4364.768035073569*coeff[0]*ql[2])*rdx2Sq; 
  out[4] += -0.0078125*(405.0*coeff[0]*qr[4]+405.0*coeff[0]*ql[4]+18090.0*coeff[0]*qc[4]+1568.558255214003*coeff[0]*qr[1]-1568.558255214003*coeff[0]*ql[1]-1609.968943799849*coeff[0]*qr[0]-1609.968943799849*coeff[0]*ql[0]+3219.937887599698*coeff[0]*qc[0])*rdx2Sq; 
  out[5] += -0.0625*(545.5960043841964*coeff[0]*qr[7]-545.5960043841964*coeff[0]*ql[7]-315.0*coeff[0]*qr[5]-315.0*coeff[0]*ql[5]+630.0*coeff[0]*qc[5])*rdx2Sq; 
  out[6] += -0.0078125*(405.0*coeff[0]*qr[6]+405.0*coeff[0]*ql[6]+18090.0*coeff[0]*qc[6]+1568.558255214004*coeff[0]*qr[3]-1568.558255214004*coeff[0]*ql[3]-1609.968943799848*coeff[0]*qr[2]-1609.968943799848*coeff[0]*ql[2]+3219.937887599697*coeff[0]*qc[2])*rdx2Sq; 
  out[7] += -0.0078125*(7245.0*coeff[0]*qr[7]+7245.0*coeff[0]*ql[7]+15750.0*coeff[0]*qc[7]-4364.768035073571*coeff[0]*qr[5]+4364.768035073571*coeff[0]*ql[5])*rdx2Sq; 

  return 0.;

}

