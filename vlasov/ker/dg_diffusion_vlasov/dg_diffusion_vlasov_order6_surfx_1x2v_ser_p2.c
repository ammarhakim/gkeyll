#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order6_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += 0.0625*(563.489130329947*coeff[0]*qr[7]+563.489130329947*coeff[0]*ql[7]-1126.978260659894*coeff[0]*qc[7]-545.5960043841961*coeff[0]*qr[1]+545.5960043841961*coeff[0]*ql[1]+315.0*coeff[0]*qr[0]+315.0*coeff[0]*ql[0]-630.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*qr[7]-6587.944671898812*coeff[0]*ql[7]-7245.0*coeff[0]*qr[1]-7245.0*coeff[0]*ql[1]-15750.0*coeff[0]*qc[1]+4364.768035073569*coeff[0]*qr[0]-4364.768035073569*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0625*(563.4891303299469*coeff[0]*qr[11]+563.4891303299469*coeff[0]*ql[11]-1126.978260659894*coeff[0]*qc[11]-545.5960043841961*coeff[0]*qr[4]+545.5960043841961*coeff[0]*ql[4]+315.0*coeff[0]*qr[2]+315.0*coeff[0]*ql[2]-630.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.0625*(563.4891303299469*coeff[0]*qr[13]+563.4891303299469*coeff[0]*ql[13]-1126.978260659894*coeff[0]*qc[13]-545.5960043841961*coeff[0]*qr[5]+545.5960043841961*coeff[0]*ql[5]+315.0*coeff[0]*qr[3]+315.0*coeff[0]*ql[3]-630.0*coeff[0]*qc[3])*Jfac; 
  out[4] += 0.0078125*(6587.944671898817*coeff[0]*qr[11]-6587.944671898817*coeff[0]*ql[11]-7245.0*coeff[0]*qr[4]-7245.0*coeff[0]*ql[4]-15750.0*coeff[0]*qc[4]+4364.768035073569*coeff[0]*qr[2]-4364.768035073569*coeff[0]*ql[2])*Jfac; 
  out[5] += 0.0078125*(6587.944671898817*coeff[0]*qr[13]-6587.944671898817*coeff[0]*ql[13]-7245.0*coeff[0]*qr[5]-7245.0*coeff[0]*ql[5]-15750.0*coeff[0]*qc[5]+4364.768035073569*coeff[0]*qr[3]-4364.768035073569*coeff[0]*ql[3])*Jfac; 
  out[6] += 0.0625*(563.489130329947*coeff[0]*qr[17]+563.489130329947*coeff[0]*ql[17]-1126.978260659894*coeff[0]*qc[17]-545.5960043841961*coeff[0]*qr[10]+545.5960043841961*coeff[0]*ql[10]+315.0*coeff[0]*qr[6]+315.0*coeff[0]*ql[6]-630.0*coeff[0]*qc[6])*Jfac; 
  out[7] += -0.0078125*(405.0*coeff[0]*qr[7]+405.0*coeff[0]*ql[7]+18090.0*coeff[0]*qc[7]+1568.558255214003*coeff[0]*qr[1]-1568.558255214003*coeff[0]*ql[1]-1609.968943799849*coeff[0]*qr[0]-1609.968943799849*coeff[0]*ql[0]+3219.937887599698*coeff[0]*qc[0])*Jfac; 
  out[8] += -0.0625*(545.5960043841964*coeff[0]*qr[12]-545.5960043841964*coeff[0]*ql[12]-315.0*coeff[0]*qr[8]-315.0*coeff[0]*ql[8]+630.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0625*(545.5960043841964*coeff[0]*qr[15]-545.5960043841964*coeff[0]*ql[15]-315.0*coeff[0]*qr[9]-315.0*coeff[0]*ql[9]+630.0*coeff[0]*qc[9])*Jfac; 
  out[10] += 0.0078125*(6587.944671898812*coeff[0]*qr[17]-6587.944671898812*coeff[0]*ql[17]-7245.0*coeff[0]*qr[10]-7245.0*coeff[0]*ql[10]-15750.0*coeff[0]*qc[10]+4364.768035073569*coeff[0]*qr[6]-4364.768035073569*coeff[0]*ql[6])*Jfac; 
  out[11] += -0.0078125*(405.0*coeff[0]*qr[11]+405.0*coeff[0]*ql[11]+18090.0*coeff[0]*qc[11]+1568.558255214004*coeff[0]*qr[4]-1568.558255214004*coeff[0]*ql[4]-1609.968943799848*coeff[0]*qr[2]-1609.968943799848*coeff[0]*ql[2]+3219.937887599697*coeff[0]*qc[2])*Jfac; 
  out[12] += -0.0078125*(7245.0*coeff[0]*qr[12]+7245.0*coeff[0]*ql[12]+15750.0*coeff[0]*qc[12]-4364.768035073571*coeff[0]*qr[8]+4364.768035073571*coeff[0]*ql[8])*Jfac; 
  out[13] += -0.0078125*(405.0*coeff[0]*qr[13]+405.0*coeff[0]*ql[13]+18090.0*coeff[0]*qc[13]+1568.558255214004*coeff[0]*qr[5]-1568.558255214004*coeff[0]*ql[5]-1609.968943799848*coeff[0]*qr[3]-1609.968943799848*coeff[0]*ql[3]+3219.937887599697*coeff[0]*qc[3])*Jfac; 
  out[14] += -0.0625*(545.5960043841964*coeff[0]*qr[18]-545.5960043841964*coeff[0]*ql[18]-315.0*coeff[0]*qr[14]-315.0*coeff[0]*ql[14]+630.0*coeff[0]*qc[14])*Jfac; 
  out[15] += -0.0078125*(7245.0*coeff[0]*qr[15]+7245.0*coeff[0]*ql[15]+15750.0*coeff[0]*qc[15]-4364.768035073571*coeff[0]*qr[9]+4364.768035073571*coeff[0]*ql[9])*Jfac; 
  out[16] += -0.0625*(545.5960043841964*coeff[0]*qr[19]-545.5960043841964*coeff[0]*ql[19]-315.0*coeff[0]*qr[16]-315.0*coeff[0]*ql[16]+630.0*coeff[0]*qc[16])*Jfac; 
  out[17] += -0.0078125*(405.0*coeff[0]*qr[17]+405.0*coeff[0]*ql[17]+18090.0*coeff[0]*qc[17]+1568.558255214003*coeff[0]*qr[10]-1568.558255214003*coeff[0]*ql[10]-1609.968943799849*coeff[0]*qr[6]-1609.968943799849*coeff[0]*ql[6]+3219.937887599698*coeff[0]*qc[6])*Jfac; 
  out[18] += -0.0078125*(7245.0*coeff[0]*qr[18]+7245.0*coeff[0]*ql[18]+15750.0*coeff[0]*qc[18]-4364.768035073571*coeff[0]*qr[14]+4364.768035073571*coeff[0]*ql[14])*Jfac; 
  out[19] += -0.0078125*(7245.0*coeff[0]*qr[19]+7245.0*coeff[0]*ql[19]+15750.0*coeff[0]*qc[19]-4364.768035073571*coeff[0]*qr[16]+4364.768035073571*coeff[0]*ql[16])*Jfac; 

  return 0.;

}

