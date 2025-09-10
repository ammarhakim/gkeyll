#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -0.0625*(25.98076211353316*coeff[0]*qr[1]-25.98076211353316*coeff[0]*ql[1]-15.0*coeff[0]*qr[0]-15.0*coeff[0]*ql[0]+30.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.0625*(33.0*coeff[0]*qr[1]+33.0*coeff[0]*ql[1]+114.0*coeff[0]*qc[1]-25.98076211353316*coeff[0]*qr[0]+25.98076211353316*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.0625*(25.98076211353316*coeff[0]*qr[4]-25.98076211353316*coeff[0]*ql[4]-15.0*coeff[0]*qr[2]-15.0*coeff[0]*ql[2]+30.0*coeff[0]*qc[2])*Jfac; 
  out[3] += -0.0625*(25.98076211353316*coeff[0]*qr[5]-25.98076211353316*coeff[0]*ql[5]-15.0*coeff[0]*qr[3]-15.0*coeff[0]*ql[3]+30.0*coeff[0]*qc[3])*Jfac; 
  out[4] += -0.0625*(33.0*coeff[0]*qr[4]+33.0*coeff[0]*ql[4]+114.0*coeff[0]*qc[4]-25.98076211353316*coeff[0]*qr[2]+25.98076211353316*coeff[0]*ql[2])*Jfac; 
  out[5] += -0.0625*(33.0*coeff[0]*qr[5]+33.0*coeff[0]*ql[5]+114.0*coeff[0]*qc[5]-25.98076211353316*coeff[0]*qr[3]+25.98076211353316*coeff[0]*ql[3])*Jfac; 
  out[6] += -0.0625*(25.98076211353316*coeff[0]*qr[7]-25.98076211353316*coeff[0]*ql[7]-15.0*coeff[0]*qr[6]-15.0*coeff[0]*ql[6]+30.0*coeff[0]*qc[6])*Jfac; 
  out[7] += -0.0625*(33.0*coeff[0]*qr[7]+33.0*coeff[0]*ql[7]+114.0*coeff[0]*qc[7]-25.98076211353316*coeff[0]*qr[6]+25.98076211353316*coeff[0]*ql[6])*Jfac; 
  out[8] += -0.0625*(25.98076211353316*coeff[0]*qr[9]-25.98076211353316*coeff[0]*ql[9]-15.0*coeff[0]*qr[8]-15.0*coeff[0]*ql[8]+30.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0625*(33.0*coeff[0]*qr[9]+33.0*coeff[0]*ql[9]+114.0*coeff[0]*qc[9]-25.98076211353316*coeff[0]*qr[8]+25.98076211353316*coeff[0]*ql[8])*Jfac; 
  out[10] += -0.0625*(25.98076211353316*coeff[0]*qr[11]-25.98076211353316*coeff[0]*ql[11]-15.0*coeff[0]*qr[10]-15.0*coeff[0]*ql[10]+30.0*coeff[0]*qc[10])*Jfac; 
  out[11] += -0.0625*(33.0*coeff[0]*qr[11]+33.0*coeff[0]*ql[11]+114.0*coeff[0]*qc[11]-25.98076211353316*coeff[0]*qr[10]+25.98076211353316*coeff[0]*ql[10])*Jfac; 
  out[12] += -0.0625*(25.98076211353316*coeff[0]*qr[13]-25.98076211353316*coeff[0]*ql[13]-15.0*coeff[0]*qr[12]-15.0*coeff[0]*ql[12]+30.0*coeff[0]*qc[12])*Jfac; 
  out[13] += -0.0625*(33.0*coeff[0]*qr[13]+33.0*coeff[0]*ql[13]+114.0*coeff[0]*qc[13]-25.98076211353316*coeff[0]*qr[12]+25.98076211353316*coeff[0]*ql[12])*Jfac; 
  out[14] += -0.0625*(25.98076211353316*coeff[0]*qr[15]-25.98076211353316*coeff[0]*ql[15]-15.0*coeff[0]*qr[14]-15.0*coeff[0]*ql[14]+30.0*coeff[0]*qc[14])*Jfac; 
  out[15] += -0.0625*(33.0*coeff[0]*qr[15]+33.0*coeff[0]*ql[15]+114.0*coeff[0]*qc[15]-25.98076211353316*coeff[0]*qr[14]+25.98076211353316*coeff[0]*ql[14])*Jfac; 

  return 0.;

}

