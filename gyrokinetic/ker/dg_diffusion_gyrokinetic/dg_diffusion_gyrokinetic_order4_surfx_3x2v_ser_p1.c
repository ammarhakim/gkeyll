#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_3x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  out[0] += -0.0625*(25.98076211353316*coeff[0]*qr[1]-25.98076211353316*coeff[0]*ql[1]-15.0*coeff[0]*qr[0]-15.0*coeff[0]*ql[0]+30.0*coeff[0]*qc[0])*rdx2Sq; 
  out[1] += -0.0625*(33.0*coeff[0]*qr[1]+33.0*coeff[0]*ql[1]+114.0*coeff[0]*qc[1]-25.98076211353316*coeff[0]*qr[0]+25.98076211353316*coeff[0]*ql[0])*rdx2Sq; 
  out[2] += -0.0625*(25.98076211353316*coeff[0]*qr[6]-25.98076211353316*coeff[0]*ql[6]-15.0*coeff[0]*qr[2]-15.0*coeff[0]*ql[2]+30.0*coeff[0]*qc[2])*rdx2Sq; 
  out[3] += -0.0625*(25.98076211353316*coeff[0]*qr[7]-25.98076211353316*coeff[0]*ql[7]-15.0*coeff[0]*qr[3]-15.0*coeff[0]*ql[3]+30.0*coeff[0]*qc[3])*rdx2Sq; 
  out[4] += -0.0625*(25.98076211353316*coeff[0]*qr[9]-25.98076211353316*coeff[0]*ql[9]-15.0*coeff[0]*qr[4]-15.0*coeff[0]*ql[4]+30.0*coeff[0]*qc[4])*rdx2Sq; 
  out[5] += -0.0625*(25.98076211353316*coeff[0]*qr[12]-25.98076211353316*coeff[0]*ql[12]-15.0*coeff[0]*qr[5]-15.0*coeff[0]*ql[5]+30.0*coeff[0]*qc[5])*rdx2Sq; 
  out[6] += -0.0625*(33.0*coeff[0]*qr[6]+33.0*coeff[0]*ql[6]+114.0*coeff[0]*qc[6]-25.98076211353316*coeff[0]*qr[2]+25.98076211353316*coeff[0]*ql[2])*rdx2Sq; 
  out[7] += -0.0625*(33.0*coeff[0]*qr[7]+33.0*coeff[0]*ql[7]+114.0*coeff[0]*qc[7]-25.98076211353316*coeff[0]*qr[3]+25.98076211353316*coeff[0]*ql[3])*rdx2Sq; 
  out[8] += -0.0625*(25.98076211353316*coeff[0]*qr[16]-25.98076211353316*coeff[0]*ql[16]-15.0*coeff[0]*qr[8]-15.0*coeff[0]*ql[8]+30.0*coeff[0]*qc[8])*rdx2Sq; 
  out[9] += -0.0625*(33.0*coeff[0]*qr[9]+33.0*coeff[0]*ql[9]+114.0*coeff[0]*qc[9]-25.98076211353316*coeff[0]*qr[4]+25.98076211353316*coeff[0]*ql[4])*rdx2Sq; 
  out[10] += -0.0625*(25.98076211353316*coeff[0]*qr[17]-25.98076211353316*coeff[0]*ql[17]-15.0*coeff[0]*qr[10]-15.0*coeff[0]*ql[10]+30.0*coeff[0]*qc[10])*rdx2Sq; 
  out[11] += -0.0625*(25.98076211353316*coeff[0]*qr[18]-25.98076211353316*coeff[0]*ql[18]-15.0*coeff[0]*qr[11]-15.0*coeff[0]*ql[11]+30.0*coeff[0]*qc[11])*rdx2Sq; 
  out[12] += -0.0625*(33.0*coeff[0]*qr[12]+33.0*coeff[0]*ql[12]+114.0*coeff[0]*qc[12]-25.98076211353316*coeff[0]*qr[5]+25.98076211353316*coeff[0]*ql[5])*rdx2Sq; 
  out[13] += -0.0625*(25.98076211353316*coeff[0]*qr[20]-25.98076211353316*coeff[0]*ql[20]-15.0*coeff[0]*qr[13]-15.0*coeff[0]*ql[13]+30.0*coeff[0]*qc[13])*rdx2Sq; 
  out[14] += -0.0625*(25.98076211353316*coeff[0]*qr[21]-25.98076211353316*coeff[0]*ql[21]-15.0*coeff[0]*qr[14]-15.0*coeff[0]*ql[14]+30.0*coeff[0]*qc[14])*rdx2Sq; 
  out[15] += -0.0625*(25.98076211353316*coeff[0]*qr[23]-25.98076211353316*coeff[0]*ql[23]-15.0*coeff[0]*qr[15]-15.0*coeff[0]*ql[15]+30.0*coeff[0]*qc[15])*rdx2Sq; 
  out[16] += -0.0625*(33.0*coeff[0]*qr[16]+33.0*coeff[0]*ql[16]+114.0*coeff[0]*qc[16]-25.98076211353316*coeff[0]*qr[8]+25.98076211353316*coeff[0]*ql[8])*rdx2Sq; 
  out[17] += -0.0625*(33.0*coeff[0]*qr[17]+33.0*coeff[0]*ql[17]+114.0*coeff[0]*qc[17]-25.98076211353316*coeff[0]*qr[10]+25.98076211353316*coeff[0]*ql[10])*rdx2Sq; 
  out[18] += -0.0625*(33.0*coeff[0]*qr[18]+33.0*coeff[0]*ql[18]+114.0*coeff[0]*qc[18]-25.98076211353316*coeff[0]*qr[11]+25.98076211353316*coeff[0]*ql[11])*rdx2Sq; 
  out[19] += -0.0625*(25.98076211353316*coeff[0]*qr[26]-25.98076211353316*coeff[0]*ql[26]-15.0*coeff[0]*qr[19]-15.0*coeff[0]*ql[19]+30.0*coeff[0]*qc[19])*rdx2Sq; 
  out[20] += -0.0625*(33.0*coeff[0]*qr[20]+33.0*coeff[0]*ql[20]+114.0*coeff[0]*qc[20]-25.98076211353316*coeff[0]*qr[13]+25.98076211353316*coeff[0]*ql[13])*rdx2Sq; 
  out[21] += -0.0625*(33.0*coeff[0]*qr[21]+33.0*coeff[0]*ql[21]+114.0*coeff[0]*qc[21]-25.98076211353316*coeff[0]*qr[14]+25.98076211353316*coeff[0]*ql[14])*rdx2Sq; 
  out[22] += -0.0625*(25.98076211353316*coeff[0]*qr[27]-25.98076211353316*coeff[0]*ql[27]-15.0*coeff[0]*qr[22]-15.0*coeff[0]*ql[22]+30.0*coeff[0]*qc[22])*rdx2Sq; 
  out[23] += -0.0625*(33.0*coeff[0]*qr[23]+33.0*coeff[0]*ql[23]+114.0*coeff[0]*qc[23]-25.98076211353316*coeff[0]*qr[15]+25.98076211353316*coeff[0]*ql[15])*rdx2Sq; 
  out[24] += -0.0625*(25.98076211353316*coeff[0]*qr[28]-25.98076211353316*coeff[0]*ql[28]-15.0*coeff[0]*qr[24]-15.0*coeff[0]*ql[24]+30.0*coeff[0]*qc[24])*rdx2Sq; 
  out[25] += -0.0625*(25.98076211353316*coeff[0]*qr[29]-25.98076211353316*coeff[0]*ql[29]-15.0*coeff[0]*qr[25]-15.0*coeff[0]*ql[25]+30.0*coeff[0]*qc[25])*rdx2Sq; 
  out[26] += -0.0625*(33.0*coeff[0]*qr[26]+33.0*coeff[0]*ql[26]+114.0*coeff[0]*qc[26]-25.98076211353316*coeff[0]*qr[19]+25.98076211353316*coeff[0]*ql[19])*rdx2Sq; 
  out[27] += -0.0625*(33.0*coeff[0]*qr[27]+33.0*coeff[0]*ql[27]+114.0*coeff[0]*qc[27]-25.98076211353316*coeff[0]*qr[22]+25.98076211353316*coeff[0]*ql[22])*rdx2Sq; 
  out[28] += -0.0625*(33.0*coeff[0]*qr[28]+33.0*coeff[0]*ql[28]+114.0*coeff[0]*qc[28]-25.98076211353316*coeff[0]*qr[24]+25.98076211353316*coeff[0]*ql[24])*rdx2Sq; 
  out[29] += -0.0625*(33.0*coeff[0]*qr[29]+33.0*coeff[0]*ql[29]+114.0*coeff[0]*qc[29]-25.98076211353316*coeff[0]*qr[25]+25.98076211353316*coeff[0]*ql[25])*rdx2Sq; 
  out[30] += -0.0625*(25.98076211353316*coeff[0]*qr[31]-25.98076211353316*coeff[0]*ql[31]-15.0*coeff[0]*qr[30]-15.0*coeff[0]*ql[30]+30.0*coeff[0]*qc[30])*rdx2Sq; 
  out[31] += -0.0625*(33.0*coeff[0]*qr[31]+33.0*coeff[0]*ql[31]+114.0*coeff[0]*qc[31]-25.98076211353316*coeff[0]*qr[30]+25.98076211353316*coeff[0]*ql[30])*rdx2Sq; 
  out[32] += -0.0625*(25.98076211353316*coeff[0]*qr[33]-25.98076211353316*coeff[0]*ql[33]-15.0*coeff[0]*qr[32]-15.0*coeff[0]*ql[32]+30.0*coeff[0]*qc[32])*rdx2Sq; 
  out[33] += -0.0625*(33.0*coeff[0]*qr[33]+33.0*coeff[0]*ql[33]+114.0*coeff[0]*qc[33]-25.98076211353316*coeff[0]*qr[32]+25.98076211353316*coeff[0]*ql[32])*rdx2Sq; 
  out[34] += -0.0625*(25.98076211353316*coeff[0]*qr[37]-25.98076211353316*coeff[0]*ql[37]-15.0*coeff[0]*qr[34]-15.0*coeff[0]*ql[34]+30.0*coeff[0]*qc[34])*rdx2Sq; 
  out[35] += -0.0625*(25.98076211353316*coeff[0]*qr[38]-25.98076211353316*coeff[0]*ql[38]-15.0*coeff[0]*qr[35]-15.0*coeff[0]*ql[35]+30.0*coeff[0]*qc[35])*rdx2Sq; 
  out[36] += -0.0625*(25.98076211353316*coeff[0]*qr[40]-25.98076211353316*coeff[0]*ql[40]-15.0*coeff[0]*qr[36]-15.0*coeff[0]*ql[36]+30.0*coeff[0]*qc[36])*rdx2Sq; 
  out[37] += -0.0625*(33.0*coeff[0]*qr[37]+33.0*coeff[0]*ql[37]+114.0*coeff[0]*qc[37]-25.98076211353316*coeff[0]*qr[34]+25.98076211353316*coeff[0]*ql[34])*rdx2Sq; 
  out[38] += -0.0625*(33.0*coeff[0]*qr[38]+33.0*coeff[0]*ql[38]+114.0*coeff[0]*qc[38]-25.98076211353316*coeff[0]*qr[35]+25.98076211353316*coeff[0]*ql[35])*rdx2Sq; 
  out[39] += -0.0625*(25.98076211353316*coeff[0]*qr[43]-25.98076211353316*coeff[0]*ql[43]-15.0*coeff[0]*qr[39]-15.0*coeff[0]*ql[39]+30.0*coeff[0]*qc[39])*rdx2Sq; 
  out[40] += -0.0625*(33.0*coeff[0]*qr[40]+33.0*coeff[0]*ql[40]+114.0*coeff[0]*qc[40]-25.98076211353316*coeff[0]*qr[36]+25.98076211353316*coeff[0]*ql[36])*rdx2Sq; 
  out[41] += -0.0625*(25.98076211353316*coeff[0]*qr[44]-25.98076211353316*coeff[0]*ql[44]-15.0*coeff[0]*qr[41]-15.0*coeff[0]*ql[41]+30.0*coeff[0]*qc[41])*rdx2Sq; 
  out[42] += -0.0625*(25.98076211353316*coeff[0]*qr[45]-25.98076211353316*coeff[0]*ql[45]-15.0*coeff[0]*qr[42]-15.0*coeff[0]*ql[42]+30.0*coeff[0]*qc[42])*rdx2Sq; 
  out[43] += -0.0625*(33.0*coeff[0]*qr[43]+33.0*coeff[0]*ql[43]+114.0*coeff[0]*qc[43]-25.98076211353316*coeff[0]*qr[39]+25.98076211353316*coeff[0]*ql[39])*rdx2Sq; 
  out[44] += -0.0625*(33.0*coeff[0]*qr[44]+33.0*coeff[0]*ql[44]+114.0*coeff[0]*qc[44]-25.98076211353316*coeff[0]*qr[41]+25.98076211353316*coeff[0]*ql[41])*rdx2Sq; 
  out[45] += -0.0625*(33.0*coeff[0]*qr[45]+33.0*coeff[0]*ql[45]+114.0*coeff[0]*qc[45]-25.98076211353316*coeff[0]*qr[42]+25.98076211353316*coeff[0]*ql[42])*rdx2Sq; 
  out[46] += -0.0625*(25.98076211353316*coeff[0]*qr[47]-25.98076211353316*coeff[0]*ql[47]-15.0*coeff[0]*qr[46]-15.0*coeff[0]*ql[46]+30.0*coeff[0]*qc[46])*rdx2Sq; 
  out[47] += -0.0625*(33.0*coeff[0]*qr[47]+33.0*coeff[0]*ql[47]+114.0*coeff[0]*qc[47]-25.98076211353316*coeff[0]*qr[46]+25.98076211353316*coeff[0]*ql[46])*rdx2Sq; 

  return 0.;

}

