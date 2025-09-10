#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[0]*qr[7]+107.3312629199899*coeff[0]*ql[7]-214.6625258399798*coeff[0]*qc[7]-129.9038105676658*coeff[0]*qr[1]+129.9038105676658*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.03125*(290.4737509655563*coeff[0]*qr[7]-290.4737509655563*coeff[0]*ql[7]-405.0*coeff[0]*qr[1]-405.0*coeff[0]*ql[1]-990.0*coeff[0]*qc[1]+259.8076211353315*coeff[0]*qr[0]-259.8076211353315*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0625*(107.3312629199899*coeff[0]*qr[11]+107.3312629199899*coeff[0]*ql[11]-214.6625258399798*coeff[0]*qc[11]-129.9038105676658*coeff[0]*qr[4]+129.9038105676658*coeff[0]*ql[4]+75.0*coeff[0]*qr[2]+75.0*coeff[0]*ql[2]-150.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.0625*(107.3312629199899*coeff[0]*qr[13]+107.3312629199899*coeff[0]*ql[13]-214.6625258399798*coeff[0]*qc[13]-129.9038105676658*coeff[0]*qr[5]+129.9038105676658*coeff[0]*ql[5]+75.0*coeff[0]*qr[3]+75.0*coeff[0]*ql[3]-150.0*coeff[0]*qc[3])*Jfac; 
  out[4] += 0.03125*(290.4737509655563*coeff[0]*qr[11]-290.4737509655563*coeff[0]*ql[11]-405.0*coeff[0]*qr[4]-405.0*coeff[0]*ql[4]-990.0*coeff[0]*qc[4]+259.8076211353315*coeff[0]*qr[2]-259.8076211353315*coeff[0]*ql[2])*Jfac; 
  out[5] += 0.03125*(290.4737509655563*coeff[0]*qr[13]-290.4737509655563*coeff[0]*ql[13]-405.0*coeff[0]*qr[5]-405.0*coeff[0]*ql[5]-990.0*coeff[0]*qc[5]+259.8076211353315*coeff[0]*qr[3]-259.8076211353315*coeff[0]*ql[3])*Jfac; 
  out[6] += 0.0625*(107.3312629199899*coeff[0]*qr[17]+107.3312629199899*coeff[0]*ql[17]-214.6625258399798*coeff[0]*qc[17]-129.9038105676658*coeff[0]*qr[10]+129.9038105676658*coeff[0]*ql[10]+75.0*coeff[0]*qr[6]+75.0*coeff[0]*ql[6]-150.0*coeff[0]*qc[6])*Jfac; 
  out[7] += 0.03125*(21.0*coeff[0]*qr[7]+21.0*coeff[0]*ql[7]-1302.0*coeff[0]*qc[7]-151.0463505020892*coeff[0]*qr[1]+151.0463505020892*coeff[0]*ql[1]+134.1640786499874*coeff[0]*qr[0]+134.1640786499874*coeff[0]*ql[0]-268.3281572999748*coeff[0]*qc[0])*Jfac; 
  out[8] += -0.0625*(129.9038105676658*coeff[0]*qr[12]-129.9038105676658*coeff[0]*ql[12]-75.0*coeff[0]*qr[8]-75.0*coeff[0]*ql[8]+150.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0625*(129.9038105676658*coeff[0]*qr[15]-129.9038105676658*coeff[0]*ql[15]-75.0*coeff[0]*qr[9]-75.0*coeff[0]*ql[9]+150.0*coeff[0]*qc[9])*Jfac; 
  out[10] += 0.03125*(290.4737509655563*coeff[0]*qr[17]-290.4737509655563*coeff[0]*ql[17]-405.0*coeff[0]*qr[10]-405.0*coeff[0]*ql[10]-990.0*coeff[0]*qc[10]+259.8076211353315*coeff[0]*qr[6]-259.8076211353315*coeff[0]*ql[6])*Jfac; 
  out[11] += 0.03125*(21.0*coeff[0]*qr[11]+21.0*coeff[0]*ql[11]-1302.0*coeff[0]*qc[11]-151.0463505020893*coeff[0]*qr[4]+151.0463505020893*coeff[0]*ql[4]+134.1640786499874*coeff[0]*qr[2]+134.1640786499874*coeff[0]*ql[2]-268.3281572999747*coeff[0]*qc[2])*Jfac; 
  out[12] += -0.03125*(405.0*coeff[0]*qr[12]+405.0*coeff[0]*ql[12]+990.0*coeff[0]*qc[12]-259.8076211353317*coeff[0]*qr[8]+259.8076211353317*coeff[0]*ql[8])*Jfac; 
  out[13] += 0.03125*(21.0*coeff[0]*qr[13]+21.0*coeff[0]*ql[13]-1302.0*coeff[0]*qc[13]-151.0463505020893*coeff[0]*qr[5]+151.0463505020893*coeff[0]*ql[5]+134.1640786499874*coeff[0]*qr[3]+134.1640786499874*coeff[0]*ql[3]-268.3281572999747*coeff[0]*qc[3])*Jfac; 
  out[14] += -0.0625*(129.9038105676658*coeff[0]*qr[18]-129.9038105676658*coeff[0]*ql[18]-75.0*coeff[0]*qr[14]-75.0*coeff[0]*ql[14]+150.0*coeff[0]*qc[14])*Jfac; 
  out[15] += -0.03125*(405.0*coeff[0]*qr[15]+405.0*coeff[0]*ql[15]+990.0*coeff[0]*qc[15]-259.8076211353317*coeff[0]*qr[9]+259.8076211353317*coeff[0]*ql[9])*Jfac; 
  out[16] += -0.0625*(129.9038105676658*coeff[0]*qr[19]-129.9038105676658*coeff[0]*ql[19]-75.0*coeff[0]*qr[16]-75.0*coeff[0]*ql[16]+150.0*coeff[0]*qc[16])*Jfac; 
  out[17] += 0.03125*(21.0*coeff[0]*qr[17]+21.0*coeff[0]*ql[17]-1302.0*coeff[0]*qc[17]-151.0463505020892*coeff[0]*qr[10]+151.0463505020892*coeff[0]*ql[10]+134.1640786499874*coeff[0]*qr[6]+134.1640786499874*coeff[0]*ql[6]-268.3281572999748*coeff[0]*qc[6])*Jfac; 
  out[18] += -0.03125*(405.0*coeff[0]*qr[18]+405.0*coeff[0]*ql[18]+990.0*coeff[0]*qc[18]-259.8076211353317*coeff[0]*qr[14]+259.8076211353317*coeff[0]*ql[14])*Jfac; 
  out[19] += -0.03125*(405.0*coeff[0]*qr[19]+405.0*coeff[0]*ql[19]+990.0*coeff[0]*qc[19]-259.8076211353317*coeff[0]*qr[16]+259.8076211353317*coeff[0]*ql[16])*Jfac; 

  return 0.;

}

