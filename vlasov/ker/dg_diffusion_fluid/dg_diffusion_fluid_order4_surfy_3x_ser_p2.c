#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[1]*qr[8]+107.3312629199899*coeff[1]*ql[8]-214.6625258399798*coeff[1]*qc[8]-129.9038105676658*coeff[1]*qr[2]+129.9038105676658*coeff[1]*ql[2]+(75.0*qr[0]+75.0*ql[0]-150.0*qc[0])*coeff[1])*Jfac; 
  out[1] += 0.0625*(107.3312629199899*coeff[1]*qr[12]+107.3312629199899*coeff[1]*ql[12]-214.6625258399798*coeff[1]*qc[12]-129.9038105676658*coeff[1]*qr[4]+129.9038105676658*coeff[1]*ql[4]+75.0*coeff[1]*qr[1]+75.0*coeff[1]*ql[1]-150.0*coeff[1]*qc[1])*Jfac; 
  out[2] += 0.03125*(290.4737509655563*coeff[1]*qr[8]-290.4737509655563*coeff[1]*ql[8]-405.0*coeff[1]*qr[2]-405.0*coeff[1]*ql[2]-990.0*coeff[1]*qc[2]+(259.8076211353315*qr[0]-259.8076211353315*ql[0])*coeff[1])*Jfac; 
  out[3] += 0.0625*(107.3312629199899*coeff[1]*qr[14]+107.3312629199899*coeff[1]*ql[14]-214.6625258399798*coeff[1]*qc[14]-129.9038105676658*coeff[1]*qr[6]+129.9038105676658*coeff[1]*ql[6]+75.0*coeff[1]*qr[3]+75.0*coeff[1]*ql[3]-150.0*coeff[1]*qc[3])*Jfac; 
  out[4] += 0.03125*(290.4737509655563*coeff[1]*qr[12]-290.4737509655563*coeff[1]*ql[12]-405.0*coeff[1]*qr[4]-405.0*coeff[1]*ql[4]-990.0*coeff[1]*qc[4]+259.8076211353315*coeff[1]*qr[1]-259.8076211353315*coeff[1]*ql[1])*Jfac; 
  out[5] += 0.0625*(107.3312629199899*coeff[1]*qr[18]+107.3312629199899*coeff[1]*ql[18]-214.6625258399798*coeff[1]*qc[18]-129.9038105676658*coeff[1]*qr[10]+129.9038105676658*coeff[1]*ql[10]+75.0*coeff[1]*qr[5]+75.0*coeff[1]*ql[5]-150.0*coeff[1]*qc[5])*Jfac; 
  out[6] += 0.03125*(290.4737509655563*coeff[1]*qr[14]-290.4737509655563*coeff[1]*ql[14]-405.0*coeff[1]*qr[6]-405.0*coeff[1]*ql[6]-990.0*coeff[1]*qc[6]+259.8076211353315*coeff[1]*qr[3]-259.8076211353315*coeff[1]*ql[3])*Jfac; 
  out[7] += -0.0625*(129.9038105676658*coeff[1]*qr[11]-129.9038105676658*coeff[1]*ql[11]-75.0*coeff[1]*qr[7]-75.0*coeff[1]*ql[7]+150.0*coeff[1]*qc[7])*Jfac; 
  out[8] += 0.03125*(21.0*coeff[1]*qr[8]+21.0*coeff[1]*ql[8]-1302.0*coeff[1]*qc[8]-151.0463505020892*coeff[1]*qr[2]+151.0463505020892*coeff[1]*ql[2]+(134.1640786499874*qr[0]+134.1640786499874*ql[0]-268.3281572999748*qc[0])*coeff[1])*Jfac; 
  out[9] += -0.0625*(129.9038105676658*coeff[1]*qr[16]-129.9038105676658*coeff[1]*ql[16]-75.0*coeff[1]*qr[9]-75.0*coeff[1]*ql[9]+150.0*coeff[1]*qc[9])*Jfac; 
  out[10] += 0.03125*(290.4737509655563*coeff[1]*qr[18]-290.4737509655563*coeff[1]*ql[18]-405.0*coeff[1]*qr[10]-405.0*coeff[1]*ql[10]-990.0*coeff[1]*qc[10]+259.8076211353315*coeff[1]*qr[5]-259.8076211353315*coeff[1]*ql[5])*Jfac; 
  out[11] += -0.03125*(405.0*coeff[1]*qr[11]+405.0*coeff[1]*ql[11]+990.0*coeff[1]*qc[11]-259.8076211353317*coeff[1]*qr[7]+259.8076211353317*coeff[1]*ql[7])*Jfac; 
  out[12] += 0.03125*(21.0*coeff[1]*qr[12]+21.0*coeff[1]*ql[12]-1302.0*coeff[1]*qc[12]-151.0463505020893*coeff[1]*qr[4]+151.0463505020893*coeff[1]*ql[4]+134.1640786499874*coeff[1]*qr[1]+134.1640786499874*coeff[1]*ql[1]-268.3281572999747*coeff[1]*qc[1])*Jfac; 
  out[13] += -0.0625*(129.9038105676658*coeff[1]*qr[17]-129.9038105676658*coeff[1]*ql[17]-75.0*coeff[1]*qr[13]-75.0*coeff[1]*ql[13]+150.0*coeff[1]*qc[13])*Jfac; 
  out[14] += 0.03125*(21.0*coeff[1]*qr[14]+21.0*coeff[1]*ql[14]-1302.0*coeff[1]*qc[14]-151.0463505020893*coeff[1]*qr[6]+151.0463505020893*coeff[1]*ql[6]+134.1640786499874*coeff[1]*qr[3]+134.1640786499874*coeff[1]*ql[3]-268.3281572999747*coeff[1]*qc[3])*Jfac; 
  out[15] += -0.0625*(129.9038105676658*coeff[1]*qr[19]-129.9038105676658*coeff[1]*ql[19]-75.0*coeff[1]*qr[15]-75.0*coeff[1]*ql[15]+150.0*coeff[1]*qc[15])*Jfac; 
  out[16] += -0.03125*(405.0*coeff[1]*qr[16]+405.0*coeff[1]*ql[16]+990.0*coeff[1]*qc[16]-259.8076211353317*coeff[1]*qr[9]+259.8076211353317*coeff[1]*ql[9])*Jfac; 
  out[17] += -0.03125*(405.0*coeff[1]*qr[17]+405.0*coeff[1]*ql[17]+990.0*coeff[1]*qc[17]-259.8076211353317*coeff[1]*qr[13]+259.8076211353317*coeff[1]*ql[13])*Jfac; 
  out[18] += 0.03125*(21.0*coeff[1]*qr[18]+21.0*coeff[1]*ql[18]-1302.0*coeff[1]*qc[18]-151.0463505020892*coeff[1]*qr[10]+151.0463505020892*coeff[1]*ql[10]+134.1640786499874*coeff[1]*qr[5]+134.1640786499874*coeff[1]*ql[5]-268.3281572999748*coeff[1]*qc[5])*Jfac; 
  out[19] += -0.03125*(405.0*coeff[1]*qr[19]+405.0*coeff[1]*ql[19]+990.0*coeff[1]*qc[19]-259.8076211353317*coeff[1]*qr[15]+259.8076211353317*coeff[1]*ql[15])*Jfac; 

  return 0.;

}

