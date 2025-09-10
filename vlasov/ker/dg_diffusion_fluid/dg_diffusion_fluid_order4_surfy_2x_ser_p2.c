#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += 0.0625*(107.3312629199899*coeff[1]*qr[5]+107.3312629199899*coeff[1]*ql[5]-214.6625258399798*coeff[1]*qc[5]-129.9038105676658*coeff[1]*qr[2]+129.9038105676658*coeff[1]*ql[2]+(75.0*qr[0]+75.0*ql[0]-150.0*qc[0])*coeff[1])*Jfac; 
  out[1] += 0.0625*(107.3312629199899*coeff[1]*qr[7]+107.3312629199899*coeff[1]*ql[7]-214.6625258399798*coeff[1]*qc[7]-129.9038105676658*coeff[1]*qr[3]+129.9038105676658*coeff[1]*ql[3]+75.0*coeff[1]*qr[1]+75.0*coeff[1]*ql[1]-150.0*coeff[1]*qc[1])*Jfac; 
  out[2] += 0.03125*(290.4737509655563*coeff[1]*qr[5]-290.4737509655563*coeff[1]*ql[5]-405.0*coeff[1]*qr[2]-405.0*coeff[1]*ql[2]-990.0*coeff[1]*qc[2]+(259.8076211353315*qr[0]-259.8076211353315*ql[0])*coeff[1])*Jfac; 
  out[3] += 0.03125*(290.4737509655563*coeff[1]*qr[7]-290.4737509655563*coeff[1]*ql[7]-405.0*coeff[1]*qr[3]-405.0*coeff[1]*ql[3]-990.0*coeff[1]*qc[3]+259.8076211353315*coeff[1]*qr[1]-259.8076211353315*coeff[1]*ql[1])*Jfac; 
  out[4] += -0.0625*(129.9038105676658*coeff[1]*qr[6]-129.9038105676658*coeff[1]*ql[6]-75.0*coeff[1]*qr[4]-75.0*coeff[1]*ql[4]+150.0*coeff[1]*qc[4])*Jfac; 
  out[5] += 0.03125*(21.0*coeff[1]*qr[5]+21.0*coeff[1]*ql[5]-1302.0*coeff[1]*qc[5]-151.0463505020892*coeff[1]*qr[2]+151.0463505020892*coeff[1]*ql[2]+(134.1640786499874*qr[0]+134.1640786499874*ql[0]-268.3281572999748*qc[0])*coeff[1])*Jfac; 
  out[6] += -0.03125*(405.0*coeff[1]*qr[6]+405.0*coeff[1]*ql[6]+990.0*coeff[1]*qc[6]-259.8076211353317*coeff[1]*qr[4]+259.8076211353317*coeff[1]*ql[4])*Jfac; 
  out[7] += 0.03125*(21.0*coeff[1]*qr[7]+21.0*coeff[1]*ql[7]-1302.0*coeff[1]*qc[7]-151.0463505020893*coeff[1]*qr[3]+151.0463505020893*coeff[1]*ql[3]+134.1640786499874*coeff[1]*qr[1]+134.1640786499874*coeff[1]*ql[1]-268.3281572999747*coeff[1]*qc[1])*Jfac; 

  return 0.;

}

