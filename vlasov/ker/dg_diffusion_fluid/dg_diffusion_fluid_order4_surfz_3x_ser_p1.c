#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfz_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],4.);

  out[0] += -0.0625*(25.98076211353316*coeff[2]*qr[3]-25.98076211353316*coeff[2]*ql[3]+((-15.0*qr[0])-15.0*ql[0]+30.0*qc[0])*coeff[2])*Jfac; 
  out[1] += -0.0625*(25.98076211353316*coeff[2]*qr[5]-25.98076211353316*coeff[2]*ql[5]+((-15.0*qr[1])-15.0*ql[1]+30.0*qc[1])*coeff[2])*Jfac; 
  out[2] += -0.0625*(25.98076211353316*coeff[2]*qr[6]-25.98076211353316*coeff[2]*ql[6]-15.0*coeff[2]*qr[2]-15.0*coeff[2]*ql[2]+30.0*coeff[2]*qc[2])*Jfac; 
  out[3] += -0.0625*(33.0*coeff[2]*qr[3]+33.0*coeff[2]*ql[3]+114.0*coeff[2]*qc[3]+(25.98076211353316*ql[0]-25.98076211353316*qr[0])*coeff[2])*Jfac; 
  out[4] += -0.0625*(25.98076211353316*coeff[2]*qr[7]-25.98076211353316*coeff[2]*ql[7]-15.0*coeff[2]*qr[4]-15.0*coeff[2]*ql[4]+30.0*coeff[2]*qc[4])*Jfac; 
  out[5] += -0.0625*(33.0*coeff[2]*qr[5]+33.0*coeff[2]*ql[5]+114.0*coeff[2]*qc[5]+(25.98076211353316*ql[1]-25.98076211353316*qr[1])*coeff[2])*Jfac; 
  out[6] += -0.0625*(33.0*coeff[2]*qr[6]+33.0*coeff[2]*ql[6]+114.0*coeff[2]*qc[6]-25.98076211353316*coeff[2]*qr[2]+25.98076211353316*coeff[2]*ql[2])*Jfac; 
  out[7] += -0.0625*(33.0*coeff[2]*qr[7]+33.0*coeff[2]*ql[7]+114.0*coeff[2]*qc[7]-25.98076211353316*coeff[2]*qr[4]+25.98076211353316*coeff[2]*ql[4])*Jfac; 

  return 0.;

}

