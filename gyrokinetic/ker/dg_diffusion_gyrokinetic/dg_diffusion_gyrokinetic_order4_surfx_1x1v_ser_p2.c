#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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

  out[0] += 0.0625*(107.3312629199899*coeff[0]*qr[4]+107.3312629199899*coeff[0]*ql[4]-214.6625258399798*coeff[0]*qc[4]-129.9038105676658*coeff[0]*qr[1]+129.9038105676658*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*rdx2Sq; 
  out[1] += 0.03125*(290.4737509655563*coeff[0]*qr[4]-290.4737509655563*coeff[0]*ql[4]-405.0*coeff[0]*qr[1]-405.0*coeff[0]*ql[1]-990.0*coeff[0]*qc[1]+259.8076211353315*coeff[0]*qr[0]-259.8076211353315*coeff[0]*ql[0])*rdx2Sq; 
  out[2] += 0.0625*(107.3312629199899*coeff[0]*qr[6]+107.3312629199899*coeff[0]*ql[6]-214.6625258399798*coeff[0]*qc[6]-129.9038105676658*coeff[0]*qr[3]+129.9038105676658*coeff[0]*ql[3]+75.0*coeff[0]*qr[2]+75.0*coeff[0]*ql[2]-150.0*coeff[0]*qc[2])*rdx2Sq; 
  out[3] += 0.03125*(290.4737509655563*coeff[0]*qr[6]-290.4737509655563*coeff[0]*ql[6]-405.0*coeff[0]*qr[3]-405.0*coeff[0]*ql[3]-990.0*coeff[0]*qc[3]+259.8076211353315*coeff[0]*qr[2]-259.8076211353315*coeff[0]*ql[2])*rdx2Sq; 
  out[4] += 0.03125*(21.0*coeff[0]*qr[4]+21.0*coeff[0]*ql[4]-1302.0*coeff[0]*qc[4]-151.0463505020892*coeff[0]*qr[1]+151.0463505020892*coeff[0]*ql[1]+134.1640786499874*coeff[0]*qr[0]+134.1640786499874*coeff[0]*ql[0]-268.3281572999748*coeff[0]*qc[0])*rdx2Sq; 
  out[5] += -0.0625*(129.9038105676658*coeff[0]*qr[7]-129.9038105676658*coeff[0]*ql[7]-75.0*coeff[0]*qr[5]-75.0*coeff[0]*ql[5]+150.0*coeff[0]*qc[5])*rdx2Sq; 
  out[6] += 0.03125*(21.0*coeff[0]*qr[6]+21.0*coeff[0]*ql[6]-1302.0*coeff[0]*qc[6]-151.0463505020893*coeff[0]*qr[3]+151.0463505020893*coeff[0]*ql[3]+134.1640786499874*coeff[0]*qr[2]+134.1640786499874*coeff[0]*ql[2]-268.3281572999747*coeff[0]*qc[2])*rdx2Sq; 
  out[7] += -0.03125*(405.0*coeff[0]*qr[7]+405.0*coeff[0]*ql[7]+990.0*coeff[0]*qc[7]-259.8076211353317*coeff[0]*qr[5]+259.8076211353317*coeff[0]*ql[5])*rdx2Sq; 

  return 0.;

}

