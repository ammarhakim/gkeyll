#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  out[0] += -0.0625*(8.660254037844386*coeff[0]*qr[1]-8.660254037844386*coeff[0]*ql[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*rdx2Sq; 
  out[1] += -0.0625*(7.0*coeff[0]*qr[1]+7.0*coeff[0]*ql[1]+46.0*coeff[0]*qc[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*rdx2Sq; 
  out[2] += -0.0625*(8.660254037844386*coeff[0]*qr[3]-8.660254037844386*coeff[0]*ql[3]-9.0*coeff[0]*qr[2]-9.0*coeff[0]*ql[2]+18.0*coeff[0]*qc[2])*rdx2Sq; 
  out[3] += -0.0625*(7.0*coeff[0]*qr[3]+7.0*coeff[0]*ql[3]+46.0*coeff[0]*qc[3]-8.660254037844386*coeff[0]*qr[2]+8.660254037844386*coeff[0]*ql[2])*rdx2Sq; 
  out[4] += -0.0625*(8.660254037844387*coeff[0]*qr[5]-8.660254037844387*coeff[0]*ql[5]-9.0*coeff[0]*qr[4]-9.0*coeff[0]*ql[4]+18.0*coeff[0]*qc[4])*rdx2Sq; 
  out[5] += -0.0625*(7.0*coeff[0]*qr[5]+7.0*coeff[0]*ql[5]+46.0*coeff[0]*qc[5]-8.660254037844387*coeff[0]*qr[4]+8.660254037844387*coeff[0]*ql[4])*rdx2Sq; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fl[6];
  fl[0] = 0.7071067811865476*(jacobgeo_inv[1]*ql[1]+jacobgeo_inv[0]*ql[0]); 
  fl[1] = 0.7071067811865476*(jacobgeo_inv[0]*ql[1]+ql[0]*jacobgeo_inv[1]); 
  fl[2] = 0.7071067811865476*(jacobgeo_inv[1]*ql[3]+jacobgeo_inv[0]*ql[2]); 
  fl[3] = 0.7071067811865476*(jacobgeo_inv[0]*ql[3]+jacobgeo_inv[1]*ql[2]); 
  fl[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*ql[5]+21.21320343559643*jacobgeo_inv[0]*ql[4]); 
  fl[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*ql[5]+21.21320343559643*jacobgeo_inv[1]*ql[4]); 

  double fc[6];
  fc[0] = 0.7071067811865476*(jacobgeo_inv[1]*qc[1]+jacobgeo_inv[0]*qc[0]); 
  fc[1] = 0.7071067811865476*(jacobgeo_inv[0]*qc[1]+qc[0]*jacobgeo_inv[1]); 
  fc[2] = 0.7071067811865476*(jacobgeo_inv[1]*qc[3]+jacobgeo_inv[0]*qc[2]); 
  fc[3] = 0.7071067811865476*(jacobgeo_inv[0]*qc[3]+jacobgeo_inv[1]*qc[2]); 
  fc[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qc[5]+21.21320343559643*jacobgeo_inv[0]*qc[4]); 
  fc[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*qc[5]+21.21320343559643*jacobgeo_inv[1]*qc[4]); 

  double fr[6];
  fr[0] = 0.7071067811865476*(jacobgeo_inv[1]*qr[1]+jacobgeo_inv[0]*qr[0]); 
  fr[1] = 0.7071067811865476*(jacobgeo_inv[0]*qr[1]+qr[0]*jacobgeo_inv[1]); 
  fr[2] = 0.7071067811865476*(jacobgeo_inv[1]*qr[3]+jacobgeo_inv[0]*qr[2]); 
  fr[3] = 0.7071067811865476*(jacobgeo_inv[0]*qr[3]+jacobgeo_inv[1]*qr[2]); 
  fr[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qr[5]+21.21320343559643*jacobgeo_inv[0]*qr[4]); 
  fr[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*qr[5]+21.21320343559643*jacobgeo_inv[1]*qr[4]); 

  out[0] += -0.03125*((21.21320343559643*coeff[1]+12.24744871391589*coeff[0])*fr[1]+(21.21320343559643*coeff[1]-12.24744871391589*coeff[0])*fl[1]+42.42640687119286*coeff[1]*fc[1]+(22.0454076850486*fl[0]-22.0454076850486*fr[0])*coeff[1]-12.72792206135786*coeff[0]*fr[0]-12.72792206135786*coeff[0]*fl[0]+25.45584412271572*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += -0.03125*((17.14642819948224*coeff[1]+9.899494936611665*coeff[0])*fr[1]+(9.899494936611665*coeff[0]-17.14642819948224*coeff[1])*fl[1]+65.05382386916239*coeff[0]*fc[1]+((-21.21320343559643*fr[0])-21.21320343559643*fl[0]+42.42640687119286*fc[0])*coeff[1]-12.24744871391589*coeff[0]*fr[0]+12.24744871391589*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += -0.03125*((21.21320343559643*coeff[1]+12.24744871391589*coeff[0])*fr[3]+(21.21320343559643*coeff[1]-12.24744871391589*coeff[0])*fl[3]+42.42640687119286*coeff[1]*fc[3]+((-22.0454076850486*coeff[1])-12.72792206135786*coeff[0])*fr[2]+(22.0454076850486*coeff[1]-12.72792206135786*coeff[0])*fl[2]+25.45584412271572*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += -0.03125*((17.14642819948224*coeff[1]+9.899494936611665*coeff[0])*fr[3]+(9.899494936611665*coeff[0]-17.14642819948224*coeff[1])*fl[3]+65.05382386916239*coeff[0]*fc[3]+((-21.21320343559643*coeff[1])-12.24744871391589*coeff[0])*fr[2]+(12.24744871391589*coeff[0]-21.21320343559643*coeff[1])*fl[2]+42.42640687119286*coeff[1]*fc[2])*rdx2Sq; 
  out[4] += -0.03125*((21.21320343559643*coeff[1]+12.24744871391589*coeff[0])*fr[5]+(21.21320343559643*coeff[1]-12.24744871391589*coeff[0])*fl[5]+42.42640687119286*coeff[1]*fc[5]+((-22.0454076850486*coeff[1])-12.72792206135786*coeff[0])*fr[4]+(22.0454076850486*coeff[1]-12.72792206135786*coeff[0])*fl[4]+25.45584412271572*coeff[0]*fc[4])*rdx2Sq; 
  out[5] += -0.03125*((17.14642819948224*coeff[1]+9.899494936611665*coeff[0])*fr[5]+(9.899494936611665*coeff[0]-17.14642819948224*coeff[1])*fl[5]+65.05382386916239*coeff[0]*fc[5]+((-21.21320343559643*coeff[1])-12.24744871391589*coeff[0])*fr[4]+(12.24744871391589*coeff[0]-21.21320343559643*coeff[1])*fl[4]+42.42640687119286*coeff[1]*fc[4])*rdx2Sq; 

  return 0.;

}

