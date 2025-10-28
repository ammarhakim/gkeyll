#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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
  out[2] += -0.0625*(8.660254037844386*coeff[0]*qr[4]-8.660254037844386*coeff[0]*ql[4]-9.0*coeff[0]*qr[2]-9.0*coeff[0]*ql[2]+18.0*coeff[0]*qc[2])*rdx2Sq; 
  out[3] += -0.0625*(8.660254037844386*coeff[0]*qr[5]-8.660254037844386*coeff[0]*ql[5]-9.0*coeff[0]*qr[3]-9.0*coeff[0]*ql[3]+18.0*coeff[0]*qc[3])*rdx2Sq; 
  out[4] += -0.0625*(7.0*coeff[0]*qr[4]+7.0*coeff[0]*ql[4]+46.0*coeff[0]*qc[4]-8.660254037844386*coeff[0]*qr[2]+8.660254037844386*coeff[0]*ql[2])*rdx2Sq; 
  out[5] += -0.0625*(7.0*coeff[0]*qr[5]+7.0*coeff[0]*ql[5]+46.0*coeff[0]*qc[5]-8.660254037844386*coeff[0]*qr[3]+8.660254037844386*coeff[0]*ql[3])*rdx2Sq; 
  out[6] += -0.0625*(8.660254037844386*coeff[0]*qr[7]-8.660254037844386*coeff[0]*ql[7]-9.0*coeff[0]*qr[6]-9.0*coeff[0]*ql[6]+18.0*coeff[0]*qc[6])*rdx2Sq; 
  out[7] += -0.0625*(7.0*coeff[0]*qr[7]+7.0*coeff[0]*ql[7]+46.0*coeff[0]*qc[7]-8.660254037844386*coeff[0]*qr[6]+8.660254037844386*coeff[0]*ql[6])*rdx2Sq; 
  out[8] += -0.0625*(8.660254037844387*coeff[0]*qr[9]-8.660254037844387*coeff[0]*ql[9]-9.0*coeff[0]*qr[8]-9.0*coeff[0]*ql[8]+18.0*coeff[0]*qc[8])*rdx2Sq; 
  out[9] += -0.0625*(7.0*coeff[0]*qr[9]+7.0*coeff[0]*ql[9]+46.0*coeff[0]*qc[9]-8.660254037844387*coeff[0]*qr[8]+8.660254037844387*coeff[0]*ql[8])*rdx2Sq; 
  out[10] += -0.0625*(8.660254037844387*coeff[0]*qr[11]-8.660254037844387*coeff[0]*ql[11]-9.0*coeff[0]*qr[10]-9.0*coeff[0]*ql[10]+18.0*coeff[0]*qc[10])*rdx2Sq; 
  out[11] += -0.0625*(7.0*coeff[0]*qr[11]+7.0*coeff[0]*ql[11]+46.0*coeff[0]*qc[11]-8.660254037844387*coeff[0]*qr[10]+8.660254037844387*coeff[0]*ql[10])*rdx2Sq; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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

  double fl[12];
  fl[0] = 0.7071067811865475*(jacobgeo_inv[1]*ql[1]+jacobgeo_inv[0]*ql[0]); 
  fl[1] = 0.7071067811865475*(jacobgeo_inv[0]*ql[1]+ql[0]*jacobgeo_inv[1]); 
  fl[2] = 0.7071067811865475*(jacobgeo_inv[1]*ql[4]+jacobgeo_inv[0]*ql[2]); 
  fl[3] = 0.7071067811865475*(jacobgeo_inv[1]*ql[5]+jacobgeo_inv[0]*ql[3]); 
  fl[4] = 0.7071067811865475*(jacobgeo_inv[0]*ql[4]+jacobgeo_inv[1]*ql[2]); 
  fl[5] = 0.7071067811865475*(jacobgeo_inv[0]*ql[5]+jacobgeo_inv[1]*ql[3]); 
  fl[6] = 0.7071067811865475*(jacobgeo_inv[1]*ql[7]+jacobgeo_inv[0]*ql[6]); 
  fl[7] = 0.7071067811865475*(jacobgeo_inv[0]*ql[7]+jacobgeo_inv[1]*ql[6]); 
  fl[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*ql[9]+15.0*jacobgeo_inv[0]*ql[8]); 
  fl[9] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*ql[9]+15.0*jacobgeo_inv[1]*ql[8]); 
  fl[10] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*ql[11]+15.0*jacobgeo_inv[0]*ql[10]); 
  fl[11] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*ql[11]+15.0*jacobgeo_inv[1]*ql[10]); 

  double fc[12];
  fc[0] = 0.7071067811865475*(jacobgeo_inv[1]*qc[1]+jacobgeo_inv[0]*qc[0]); 
  fc[1] = 0.7071067811865475*(jacobgeo_inv[0]*qc[1]+qc[0]*jacobgeo_inv[1]); 
  fc[2] = 0.7071067811865475*(jacobgeo_inv[1]*qc[4]+jacobgeo_inv[0]*qc[2]); 
  fc[3] = 0.7071067811865475*(jacobgeo_inv[1]*qc[5]+jacobgeo_inv[0]*qc[3]); 
  fc[4] = 0.7071067811865475*(jacobgeo_inv[0]*qc[4]+jacobgeo_inv[1]*qc[2]); 
  fc[5] = 0.7071067811865475*(jacobgeo_inv[0]*qc[5]+jacobgeo_inv[1]*qc[3]); 
  fc[6] = 0.7071067811865475*(jacobgeo_inv[1]*qc[7]+jacobgeo_inv[0]*qc[6]); 
  fc[7] = 0.7071067811865475*(jacobgeo_inv[0]*qc[7]+jacobgeo_inv[1]*qc[6]); 
  fc[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qc[9]+15.0*jacobgeo_inv[0]*qc[8]); 
  fc[9] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qc[9]+15.0*jacobgeo_inv[1]*qc[8]); 
  fc[10] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qc[11]+15.0*jacobgeo_inv[0]*qc[10]); 
  fc[11] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qc[11]+15.0*jacobgeo_inv[1]*qc[10]); 

  double fr[12];
  fr[0] = 0.7071067811865475*(jacobgeo_inv[1]*qr[1]+jacobgeo_inv[0]*qr[0]); 
  fr[1] = 0.7071067811865475*(jacobgeo_inv[0]*qr[1]+qr[0]*jacobgeo_inv[1]); 
  fr[2] = 0.7071067811865475*(jacobgeo_inv[1]*qr[4]+jacobgeo_inv[0]*qr[2]); 
  fr[3] = 0.7071067811865475*(jacobgeo_inv[1]*qr[5]+jacobgeo_inv[0]*qr[3]); 
  fr[4] = 0.7071067811865475*(jacobgeo_inv[0]*qr[4]+jacobgeo_inv[1]*qr[2]); 
  fr[5] = 0.7071067811865475*(jacobgeo_inv[0]*qr[5]+jacobgeo_inv[1]*qr[3]); 
  fr[6] = 0.7071067811865475*(jacobgeo_inv[1]*qr[7]+jacobgeo_inv[0]*qr[6]); 
  fr[7] = 0.7071067811865475*(jacobgeo_inv[0]*qr[7]+jacobgeo_inv[1]*qr[6]); 
  fr[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qr[9]+15.0*jacobgeo_inv[0]*qr[8]); 
  fr[9] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qr[9]+15.0*jacobgeo_inv[1]*qr[8]); 
  fr[10] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qr[11]+15.0*jacobgeo_inv[0]*qr[10]); 
  fr[11] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qr[11]+15.0*jacobgeo_inv[1]*qr[10]); 

  out[0] += -0.0441941738241592*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[1]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[1]+30.0*coeff[1]*fc[1]+(15.58845726811989*fl[0]-15.58845726811989*fr[0])*coeff[1]-9.0*coeff[0]*fr[0]-9.0*coeff[0]*fl[0]+18.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += -0.0441941738241592*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[1]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[1]+46.0*coeff[0]*fc[1]+((-15.0*fr[0])-15.0*fl[0]+30.0*fc[0])*coeff[1]-8.660254037844386*coeff[0]*fr[0]+8.660254037844386*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += -0.0441941738241592*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[4]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[4]+30.0*coeff[1]*fc[4]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[2]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[2]+18.0*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += -0.0441941738241592*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[5]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[5]+30.0*coeff[1]*fc[5]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[3]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[3]+18.0*coeff[0]*fc[3])*rdx2Sq; 
  out[4] += -0.0441941738241592*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[4]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[4]+46.0*coeff[0]*fc[4]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[2]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[2]+30.0*coeff[1]*fc[2])*rdx2Sq; 
  out[5] += -0.0441941738241592*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[5]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[5]+46.0*coeff[0]*fc[5]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[3]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[3]+30.0*coeff[1]*fc[3])*rdx2Sq; 
  out[6] += -0.0441941738241592*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[7]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[7]+30.0*coeff[1]*fc[7]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[6]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[6]+18.0*coeff[0]*fc[6])*rdx2Sq; 
  out[7] += -0.0441941738241592*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[7]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[7]+46.0*coeff[0]*fc[7]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[6]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[6]+30.0*coeff[1]*fc[6])*rdx2Sq; 
  out[8] += -0.0441941738241592*((15.0*coeff[1]+8.660254037844387*coeff[0])*fr[9]+(15.0*coeff[1]-8.660254037844387*coeff[0])*fl[9]+30.0*coeff[1]*fc[9]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[8]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[8]+18.0*coeff[0]*fc[8])*rdx2Sq; 
  out[9] += -0.0441941738241592*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[9]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[9]+46.0*coeff[0]*fc[9]+((-15.0*coeff[1])-8.660254037844387*coeff[0])*fr[8]+(8.660254037844387*coeff[0]-15.0*coeff[1])*fl[8]+30.0*coeff[1]*fc[8])*rdx2Sq; 
  out[10] += -0.0441941738241592*((15.0*coeff[1]+8.660254037844387*coeff[0])*fr[11]+(15.0*coeff[1]-8.660254037844387*coeff[0])*fl[11]+30.0*coeff[1]*fc[11]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[10]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[10]+18.0*coeff[0]*fc[10])*rdx2Sq; 
  out[11] += -0.0441941738241592*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[11]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[11]+46.0*coeff[0]*fc[11]+((-15.0*coeff[1])-8.660254037844387*coeff[0])*fr[10]+(8.660254037844387*coeff[0]-15.0*coeff[1])*fl[10]+30.0*coeff[1]*fc[10])*rdx2Sq; 

  return 0.;

}

