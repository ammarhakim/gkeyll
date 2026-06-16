#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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
  out[2] += -0.0625*(8.660254037844386*coeff[0]*qr[5]-8.660254037844386*coeff[0]*ql[5]-9.0*coeff[0]*qr[2]-9.0*coeff[0]*ql[2]+18.0*coeff[0]*qc[2])*rdx2Sq; 
  out[3] += -0.0625*(8.660254037844386*coeff[0]*qr[6]-8.660254037844386*coeff[0]*ql[6]-9.0*coeff[0]*qr[3]-9.0*coeff[0]*ql[3]+18.0*coeff[0]*qc[3])*rdx2Sq; 
  out[4] += -0.0625*(8.660254037844386*coeff[0]*qr[8]-8.660254037844386*coeff[0]*ql[8]-9.0*coeff[0]*qr[4]-9.0*coeff[0]*ql[4]+18.0*coeff[0]*qc[4])*rdx2Sq; 
  out[5] += -0.0625*(7.0*coeff[0]*qr[5]+7.0*coeff[0]*ql[5]+46.0*coeff[0]*qc[5]-8.660254037844386*coeff[0]*qr[2]+8.660254037844386*coeff[0]*ql[2])*rdx2Sq; 
  out[6] += -0.0625*(7.0*coeff[0]*qr[6]+7.0*coeff[0]*ql[6]+46.0*coeff[0]*qc[6]-8.660254037844386*coeff[0]*qr[3]+8.660254037844386*coeff[0]*ql[3])*rdx2Sq; 
  out[7] += -0.0625*(8.660254037844386*coeff[0]*qr[11]-8.660254037844386*coeff[0]*ql[11]-9.0*coeff[0]*qr[7]-9.0*coeff[0]*ql[7]+18.0*coeff[0]*qc[7])*rdx2Sq; 
  out[8] += -0.0625*(7.0*coeff[0]*qr[8]+7.0*coeff[0]*ql[8]+46.0*coeff[0]*qc[8]-8.660254037844386*coeff[0]*qr[4]+8.660254037844386*coeff[0]*ql[4])*rdx2Sq; 
  out[9] += -0.0625*(8.660254037844386*coeff[0]*qr[12]-8.660254037844386*coeff[0]*ql[12]-9.0*coeff[0]*qr[9]-9.0*coeff[0]*ql[9]+18.0*coeff[0]*qc[9])*rdx2Sq; 
  out[10] += -0.0625*(8.660254037844386*coeff[0]*qr[13]-8.660254037844386*coeff[0]*ql[13]-9.0*coeff[0]*qr[10]-9.0*coeff[0]*ql[10]+18.0*coeff[0]*qc[10])*rdx2Sq; 
  out[11] += -0.0625*(7.0*coeff[0]*qr[11]+7.0*coeff[0]*ql[11]+46.0*coeff[0]*qc[11]-8.660254037844386*coeff[0]*qr[7]+8.660254037844386*coeff[0]*ql[7])*rdx2Sq; 
  out[12] += -0.0625*(7.0*coeff[0]*qr[12]+7.0*coeff[0]*ql[12]+46.0*coeff[0]*qc[12]-8.660254037844386*coeff[0]*qr[9]+8.660254037844386*coeff[0]*ql[9])*rdx2Sq; 
  out[13] += -0.0625*(7.0*coeff[0]*qr[13]+7.0*coeff[0]*ql[13]+46.0*coeff[0]*qc[13]-8.660254037844386*coeff[0]*qr[10]+8.660254037844386*coeff[0]*ql[10])*rdx2Sq; 
  out[14] += -0.0625*(8.660254037844386*coeff[0]*qr[15]-8.660254037844386*coeff[0]*ql[15]-9.0*coeff[0]*qr[14]-9.0*coeff[0]*ql[14]+18.0*coeff[0]*qc[14])*rdx2Sq; 
  out[15] += -0.0625*(7.0*coeff[0]*qr[15]+7.0*coeff[0]*ql[15]+46.0*coeff[0]*qc[15]-8.660254037844386*coeff[0]*qr[14]+8.660254037844386*coeff[0]*ql[14])*rdx2Sq; 
  out[16] += -0.0625*(8.660254037844387*coeff[0]*qr[17]-8.660254037844387*coeff[0]*ql[17]-9.0*coeff[0]*qr[16]-9.0*coeff[0]*ql[16]+18.0*coeff[0]*qc[16])*rdx2Sq; 
  out[17] += -0.0625*(7.0*coeff[0]*qr[17]+7.0*coeff[0]*ql[17]+46.0*coeff[0]*qc[17]-8.660254037844387*coeff[0]*qr[16]+8.660254037844387*coeff[0]*ql[16])*rdx2Sq; 
  out[18] += -0.0625*(8.660254037844387*coeff[0]*qr[20]-8.660254037844387*coeff[0]*ql[20]-9.0*coeff[0]*qr[18]-9.0*coeff[0]*ql[18]+18.0*coeff[0]*qc[18])*rdx2Sq; 
  out[19] += -0.0625*(8.660254037844387*coeff[0]*qr[21]-8.660254037844387*coeff[0]*ql[21]-9.0*coeff[0]*qr[19]-9.0*coeff[0]*ql[19]+18.0*coeff[0]*qc[19])*rdx2Sq; 
  out[20] += -0.0625*(7.0*coeff[0]*qr[20]+7.0*coeff[0]*ql[20]+46.0*coeff[0]*qc[20]-8.660254037844387*coeff[0]*qr[18]+8.660254037844387*coeff[0]*ql[18])*rdx2Sq; 
  out[21] += -0.0625*(7.0*coeff[0]*qr[21]+7.0*coeff[0]*ql[21]+46.0*coeff[0]*qc[21]-8.660254037844387*coeff[0]*qr[19]+8.660254037844387*coeff[0]*ql[19])*rdx2Sq; 
  out[22] += -0.0625*(8.660254037844387*coeff[0]*qr[23]-8.660254037844387*coeff[0]*ql[23]-9.0*coeff[0]*qr[22]-9.0*coeff[0]*ql[22]+18.0*coeff[0]*qc[22])*rdx2Sq; 
  out[23] += -0.0625*(7.0*coeff[0]*qr[23]+7.0*coeff[0]*ql[23]+46.0*coeff[0]*qc[23]-8.660254037844387*coeff[0]*qr[22]+8.660254037844387*coeff[0]*ql[22])*rdx2Sq; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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

  double fl[24];
  fl[0] = 0.5*(jacobgeo_inv[3]*ql[5]+jacobgeo_inv[2]*ql[2]+jacobgeo_inv[1]*ql[1]+jacobgeo_inv[0]*ql[0]); 
  fl[1] = 0.5*(jacobgeo_inv[2]*ql[5]+ql[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*ql[1]+ql[0]*jacobgeo_inv[1]); 
  fl[2] = 0.5*(jacobgeo_inv[1]*ql[5]+ql[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*ql[2]+ql[0]*jacobgeo_inv[2]); 
  fl[3] = 0.5*(jacobgeo_inv[3]*ql[11]+jacobgeo_inv[2]*ql[7]+jacobgeo_inv[1]*ql[6]+jacobgeo_inv[0]*ql[3]); 
  fl[4] = 0.5*(jacobgeo_inv[3]*ql[12]+jacobgeo_inv[2]*ql[9]+jacobgeo_inv[1]*ql[8]+jacobgeo_inv[0]*ql[4]); 
  fl[5] = 0.5*(jacobgeo_inv[0]*ql[5]+ql[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*ql[2]+ql[1]*jacobgeo_inv[2]); 
  fl[6] = 0.5*(jacobgeo_inv[2]*ql[11]+jacobgeo_inv[3]*ql[7]+jacobgeo_inv[0]*ql[6]+jacobgeo_inv[1]*ql[3]); 
  fl[7] = 0.5*(jacobgeo_inv[1]*ql[11]+jacobgeo_inv[0]*ql[7]+jacobgeo_inv[3]*ql[6]+jacobgeo_inv[2]*ql[3]); 
  fl[8] = 0.5*(jacobgeo_inv[2]*ql[12]+jacobgeo_inv[3]*ql[9]+jacobgeo_inv[0]*ql[8]+jacobgeo_inv[1]*ql[4]); 
  fl[9] = 0.5*(jacobgeo_inv[1]*ql[12]+jacobgeo_inv[0]*ql[9]+jacobgeo_inv[3]*ql[8]+jacobgeo_inv[2]*ql[4]); 
  fl[10] = 0.5*(jacobgeo_inv[3]*ql[15]+jacobgeo_inv[2]*ql[14]+jacobgeo_inv[1]*ql[13]+jacobgeo_inv[0]*ql[10]); 
  fl[11] = 0.5*(jacobgeo_inv[0]*ql[11]+jacobgeo_inv[1]*ql[7]+jacobgeo_inv[2]*ql[6]+jacobgeo_inv[3]*ql[3]); 
  fl[12] = 0.5*(jacobgeo_inv[0]*ql[12]+jacobgeo_inv[1]*ql[9]+jacobgeo_inv[2]*ql[8]+jacobgeo_inv[3]*ql[4]); 
  fl[13] = 0.5*(jacobgeo_inv[2]*ql[15]+jacobgeo_inv[3]*ql[14]+jacobgeo_inv[0]*ql[13]+jacobgeo_inv[1]*ql[10]); 
  fl[14] = 0.5*(jacobgeo_inv[1]*ql[15]+jacobgeo_inv[0]*ql[14]+jacobgeo_inv[3]*ql[13]+jacobgeo_inv[2]*ql[10]); 
  fl[15] = 0.5*(jacobgeo_inv[0]*ql[15]+jacobgeo_inv[1]*ql[14]+jacobgeo_inv[2]*ql[13]+jacobgeo_inv[3]*ql[10]); 
  fl[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*ql[20]+15.0*(jacobgeo_inv[2]*ql[18]+jacobgeo_inv[1]*ql[17])+15.0*jacobgeo_inv[0]*ql[16]); 
  fl[17] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*ql[20]+15.0*(jacobgeo_inv[3]*ql[18]+jacobgeo_inv[0]*ql[17])+15.0*jacobgeo_inv[1]*ql[16]); 
  fl[18] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*ql[20]+15.0*(jacobgeo_inv[0]*ql[18]+jacobgeo_inv[3]*ql[17])+15.0*jacobgeo_inv[2]*ql[16]); 
  fl[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*ql[23]+15.0*(jacobgeo_inv[2]*ql[22]+jacobgeo_inv[1]*ql[21])+15.0*jacobgeo_inv[0]*ql[19]); 
  fl[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*ql[20]+15.0*(jacobgeo_inv[1]*ql[18]+jacobgeo_inv[2]*ql[17])+15.0*jacobgeo_inv[3]*ql[16]); 
  fl[21] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*ql[23]+15.0*(jacobgeo_inv[3]*ql[22]+jacobgeo_inv[0]*ql[21])+15.0*jacobgeo_inv[1]*ql[19]); 
  fl[22] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*ql[23]+15.0*(jacobgeo_inv[0]*ql[22]+jacobgeo_inv[3]*ql[21])+15.0*jacobgeo_inv[2]*ql[19]); 
  fl[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*ql[23]+15.0*(jacobgeo_inv[1]*ql[22]+jacobgeo_inv[2]*ql[21])+15.0*jacobgeo_inv[3]*ql[19]); 

  double fc[24];
  fc[0] = 0.5*(jacobgeo_inv[3]*qc[5]+jacobgeo_inv[2]*qc[2]+jacobgeo_inv[1]*qc[1]+jacobgeo_inv[0]*qc[0]); 
  fc[1] = 0.5*(jacobgeo_inv[2]*qc[5]+qc[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*qc[1]+qc[0]*jacobgeo_inv[1]); 
  fc[2] = 0.5*(jacobgeo_inv[1]*qc[5]+qc[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*qc[2]+qc[0]*jacobgeo_inv[2]); 
  fc[3] = 0.5*(jacobgeo_inv[3]*qc[11]+jacobgeo_inv[2]*qc[7]+jacobgeo_inv[1]*qc[6]+jacobgeo_inv[0]*qc[3]); 
  fc[4] = 0.5*(jacobgeo_inv[3]*qc[12]+jacobgeo_inv[2]*qc[9]+jacobgeo_inv[1]*qc[8]+jacobgeo_inv[0]*qc[4]); 
  fc[5] = 0.5*(jacobgeo_inv[0]*qc[5]+qc[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*qc[2]+qc[1]*jacobgeo_inv[2]); 
  fc[6] = 0.5*(jacobgeo_inv[2]*qc[11]+jacobgeo_inv[3]*qc[7]+jacobgeo_inv[0]*qc[6]+jacobgeo_inv[1]*qc[3]); 
  fc[7] = 0.5*(jacobgeo_inv[1]*qc[11]+jacobgeo_inv[0]*qc[7]+jacobgeo_inv[3]*qc[6]+jacobgeo_inv[2]*qc[3]); 
  fc[8] = 0.5*(jacobgeo_inv[2]*qc[12]+jacobgeo_inv[3]*qc[9]+jacobgeo_inv[0]*qc[8]+jacobgeo_inv[1]*qc[4]); 
  fc[9] = 0.5*(jacobgeo_inv[1]*qc[12]+jacobgeo_inv[0]*qc[9]+jacobgeo_inv[3]*qc[8]+jacobgeo_inv[2]*qc[4]); 
  fc[10] = 0.5*(jacobgeo_inv[3]*qc[15]+jacobgeo_inv[2]*qc[14]+jacobgeo_inv[1]*qc[13]+jacobgeo_inv[0]*qc[10]); 
  fc[11] = 0.5*(jacobgeo_inv[0]*qc[11]+jacobgeo_inv[1]*qc[7]+jacobgeo_inv[2]*qc[6]+jacobgeo_inv[3]*qc[3]); 
  fc[12] = 0.5*(jacobgeo_inv[0]*qc[12]+jacobgeo_inv[1]*qc[9]+jacobgeo_inv[2]*qc[8]+jacobgeo_inv[3]*qc[4]); 
  fc[13] = 0.5*(jacobgeo_inv[2]*qc[15]+jacobgeo_inv[3]*qc[14]+jacobgeo_inv[0]*qc[13]+jacobgeo_inv[1]*qc[10]); 
  fc[14] = 0.5*(jacobgeo_inv[1]*qc[15]+jacobgeo_inv[0]*qc[14]+jacobgeo_inv[3]*qc[13]+jacobgeo_inv[2]*qc[10]); 
  fc[15] = 0.5*(jacobgeo_inv[0]*qc[15]+jacobgeo_inv[1]*qc[14]+jacobgeo_inv[2]*qc[13]+jacobgeo_inv[3]*qc[10]); 
  fc[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qc[20]+15.0*(jacobgeo_inv[2]*qc[18]+jacobgeo_inv[1]*qc[17])+15.0*jacobgeo_inv[0]*qc[16]); 
  fc[17] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qc[20]+15.0*(jacobgeo_inv[3]*qc[18]+jacobgeo_inv[0]*qc[17])+15.0*jacobgeo_inv[1]*qc[16]); 
  fc[18] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qc[20]+15.0*(jacobgeo_inv[0]*qc[18]+jacobgeo_inv[3]*qc[17])+15.0*jacobgeo_inv[2]*qc[16]); 
  fc[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qc[23]+15.0*(jacobgeo_inv[2]*qc[22]+jacobgeo_inv[1]*qc[21])+15.0*jacobgeo_inv[0]*qc[19]); 
  fc[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qc[20]+15.0*(jacobgeo_inv[1]*qc[18]+jacobgeo_inv[2]*qc[17])+15.0*jacobgeo_inv[3]*qc[16]); 
  fc[21] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qc[23]+15.0*(jacobgeo_inv[3]*qc[22]+jacobgeo_inv[0]*qc[21])+15.0*jacobgeo_inv[1]*qc[19]); 
  fc[22] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qc[23]+15.0*(jacobgeo_inv[0]*qc[22]+jacobgeo_inv[3]*qc[21])+15.0*jacobgeo_inv[2]*qc[19]); 
  fc[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qc[23]+15.0*(jacobgeo_inv[1]*qc[22]+jacobgeo_inv[2]*qc[21])+15.0*jacobgeo_inv[3]*qc[19]); 

  double fr[24];
  fr[0] = 0.5*(jacobgeo_inv[3]*qr[5]+jacobgeo_inv[2]*qr[2]+jacobgeo_inv[1]*qr[1]+jacobgeo_inv[0]*qr[0]); 
  fr[1] = 0.5*(jacobgeo_inv[2]*qr[5]+qr[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*qr[1]+qr[0]*jacobgeo_inv[1]); 
  fr[2] = 0.5*(jacobgeo_inv[1]*qr[5]+qr[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*qr[2]+qr[0]*jacobgeo_inv[2]); 
  fr[3] = 0.5*(jacobgeo_inv[3]*qr[11]+jacobgeo_inv[2]*qr[7]+jacobgeo_inv[1]*qr[6]+jacobgeo_inv[0]*qr[3]); 
  fr[4] = 0.5*(jacobgeo_inv[3]*qr[12]+jacobgeo_inv[2]*qr[9]+jacobgeo_inv[1]*qr[8]+jacobgeo_inv[0]*qr[4]); 
  fr[5] = 0.5*(jacobgeo_inv[0]*qr[5]+qr[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*qr[2]+qr[1]*jacobgeo_inv[2]); 
  fr[6] = 0.5*(jacobgeo_inv[2]*qr[11]+jacobgeo_inv[3]*qr[7]+jacobgeo_inv[0]*qr[6]+jacobgeo_inv[1]*qr[3]); 
  fr[7] = 0.5*(jacobgeo_inv[1]*qr[11]+jacobgeo_inv[0]*qr[7]+jacobgeo_inv[3]*qr[6]+jacobgeo_inv[2]*qr[3]); 
  fr[8] = 0.5*(jacobgeo_inv[2]*qr[12]+jacobgeo_inv[3]*qr[9]+jacobgeo_inv[0]*qr[8]+jacobgeo_inv[1]*qr[4]); 
  fr[9] = 0.5*(jacobgeo_inv[1]*qr[12]+jacobgeo_inv[0]*qr[9]+jacobgeo_inv[3]*qr[8]+jacobgeo_inv[2]*qr[4]); 
  fr[10] = 0.5*(jacobgeo_inv[3]*qr[15]+jacobgeo_inv[2]*qr[14]+jacobgeo_inv[1]*qr[13]+jacobgeo_inv[0]*qr[10]); 
  fr[11] = 0.5*(jacobgeo_inv[0]*qr[11]+jacobgeo_inv[1]*qr[7]+jacobgeo_inv[2]*qr[6]+jacobgeo_inv[3]*qr[3]); 
  fr[12] = 0.5*(jacobgeo_inv[0]*qr[12]+jacobgeo_inv[1]*qr[9]+jacobgeo_inv[2]*qr[8]+jacobgeo_inv[3]*qr[4]); 
  fr[13] = 0.5*(jacobgeo_inv[2]*qr[15]+jacobgeo_inv[3]*qr[14]+jacobgeo_inv[0]*qr[13]+jacobgeo_inv[1]*qr[10]); 
  fr[14] = 0.5*(jacobgeo_inv[1]*qr[15]+jacobgeo_inv[0]*qr[14]+jacobgeo_inv[3]*qr[13]+jacobgeo_inv[2]*qr[10]); 
  fr[15] = 0.5*(jacobgeo_inv[0]*qr[15]+jacobgeo_inv[1]*qr[14]+jacobgeo_inv[2]*qr[13]+jacobgeo_inv[3]*qr[10]); 
  fr[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qr[20]+15.0*(jacobgeo_inv[2]*qr[18]+jacobgeo_inv[1]*qr[17])+15.0*jacobgeo_inv[0]*qr[16]); 
  fr[17] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qr[20]+15.0*(jacobgeo_inv[3]*qr[18]+jacobgeo_inv[0]*qr[17])+15.0*jacobgeo_inv[1]*qr[16]); 
  fr[18] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qr[20]+15.0*(jacobgeo_inv[0]*qr[18]+jacobgeo_inv[3]*qr[17])+15.0*jacobgeo_inv[2]*qr[16]); 
  fr[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qr[23]+15.0*(jacobgeo_inv[2]*qr[22]+jacobgeo_inv[1]*qr[21])+15.0*jacobgeo_inv[0]*qr[19]); 
  fr[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qr[20]+15.0*(jacobgeo_inv[1]*qr[18]+jacobgeo_inv[2]*qr[17])+15.0*jacobgeo_inv[3]*qr[16]); 
  fr[21] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qr[23]+15.0*(jacobgeo_inv[3]*qr[22]+jacobgeo_inv[0]*qr[21])+15.0*jacobgeo_inv[1]*qr[19]); 
  fr[22] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qr[23]+15.0*(jacobgeo_inv[0]*qr[22]+jacobgeo_inv[3]*qr[21])+15.0*jacobgeo_inv[2]*qr[19]); 
  fr[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qr[23]+15.0*(jacobgeo_inv[1]*qr[22]+jacobgeo_inv[2]*qr[21])+15.0*jacobgeo_inv[3]*qr[19]); 

  out[0] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*fr[5]+(15.0*coeff[3]-8.660254037844386*coeff[2])*fl[5]+30.0*coeff[3]*fc[5]+(15.58845726811989*fl[2]-15.58845726811989*fr[2])*coeff[3]-9.0*coeff[2]*fr[2]-9.0*coeff[2]*fl[2]+18.0*coeff[2]*fc[2]+(15.0*coeff[1]+8.660254037844386*coeff[0])*fr[1]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[1]+30.0*coeff[1]*fc[1]+(15.58845726811989*fl[0]-15.58845726811989*fr[0])*coeff[1]-9.0*coeff[0]*fr[0]-9.0*coeff[0]*fl[0]+18.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*fr[5]+(7.0*coeff[2]-12.12435565298214*coeff[3])*fl[5]+46.0*coeff[2]*fc[5]+((-15.0*fr[2])-15.0*fl[2]+30.0*fc[2])*coeff[3]-8.660254037844386*coeff[2]*fr[2]+8.660254037844386*coeff[2]*fl[2]+(12.12435565298214*coeff[1]+7.0*coeff[0])*fr[1]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[1]+46.0*coeff[0]*fc[1]+((-15.0*fr[0])-15.0*fl[0]+30.0*fc[0])*coeff[1]-8.660254037844386*coeff[0]*fr[0]+8.660254037844386*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[5]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[5]+30.0*coeff[1]*fc[5]+(15.0*fr[1]+15.0*fl[1]+30.0*fc[1]-15.58845726811989*fr[0]+15.58845726811989*fl[0])*coeff[3]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[2]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[2]+18.0*coeff[0]*fc[2]+(8.660254037844386*fr[1]-8.660254037844386*fl[1]-9.0*fr[0]-9.0*fl[0]+18.0*fc[0])*coeff[2])*rdx2Sq; 
  out[3] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*fr[11]+(15.0*coeff[3]-8.660254037844386*coeff[2])*fl[11]+30.0*coeff[3]*fc[11]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*fr[7]+(15.58845726811989*coeff[3]-9.0*coeff[2])*fl[7]+18.0*coeff[2]*fc[7]+(15.0*coeff[1]+8.660254037844386*coeff[0])*fr[6]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[6]+30.0*coeff[1]*fc[6]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[3]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[3]+18.0*coeff[0]*fc[3])*rdx2Sq; 
  out[4] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*fr[12]+(15.0*coeff[3]-8.660254037844386*coeff[2])*fl[12]+30.0*coeff[3]*fc[12]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*fr[9]+(15.58845726811989*coeff[3]-9.0*coeff[2])*fl[9]+18.0*coeff[2]*fc[9]+(15.0*coeff[1]+8.660254037844386*coeff[0])*fr[8]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[8]+30.0*coeff[1]*fc[8]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[4]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[4]+18.0*coeff[0]*fc[4])*rdx2Sq; 
  out[5] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[5]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[5]+46.0*coeff[0]*fc[5]+(12.12435565298214*fr[1]-12.12435565298214*fl[1]-15.0*fr[0]-15.0*fl[0]+30.0*fc[0])*coeff[3]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[2]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[2]+30.0*coeff[1]*fc[2]+(7.0*fr[1]+7.0*fl[1]+46.0*fc[1]-8.660254037844386*fr[0]+8.660254037844386*fl[0])*coeff[2])*rdx2Sq; 
  out[6] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*fr[11]+(7.0*coeff[2]-12.12435565298214*coeff[3])*fl[11]+46.0*coeff[2]*fc[11]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*fr[7]+(8.660254037844386*coeff[2]-15.0*coeff[3])*fl[7]+30.0*coeff[3]*fc[7]+(12.12435565298214*coeff[1]+7.0*coeff[0])*fr[6]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[6]+46.0*coeff[0]*fc[6]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[3]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[3]+30.0*coeff[1]*fc[3])*rdx2Sq; 
  out[7] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[11]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[11]+30.0*coeff[1]*fc[11]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[7]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[7]+18.0*coeff[0]*fc[7]+(15.0*coeff[3]+8.660254037844386*coeff[2])*fr[6]+(15.0*coeff[3]-8.660254037844386*coeff[2])*fl[6]+30.0*coeff[3]*fc[6]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*fr[3]+(15.58845726811989*coeff[3]-9.0*coeff[2])*fl[3]+18.0*coeff[2]*fc[3])*rdx2Sq; 
  out[8] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*fr[12]+(7.0*coeff[2]-12.12435565298214*coeff[3])*fl[12]+46.0*coeff[2]*fc[12]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*fr[9]+(8.660254037844386*coeff[2]-15.0*coeff[3])*fl[9]+30.0*coeff[3]*fc[9]+(12.12435565298214*coeff[1]+7.0*coeff[0])*fr[8]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[8]+46.0*coeff[0]*fc[8]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[4]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[4]+30.0*coeff[1]*fc[4])*rdx2Sq; 
  out[9] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[12]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[12]+30.0*coeff[1]*fc[12]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[9]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[9]+18.0*coeff[0]*fc[9]+(15.0*coeff[3]+8.660254037844386*coeff[2])*fr[8]+(15.0*coeff[3]-8.660254037844386*coeff[2])*fl[8]+30.0*coeff[3]*fc[8]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*fr[4]+(15.58845726811989*coeff[3]-9.0*coeff[2])*fl[4]+18.0*coeff[2]*fc[4])*rdx2Sq; 
  out[10] += -0.03125*((15.0*coeff[3]+8.660254037844386*coeff[2])*fr[15]+(15.0*coeff[3]-8.660254037844386*coeff[2])*fl[15]+30.0*coeff[3]*fc[15]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*fr[14]+(15.58845726811989*coeff[3]-9.0*coeff[2])*fl[14]+18.0*coeff[2]*fc[14]+(15.0*coeff[1]+8.660254037844386*coeff[0])*fr[13]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[13]+30.0*coeff[1]*fc[13]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[10]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[10]+18.0*coeff[0]*fc[10])*rdx2Sq; 
  out[11] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[11]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[11]+46.0*coeff[0]*fc[11]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[7]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[7]+30.0*coeff[1]*fc[7]+(12.12435565298214*coeff[3]+7.0*coeff[2])*fr[6]+(7.0*coeff[2]-12.12435565298214*coeff[3])*fl[6]+46.0*coeff[2]*fc[6]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*fr[3]+(8.660254037844386*coeff[2]-15.0*coeff[3])*fl[3]+30.0*coeff[3]*fc[3])*rdx2Sq; 
  out[12] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[12]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[12]+46.0*coeff[0]*fc[12]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[9]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[9]+30.0*coeff[1]*fc[9]+(12.12435565298214*coeff[3]+7.0*coeff[2])*fr[8]+(7.0*coeff[2]-12.12435565298214*coeff[3])*fl[8]+46.0*coeff[2]*fc[8]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*fr[4]+(8.660254037844386*coeff[2]-15.0*coeff[3])*fl[4]+30.0*coeff[3]*fc[4])*rdx2Sq; 
  out[13] += -0.03125*((12.12435565298214*coeff[3]+7.0*coeff[2])*fr[15]+(7.0*coeff[2]-12.12435565298214*coeff[3])*fl[15]+46.0*coeff[2]*fc[15]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*fr[14]+(8.660254037844386*coeff[2]-15.0*coeff[3])*fl[14]+30.0*coeff[3]*fc[14]+(12.12435565298214*coeff[1]+7.0*coeff[0])*fr[13]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[13]+46.0*coeff[0]*fc[13]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[10]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[10]+30.0*coeff[1]*fc[10])*rdx2Sq; 
  out[14] += -0.03125*((15.0*coeff[1]+8.660254037844386*coeff[0])*fr[15]+(15.0*coeff[1]-8.660254037844386*coeff[0])*fl[15]+30.0*coeff[1]*fc[15]+((-15.58845726811989*coeff[1])-9.0*coeff[0])*fr[14]+(15.58845726811989*coeff[1]-9.0*coeff[0])*fl[14]+18.0*coeff[0]*fc[14]+(15.0*coeff[3]+8.660254037844386*coeff[2])*fr[13]+(15.0*coeff[3]-8.660254037844386*coeff[2])*fl[13]+30.0*coeff[3]*fc[13]+((-15.58845726811989*coeff[3])-9.0*coeff[2])*fr[10]+(15.58845726811989*coeff[3]-9.0*coeff[2])*fl[10]+18.0*coeff[2]*fc[10])*rdx2Sq; 
  out[15] += -0.03125*((12.12435565298214*coeff[1]+7.0*coeff[0])*fr[15]+(7.0*coeff[0]-12.12435565298214*coeff[1])*fl[15]+46.0*coeff[0]*fc[15]+((-15.0*coeff[1])-8.660254037844386*coeff[0])*fr[14]+(8.660254037844386*coeff[0]-15.0*coeff[1])*fl[14]+30.0*coeff[1]*fc[14]+(12.12435565298214*coeff[3]+7.0*coeff[2])*fr[13]+(7.0*coeff[2]-12.12435565298214*coeff[3])*fl[13]+46.0*coeff[2]*fc[13]+((-15.0*coeff[3])-8.660254037844386*coeff[2])*fr[10]+(8.660254037844386*coeff[2]-15.0*coeff[3])*fl[10]+30.0*coeff[3]*fc[10])*rdx2Sq; 
  out[16] += -0.00625*((75.0*coeff[3]+43.30127018922193*coeff[2])*fr[20]+(75.0*coeff[3]-43.30127018922193*coeff[2])*fl[20]+150.0*coeff[3]*fc[20]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*fr[18]+(77.94228634059948*coeff[3]-45.0*coeff[2])*fl[18]+90.0*coeff[2]*fc[18]+(75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*fr[17]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*fl[17]+150.0*coeff[1]*fc[17]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*fr[16]+(77.94228634059945*coeff[1]-45.0*coeff[0])*fl[16]+90.0*coeff[0]*fc[16])*rdx2Sq; 
  out[17] += -0.002083333333333333*((181.8653347947321*coeff[3]+105.0*coeff[2])*fr[20]+(105.0*coeff[2]-181.8653347947321*coeff[3])*fl[20]+690.0*coeff[2]*fc[20]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*fr[18]+(129.9038105676658*coeff[2]-225.0*coeff[3])*fl[18]+450.0*coeff[3]*fc[18]+(181.8653347947321*coeff[1]+105.0*coeff[0])*fr[17]+(105.0*coeff[0]-181.8653347947321*coeff[1])*fl[17]+690.0*coeff[0]*fc[17]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*fr[16]+(129.9038105676658*coeff[0]-225.0*coeff[1])*fl[16]+450.0000000000001*coeff[1]*fc[16])*rdx2Sq; 
  out[18] += -0.00625*((75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*fr[20]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*fl[20]+150.0*coeff[1]*fc[20]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*fr[18]+(77.94228634059945*coeff[1]-45.0*coeff[0])*fl[18]+90.0*coeff[0]*fc[18]+(75.0*coeff[3]+43.30127018922193*coeff[2])*fr[17]+(75.0*coeff[3]-43.30127018922193*coeff[2])*fl[17]+150.0*coeff[3]*fc[17]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*fr[16]+(77.94228634059948*coeff[3]-45.0*coeff[2])*fl[16]+90.0*coeff[2]*fc[16])*rdx2Sq; 
  out[19] += -0.00625*((75.0*coeff[3]+43.30127018922193*coeff[2])*fr[23]+(75.0*coeff[3]-43.30127018922193*coeff[2])*fl[23]+150.0*coeff[3]*fc[23]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*fr[22]+(77.94228634059948*coeff[3]-45.0*coeff[2])*fl[22]+90.0*coeff[2]*fc[22]+(75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*fr[21]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*fl[21]+150.0*coeff[1]*fc[21]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*fr[19]+(77.94228634059945*coeff[1]-45.0*coeff[0])*fl[19]+90.0*coeff[0]*fc[19])*rdx2Sq; 
  out[20] += -0.002083333333333333*((181.8653347947321*coeff[1]+105.0*coeff[0])*fr[20]+(105.0*coeff[0]-181.8653347947321*coeff[1])*fl[20]+690.0*coeff[0]*fc[20]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*fr[18]+(129.9038105676658*coeff[0]-225.0*coeff[1])*fl[18]+450.0000000000001*coeff[1]*fc[18]+(181.8653347947321*coeff[3]+105.0*coeff[2])*fr[17]+(105.0*coeff[2]-181.8653347947321*coeff[3])*fl[17]+690.0*coeff[2]*fc[17]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*fr[16]+(129.9038105676658*coeff[2]-225.0*coeff[3])*fl[16]+450.0*coeff[3]*fc[16])*rdx2Sq; 
  out[21] += -0.002083333333333333*((181.8653347947321*coeff[3]+105.0*coeff[2])*fr[23]+(105.0*coeff[2]-181.8653347947321*coeff[3])*fl[23]+690.0*coeff[2]*fc[23]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*fr[22]+(129.9038105676658*coeff[2]-225.0*coeff[3])*fl[22]+450.0*coeff[3]*fc[22]+(181.8653347947321*coeff[1]+105.0*coeff[0])*fr[21]+(105.0*coeff[0]-181.8653347947321*coeff[1])*fl[21]+690.0*coeff[0]*fc[21]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*fr[19]+(129.9038105676658*coeff[0]-225.0*coeff[1])*fl[19]+450.0000000000001*coeff[1]*fc[19])*rdx2Sq; 
  out[22] += -0.00625*((75.00000000000001*coeff[1]+43.30127018922195*coeff[0])*fr[23]+(75.00000000000001*coeff[1]-43.30127018922195*coeff[0])*fl[23]+150.0*coeff[1]*fc[23]+((-77.94228634059945*coeff[1])-45.0*coeff[0])*fr[22]+(77.94228634059945*coeff[1]-45.0*coeff[0])*fl[22]+90.0*coeff[0]*fc[22]+(75.0*coeff[3]+43.30127018922193*coeff[2])*fr[21]+(75.0*coeff[3]-43.30127018922193*coeff[2])*fl[21]+150.0*coeff[3]*fc[21]+((-77.94228634059948*coeff[3])-45.0*coeff[2])*fr[19]+(77.94228634059948*coeff[3]-45.0*coeff[2])*fl[19]+90.0*coeff[2]*fc[19])*rdx2Sq; 
  out[23] += -0.002083333333333333*((181.8653347947321*coeff[1]+105.0*coeff[0])*fr[23]+(105.0*coeff[0]-181.8653347947321*coeff[1])*fl[23]+690.0*coeff[0]*fc[23]+((-225.0*coeff[1])-129.9038105676658*coeff[0])*fr[22]+(129.9038105676658*coeff[0]-225.0*coeff[1])*fl[22]+450.0000000000001*coeff[1]*fc[22]+(181.8653347947321*coeff[3]+105.0*coeff[2])*fr[21]+(105.0*coeff[2]-181.8653347947321*coeff[3])*fl[21]+690.0*coeff[2]*fc[21]+((-225.0*coeff[3])-129.9038105676658*coeff[2])*fr[19]+(129.9038105676658*coeff[2]-225.0*coeff[3])*fl[19]+450.0*coeff[3]*fc[19])*rdx2Sq; 

  return 0.;

}

