#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],2.);

  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[2] = 0.8660254037844386*qSkin[0]*coeff[1]-1.0*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[1]*qSkin[1]-1.0*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[1]*qSkin[3]-1.0*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = 0.8660254037844386*coeff[1]*qSkin[4]-1.0*coeff[1]*qSkin[9]; 
  boundSurf_incr[11] = 0.8660254037844386*coeff[1]*qSkin[6]-1.0*coeff[1]*qSkin[11]; 
  boundSurf_incr[12] = 0.8660254037844386*coeff[1]*qSkin[8]-1.0*coeff[1]*qSkin[12]; 
  boundSurf_incr[14] = 0.8660254037844386*coeff[1]*qSkin[10]-1.0*coeff[1]*qSkin[14]; 
  boundSurf_incr[15] = 0.8660254037844386*coeff[1]*qSkin[13]-1.0*coeff[1]*qSkin[15]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[1]*qSkin[16]-1.0*coeff[1]*qSkin[18]; 
  boundSurf_incr[20] = 0.8660254037844387*coeff[1]*qSkin[17]-1.0*coeff[1]*qSkin[20]; 
  boundSurf_incr[22] = 0.8660254037844387*coeff[1]*qSkin[19]-1.0*coeff[1]*qSkin[22]; 
  boundSurf_incr[23] = 0.8660254037844387*coeff[1]*qSkin[21]-1.0*coeff[1]*qSkin[23]; 

  } else { 

  boundSurf_incr[2] = -(1.0*coeff[1]*qSkin[2])-0.8660254037844386*qSkin[0]*coeff[1]; 
  boundSurf_incr[5] = -(1.0*coeff[1]*qSkin[5])-0.8660254037844386*coeff[1]*qSkin[1]; 
  boundSurf_incr[7] = -(1.0*coeff[1]*qSkin[7])-0.8660254037844386*coeff[1]*qSkin[3]; 
  boundSurf_incr[9] = -(1.0*coeff[1]*qSkin[9])-0.8660254037844386*coeff[1]*qSkin[4]; 
  boundSurf_incr[11] = -(1.0*coeff[1]*qSkin[11])-0.8660254037844386*coeff[1]*qSkin[6]; 
  boundSurf_incr[12] = -(1.0*coeff[1]*qSkin[12])-0.8660254037844386*coeff[1]*qSkin[8]; 
  boundSurf_incr[14] = -(1.0*coeff[1]*qSkin[14])-0.8660254037844386*coeff[1]*qSkin[10]; 
  boundSurf_incr[15] = -(1.0*coeff[1]*qSkin[15])-0.8660254037844386*coeff[1]*qSkin[13]; 
  boundSurf_incr[18] = -(1.0*coeff[1]*qSkin[18])-0.8660254037844387*coeff[1]*qSkin[16]; 
  boundSurf_incr[20] = -(1.0*coeff[1]*qSkin[20])-0.8660254037844387*coeff[1]*qSkin[17]; 
  boundSurf_incr[22] = -(1.0*coeff[1]*qSkin[22])-0.8660254037844387*coeff[1]*qSkin[19]; 
  boundSurf_incr[23] = -(1.0*coeff[1]*qSkin[23])-0.8660254037844387*coeff[1]*qSkin[21]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 
  out[8] += boundSurf_incr[8]*rdx2Sq; 
  out[9] += boundSurf_incr[9]*rdx2Sq; 
  out[10] += boundSurf_incr[10]*rdx2Sq; 
  out[11] += boundSurf_incr[11]*rdx2Sq; 
  out[12] += boundSurf_incr[12]*rdx2Sq; 
  out[13] += boundSurf_incr[13]*rdx2Sq; 
  out[14] += boundSurf_incr[14]*rdx2Sq; 
  out[15] += boundSurf_incr[15]*rdx2Sq; 
  out[16] += boundSurf_incr[16]*rdx2Sq; 
  out[17] += boundSurf_incr[17]*rdx2Sq; 
  out[18] += boundSurf_incr[18]*rdx2Sq; 
  out[19] += boundSurf_incr[19]*rdx2Sq; 
  out[20] += boundSurf_incr[20]*rdx2Sq; 
  out[21] += boundSurf_incr[21]*rdx2Sq; 
  out[22] += boundSurf_incr[22]*rdx2Sq; 
  out[23] += boundSurf_incr[23]*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],2.);

  double fGhost[24];
  fGhost[0] = 0.5*(jacobgeo_inv[3]*qGhost[5]+jacobgeo_inv[2]*qGhost[2]+jacobgeo_inv[1]*qGhost[1]+jacobgeo_inv[0]*qGhost[0]); 
  fGhost[1] = 0.5*(jacobgeo_inv[2]*qGhost[5]+qGhost[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*qGhost[1]+qGhost[0]*jacobgeo_inv[1]); 
  fGhost[2] = 0.5*(jacobgeo_inv[1]*qGhost[5]+qGhost[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*qGhost[2]+qGhost[0]*jacobgeo_inv[2]); 
  fGhost[3] = 0.5*(jacobgeo_inv[3]*qGhost[11]+jacobgeo_inv[2]*qGhost[7]+jacobgeo_inv[1]*qGhost[6]+jacobgeo_inv[0]*qGhost[3]); 
  fGhost[4] = 0.5*(jacobgeo_inv[3]*qGhost[12]+jacobgeo_inv[2]*qGhost[9]+jacobgeo_inv[1]*qGhost[8]+jacobgeo_inv[0]*qGhost[4]); 
  fGhost[5] = 0.5*(jacobgeo_inv[0]*qGhost[5]+qGhost[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*qGhost[2]+qGhost[1]*jacobgeo_inv[2]); 
  fGhost[6] = 0.5*(jacobgeo_inv[2]*qGhost[11]+jacobgeo_inv[3]*qGhost[7]+jacobgeo_inv[0]*qGhost[6]+jacobgeo_inv[1]*qGhost[3]); 
  fGhost[7] = 0.5*(jacobgeo_inv[1]*qGhost[11]+jacobgeo_inv[0]*qGhost[7]+jacobgeo_inv[3]*qGhost[6]+jacobgeo_inv[2]*qGhost[3]); 
  fGhost[8] = 0.5*(jacobgeo_inv[2]*qGhost[12]+jacobgeo_inv[3]*qGhost[9]+jacobgeo_inv[0]*qGhost[8]+jacobgeo_inv[1]*qGhost[4]); 
  fGhost[9] = 0.5*(jacobgeo_inv[1]*qGhost[12]+jacobgeo_inv[0]*qGhost[9]+jacobgeo_inv[3]*qGhost[8]+jacobgeo_inv[2]*qGhost[4]); 
  fGhost[10] = 0.5*(jacobgeo_inv[3]*qGhost[15]+jacobgeo_inv[2]*qGhost[14]+jacobgeo_inv[1]*qGhost[13]+jacobgeo_inv[0]*qGhost[10]); 
  fGhost[11] = 0.5*(jacobgeo_inv[0]*qGhost[11]+jacobgeo_inv[1]*qGhost[7]+jacobgeo_inv[2]*qGhost[6]+jacobgeo_inv[3]*qGhost[3]); 
  fGhost[12] = 0.5*(jacobgeo_inv[0]*qGhost[12]+jacobgeo_inv[1]*qGhost[9]+jacobgeo_inv[2]*qGhost[8]+jacobgeo_inv[3]*qGhost[4]); 
  fGhost[13] = 0.5*(jacobgeo_inv[2]*qGhost[15]+jacobgeo_inv[3]*qGhost[14]+jacobgeo_inv[0]*qGhost[13]+jacobgeo_inv[1]*qGhost[10]); 
  fGhost[14] = 0.5*(jacobgeo_inv[1]*qGhost[15]+jacobgeo_inv[0]*qGhost[14]+jacobgeo_inv[3]*qGhost[13]+jacobgeo_inv[2]*qGhost[10]); 
  fGhost[15] = 0.5*(jacobgeo_inv[0]*qGhost[15]+jacobgeo_inv[1]*qGhost[14]+jacobgeo_inv[2]*qGhost[13]+jacobgeo_inv[3]*qGhost[10]); 
  fGhost[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qGhost[20]+15.000000000000002*(jacobgeo_inv[2]*qGhost[18]+jacobgeo_inv[1]*qGhost[17])+15.0*jacobgeo_inv[0]*qGhost[16]); 
  fGhost[17] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*qGhost[20]+15.0*(jacobgeo_inv[3]*qGhost[18]+jacobgeo_inv[0]*qGhost[17])+15.000000000000002*jacobgeo_inv[1]*qGhost[16]); 
  fGhost[18] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*qGhost[20]+15.0*(jacobgeo_inv[0]*qGhost[18]+jacobgeo_inv[3]*qGhost[17])+15.000000000000002*jacobgeo_inv[2]*qGhost[16]); 
  fGhost[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qGhost[23]+15.000000000000002*(jacobgeo_inv[2]*qGhost[22]+jacobgeo_inv[1]*qGhost[21])+15.0*jacobgeo_inv[0]*qGhost[19]); 
  fGhost[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qGhost[20]+15.000000000000002*(jacobgeo_inv[1]*qGhost[18]+jacobgeo_inv[2]*qGhost[17])+15.0*jacobgeo_inv[3]*qGhost[16]); 
  fGhost[21] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*qGhost[23]+15.0*(jacobgeo_inv[3]*qGhost[22]+jacobgeo_inv[0]*qGhost[21])+15.000000000000002*jacobgeo_inv[1]*qGhost[19]); 
  fGhost[22] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*qGhost[23]+15.0*(jacobgeo_inv[0]*qGhost[22]+jacobgeo_inv[3]*qGhost[21])+15.000000000000002*jacobgeo_inv[2]*qGhost[19]); 
  fGhost[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qGhost[23]+15.000000000000002*(jacobgeo_inv[1]*qGhost[22]+jacobgeo_inv[2]*qGhost[21])+15.0*jacobgeo_inv[3]*qGhost[19]); 

  double fSkin[24];
  fSkin[0] = 0.5*(jacobgeo_inv[3]*qSkin[5]+jacobgeo_inv[2]*qSkin[2]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.5*(jacobgeo_inv[2]*qSkin[5]+qSkin[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1]); 
  fSkin[2] = 0.5*(jacobgeo_inv[1]*qSkin[5]+qSkin[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*qSkin[2]+qSkin[0]*jacobgeo_inv[2]); 
  fSkin[3] = 0.5*(jacobgeo_inv[3]*qSkin[11]+jacobgeo_inv[2]*qSkin[7]+jacobgeo_inv[1]*qSkin[6]+jacobgeo_inv[0]*qSkin[3]); 
  fSkin[4] = 0.5*(jacobgeo_inv[3]*qSkin[12]+jacobgeo_inv[2]*qSkin[9]+jacobgeo_inv[1]*qSkin[8]+jacobgeo_inv[0]*qSkin[4]); 
  fSkin[5] = 0.5*(jacobgeo_inv[0]*qSkin[5]+qSkin[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*qSkin[2]+qSkin[1]*jacobgeo_inv[2]); 
  fSkin[6] = 0.5*(jacobgeo_inv[2]*qSkin[11]+jacobgeo_inv[3]*qSkin[7]+jacobgeo_inv[0]*qSkin[6]+jacobgeo_inv[1]*qSkin[3]); 
  fSkin[7] = 0.5*(jacobgeo_inv[1]*qSkin[11]+jacobgeo_inv[0]*qSkin[7]+jacobgeo_inv[3]*qSkin[6]+jacobgeo_inv[2]*qSkin[3]); 
  fSkin[8] = 0.5*(jacobgeo_inv[2]*qSkin[12]+jacobgeo_inv[3]*qSkin[9]+jacobgeo_inv[0]*qSkin[8]+jacobgeo_inv[1]*qSkin[4]); 
  fSkin[9] = 0.5*(jacobgeo_inv[1]*qSkin[12]+jacobgeo_inv[0]*qSkin[9]+jacobgeo_inv[3]*qSkin[8]+jacobgeo_inv[2]*qSkin[4]); 
  fSkin[10] = 0.5*(jacobgeo_inv[3]*qSkin[15]+jacobgeo_inv[2]*qSkin[14]+jacobgeo_inv[1]*qSkin[13]+jacobgeo_inv[0]*qSkin[10]); 
  fSkin[11] = 0.5*(jacobgeo_inv[0]*qSkin[11]+jacobgeo_inv[1]*qSkin[7]+jacobgeo_inv[2]*qSkin[6]+jacobgeo_inv[3]*qSkin[3]); 
  fSkin[12] = 0.5*(jacobgeo_inv[0]*qSkin[12]+jacobgeo_inv[1]*qSkin[9]+jacobgeo_inv[2]*qSkin[8]+jacobgeo_inv[3]*qSkin[4]); 
  fSkin[13] = 0.5*(jacobgeo_inv[2]*qSkin[15]+jacobgeo_inv[3]*qSkin[14]+jacobgeo_inv[0]*qSkin[13]+jacobgeo_inv[1]*qSkin[10]); 
  fSkin[14] = 0.5*(jacobgeo_inv[1]*qSkin[15]+jacobgeo_inv[0]*qSkin[14]+jacobgeo_inv[3]*qSkin[13]+jacobgeo_inv[2]*qSkin[10]); 
  fSkin[15] = 0.5*(jacobgeo_inv[0]*qSkin[15]+jacobgeo_inv[1]*qSkin[14]+jacobgeo_inv[2]*qSkin[13]+jacobgeo_inv[3]*qSkin[10]); 
  fSkin[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qSkin[20]+15.000000000000002*(jacobgeo_inv[2]*qSkin[18]+jacobgeo_inv[1]*qSkin[17])+15.0*jacobgeo_inv[0]*qSkin[16]); 
  fSkin[17] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*qSkin[20]+15.0*(jacobgeo_inv[3]*qSkin[18]+jacobgeo_inv[0]*qSkin[17])+15.000000000000002*jacobgeo_inv[1]*qSkin[16]); 
  fSkin[18] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*qSkin[20]+15.0*(jacobgeo_inv[0]*qSkin[18]+jacobgeo_inv[3]*qSkin[17])+15.000000000000002*jacobgeo_inv[2]*qSkin[16]); 
  fSkin[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qSkin[23]+15.000000000000002*(jacobgeo_inv[2]*qSkin[22]+jacobgeo_inv[1]*qSkin[21])+15.0*jacobgeo_inv[0]*qSkin[19]); 
  fSkin[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qSkin[20]+15.000000000000002*(jacobgeo_inv[1]*qSkin[18]+jacobgeo_inv[2]*qSkin[17])+15.0*jacobgeo_inv[3]*qSkin[16]); 
  fSkin[21] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*qSkin[23]+15.0*(jacobgeo_inv[3]*qSkin[22]+jacobgeo_inv[0]*qSkin[21])+15.000000000000002*jacobgeo_inv[1]*qSkin[19]); 
  fSkin[22] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*qSkin[23]+15.0*(jacobgeo_inv[0]*qSkin[22]+jacobgeo_inv[3]*qSkin[21])+15.000000000000002*jacobgeo_inv[2]*qSkin[19]); 
  fSkin[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qSkin[23]+15.000000000000002*(jacobgeo_inv[1]*qSkin[22]+jacobgeo_inv[2]*qSkin[21])+15.0*jacobgeo_inv[3]*qSkin[19]); 

  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[2] = 0.8660254037844386*fSkin[5]*coeff[7]-0.75*fSkin[1]*coeff[7]+0.8660254037844386*fSkin[2]*coeff[6]-0.75*fSkin[0]*coeff[6]-0.5*coeff[5]*fSkin[5]+0.4330127018922193*fSkin[1]*coeff[5]-0.5*fSkin[2]*coeff[4]+0.4330127018922193*fSkin[0]*coeff[4]; 
  boundSurf_incr[5] = 0.8660254037844386*fSkin[2]*coeff[7]-0.75*fSkin[0]*coeff[7]+0.8660254037844386*fSkin[5]*coeff[6]-0.75*fSkin[1]*coeff[6]-0.5*coeff[4]*fSkin[5]-0.5*fSkin[2]*coeff[5]+0.4330127018922193*fSkin[0]*coeff[5]+0.4330127018922193*fSkin[1]*coeff[4]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[7]*fSkin[11]-0.5*coeff[5]*fSkin[11]+0.8660254037844386*coeff[6]*fSkin[7]-0.5*coeff[4]*fSkin[7]-0.75*fSkin[6]*coeff[7]+0.4330127018922193*coeff[5]*fSkin[6]-0.75*fSkin[3]*coeff[6]+0.4330127018922193*fSkin[3]*coeff[4]; 
  boundSurf_incr[9] = 0.8660254037844386*coeff[7]*fSkin[12]-0.5*coeff[5]*fSkin[12]+0.8660254037844386*coeff[6]*fSkin[9]-0.5*coeff[4]*fSkin[9]-0.75*coeff[7]*fSkin[8]+0.4330127018922193*coeff[5]*fSkin[8]-0.75*fSkin[4]*coeff[6]+0.4330127018922193*coeff[4]*fSkin[4]; 
  boundSurf_incr[11] = 0.8660254037844386*coeff[6]*fSkin[11]-0.5*coeff[4]*fSkin[11]+0.8660254037844386*coeff[7]*fSkin[7]-0.5*coeff[5]*fSkin[7]-0.75*fSkin[3]*coeff[7]-0.75*coeff[6]*fSkin[6]+0.4330127018922193*coeff[4]*fSkin[6]+0.4330127018922193*fSkin[3]*coeff[5]; 
  boundSurf_incr[12] = 0.8660254037844386*coeff[6]*fSkin[12]-0.5*coeff[4]*fSkin[12]+0.8660254037844386*coeff[7]*fSkin[9]-0.5*coeff[5]*fSkin[9]-0.75*coeff[6]*fSkin[8]+0.4330127018922193*coeff[4]*fSkin[8]-0.75*fSkin[4]*coeff[7]+0.4330127018922193*fSkin[4]*coeff[5]; 
  boundSurf_incr[14] = 0.8660254037844386*coeff[7]*fSkin[15]-0.5*coeff[5]*fSkin[15]+0.8660254037844386*coeff[6]*fSkin[14]-0.5*coeff[4]*fSkin[14]-0.75*coeff[7]*fSkin[13]+0.4330127018922193*coeff[5]*fSkin[13]-0.75*coeff[6]*fSkin[10]+0.4330127018922193*coeff[4]*fSkin[10]; 
  boundSurf_incr[15] = 0.8660254037844386*coeff[6]*fSkin[15]-0.5*coeff[4]*fSkin[15]+0.8660254037844386*coeff[7]*fSkin[14]-0.5*coeff[5]*fSkin[14]-0.75*coeff[6]*fSkin[13]+0.4330127018922193*coeff[4]*fSkin[13]-0.75*coeff[7]*fSkin[10]+0.4330127018922193*coeff[5]*fSkin[10]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[7]*fSkin[20]-0.5000000000000001*coeff[5]*fSkin[20]+0.8660254037844386*coeff[6]*fSkin[18]-0.5*coeff[4]*fSkin[18]-0.75*coeff[7]*fSkin[17]+0.4330127018922193*coeff[5]*fSkin[17]-0.75*coeff[6]*fSkin[16]+0.43301270189221935*coeff[4]*fSkin[16]; 
  boundSurf_incr[20] = 0.8660254037844386*coeff[6]*fSkin[20]-0.5*coeff[4]*fSkin[20]+0.8660254037844387*coeff[7]*fSkin[18]-0.5000000000000001*coeff[5]*fSkin[18]-0.75*coeff[6]*fSkin[17]+0.43301270189221935*coeff[4]*fSkin[17]-0.75*coeff[7]*fSkin[16]+0.4330127018922193*coeff[5]*fSkin[16]; 
  boundSurf_incr[22] = 0.8660254037844387*coeff[7]*fSkin[23]-0.5000000000000001*coeff[5]*fSkin[23]+0.8660254037844386*coeff[6]*fSkin[22]-0.5*coeff[4]*fSkin[22]-0.75*coeff[7]*fSkin[21]+0.4330127018922193*coeff[5]*fSkin[21]-0.75*coeff[6]*fSkin[19]+0.43301270189221935*coeff[4]*fSkin[19]; 
  boundSurf_incr[23] = 0.8660254037844386*coeff[6]*fSkin[23]-0.5*coeff[4]*fSkin[23]+0.8660254037844387*coeff[7]*fSkin[22]-0.5000000000000001*coeff[5]*fSkin[22]-0.75*coeff[6]*fSkin[21]+0.43301270189221935*coeff[4]*fSkin[21]-0.75*coeff[7]*fSkin[19]+0.4330127018922193*coeff[5]*fSkin[19]; 

  } else { 

  boundSurf_incr[2] = -(0.8660254037844386*fSkin[5]*coeff[7])-0.75*fSkin[1]*coeff[7]-0.8660254037844386*fSkin[2]*coeff[6]-0.75*fSkin[0]*coeff[6]-0.5*coeff[5]*fSkin[5]-0.4330127018922193*fSkin[1]*coeff[5]-0.5*fSkin[2]*coeff[4]-0.4330127018922193*fSkin[0]*coeff[4]; 
  boundSurf_incr[5] = -(0.8660254037844386*fSkin[2]*coeff[7])-0.75*fSkin[0]*coeff[7]-0.8660254037844386*fSkin[5]*coeff[6]-0.75*fSkin[1]*coeff[6]-0.5*coeff[4]*fSkin[5]-0.5*fSkin[2]*coeff[5]-0.4330127018922193*fSkin[0]*coeff[5]-0.4330127018922193*fSkin[1]*coeff[4]; 
  boundSurf_incr[7] = -(0.8660254037844386*coeff[7]*fSkin[11])-0.5*coeff[5]*fSkin[11]-0.8660254037844386*coeff[6]*fSkin[7]-0.5*coeff[4]*fSkin[7]-0.75*fSkin[6]*coeff[7]-0.4330127018922193*coeff[5]*fSkin[6]-0.75*fSkin[3]*coeff[6]-0.4330127018922193*fSkin[3]*coeff[4]; 
  boundSurf_incr[9] = -(0.8660254037844386*coeff[7]*fSkin[12])-0.5*coeff[5]*fSkin[12]-0.8660254037844386*coeff[6]*fSkin[9]-0.5*coeff[4]*fSkin[9]-0.75*coeff[7]*fSkin[8]-0.4330127018922193*coeff[5]*fSkin[8]-0.75*fSkin[4]*coeff[6]-0.4330127018922193*coeff[4]*fSkin[4]; 
  boundSurf_incr[11] = -(0.8660254037844386*coeff[6]*fSkin[11])-0.5*coeff[4]*fSkin[11]-0.8660254037844386*coeff[7]*fSkin[7]-0.5*coeff[5]*fSkin[7]-0.75*fSkin[3]*coeff[7]-0.75*coeff[6]*fSkin[6]-0.4330127018922193*coeff[4]*fSkin[6]-0.4330127018922193*fSkin[3]*coeff[5]; 
  boundSurf_incr[12] = -(0.8660254037844386*coeff[6]*fSkin[12])-0.5*coeff[4]*fSkin[12]-0.8660254037844386*coeff[7]*fSkin[9]-0.5*coeff[5]*fSkin[9]-0.75*coeff[6]*fSkin[8]-0.4330127018922193*coeff[4]*fSkin[8]-0.75*fSkin[4]*coeff[7]-0.4330127018922193*fSkin[4]*coeff[5]; 
  boundSurf_incr[14] = -(0.8660254037844386*coeff[7]*fSkin[15])-0.5*coeff[5]*fSkin[15]-0.8660254037844386*coeff[6]*fSkin[14]-0.5*coeff[4]*fSkin[14]-0.75*coeff[7]*fSkin[13]-0.4330127018922193*coeff[5]*fSkin[13]-0.75*coeff[6]*fSkin[10]-0.4330127018922193*coeff[4]*fSkin[10]; 
  boundSurf_incr[15] = -(0.8660254037844386*coeff[6]*fSkin[15])-0.5*coeff[4]*fSkin[15]-0.8660254037844386*coeff[7]*fSkin[14]-0.5*coeff[5]*fSkin[14]-0.75*coeff[6]*fSkin[13]-0.4330127018922193*coeff[4]*fSkin[13]-0.75*coeff[7]*fSkin[10]-0.4330127018922193*coeff[5]*fSkin[10]; 
  boundSurf_incr[18] = -(0.8660254037844387*coeff[7]*fSkin[20])-0.5000000000000001*coeff[5]*fSkin[20]-0.8660254037844386*coeff[6]*fSkin[18]-0.5*coeff[4]*fSkin[18]-0.75*coeff[7]*fSkin[17]-0.4330127018922193*coeff[5]*fSkin[17]-0.75*coeff[6]*fSkin[16]-0.43301270189221935*coeff[4]*fSkin[16]; 
  boundSurf_incr[20] = -(0.8660254037844386*coeff[6]*fSkin[20])-0.5*coeff[4]*fSkin[20]-0.8660254037844387*coeff[7]*fSkin[18]-0.5000000000000001*coeff[5]*fSkin[18]-0.75*coeff[6]*fSkin[17]-0.43301270189221935*coeff[4]*fSkin[17]-0.75*coeff[7]*fSkin[16]-0.4330127018922193*coeff[5]*fSkin[16]; 
  boundSurf_incr[22] = -(0.8660254037844387*coeff[7]*fSkin[23])-0.5000000000000001*coeff[5]*fSkin[23]-0.8660254037844386*coeff[6]*fSkin[22]-0.5*coeff[4]*fSkin[22]-0.75*coeff[7]*fSkin[21]-0.4330127018922193*coeff[5]*fSkin[21]-0.75*coeff[6]*fSkin[19]-0.43301270189221935*coeff[4]*fSkin[19]; 
  boundSurf_incr[23] = -(0.8660254037844386*coeff[6]*fSkin[23])-0.5*coeff[4]*fSkin[23]-0.8660254037844387*coeff[7]*fSkin[22]-0.5000000000000001*coeff[5]*fSkin[22]-0.75*coeff[6]*fSkin[21]-0.43301270189221935*coeff[4]*fSkin[21]-0.75*coeff[7]*fSkin[19]-0.4330127018922193*coeff[5]*fSkin[19]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 
  out[8] += boundSurf_incr[8]*rdx2Sq; 
  out[9] += boundSurf_incr[9]*rdx2Sq; 
  out[10] += boundSurf_incr[10]*rdx2Sq; 
  out[11] += boundSurf_incr[11]*rdx2Sq; 
  out[12] += boundSurf_incr[12]*rdx2Sq; 
  out[13] += boundSurf_incr[13]*rdx2Sq; 
  out[14] += boundSurf_incr[14]*rdx2Sq; 
  out[15] += boundSurf_incr[15]*rdx2Sq; 
  out[16] += boundSurf_incr[16]*rdx2Sq; 
  out[17] += boundSurf_incr[17]*rdx2Sq; 
  out[18] += boundSurf_incr[18]*rdx2Sq; 
  out[19] += boundSurf_incr[19]*rdx2Sq; 
  out[20] += boundSurf_incr[20]*rdx2Sq; 
  out[21] += boundSurf_incr[21]*rdx2Sq; 
  out[22] += boundSurf_incr[22]*rdx2Sq; 
  out[23] += boundSurf_incr[23]*rdx2Sq; 

  return 0.;
}

