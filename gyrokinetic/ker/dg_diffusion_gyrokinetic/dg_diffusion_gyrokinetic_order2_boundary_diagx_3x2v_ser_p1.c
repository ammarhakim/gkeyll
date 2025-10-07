#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_3x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*qSkin[0]-1.0*coeff[0]*qSkin[1]; 
  boundSurf_incr[6] = 0.8660254037844386*coeff[0]*qSkin[2]-1.0*coeff[0]*qSkin[6]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[0]*qSkin[3]-1.0*coeff[0]*qSkin[7]; 
  boundSurf_incr[9] = 0.8660254037844386*coeff[0]*qSkin[4]-1.0*coeff[0]*qSkin[9]; 
  boundSurf_incr[12] = 0.8660254037844386*coeff[0]*qSkin[5]-1.0*coeff[0]*qSkin[12]; 
  boundSurf_incr[16] = 0.8660254037844386*coeff[0]*qSkin[8]-1.0*coeff[0]*qSkin[16]; 
  boundSurf_incr[17] = 0.8660254037844386*coeff[0]*qSkin[10]-1.0*coeff[0]*qSkin[17]; 
  boundSurf_incr[18] = 0.8660254037844386*coeff[0]*qSkin[11]-1.0*coeff[0]*qSkin[18]; 
  boundSurf_incr[20] = 0.8660254037844386*coeff[0]*qSkin[13]-1.0*coeff[0]*qSkin[20]; 
  boundSurf_incr[21] = 0.8660254037844386*coeff[0]*qSkin[14]-1.0*coeff[0]*qSkin[21]; 
  boundSurf_incr[23] = 0.8660254037844386*coeff[0]*qSkin[15]-1.0*coeff[0]*qSkin[23]; 
  boundSurf_incr[26] = 0.8660254037844386*coeff[0]*qSkin[19]-1.0*coeff[0]*qSkin[26]; 
  boundSurf_incr[27] = 0.8660254037844386*coeff[0]*qSkin[22]-1.0*coeff[0]*qSkin[27]; 
  boundSurf_incr[28] = 0.8660254037844386*coeff[0]*qSkin[24]-1.0*coeff[0]*qSkin[28]; 
  boundSurf_incr[29] = 0.8660254037844386*coeff[0]*qSkin[25]-1.0*coeff[0]*qSkin[29]; 
  boundSurf_incr[31] = 0.8660254037844386*coeff[0]*qSkin[30]-1.0*coeff[0]*qSkin[31]; 
  boundSurf_incr[33] = 0.8660254037844387*coeff[0]*qSkin[32]-1.0*coeff[0]*qSkin[33]; 
  boundSurf_incr[37] = 0.8660254037844387*coeff[0]*qSkin[34]-1.0*coeff[0]*qSkin[37]; 
  boundSurf_incr[38] = 0.8660254037844387*coeff[0]*qSkin[35]-1.0*coeff[0]*qSkin[38]; 
  boundSurf_incr[40] = 0.8660254037844387*coeff[0]*qSkin[36]-1.0*coeff[0]*qSkin[40]; 
  boundSurf_incr[43] = 0.8660254037844387*coeff[0]*qSkin[39]-1.0*coeff[0]*qSkin[43]; 
  boundSurf_incr[44] = 0.8660254037844387*coeff[0]*qSkin[41]-1.0*coeff[0]*qSkin[44]; 
  boundSurf_incr[45] = 0.8660254037844387*coeff[0]*qSkin[42]-1.0*coeff[0]*qSkin[45]; 
  boundSurf_incr[47] = 0.8660254037844387*coeff[0]*qSkin[46]-1.0*coeff[0]*qSkin[47]; 

  } else { 

  boundSurf_incr[1] = -(1.0*coeff[0]*qSkin[1])-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[6] = -(1.0*coeff[0]*qSkin[6])-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[7] = -(1.0*coeff[0]*qSkin[7])-0.8660254037844386*coeff[0]*qSkin[3]; 
  boundSurf_incr[9] = -(1.0*coeff[0]*qSkin[9])-0.8660254037844386*coeff[0]*qSkin[4]; 
  boundSurf_incr[12] = -(1.0*coeff[0]*qSkin[12])-0.8660254037844386*coeff[0]*qSkin[5]; 
  boundSurf_incr[16] = -(1.0*coeff[0]*qSkin[16])-0.8660254037844386*coeff[0]*qSkin[8]; 
  boundSurf_incr[17] = -(1.0*coeff[0]*qSkin[17])-0.8660254037844386*coeff[0]*qSkin[10]; 
  boundSurf_incr[18] = -(1.0*coeff[0]*qSkin[18])-0.8660254037844386*coeff[0]*qSkin[11]; 
  boundSurf_incr[20] = -(1.0*coeff[0]*qSkin[20])-0.8660254037844386*coeff[0]*qSkin[13]; 
  boundSurf_incr[21] = -(1.0*coeff[0]*qSkin[21])-0.8660254037844386*coeff[0]*qSkin[14]; 
  boundSurf_incr[23] = -(1.0*coeff[0]*qSkin[23])-0.8660254037844386*coeff[0]*qSkin[15]; 
  boundSurf_incr[26] = -(1.0*coeff[0]*qSkin[26])-0.8660254037844386*coeff[0]*qSkin[19]; 
  boundSurf_incr[27] = -(1.0*coeff[0]*qSkin[27])-0.8660254037844386*coeff[0]*qSkin[22]; 
  boundSurf_incr[28] = -(1.0*coeff[0]*qSkin[28])-0.8660254037844386*coeff[0]*qSkin[24]; 
  boundSurf_incr[29] = -(1.0*coeff[0]*qSkin[29])-0.8660254037844386*coeff[0]*qSkin[25]; 
  boundSurf_incr[31] = -(1.0*coeff[0]*qSkin[31])-0.8660254037844386*coeff[0]*qSkin[30]; 
  boundSurf_incr[33] = -(1.0*coeff[0]*qSkin[33])-0.8660254037844387*coeff[0]*qSkin[32]; 
  boundSurf_incr[37] = -(1.0*coeff[0]*qSkin[37])-0.8660254037844387*coeff[0]*qSkin[34]; 
  boundSurf_incr[38] = -(1.0*coeff[0]*qSkin[38])-0.8660254037844387*coeff[0]*qSkin[35]; 
  boundSurf_incr[40] = -(1.0*coeff[0]*qSkin[40])-0.8660254037844387*coeff[0]*qSkin[36]; 
  boundSurf_incr[43] = -(1.0*coeff[0]*qSkin[43])-0.8660254037844387*coeff[0]*qSkin[39]; 
  boundSurf_incr[44] = -(1.0*coeff[0]*qSkin[44])-0.8660254037844387*coeff[0]*qSkin[41]; 
  boundSurf_incr[45] = -(1.0*coeff[0]*qSkin[45])-0.8660254037844387*coeff[0]*qSkin[42]; 
  boundSurf_incr[47] = -(1.0*coeff[0]*qSkin[47])-0.8660254037844387*coeff[0]*qSkin[46]; 

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
  out[24] += boundSurf_incr[24]*rdx2Sq; 
  out[25] += boundSurf_incr[25]*rdx2Sq; 
  out[26] += boundSurf_incr[26]*rdx2Sq; 
  out[27] += boundSurf_incr[27]*rdx2Sq; 
  out[28] += boundSurf_incr[28]*rdx2Sq; 
  out[29] += boundSurf_incr[29]*rdx2Sq; 
  out[30] += boundSurf_incr[30]*rdx2Sq; 
  out[31] += boundSurf_incr[31]*rdx2Sq; 
  out[32] += boundSurf_incr[32]*rdx2Sq; 
  out[33] += boundSurf_incr[33]*rdx2Sq; 
  out[34] += boundSurf_incr[34]*rdx2Sq; 
  out[35] += boundSurf_incr[35]*rdx2Sq; 
  out[36] += boundSurf_incr[36]*rdx2Sq; 
  out[37] += boundSurf_incr[37]*rdx2Sq; 
  out[38] += boundSurf_incr[38]*rdx2Sq; 
  out[39] += boundSurf_incr[39]*rdx2Sq; 
  out[40] += boundSurf_incr[40]*rdx2Sq; 
  out[41] += boundSurf_incr[41]*rdx2Sq; 
  out[42] += boundSurf_incr[42]*rdx2Sq; 
  out[43] += boundSurf_incr[43]*rdx2Sq; 
  out[44] += boundSurf_incr[44]*rdx2Sq; 
  out[45] += boundSurf_incr[45]*rdx2Sq; 
  out[46] += boundSurf_incr[46]*rdx2Sq; 
  out[47] += boundSurf_incr[47]*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_3x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fGhost[48];
  fGhost[0] = 0.3535533905932737*(jacobgeo_inv[7]*qGhost[16]+jacobgeo_inv[6]*qGhost[8]+jacobgeo_inv[5]*qGhost[7]+jacobgeo_inv[4]*qGhost[6]+jacobgeo_inv[3]*qGhost[3]+jacobgeo_inv[2]*qGhost[2]+jacobgeo_inv[1]*qGhost[1]+jacobgeo_inv[0]*qGhost[0]); 
  fGhost[1] = 0.3535533905932737*(jacobgeo_inv[6]*qGhost[16]+jacobgeo_inv[7]*qGhost[8]+jacobgeo_inv[3]*qGhost[7]+jacobgeo_inv[2]*qGhost[6]+qGhost[3]*jacobgeo_inv[5]+qGhost[2]*jacobgeo_inv[4]+jacobgeo_inv[0]*qGhost[1]+qGhost[0]*jacobgeo_inv[1]); 
  fGhost[2] = 0.3535533905932737*(jacobgeo_inv[5]*qGhost[16]+jacobgeo_inv[3]*qGhost[8]+jacobgeo_inv[7]*qGhost[7]+jacobgeo_inv[1]*qGhost[6]+qGhost[3]*jacobgeo_inv[6]+qGhost[1]*jacobgeo_inv[4]+jacobgeo_inv[0]*qGhost[2]+qGhost[0]*jacobgeo_inv[2]); 
  fGhost[3] = 0.3535533905932737*(jacobgeo_inv[4]*qGhost[16]+jacobgeo_inv[2]*qGhost[8]+jacobgeo_inv[1]*qGhost[7]+qGhost[6]*jacobgeo_inv[7]+qGhost[2]*jacobgeo_inv[6]+qGhost[1]*jacobgeo_inv[5]+jacobgeo_inv[0]*qGhost[3]+qGhost[0]*jacobgeo_inv[3]); 
  fGhost[4] = 0.3535533905932737*(jacobgeo_inv[7]*qGhost[26]+jacobgeo_inv[6]*qGhost[19]+jacobgeo_inv[5]*qGhost[18]+jacobgeo_inv[4]*qGhost[17]+jacobgeo_inv[3]*qGhost[11]+jacobgeo_inv[2]*qGhost[10]+jacobgeo_inv[1]*qGhost[9]+jacobgeo_inv[0]*qGhost[4]); 
  fGhost[5] = 0.3535533905932737*(jacobgeo_inv[7]*qGhost[27]+jacobgeo_inv[6]*qGhost[22]+jacobgeo_inv[5]*qGhost[21]+jacobgeo_inv[4]*qGhost[20]+jacobgeo_inv[3]*qGhost[14]+jacobgeo_inv[2]*qGhost[13]+jacobgeo_inv[1]*qGhost[12]+jacobgeo_inv[0]*qGhost[5]); 
  fGhost[6] = 0.3535533905932737*(jacobgeo_inv[3]*qGhost[16]+jacobgeo_inv[5]*qGhost[8]+jacobgeo_inv[6]*qGhost[7]+qGhost[3]*jacobgeo_inv[7]+jacobgeo_inv[0]*qGhost[6]+qGhost[0]*jacobgeo_inv[4]+jacobgeo_inv[1]*qGhost[2]+qGhost[1]*jacobgeo_inv[2]); 
  fGhost[7] = 0.3535533905932737*(jacobgeo_inv[2]*qGhost[16]+jacobgeo_inv[4]*qGhost[8]+jacobgeo_inv[0]*qGhost[7]+qGhost[2]*jacobgeo_inv[7]+jacobgeo_inv[6]*qGhost[6]+qGhost[0]*jacobgeo_inv[5]+jacobgeo_inv[1]*qGhost[3]+qGhost[1]*jacobgeo_inv[3]); 
  fGhost[8] = 0.3535533905932737*(jacobgeo_inv[1]*qGhost[16]+jacobgeo_inv[0]*qGhost[8]+jacobgeo_inv[4]*qGhost[7]+qGhost[1]*jacobgeo_inv[7]+jacobgeo_inv[5]*qGhost[6]+qGhost[0]*jacobgeo_inv[6]+jacobgeo_inv[2]*qGhost[3]+qGhost[2]*jacobgeo_inv[3]); 
  fGhost[9] = 0.3535533905932737*(jacobgeo_inv[6]*qGhost[26]+jacobgeo_inv[7]*qGhost[19]+jacobgeo_inv[3]*qGhost[18]+jacobgeo_inv[2]*qGhost[17]+jacobgeo_inv[5]*qGhost[11]+jacobgeo_inv[4]*qGhost[10]+jacobgeo_inv[0]*qGhost[9]+jacobgeo_inv[1]*qGhost[4]); 
  fGhost[10] = 0.3535533905932737*(jacobgeo_inv[5]*qGhost[26]+jacobgeo_inv[3]*qGhost[19]+jacobgeo_inv[7]*qGhost[18]+jacobgeo_inv[1]*qGhost[17]+jacobgeo_inv[6]*qGhost[11]+jacobgeo_inv[0]*qGhost[10]+jacobgeo_inv[4]*qGhost[9]+jacobgeo_inv[2]*qGhost[4]); 
  fGhost[11] = 0.3535533905932737*(jacobgeo_inv[4]*qGhost[26]+jacobgeo_inv[2]*qGhost[19]+jacobgeo_inv[1]*qGhost[18]+jacobgeo_inv[7]*qGhost[17]+jacobgeo_inv[0]*qGhost[11]+jacobgeo_inv[6]*qGhost[10]+jacobgeo_inv[5]*qGhost[9]+jacobgeo_inv[3]*qGhost[4]); 
  fGhost[12] = 0.3535533905932737*(jacobgeo_inv[6]*qGhost[27]+jacobgeo_inv[7]*qGhost[22]+jacobgeo_inv[3]*qGhost[21]+jacobgeo_inv[2]*qGhost[20]+jacobgeo_inv[5]*qGhost[14]+jacobgeo_inv[4]*qGhost[13]+jacobgeo_inv[0]*qGhost[12]+jacobgeo_inv[1]*qGhost[5]); 
  fGhost[13] = 0.3535533905932737*(jacobgeo_inv[5]*qGhost[27]+jacobgeo_inv[3]*qGhost[22]+jacobgeo_inv[7]*qGhost[21]+jacobgeo_inv[1]*qGhost[20]+jacobgeo_inv[6]*qGhost[14]+jacobgeo_inv[0]*qGhost[13]+jacobgeo_inv[4]*qGhost[12]+jacobgeo_inv[2]*qGhost[5]); 
  fGhost[14] = 0.3535533905932737*(jacobgeo_inv[4]*qGhost[27]+jacobgeo_inv[2]*qGhost[22]+jacobgeo_inv[1]*qGhost[21]+jacobgeo_inv[7]*qGhost[20]+jacobgeo_inv[0]*qGhost[14]+jacobgeo_inv[6]*qGhost[13]+jacobgeo_inv[5]*qGhost[12]+jacobgeo_inv[3]*qGhost[5]); 
  fGhost[15] = 0.3535533905932737*(jacobgeo_inv[7]*qGhost[31]+jacobgeo_inv[6]*qGhost[30]+jacobgeo_inv[5]*qGhost[29]+jacobgeo_inv[4]*qGhost[28]+jacobgeo_inv[3]*qGhost[25]+jacobgeo_inv[2]*qGhost[24]+jacobgeo_inv[1]*qGhost[23]+jacobgeo_inv[0]*qGhost[15]); 
  fGhost[16] = 0.3535533905932737*(jacobgeo_inv[0]*qGhost[16]+jacobgeo_inv[1]*qGhost[8]+jacobgeo_inv[2]*qGhost[7]+qGhost[0]*jacobgeo_inv[7]+jacobgeo_inv[3]*qGhost[6]+qGhost[1]*jacobgeo_inv[6]+qGhost[2]*jacobgeo_inv[5]+qGhost[3]*jacobgeo_inv[4]); 
  fGhost[17] = 0.3535533905932737*(jacobgeo_inv[3]*qGhost[26]+jacobgeo_inv[5]*qGhost[19]+jacobgeo_inv[6]*qGhost[18]+jacobgeo_inv[0]*qGhost[17]+jacobgeo_inv[7]*qGhost[11]+jacobgeo_inv[1]*qGhost[10]+jacobgeo_inv[2]*qGhost[9]+jacobgeo_inv[4]*qGhost[4]); 
  fGhost[18] = 0.3535533905932737*(jacobgeo_inv[2]*qGhost[26]+jacobgeo_inv[4]*qGhost[19]+jacobgeo_inv[0]*qGhost[18]+jacobgeo_inv[6]*qGhost[17]+jacobgeo_inv[1]*qGhost[11]+jacobgeo_inv[7]*qGhost[10]+jacobgeo_inv[3]*qGhost[9]+qGhost[4]*jacobgeo_inv[5]); 
  fGhost[19] = 0.3535533905932737*(jacobgeo_inv[1]*qGhost[26]+jacobgeo_inv[0]*qGhost[19]+jacobgeo_inv[4]*qGhost[18]+jacobgeo_inv[5]*qGhost[17]+jacobgeo_inv[2]*qGhost[11]+jacobgeo_inv[3]*qGhost[10]+jacobgeo_inv[7]*qGhost[9]+qGhost[4]*jacobgeo_inv[6]); 
  fGhost[20] = 0.3535533905932737*(jacobgeo_inv[3]*qGhost[27]+jacobgeo_inv[5]*qGhost[22]+jacobgeo_inv[6]*qGhost[21]+jacobgeo_inv[0]*qGhost[20]+jacobgeo_inv[7]*qGhost[14]+jacobgeo_inv[1]*qGhost[13]+jacobgeo_inv[2]*qGhost[12]+jacobgeo_inv[4]*qGhost[5]); 
  fGhost[21] = 0.3535533905932737*(jacobgeo_inv[2]*qGhost[27]+jacobgeo_inv[4]*qGhost[22]+jacobgeo_inv[0]*qGhost[21]+jacobgeo_inv[6]*qGhost[20]+jacobgeo_inv[1]*qGhost[14]+jacobgeo_inv[7]*qGhost[13]+jacobgeo_inv[3]*qGhost[12]+jacobgeo_inv[5]*qGhost[5]); 
  fGhost[22] = 0.3535533905932737*(jacobgeo_inv[1]*qGhost[27]+jacobgeo_inv[0]*qGhost[22]+jacobgeo_inv[4]*qGhost[21]+jacobgeo_inv[5]*qGhost[20]+jacobgeo_inv[2]*qGhost[14]+jacobgeo_inv[3]*qGhost[13]+jacobgeo_inv[7]*qGhost[12]+qGhost[5]*jacobgeo_inv[6]); 
  fGhost[23] = 0.3535533905932737*(jacobgeo_inv[6]*qGhost[31]+jacobgeo_inv[7]*qGhost[30]+jacobgeo_inv[3]*qGhost[29]+jacobgeo_inv[2]*qGhost[28]+jacobgeo_inv[5]*qGhost[25]+jacobgeo_inv[4]*qGhost[24]+jacobgeo_inv[0]*qGhost[23]+jacobgeo_inv[1]*qGhost[15]); 
  fGhost[24] = 0.3535533905932737*(jacobgeo_inv[5]*qGhost[31]+jacobgeo_inv[3]*qGhost[30]+jacobgeo_inv[7]*qGhost[29]+jacobgeo_inv[1]*qGhost[28]+jacobgeo_inv[6]*qGhost[25]+jacobgeo_inv[0]*qGhost[24]+jacobgeo_inv[4]*qGhost[23]+jacobgeo_inv[2]*qGhost[15]); 
  fGhost[25] = 0.3535533905932737*(jacobgeo_inv[4]*qGhost[31]+jacobgeo_inv[2]*qGhost[30]+jacobgeo_inv[1]*qGhost[29]+jacobgeo_inv[7]*qGhost[28]+jacobgeo_inv[0]*qGhost[25]+jacobgeo_inv[6]*qGhost[24]+jacobgeo_inv[5]*qGhost[23]+jacobgeo_inv[3]*qGhost[15]); 
  fGhost[26] = 0.3535533905932737*(jacobgeo_inv[0]*qGhost[26]+jacobgeo_inv[1]*qGhost[19]+jacobgeo_inv[2]*qGhost[18]+jacobgeo_inv[3]*qGhost[17]+jacobgeo_inv[4]*qGhost[11]+jacobgeo_inv[5]*qGhost[10]+jacobgeo_inv[6]*qGhost[9]+qGhost[4]*jacobgeo_inv[7]); 
  fGhost[27] = 0.3535533905932737*(jacobgeo_inv[0]*qGhost[27]+jacobgeo_inv[1]*qGhost[22]+jacobgeo_inv[2]*qGhost[21]+jacobgeo_inv[3]*qGhost[20]+jacobgeo_inv[4]*qGhost[14]+jacobgeo_inv[5]*qGhost[13]+jacobgeo_inv[6]*qGhost[12]+qGhost[5]*jacobgeo_inv[7]); 
  fGhost[28] = 0.3535533905932737*(jacobgeo_inv[3]*qGhost[31]+jacobgeo_inv[5]*qGhost[30]+jacobgeo_inv[6]*qGhost[29]+jacobgeo_inv[0]*qGhost[28]+jacobgeo_inv[7]*qGhost[25]+jacobgeo_inv[1]*qGhost[24]+jacobgeo_inv[2]*qGhost[23]+jacobgeo_inv[4]*qGhost[15]); 
  fGhost[29] = 0.3535533905932737*(jacobgeo_inv[2]*qGhost[31]+jacobgeo_inv[4]*qGhost[30]+jacobgeo_inv[0]*qGhost[29]+jacobgeo_inv[6]*qGhost[28]+jacobgeo_inv[1]*qGhost[25]+jacobgeo_inv[7]*qGhost[24]+jacobgeo_inv[3]*qGhost[23]+jacobgeo_inv[5]*qGhost[15]); 
  fGhost[30] = 0.3535533905932737*(jacobgeo_inv[1]*qGhost[31]+jacobgeo_inv[0]*qGhost[30]+jacobgeo_inv[4]*qGhost[29]+jacobgeo_inv[5]*qGhost[28]+jacobgeo_inv[2]*qGhost[25]+jacobgeo_inv[3]*qGhost[24]+jacobgeo_inv[7]*qGhost[23]+jacobgeo_inv[6]*qGhost[15]); 
  fGhost[31] = 0.3535533905932737*(jacobgeo_inv[0]*qGhost[31]+jacobgeo_inv[1]*qGhost[30]+jacobgeo_inv[2]*qGhost[29]+jacobgeo_inv[3]*qGhost[28]+jacobgeo_inv[4]*qGhost[25]+jacobgeo_inv[5]*qGhost[24]+jacobgeo_inv[6]*qGhost[23]+jacobgeo_inv[7]*qGhost[15]); 
  fGhost[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*qGhost[43]+15.0*(jacobgeo_inv[6]*qGhost[39]+jacobgeo_inv[5]*qGhost[38]+jacobgeo_inv[4]*qGhost[37])+15.000000000000002*(jacobgeo_inv[3]*qGhost[35]+jacobgeo_inv[2]*qGhost[34]+jacobgeo_inv[1]*qGhost[33])+15.0*jacobgeo_inv[0]*qGhost[32]); 
  fGhost[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*qGhost[43]+15.000000000000002*(jacobgeo_inv[7]*qGhost[39]+jacobgeo_inv[3]*qGhost[38]+jacobgeo_inv[2]*qGhost[37])+15.0*(jacobgeo_inv[5]*qGhost[35]+jacobgeo_inv[4]*qGhost[34]+jacobgeo_inv[0]*qGhost[33])+15.000000000000002*jacobgeo_inv[1]*qGhost[32]); 
  fGhost[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*qGhost[43]+15.000000000000002*(jacobgeo_inv[3]*qGhost[39]+jacobgeo_inv[7]*qGhost[38]+jacobgeo_inv[1]*qGhost[37])+15.0*(jacobgeo_inv[6]*qGhost[35]+jacobgeo_inv[0]*qGhost[34]+jacobgeo_inv[4]*qGhost[33])+15.000000000000002*jacobgeo_inv[2]*qGhost[32]); 
  fGhost[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*qGhost[43]+15.000000000000002*(jacobgeo_inv[2]*qGhost[39]+jacobgeo_inv[1]*qGhost[38]+jacobgeo_inv[7]*qGhost[37])+15.0*(jacobgeo_inv[0]*qGhost[35]+jacobgeo_inv[6]*qGhost[34]+jacobgeo_inv[5]*qGhost[33])+15.000000000000002*jacobgeo_inv[3]*qGhost[32]); 
  fGhost[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*qGhost[47]+15.0*(jacobgeo_inv[6]*qGhost[46]+jacobgeo_inv[5]*qGhost[45]+jacobgeo_inv[4]*qGhost[44])+15.000000000000002*(jacobgeo_inv[3]*qGhost[42]+jacobgeo_inv[2]*qGhost[41]+jacobgeo_inv[1]*qGhost[40])+15.0*jacobgeo_inv[0]*qGhost[36]); 
  fGhost[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*qGhost[43]+15.0*(jacobgeo_inv[5]*qGhost[39]+jacobgeo_inv[6]*qGhost[38]+jacobgeo_inv[0]*qGhost[37])+15.000000000000002*(jacobgeo_inv[7]*qGhost[35]+jacobgeo_inv[1]*qGhost[34]+jacobgeo_inv[2]*qGhost[33])+15.0*jacobgeo_inv[4]*qGhost[32]); 
  fGhost[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*qGhost[43]+15.0*(jacobgeo_inv[4]*qGhost[39]+jacobgeo_inv[0]*qGhost[38]+jacobgeo_inv[6]*qGhost[37])+15.000000000000002*(jacobgeo_inv[1]*qGhost[35]+jacobgeo_inv[7]*qGhost[34]+jacobgeo_inv[3]*qGhost[33])+15.0*jacobgeo_inv[5]*qGhost[32]); 
  fGhost[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*qGhost[43]+15.0*(jacobgeo_inv[0]*qGhost[39]+jacobgeo_inv[4]*qGhost[38]+jacobgeo_inv[5]*qGhost[37])+15.000000000000002*(jacobgeo_inv[2]*qGhost[35]+jacobgeo_inv[3]*qGhost[34]+jacobgeo_inv[7]*qGhost[33])+15.0*jacobgeo_inv[6]*qGhost[32]); 
  fGhost[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*qGhost[47]+15.000000000000002*(jacobgeo_inv[7]*qGhost[46]+jacobgeo_inv[3]*qGhost[45]+jacobgeo_inv[2]*qGhost[44])+15.0*(jacobgeo_inv[5]*qGhost[42]+jacobgeo_inv[4]*qGhost[41]+jacobgeo_inv[0]*qGhost[40])+15.000000000000002*jacobgeo_inv[1]*qGhost[36]); 
  fGhost[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*qGhost[47]+15.000000000000002*(jacobgeo_inv[3]*qGhost[46]+jacobgeo_inv[7]*qGhost[45]+jacobgeo_inv[1]*qGhost[44])+15.0*(jacobgeo_inv[6]*qGhost[42]+jacobgeo_inv[0]*qGhost[41]+jacobgeo_inv[4]*qGhost[40])+15.000000000000002*jacobgeo_inv[2]*qGhost[36]); 
  fGhost[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*qGhost[47]+15.000000000000002*(jacobgeo_inv[2]*qGhost[46]+jacobgeo_inv[1]*qGhost[45]+jacobgeo_inv[7]*qGhost[44])+15.0*(jacobgeo_inv[0]*qGhost[42]+jacobgeo_inv[6]*qGhost[41]+jacobgeo_inv[5]*qGhost[40])+15.000000000000002*jacobgeo_inv[3]*qGhost[36]); 
  fGhost[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*qGhost[43]+15.000000000000002*(jacobgeo_inv[1]*qGhost[39]+jacobgeo_inv[2]*qGhost[38]+jacobgeo_inv[3]*qGhost[37])+15.0*(jacobgeo_inv[4]*qGhost[35]+jacobgeo_inv[5]*qGhost[34]+jacobgeo_inv[6]*qGhost[33])+15.000000000000002*jacobgeo_inv[7]*qGhost[32]); 
  fGhost[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*qGhost[47]+15.0*(jacobgeo_inv[5]*qGhost[46]+jacobgeo_inv[6]*qGhost[45]+jacobgeo_inv[0]*qGhost[44])+15.000000000000002*(jacobgeo_inv[7]*qGhost[42]+jacobgeo_inv[1]*qGhost[41]+jacobgeo_inv[2]*qGhost[40])+15.0*jacobgeo_inv[4]*qGhost[36]); 
  fGhost[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*qGhost[47]+15.0*(jacobgeo_inv[4]*qGhost[46]+jacobgeo_inv[0]*qGhost[45]+jacobgeo_inv[6]*qGhost[44])+15.000000000000002*(jacobgeo_inv[1]*qGhost[42]+jacobgeo_inv[7]*qGhost[41]+jacobgeo_inv[3]*qGhost[40])+15.0*jacobgeo_inv[5]*qGhost[36]); 
  fGhost[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*qGhost[47]+15.0*(jacobgeo_inv[0]*qGhost[46]+jacobgeo_inv[4]*qGhost[45]+jacobgeo_inv[5]*qGhost[44])+15.000000000000002*(jacobgeo_inv[2]*qGhost[42]+jacobgeo_inv[3]*qGhost[41]+jacobgeo_inv[7]*qGhost[40])+15.0*jacobgeo_inv[6]*qGhost[36]); 
  fGhost[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*qGhost[47]+15.000000000000002*(jacobgeo_inv[1]*qGhost[46]+jacobgeo_inv[2]*qGhost[45]+jacobgeo_inv[3]*qGhost[44])+15.0*(jacobgeo_inv[4]*qGhost[42]+jacobgeo_inv[5]*qGhost[41]+jacobgeo_inv[6]*qGhost[40])+15.000000000000002*jacobgeo_inv[7]*qGhost[36]); 

  double fSkin[48];
  fSkin[0] = 0.3535533905932737*(jacobgeo_inv[7]*qSkin[16]+jacobgeo_inv[6]*qSkin[8]+jacobgeo_inv[5]*qSkin[7]+jacobgeo_inv[4]*qSkin[6]+jacobgeo_inv[3]*qSkin[3]+jacobgeo_inv[2]*qSkin[2]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.3535533905932737*(jacobgeo_inv[6]*qSkin[16]+jacobgeo_inv[7]*qSkin[8]+jacobgeo_inv[3]*qSkin[7]+jacobgeo_inv[2]*qSkin[6]+qSkin[3]*jacobgeo_inv[5]+qSkin[2]*jacobgeo_inv[4]+jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1]); 
  fSkin[2] = 0.3535533905932737*(jacobgeo_inv[5]*qSkin[16]+jacobgeo_inv[3]*qSkin[8]+jacobgeo_inv[7]*qSkin[7]+jacobgeo_inv[1]*qSkin[6]+qSkin[3]*jacobgeo_inv[6]+qSkin[1]*jacobgeo_inv[4]+jacobgeo_inv[0]*qSkin[2]+qSkin[0]*jacobgeo_inv[2]); 
  fSkin[3] = 0.3535533905932737*(jacobgeo_inv[4]*qSkin[16]+jacobgeo_inv[2]*qSkin[8]+jacobgeo_inv[1]*qSkin[7]+qSkin[6]*jacobgeo_inv[7]+qSkin[2]*jacobgeo_inv[6]+qSkin[1]*jacobgeo_inv[5]+jacobgeo_inv[0]*qSkin[3]+qSkin[0]*jacobgeo_inv[3]); 
  fSkin[4] = 0.3535533905932737*(jacobgeo_inv[7]*qSkin[26]+jacobgeo_inv[6]*qSkin[19]+jacobgeo_inv[5]*qSkin[18]+jacobgeo_inv[4]*qSkin[17]+jacobgeo_inv[3]*qSkin[11]+jacobgeo_inv[2]*qSkin[10]+jacobgeo_inv[1]*qSkin[9]+jacobgeo_inv[0]*qSkin[4]); 
  fSkin[5] = 0.3535533905932737*(jacobgeo_inv[7]*qSkin[27]+jacobgeo_inv[6]*qSkin[22]+jacobgeo_inv[5]*qSkin[21]+jacobgeo_inv[4]*qSkin[20]+jacobgeo_inv[3]*qSkin[14]+jacobgeo_inv[2]*qSkin[13]+jacobgeo_inv[1]*qSkin[12]+jacobgeo_inv[0]*qSkin[5]); 
  fSkin[6] = 0.3535533905932737*(jacobgeo_inv[3]*qSkin[16]+jacobgeo_inv[5]*qSkin[8]+jacobgeo_inv[6]*qSkin[7]+qSkin[3]*jacobgeo_inv[7]+jacobgeo_inv[0]*qSkin[6]+qSkin[0]*jacobgeo_inv[4]+jacobgeo_inv[1]*qSkin[2]+qSkin[1]*jacobgeo_inv[2]); 
  fSkin[7] = 0.3535533905932737*(jacobgeo_inv[2]*qSkin[16]+jacobgeo_inv[4]*qSkin[8]+jacobgeo_inv[0]*qSkin[7]+qSkin[2]*jacobgeo_inv[7]+jacobgeo_inv[6]*qSkin[6]+qSkin[0]*jacobgeo_inv[5]+jacobgeo_inv[1]*qSkin[3]+qSkin[1]*jacobgeo_inv[3]); 
  fSkin[8] = 0.3535533905932737*(jacobgeo_inv[1]*qSkin[16]+jacobgeo_inv[0]*qSkin[8]+jacobgeo_inv[4]*qSkin[7]+qSkin[1]*jacobgeo_inv[7]+jacobgeo_inv[5]*qSkin[6]+qSkin[0]*jacobgeo_inv[6]+jacobgeo_inv[2]*qSkin[3]+qSkin[2]*jacobgeo_inv[3]); 
  fSkin[9] = 0.3535533905932737*(jacobgeo_inv[6]*qSkin[26]+jacobgeo_inv[7]*qSkin[19]+jacobgeo_inv[3]*qSkin[18]+jacobgeo_inv[2]*qSkin[17]+jacobgeo_inv[5]*qSkin[11]+jacobgeo_inv[4]*qSkin[10]+jacobgeo_inv[0]*qSkin[9]+jacobgeo_inv[1]*qSkin[4]); 
  fSkin[10] = 0.3535533905932737*(jacobgeo_inv[5]*qSkin[26]+jacobgeo_inv[3]*qSkin[19]+jacobgeo_inv[7]*qSkin[18]+jacobgeo_inv[1]*qSkin[17]+jacobgeo_inv[6]*qSkin[11]+jacobgeo_inv[0]*qSkin[10]+jacobgeo_inv[4]*qSkin[9]+jacobgeo_inv[2]*qSkin[4]); 
  fSkin[11] = 0.3535533905932737*(jacobgeo_inv[4]*qSkin[26]+jacobgeo_inv[2]*qSkin[19]+jacobgeo_inv[1]*qSkin[18]+jacobgeo_inv[7]*qSkin[17]+jacobgeo_inv[0]*qSkin[11]+jacobgeo_inv[6]*qSkin[10]+jacobgeo_inv[5]*qSkin[9]+jacobgeo_inv[3]*qSkin[4]); 
  fSkin[12] = 0.3535533905932737*(jacobgeo_inv[6]*qSkin[27]+jacobgeo_inv[7]*qSkin[22]+jacobgeo_inv[3]*qSkin[21]+jacobgeo_inv[2]*qSkin[20]+jacobgeo_inv[5]*qSkin[14]+jacobgeo_inv[4]*qSkin[13]+jacobgeo_inv[0]*qSkin[12]+jacobgeo_inv[1]*qSkin[5]); 
  fSkin[13] = 0.3535533905932737*(jacobgeo_inv[5]*qSkin[27]+jacobgeo_inv[3]*qSkin[22]+jacobgeo_inv[7]*qSkin[21]+jacobgeo_inv[1]*qSkin[20]+jacobgeo_inv[6]*qSkin[14]+jacobgeo_inv[0]*qSkin[13]+jacobgeo_inv[4]*qSkin[12]+jacobgeo_inv[2]*qSkin[5]); 
  fSkin[14] = 0.3535533905932737*(jacobgeo_inv[4]*qSkin[27]+jacobgeo_inv[2]*qSkin[22]+jacobgeo_inv[1]*qSkin[21]+jacobgeo_inv[7]*qSkin[20]+jacobgeo_inv[0]*qSkin[14]+jacobgeo_inv[6]*qSkin[13]+jacobgeo_inv[5]*qSkin[12]+jacobgeo_inv[3]*qSkin[5]); 
  fSkin[15] = 0.3535533905932737*(jacobgeo_inv[7]*qSkin[31]+jacobgeo_inv[6]*qSkin[30]+jacobgeo_inv[5]*qSkin[29]+jacobgeo_inv[4]*qSkin[28]+jacobgeo_inv[3]*qSkin[25]+jacobgeo_inv[2]*qSkin[24]+jacobgeo_inv[1]*qSkin[23]+jacobgeo_inv[0]*qSkin[15]); 
  fSkin[16] = 0.3535533905932737*(jacobgeo_inv[0]*qSkin[16]+jacobgeo_inv[1]*qSkin[8]+jacobgeo_inv[2]*qSkin[7]+qSkin[0]*jacobgeo_inv[7]+jacobgeo_inv[3]*qSkin[6]+qSkin[1]*jacobgeo_inv[6]+qSkin[2]*jacobgeo_inv[5]+qSkin[3]*jacobgeo_inv[4]); 
  fSkin[17] = 0.3535533905932737*(jacobgeo_inv[3]*qSkin[26]+jacobgeo_inv[5]*qSkin[19]+jacobgeo_inv[6]*qSkin[18]+jacobgeo_inv[0]*qSkin[17]+jacobgeo_inv[7]*qSkin[11]+jacobgeo_inv[1]*qSkin[10]+jacobgeo_inv[2]*qSkin[9]+jacobgeo_inv[4]*qSkin[4]); 
  fSkin[18] = 0.3535533905932737*(jacobgeo_inv[2]*qSkin[26]+jacobgeo_inv[4]*qSkin[19]+jacobgeo_inv[0]*qSkin[18]+jacobgeo_inv[6]*qSkin[17]+jacobgeo_inv[1]*qSkin[11]+jacobgeo_inv[7]*qSkin[10]+jacobgeo_inv[3]*qSkin[9]+qSkin[4]*jacobgeo_inv[5]); 
  fSkin[19] = 0.3535533905932737*(jacobgeo_inv[1]*qSkin[26]+jacobgeo_inv[0]*qSkin[19]+jacobgeo_inv[4]*qSkin[18]+jacobgeo_inv[5]*qSkin[17]+jacobgeo_inv[2]*qSkin[11]+jacobgeo_inv[3]*qSkin[10]+jacobgeo_inv[7]*qSkin[9]+qSkin[4]*jacobgeo_inv[6]); 
  fSkin[20] = 0.3535533905932737*(jacobgeo_inv[3]*qSkin[27]+jacobgeo_inv[5]*qSkin[22]+jacobgeo_inv[6]*qSkin[21]+jacobgeo_inv[0]*qSkin[20]+jacobgeo_inv[7]*qSkin[14]+jacobgeo_inv[1]*qSkin[13]+jacobgeo_inv[2]*qSkin[12]+jacobgeo_inv[4]*qSkin[5]); 
  fSkin[21] = 0.3535533905932737*(jacobgeo_inv[2]*qSkin[27]+jacobgeo_inv[4]*qSkin[22]+jacobgeo_inv[0]*qSkin[21]+jacobgeo_inv[6]*qSkin[20]+jacobgeo_inv[1]*qSkin[14]+jacobgeo_inv[7]*qSkin[13]+jacobgeo_inv[3]*qSkin[12]+jacobgeo_inv[5]*qSkin[5]); 
  fSkin[22] = 0.3535533905932737*(jacobgeo_inv[1]*qSkin[27]+jacobgeo_inv[0]*qSkin[22]+jacobgeo_inv[4]*qSkin[21]+jacobgeo_inv[5]*qSkin[20]+jacobgeo_inv[2]*qSkin[14]+jacobgeo_inv[3]*qSkin[13]+jacobgeo_inv[7]*qSkin[12]+qSkin[5]*jacobgeo_inv[6]); 
  fSkin[23] = 0.3535533905932737*(jacobgeo_inv[6]*qSkin[31]+jacobgeo_inv[7]*qSkin[30]+jacobgeo_inv[3]*qSkin[29]+jacobgeo_inv[2]*qSkin[28]+jacobgeo_inv[5]*qSkin[25]+jacobgeo_inv[4]*qSkin[24]+jacobgeo_inv[0]*qSkin[23]+jacobgeo_inv[1]*qSkin[15]); 
  fSkin[24] = 0.3535533905932737*(jacobgeo_inv[5]*qSkin[31]+jacobgeo_inv[3]*qSkin[30]+jacobgeo_inv[7]*qSkin[29]+jacobgeo_inv[1]*qSkin[28]+jacobgeo_inv[6]*qSkin[25]+jacobgeo_inv[0]*qSkin[24]+jacobgeo_inv[4]*qSkin[23]+jacobgeo_inv[2]*qSkin[15]); 
  fSkin[25] = 0.3535533905932737*(jacobgeo_inv[4]*qSkin[31]+jacobgeo_inv[2]*qSkin[30]+jacobgeo_inv[1]*qSkin[29]+jacobgeo_inv[7]*qSkin[28]+jacobgeo_inv[0]*qSkin[25]+jacobgeo_inv[6]*qSkin[24]+jacobgeo_inv[5]*qSkin[23]+jacobgeo_inv[3]*qSkin[15]); 
  fSkin[26] = 0.3535533905932737*(jacobgeo_inv[0]*qSkin[26]+jacobgeo_inv[1]*qSkin[19]+jacobgeo_inv[2]*qSkin[18]+jacobgeo_inv[3]*qSkin[17]+jacobgeo_inv[4]*qSkin[11]+jacobgeo_inv[5]*qSkin[10]+jacobgeo_inv[6]*qSkin[9]+qSkin[4]*jacobgeo_inv[7]); 
  fSkin[27] = 0.3535533905932737*(jacobgeo_inv[0]*qSkin[27]+jacobgeo_inv[1]*qSkin[22]+jacobgeo_inv[2]*qSkin[21]+jacobgeo_inv[3]*qSkin[20]+jacobgeo_inv[4]*qSkin[14]+jacobgeo_inv[5]*qSkin[13]+jacobgeo_inv[6]*qSkin[12]+qSkin[5]*jacobgeo_inv[7]); 
  fSkin[28] = 0.3535533905932737*(jacobgeo_inv[3]*qSkin[31]+jacobgeo_inv[5]*qSkin[30]+jacobgeo_inv[6]*qSkin[29]+jacobgeo_inv[0]*qSkin[28]+jacobgeo_inv[7]*qSkin[25]+jacobgeo_inv[1]*qSkin[24]+jacobgeo_inv[2]*qSkin[23]+jacobgeo_inv[4]*qSkin[15]); 
  fSkin[29] = 0.3535533905932737*(jacobgeo_inv[2]*qSkin[31]+jacobgeo_inv[4]*qSkin[30]+jacobgeo_inv[0]*qSkin[29]+jacobgeo_inv[6]*qSkin[28]+jacobgeo_inv[1]*qSkin[25]+jacobgeo_inv[7]*qSkin[24]+jacobgeo_inv[3]*qSkin[23]+jacobgeo_inv[5]*qSkin[15]); 
  fSkin[30] = 0.3535533905932737*(jacobgeo_inv[1]*qSkin[31]+jacobgeo_inv[0]*qSkin[30]+jacobgeo_inv[4]*qSkin[29]+jacobgeo_inv[5]*qSkin[28]+jacobgeo_inv[2]*qSkin[25]+jacobgeo_inv[3]*qSkin[24]+jacobgeo_inv[7]*qSkin[23]+jacobgeo_inv[6]*qSkin[15]); 
  fSkin[31] = 0.3535533905932737*(jacobgeo_inv[0]*qSkin[31]+jacobgeo_inv[1]*qSkin[30]+jacobgeo_inv[2]*qSkin[29]+jacobgeo_inv[3]*qSkin[28]+jacobgeo_inv[4]*qSkin[25]+jacobgeo_inv[5]*qSkin[24]+jacobgeo_inv[6]*qSkin[23]+jacobgeo_inv[7]*qSkin[15]); 
  fSkin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*qSkin[43]+15.0*(jacobgeo_inv[6]*qSkin[39]+jacobgeo_inv[5]*qSkin[38]+jacobgeo_inv[4]*qSkin[37])+15.000000000000002*(jacobgeo_inv[3]*qSkin[35]+jacobgeo_inv[2]*qSkin[34]+jacobgeo_inv[1]*qSkin[33])+15.0*jacobgeo_inv[0]*qSkin[32]); 
  fSkin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*qSkin[43]+15.000000000000002*(jacobgeo_inv[7]*qSkin[39]+jacobgeo_inv[3]*qSkin[38]+jacobgeo_inv[2]*qSkin[37])+15.0*(jacobgeo_inv[5]*qSkin[35]+jacobgeo_inv[4]*qSkin[34]+jacobgeo_inv[0]*qSkin[33])+15.000000000000002*jacobgeo_inv[1]*qSkin[32]); 
  fSkin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*qSkin[43]+15.000000000000002*(jacobgeo_inv[3]*qSkin[39]+jacobgeo_inv[7]*qSkin[38]+jacobgeo_inv[1]*qSkin[37])+15.0*(jacobgeo_inv[6]*qSkin[35]+jacobgeo_inv[0]*qSkin[34]+jacobgeo_inv[4]*qSkin[33])+15.000000000000002*jacobgeo_inv[2]*qSkin[32]); 
  fSkin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*qSkin[43]+15.000000000000002*(jacobgeo_inv[2]*qSkin[39]+jacobgeo_inv[1]*qSkin[38]+jacobgeo_inv[7]*qSkin[37])+15.0*(jacobgeo_inv[0]*qSkin[35]+jacobgeo_inv[6]*qSkin[34]+jacobgeo_inv[5]*qSkin[33])+15.000000000000002*jacobgeo_inv[3]*qSkin[32]); 
  fSkin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*qSkin[47]+15.0*(jacobgeo_inv[6]*qSkin[46]+jacobgeo_inv[5]*qSkin[45]+jacobgeo_inv[4]*qSkin[44])+15.000000000000002*(jacobgeo_inv[3]*qSkin[42]+jacobgeo_inv[2]*qSkin[41]+jacobgeo_inv[1]*qSkin[40])+15.0*jacobgeo_inv[0]*qSkin[36]); 
  fSkin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*qSkin[43]+15.0*(jacobgeo_inv[5]*qSkin[39]+jacobgeo_inv[6]*qSkin[38]+jacobgeo_inv[0]*qSkin[37])+15.000000000000002*(jacobgeo_inv[7]*qSkin[35]+jacobgeo_inv[1]*qSkin[34]+jacobgeo_inv[2]*qSkin[33])+15.0*jacobgeo_inv[4]*qSkin[32]); 
  fSkin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*qSkin[43]+15.0*(jacobgeo_inv[4]*qSkin[39]+jacobgeo_inv[0]*qSkin[38]+jacobgeo_inv[6]*qSkin[37])+15.000000000000002*(jacobgeo_inv[1]*qSkin[35]+jacobgeo_inv[7]*qSkin[34]+jacobgeo_inv[3]*qSkin[33])+15.0*jacobgeo_inv[5]*qSkin[32]); 
  fSkin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*qSkin[43]+15.0*(jacobgeo_inv[0]*qSkin[39]+jacobgeo_inv[4]*qSkin[38]+jacobgeo_inv[5]*qSkin[37])+15.000000000000002*(jacobgeo_inv[2]*qSkin[35]+jacobgeo_inv[3]*qSkin[34]+jacobgeo_inv[7]*qSkin[33])+15.0*jacobgeo_inv[6]*qSkin[32]); 
  fSkin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*qSkin[47]+15.000000000000002*(jacobgeo_inv[7]*qSkin[46]+jacobgeo_inv[3]*qSkin[45]+jacobgeo_inv[2]*qSkin[44])+15.0*(jacobgeo_inv[5]*qSkin[42]+jacobgeo_inv[4]*qSkin[41]+jacobgeo_inv[0]*qSkin[40])+15.000000000000002*jacobgeo_inv[1]*qSkin[36]); 
  fSkin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*qSkin[47]+15.000000000000002*(jacobgeo_inv[3]*qSkin[46]+jacobgeo_inv[7]*qSkin[45]+jacobgeo_inv[1]*qSkin[44])+15.0*(jacobgeo_inv[6]*qSkin[42]+jacobgeo_inv[0]*qSkin[41]+jacobgeo_inv[4]*qSkin[40])+15.000000000000002*jacobgeo_inv[2]*qSkin[36]); 
  fSkin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*qSkin[47]+15.000000000000002*(jacobgeo_inv[2]*qSkin[46]+jacobgeo_inv[1]*qSkin[45]+jacobgeo_inv[7]*qSkin[44])+15.0*(jacobgeo_inv[0]*qSkin[42]+jacobgeo_inv[6]*qSkin[41]+jacobgeo_inv[5]*qSkin[40])+15.000000000000002*jacobgeo_inv[3]*qSkin[36]); 
  fSkin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*qSkin[43]+15.000000000000002*(jacobgeo_inv[1]*qSkin[39]+jacobgeo_inv[2]*qSkin[38]+jacobgeo_inv[3]*qSkin[37])+15.0*(jacobgeo_inv[4]*qSkin[35]+jacobgeo_inv[5]*qSkin[34]+jacobgeo_inv[6]*qSkin[33])+15.000000000000002*jacobgeo_inv[7]*qSkin[32]); 
  fSkin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*qSkin[47]+15.0*(jacobgeo_inv[5]*qSkin[46]+jacobgeo_inv[6]*qSkin[45]+jacobgeo_inv[0]*qSkin[44])+15.000000000000002*(jacobgeo_inv[7]*qSkin[42]+jacobgeo_inv[1]*qSkin[41]+jacobgeo_inv[2]*qSkin[40])+15.0*jacobgeo_inv[4]*qSkin[36]); 
  fSkin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*qSkin[47]+15.0*(jacobgeo_inv[4]*qSkin[46]+jacobgeo_inv[0]*qSkin[45]+jacobgeo_inv[6]*qSkin[44])+15.000000000000002*(jacobgeo_inv[1]*qSkin[42]+jacobgeo_inv[7]*qSkin[41]+jacobgeo_inv[3]*qSkin[40])+15.0*jacobgeo_inv[5]*qSkin[36]); 
  fSkin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*qSkin[47]+15.0*(jacobgeo_inv[0]*qSkin[46]+jacobgeo_inv[4]*qSkin[45]+jacobgeo_inv[5]*qSkin[44])+15.000000000000002*(jacobgeo_inv[2]*qSkin[42]+jacobgeo_inv[3]*qSkin[41]+jacobgeo_inv[7]*qSkin[40])+15.0*jacobgeo_inv[6]*qSkin[36]); 
  fSkin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*qSkin[47]+15.000000000000002*(jacobgeo_inv[1]*qSkin[46]+jacobgeo_inv[2]*qSkin[45]+jacobgeo_inv[3]*qSkin[44])+15.0*(jacobgeo_inv[4]*qSkin[42]+jacobgeo_inv[5]*qSkin[41]+jacobgeo_inv[6]*qSkin[40])+15.000000000000002*jacobgeo_inv[7]*qSkin[36]); 

  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 0.6123724356957944*coeff[7]*fSkin[16]-0.3535533905932737*coeff[6]*fSkin[16]-0.5303300858899105*coeff[7]*fSkin[8]+0.3061862178478971*coeff[6]*fSkin[8]+0.6123724356957944*coeff[5]*fSkin[7]-0.3535533905932737*coeff[3]*fSkin[7]+0.6123724356957944*coeff[4]*fSkin[6]-0.3535533905932737*coeff[2]*fSkin[6]-0.5303300858899105*fSkin[3]*coeff[5]-0.5303300858899105*fSkin[2]*coeff[4]+0.3061862178478971*coeff[3]*fSkin[3]+0.3061862178478971*coeff[2]*fSkin[2]+0.6123724356957944*coeff[1]*fSkin[1]-0.3535533905932737*coeff[0]*fSkin[1]-0.5303300858899105*fSkin[0]*coeff[1]+0.3061862178478971*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = 0.6123724356957944*coeff[5]*fSkin[16]-0.3535533905932737*coeff[3]*fSkin[16]-0.5303300858899105*coeff[5]*fSkin[8]+0.3061862178478971*coeff[3]*fSkin[8]+0.6123724356957944*coeff[7]*fSkin[7]-0.3535533905932737*coeff[6]*fSkin[7]-0.5303300858899105*fSkin[3]*coeff[7]+0.6123724356957944*coeff[1]*fSkin[6]-0.3535533905932737*coeff[0]*fSkin[6]+0.3061862178478971*fSkin[3]*coeff[6]+0.6123724356957944*fSkin[1]*coeff[4]-0.5303300858899105*fSkin[0]*coeff[4]-0.5303300858899105*coeff[1]*fSkin[2]+0.3061862178478971*coeff[0]*fSkin[2]-0.3535533905932737*fSkin[1]*coeff[2]+0.3061862178478971*fSkin[0]*coeff[2]; 
  boundSurf_incr[7] = 0.6123724356957944*coeff[4]*fSkin[16]-0.3535533905932737*coeff[2]*fSkin[16]-0.5303300858899105*coeff[4]*fSkin[8]+0.3061862178478971*coeff[2]*fSkin[8]+0.6123724356957944*coeff[1]*fSkin[7]-0.3535533905932737*coeff[0]*fSkin[7]+0.6123724356957944*fSkin[6]*coeff[7]-0.5303300858899105*fSkin[2]*coeff[7]-0.3535533905932737*coeff[6]*fSkin[6]+0.3061862178478971*fSkin[2]*coeff[6]+0.6123724356957944*fSkin[1]*coeff[5]-0.5303300858899105*fSkin[0]*coeff[5]-0.5303300858899105*coeff[1]*fSkin[3]+0.3061862178478971*coeff[0]*fSkin[3]-0.3535533905932737*fSkin[1]*coeff[3]+0.3061862178478971*fSkin[0]*coeff[3]; 
  boundSurf_incr[9] = 0.6123724356957944*coeff[7]*fSkin[26]-0.3535533905932737*coeff[6]*fSkin[26]-0.5303300858899105*coeff[7]*fSkin[19]+0.3061862178478971*coeff[6]*fSkin[19]+0.6123724356957944*coeff[5]*fSkin[18]-0.3535533905932737*coeff[3]*fSkin[18]+0.6123724356957944*coeff[4]*fSkin[17]-0.3535533905932737*coeff[2]*fSkin[17]-0.5303300858899105*coeff[5]*fSkin[11]+0.3061862178478971*coeff[3]*fSkin[11]-0.5303300858899105*coeff[4]*fSkin[10]+0.3061862178478971*coeff[2]*fSkin[10]+0.6123724356957944*coeff[1]*fSkin[9]-0.3535533905932737*coeff[0]*fSkin[9]-0.5303300858899105*coeff[1]*fSkin[4]+0.3061862178478971*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 0.6123724356957944*coeff[7]*fSkin[27]-0.3535533905932737*coeff[6]*fSkin[27]-0.5303300858899105*coeff[7]*fSkin[22]+0.3061862178478971*coeff[6]*fSkin[22]+0.6123724356957944*coeff[5]*fSkin[21]-0.3535533905932737*coeff[3]*fSkin[21]+0.6123724356957944*coeff[4]*fSkin[20]-0.3535533905932737*coeff[2]*fSkin[20]-0.5303300858899105*coeff[5]*fSkin[14]+0.3061862178478971*coeff[3]*fSkin[14]-0.5303300858899105*coeff[4]*fSkin[13]+0.3061862178478971*coeff[2]*fSkin[13]+0.6123724356957944*coeff[1]*fSkin[12]-0.3535533905932737*coeff[0]*fSkin[12]-0.5303300858899105*coeff[1]*fSkin[5]+0.3061862178478971*coeff[0]*fSkin[5]; 
  boundSurf_incr[16] = 0.6123724356957944*coeff[1]*fSkin[16]-0.3535533905932737*coeff[0]*fSkin[16]-0.5303300858899105*coeff[1]*fSkin[8]+0.3061862178478971*coeff[0]*fSkin[8]+0.6123724356957944*coeff[4]*fSkin[7]-0.3535533905932737*coeff[2]*fSkin[7]+0.6123724356957944*fSkin[1]*coeff[7]-0.5303300858899105*fSkin[0]*coeff[7]+0.6123724356957944*coeff[5]*fSkin[6]-0.3535533905932737*coeff[3]*fSkin[6]-0.3535533905932737*fSkin[1]*coeff[6]+0.3061862178478971*fSkin[0]*coeff[6]-0.5303300858899105*fSkin[2]*coeff[5]-0.5303300858899105*fSkin[3]*coeff[4]+0.3061862178478971*coeff[2]*fSkin[3]+0.3061862178478971*fSkin[2]*coeff[3]; 
  boundSurf_incr[17] = 0.6123724356957944*coeff[5]*fSkin[26]-0.3535533905932737*coeff[3]*fSkin[26]-0.5303300858899105*coeff[5]*fSkin[19]+0.3061862178478971*coeff[3]*fSkin[19]+0.6123724356957944*coeff[7]*fSkin[18]-0.3535533905932737*coeff[6]*fSkin[18]+0.6123724356957944*coeff[1]*fSkin[17]-0.3535533905932737*coeff[0]*fSkin[17]-0.5303300858899105*coeff[7]*fSkin[11]+0.3061862178478971*coeff[6]*fSkin[11]-0.5303300858899105*coeff[1]*fSkin[10]+0.3061862178478971*coeff[0]*fSkin[10]+0.6123724356957944*coeff[4]*fSkin[9]-0.3535533905932737*coeff[2]*fSkin[9]-0.5303300858899105*coeff[4]*fSkin[4]+0.3061862178478971*coeff[2]*fSkin[4]; 
  boundSurf_incr[18] = 0.6123724356957944*coeff[4]*fSkin[26]-0.3535533905932737*coeff[2]*fSkin[26]-0.5303300858899105*coeff[4]*fSkin[19]+0.3061862178478971*coeff[2]*fSkin[19]+0.6123724356957944*coeff[1]*fSkin[18]-0.3535533905932737*coeff[0]*fSkin[18]+0.6123724356957944*coeff[7]*fSkin[17]-0.3535533905932737*coeff[6]*fSkin[17]-0.5303300858899105*coeff[1]*fSkin[11]+0.3061862178478971*coeff[0]*fSkin[11]-0.5303300858899105*coeff[7]*fSkin[10]+0.3061862178478971*coeff[6]*fSkin[10]+0.6123724356957944*coeff[5]*fSkin[9]-0.3535533905932737*coeff[3]*fSkin[9]-0.5303300858899105*fSkin[4]*coeff[5]+0.3061862178478971*coeff[3]*fSkin[4]; 
  boundSurf_incr[20] = 0.6123724356957944*coeff[5]*fSkin[27]-0.3535533905932737*coeff[3]*fSkin[27]-0.5303300858899105*coeff[5]*fSkin[22]+0.3061862178478971*coeff[3]*fSkin[22]+0.6123724356957944*coeff[7]*fSkin[21]-0.3535533905932737*coeff[6]*fSkin[21]+0.6123724356957944*coeff[1]*fSkin[20]-0.3535533905932737*coeff[0]*fSkin[20]-0.5303300858899105*coeff[7]*fSkin[14]+0.3061862178478971*coeff[6]*fSkin[14]-0.5303300858899105*coeff[1]*fSkin[13]+0.3061862178478971*coeff[0]*fSkin[13]+0.6123724356957944*coeff[4]*fSkin[12]-0.3535533905932737*coeff[2]*fSkin[12]-0.5303300858899105*coeff[4]*fSkin[5]+0.3061862178478971*coeff[2]*fSkin[5]; 
  boundSurf_incr[21] = 0.6123724356957944*coeff[4]*fSkin[27]-0.3535533905932737*coeff[2]*fSkin[27]-0.5303300858899105*coeff[4]*fSkin[22]+0.3061862178478971*coeff[2]*fSkin[22]+0.6123724356957944*coeff[1]*fSkin[21]-0.3535533905932737*coeff[0]*fSkin[21]+0.6123724356957944*coeff[7]*fSkin[20]-0.3535533905932737*coeff[6]*fSkin[20]-0.5303300858899105*coeff[1]*fSkin[14]+0.3061862178478971*coeff[0]*fSkin[14]-0.5303300858899105*coeff[7]*fSkin[13]+0.3061862178478971*coeff[6]*fSkin[13]+0.6123724356957944*coeff[5]*fSkin[12]-0.3535533905932737*coeff[3]*fSkin[12]-0.5303300858899105*coeff[5]*fSkin[5]+0.3061862178478971*coeff[3]*fSkin[5]; 
  boundSurf_incr[23] = 0.6123724356957944*coeff[7]*fSkin[31]-0.3535533905932737*coeff[6]*fSkin[31]-0.5303300858899105*coeff[7]*fSkin[30]+0.3061862178478971*coeff[6]*fSkin[30]+0.6123724356957944*coeff[5]*fSkin[29]-0.3535533905932737*coeff[3]*fSkin[29]+0.6123724356957944*coeff[4]*fSkin[28]-0.3535533905932737*coeff[2]*fSkin[28]-0.5303300858899105*coeff[5]*fSkin[25]+0.3061862178478971*coeff[3]*fSkin[25]-0.5303300858899105*coeff[4]*fSkin[24]+0.3061862178478971*coeff[2]*fSkin[24]+0.6123724356957944*coeff[1]*fSkin[23]-0.3535533905932737*coeff[0]*fSkin[23]-0.5303300858899105*coeff[1]*fSkin[15]+0.3061862178478971*coeff[0]*fSkin[15]; 
  boundSurf_incr[26] = 0.6123724356957944*coeff[1]*fSkin[26]-0.3535533905932737*coeff[0]*fSkin[26]-0.5303300858899105*coeff[1]*fSkin[19]+0.3061862178478971*coeff[0]*fSkin[19]+0.6123724356957944*coeff[4]*fSkin[18]-0.3535533905932737*coeff[2]*fSkin[18]+0.6123724356957944*coeff[5]*fSkin[17]-0.3535533905932737*coeff[3]*fSkin[17]-0.5303300858899105*coeff[4]*fSkin[11]+0.3061862178478971*coeff[2]*fSkin[11]-0.5303300858899105*coeff[5]*fSkin[10]+0.3061862178478971*coeff[3]*fSkin[10]+0.6123724356957944*coeff[7]*fSkin[9]-0.3535533905932737*coeff[6]*fSkin[9]-0.5303300858899105*fSkin[4]*coeff[7]+0.3061862178478971*fSkin[4]*coeff[6]; 
  boundSurf_incr[27] = 0.6123724356957944*coeff[1]*fSkin[27]-0.3535533905932737*coeff[0]*fSkin[27]-0.5303300858899105*coeff[1]*fSkin[22]+0.3061862178478971*coeff[0]*fSkin[22]+0.6123724356957944*coeff[4]*fSkin[21]-0.3535533905932737*coeff[2]*fSkin[21]+0.6123724356957944*coeff[5]*fSkin[20]-0.3535533905932737*coeff[3]*fSkin[20]-0.5303300858899105*coeff[4]*fSkin[14]+0.3061862178478971*coeff[2]*fSkin[14]-0.5303300858899105*coeff[5]*fSkin[13]+0.3061862178478971*coeff[3]*fSkin[13]+0.6123724356957944*coeff[7]*fSkin[12]-0.3535533905932737*coeff[6]*fSkin[12]-0.5303300858899105*fSkin[5]*coeff[7]+0.3061862178478971*fSkin[5]*coeff[6]; 
  boundSurf_incr[28] = 0.6123724356957944*coeff[5]*fSkin[31]-0.3535533905932737*coeff[3]*fSkin[31]-0.5303300858899105*coeff[5]*fSkin[30]+0.3061862178478971*coeff[3]*fSkin[30]+0.6123724356957944*coeff[7]*fSkin[29]-0.3535533905932737*coeff[6]*fSkin[29]+0.6123724356957944*coeff[1]*fSkin[28]-0.3535533905932737*coeff[0]*fSkin[28]-0.5303300858899105*coeff[7]*fSkin[25]+0.3061862178478971*coeff[6]*fSkin[25]-0.5303300858899105*coeff[1]*fSkin[24]+0.3061862178478971*coeff[0]*fSkin[24]+0.6123724356957944*coeff[4]*fSkin[23]-0.3535533905932737*coeff[2]*fSkin[23]-0.5303300858899105*coeff[4]*fSkin[15]+0.3061862178478971*coeff[2]*fSkin[15]; 
  boundSurf_incr[29] = 0.6123724356957944*coeff[4]*fSkin[31]-0.3535533905932737*coeff[2]*fSkin[31]-0.5303300858899105*coeff[4]*fSkin[30]+0.3061862178478971*coeff[2]*fSkin[30]+0.6123724356957944*coeff[1]*fSkin[29]-0.3535533905932737*coeff[0]*fSkin[29]+0.6123724356957944*coeff[7]*fSkin[28]-0.3535533905932737*coeff[6]*fSkin[28]-0.5303300858899105*coeff[1]*fSkin[25]+0.3061862178478971*coeff[0]*fSkin[25]-0.5303300858899105*coeff[7]*fSkin[24]+0.3061862178478971*coeff[6]*fSkin[24]+0.6123724356957944*coeff[5]*fSkin[23]-0.3535533905932737*coeff[3]*fSkin[23]-0.5303300858899105*coeff[5]*fSkin[15]+0.3061862178478971*coeff[3]*fSkin[15]; 
  boundSurf_incr[31] = 0.6123724356957944*coeff[1]*fSkin[31]-0.3535533905932737*coeff[0]*fSkin[31]-0.5303300858899105*coeff[1]*fSkin[30]+0.3061862178478971*coeff[0]*fSkin[30]+0.6123724356957944*coeff[4]*fSkin[29]-0.3535533905932737*coeff[2]*fSkin[29]+0.6123724356957944*coeff[5]*fSkin[28]-0.3535533905932737*coeff[3]*fSkin[28]-0.5303300858899105*coeff[4]*fSkin[25]+0.3061862178478971*coeff[2]*fSkin[25]-0.5303300858899105*coeff[5]*fSkin[24]+0.3061862178478971*coeff[3]*fSkin[24]+0.6123724356957944*coeff[7]*fSkin[23]-0.3535533905932737*coeff[6]*fSkin[23]-0.5303300858899105*coeff[7]*fSkin[15]+0.3061862178478971*coeff[6]*fSkin[15]; 
  boundSurf_incr[33] = 0.6123724356957944*coeff[7]*fSkin[43]-0.3535533905932737*coeff[6]*fSkin[43]-0.5303300858899104*coeff[7]*fSkin[39]+0.3061862178478971*coeff[6]*fSkin[39]+0.6123724356957944*coeff[5]*fSkin[38]-0.3535533905932737*coeff[3]*fSkin[38]+0.6123724356957944*coeff[4]*fSkin[37]-0.3535533905932737*coeff[2]*fSkin[37]-0.5303300858899105*coeff[5]*fSkin[35]+0.3061862178478971*coeff[3]*fSkin[35]-0.5303300858899105*coeff[4]*fSkin[34]+0.3061862178478971*coeff[2]*fSkin[34]+0.6123724356957944*coeff[1]*fSkin[33]-0.3535533905932737*coeff[0]*fSkin[33]-0.5303300858899104*coeff[1]*fSkin[32]+0.3061862178478971*coeff[0]*fSkin[32]; 
  boundSurf_incr[37] = 0.6123724356957944*coeff[5]*fSkin[43]-0.3535533905932737*coeff[3]*fSkin[43]-0.5303300858899105*coeff[5]*fSkin[39]+0.3061862178478971*coeff[3]*fSkin[39]+0.6123724356957944*coeff[7]*fSkin[38]-0.3535533905932737*coeff[6]*fSkin[38]+0.6123724356957944*coeff[1]*fSkin[37]-0.3535533905932737*coeff[0]*fSkin[37]-0.5303300858899104*coeff[7]*fSkin[35]+0.3061862178478971*coeff[6]*fSkin[35]-0.5303300858899104*coeff[1]*fSkin[34]+0.3061862178478971*coeff[0]*fSkin[34]+0.6123724356957944*coeff[4]*fSkin[33]-0.3535533905932737*coeff[2]*fSkin[33]-0.5303300858899105*coeff[4]*fSkin[32]+0.3061862178478971*coeff[2]*fSkin[32]; 
  boundSurf_incr[38] = 0.6123724356957944*coeff[4]*fSkin[43]-0.3535533905932737*coeff[2]*fSkin[43]-0.5303300858899105*coeff[4]*fSkin[39]+0.3061862178478971*coeff[2]*fSkin[39]+0.6123724356957944*coeff[1]*fSkin[38]-0.3535533905932737*coeff[0]*fSkin[38]+0.6123724356957944*coeff[7]*fSkin[37]-0.3535533905932737*coeff[6]*fSkin[37]-0.5303300858899104*coeff[1]*fSkin[35]+0.3061862178478971*coeff[0]*fSkin[35]-0.5303300858899104*coeff[7]*fSkin[34]+0.3061862178478971*coeff[6]*fSkin[34]+0.6123724356957944*coeff[5]*fSkin[33]-0.3535533905932737*coeff[3]*fSkin[33]-0.5303300858899105*coeff[5]*fSkin[32]+0.3061862178478971*coeff[3]*fSkin[32]; 
  boundSurf_incr[40] = 0.6123724356957944*coeff[7]*fSkin[47]-0.3535533905932737*coeff[6]*fSkin[47]-0.5303300858899104*coeff[7]*fSkin[46]+0.3061862178478971*coeff[6]*fSkin[46]+0.6123724356957944*coeff[5]*fSkin[45]-0.3535533905932737*coeff[3]*fSkin[45]+0.6123724356957944*coeff[4]*fSkin[44]-0.3535533905932737*coeff[2]*fSkin[44]-0.5303300858899105*coeff[5]*fSkin[42]+0.3061862178478971*coeff[3]*fSkin[42]-0.5303300858899105*coeff[4]*fSkin[41]+0.3061862178478971*coeff[2]*fSkin[41]+0.6123724356957944*coeff[1]*fSkin[40]-0.3535533905932737*coeff[0]*fSkin[40]-0.5303300858899104*coeff[1]*fSkin[36]+0.3061862178478971*coeff[0]*fSkin[36]; 
  boundSurf_incr[43] = 0.6123724356957944*coeff[1]*fSkin[43]-0.3535533905932737*coeff[0]*fSkin[43]-0.5303300858899104*coeff[1]*fSkin[39]+0.3061862178478971*coeff[0]*fSkin[39]+0.6123724356957944*coeff[4]*fSkin[38]-0.3535533905932737*coeff[2]*fSkin[38]+0.6123724356957944*coeff[5]*fSkin[37]-0.3535533905932737*coeff[3]*fSkin[37]-0.5303300858899105*coeff[4]*fSkin[35]+0.3061862178478971*coeff[2]*fSkin[35]-0.5303300858899105*coeff[5]*fSkin[34]+0.3061862178478971*coeff[3]*fSkin[34]+0.6123724356957944*coeff[7]*fSkin[33]-0.3535533905932737*coeff[6]*fSkin[33]-0.5303300858899104*coeff[7]*fSkin[32]+0.3061862178478971*coeff[6]*fSkin[32]; 
  boundSurf_incr[44] = 0.6123724356957944*coeff[5]*fSkin[47]-0.3535533905932737*coeff[3]*fSkin[47]-0.5303300858899105*coeff[5]*fSkin[46]+0.3061862178478971*coeff[3]*fSkin[46]+0.6123724356957944*coeff[7]*fSkin[45]-0.3535533905932737*coeff[6]*fSkin[45]+0.6123724356957944*coeff[1]*fSkin[44]-0.3535533905932737*coeff[0]*fSkin[44]-0.5303300858899104*coeff[7]*fSkin[42]+0.3061862178478971*coeff[6]*fSkin[42]-0.5303300858899104*coeff[1]*fSkin[41]+0.3061862178478971*coeff[0]*fSkin[41]+0.6123724356957944*coeff[4]*fSkin[40]-0.3535533905932737*coeff[2]*fSkin[40]-0.5303300858899105*coeff[4]*fSkin[36]+0.3061862178478971*coeff[2]*fSkin[36]; 
  boundSurf_incr[45] = 0.6123724356957944*coeff[4]*fSkin[47]-0.3535533905932737*coeff[2]*fSkin[47]-0.5303300858899105*coeff[4]*fSkin[46]+0.3061862178478971*coeff[2]*fSkin[46]+0.6123724356957944*coeff[1]*fSkin[45]-0.3535533905932737*coeff[0]*fSkin[45]+0.6123724356957944*coeff[7]*fSkin[44]-0.3535533905932737*coeff[6]*fSkin[44]-0.5303300858899104*coeff[1]*fSkin[42]+0.3061862178478971*coeff[0]*fSkin[42]-0.5303300858899104*coeff[7]*fSkin[41]+0.3061862178478971*coeff[6]*fSkin[41]+0.6123724356957944*coeff[5]*fSkin[40]-0.3535533905932737*coeff[3]*fSkin[40]-0.5303300858899105*coeff[5]*fSkin[36]+0.3061862178478971*coeff[3]*fSkin[36]; 
  boundSurf_incr[47] = 0.6123724356957944*coeff[1]*fSkin[47]-0.3535533905932737*coeff[0]*fSkin[47]-0.5303300858899104*coeff[1]*fSkin[46]+0.3061862178478971*coeff[0]*fSkin[46]+0.6123724356957944*coeff[4]*fSkin[45]-0.3535533905932737*coeff[2]*fSkin[45]+0.6123724356957944*coeff[5]*fSkin[44]-0.3535533905932737*coeff[3]*fSkin[44]-0.5303300858899105*coeff[4]*fSkin[42]+0.3061862178478971*coeff[2]*fSkin[42]-0.5303300858899105*coeff[5]*fSkin[41]+0.3061862178478971*coeff[3]*fSkin[41]+0.6123724356957944*coeff[7]*fSkin[40]-0.3535533905932737*coeff[6]*fSkin[40]-0.5303300858899104*coeff[7]*fSkin[36]+0.3061862178478971*coeff[6]*fSkin[36]; 

  } else { 

  boundSurf_incr[1] = -(0.6123724356957944*coeff[7]*fSkin[16])-0.3535533905932737*coeff[6]*fSkin[16]-0.5303300858899105*coeff[7]*fSkin[8]-0.3061862178478971*coeff[6]*fSkin[8]-0.6123724356957944*coeff[5]*fSkin[7]-0.3535533905932737*coeff[3]*fSkin[7]-0.6123724356957944*coeff[4]*fSkin[6]-0.3535533905932737*coeff[2]*fSkin[6]-0.5303300858899105*fSkin[3]*coeff[5]-0.5303300858899105*fSkin[2]*coeff[4]-0.3061862178478971*coeff[3]*fSkin[3]-0.3061862178478971*coeff[2]*fSkin[2]-0.6123724356957944*coeff[1]*fSkin[1]-0.3535533905932737*coeff[0]*fSkin[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.3061862178478971*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = -(0.6123724356957944*coeff[5]*fSkin[16])-0.3535533905932737*coeff[3]*fSkin[16]-0.5303300858899105*coeff[5]*fSkin[8]-0.3061862178478971*coeff[3]*fSkin[8]-0.6123724356957944*coeff[7]*fSkin[7]-0.3535533905932737*coeff[6]*fSkin[7]-0.5303300858899105*fSkin[3]*coeff[7]-0.6123724356957944*coeff[1]*fSkin[6]-0.3535533905932737*coeff[0]*fSkin[6]-0.3061862178478971*fSkin[3]*coeff[6]-0.6123724356957944*fSkin[1]*coeff[4]-0.5303300858899105*fSkin[0]*coeff[4]-0.5303300858899105*coeff[1]*fSkin[2]-0.3061862178478971*coeff[0]*fSkin[2]-0.3535533905932737*fSkin[1]*coeff[2]-0.3061862178478971*fSkin[0]*coeff[2]; 
  boundSurf_incr[7] = -(0.6123724356957944*coeff[4]*fSkin[16])-0.3535533905932737*coeff[2]*fSkin[16]-0.5303300858899105*coeff[4]*fSkin[8]-0.3061862178478971*coeff[2]*fSkin[8]-0.6123724356957944*coeff[1]*fSkin[7]-0.3535533905932737*coeff[0]*fSkin[7]-0.6123724356957944*fSkin[6]*coeff[7]-0.5303300858899105*fSkin[2]*coeff[7]-0.3535533905932737*coeff[6]*fSkin[6]-0.3061862178478971*fSkin[2]*coeff[6]-0.6123724356957944*fSkin[1]*coeff[5]-0.5303300858899105*fSkin[0]*coeff[5]-0.5303300858899105*coeff[1]*fSkin[3]-0.3061862178478971*coeff[0]*fSkin[3]-0.3535533905932737*fSkin[1]*coeff[3]-0.3061862178478971*fSkin[0]*coeff[3]; 
  boundSurf_incr[9] = -(0.6123724356957944*coeff[7]*fSkin[26])-0.3535533905932737*coeff[6]*fSkin[26]-0.5303300858899105*coeff[7]*fSkin[19]-0.3061862178478971*coeff[6]*fSkin[19]-0.6123724356957944*coeff[5]*fSkin[18]-0.3535533905932737*coeff[3]*fSkin[18]-0.6123724356957944*coeff[4]*fSkin[17]-0.3535533905932737*coeff[2]*fSkin[17]-0.5303300858899105*coeff[5]*fSkin[11]-0.3061862178478971*coeff[3]*fSkin[11]-0.5303300858899105*coeff[4]*fSkin[10]-0.3061862178478971*coeff[2]*fSkin[10]-0.6123724356957944*coeff[1]*fSkin[9]-0.3535533905932737*coeff[0]*fSkin[9]-0.5303300858899105*coeff[1]*fSkin[4]-0.3061862178478971*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = -(0.6123724356957944*coeff[7]*fSkin[27])-0.3535533905932737*coeff[6]*fSkin[27]-0.5303300858899105*coeff[7]*fSkin[22]-0.3061862178478971*coeff[6]*fSkin[22]-0.6123724356957944*coeff[5]*fSkin[21]-0.3535533905932737*coeff[3]*fSkin[21]-0.6123724356957944*coeff[4]*fSkin[20]-0.3535533905932737*coeff[2]*fSkin[20]-0.5303300858899105*coeff[5]*fSkin[14]-0.3061862178478971*coeff[3]*fSkin[14]-0.5303300858899105*coeff[4]*fSkin[13]-0.3061862178478971*coeff[2]*fSkin[13]-0.6123724356957944*coeff[1]*fSkin[12]-0.3535533905932737*coeff[0]*fSkin[12]-0.5303300858899105*coeff[1]*fSkin[5]-0.3061862178478971*coeff[0]*fSkin[5]; 
  boundSurf_incr[16] = -(0.6123724356957944*coeff[1]*fSkin[16])-0.3535533905932737*coeff[0]*fSkin[16]-0.5303300858899105*coeff[1]*fSkin[8]-0.3061862178478971*coeff[0]*fSkin[8]-0.6123724356957944*coeff[4]*fSkin[7]-0.3535533905932737*coeff[2]*fSkin[7]-0.6123724356957944*fSkin[1]*coeff[7]-0.5303300858899105*fSkin[0]*coeff[7]-0.6123724356957944*coeff[5]*fSkin[6]-0.3535533905932737*coeff[3]*fSkin[6]-0.3535533905932737*fSkin[1]*coeff[6]-0.3061862178478971*fSkin[0]*coeff[6]-0.5303300858899105*fSkin[2]*coeff[5]-0.5303300858899105*fSkin[3]*coeff[4]-0.3061862178478971*coeff[2]*fSkin[3]-0.3061862178478971*fSkin[2]*coeff[3]; 
  boundSurf_incr[17] = -(0.6123724356957944*coeff[5]*fSkin[26])-0.3535533905932737*coeff[3]*fSkin[26]-0.5303300858899105*coeff[5]*fSkin[19]-0.3061862178478971*coeff[3]*fSkin[19]-0.6123724356957944*coeff[7]*fSkin[18]-0.3535533905932737*coeff[6]*fSkin[18]-0.6123724356957944*coeff[1]*fSkin[17]-0.3535533905932737*coeff[0]*fSkin[17]-0.5303300858899105*coeff[7]*fSkin[11]-0.3061862178478971*coeff[6]*fSkin[11]-0.5303300858899105*coeff[1]*fSkin[10]-0.3061862178478971*coeff[0]*fSkin[10]-0.6123724356957944*coeff[4]*fSkin[9]-0.3535533905932737*coeff[2]*fSkin[9]-0.5303300858899105*coeff[4]*fSkin[4]-0.3061862178478971*coeff[2]*fSkin[4]; 
  boundSurf_incr[18] = -(0.6123724356957944*coeff[4]*fSkin[26])-0.3535533905932737*coeff[2]*fSkin[26]-0.5303300858899105*coeff[4]*fSkin[19]-0.3061862178478971*coeff[2]*fSkin[19]-0.6123724356957944*coeff[1]*fSkin[18]-0.3535533905932737*coeff[0]*fSkin[18]-0.6123724356957944*coeff[7]*fSkin[17]-0.3535533905932737*coeff[6]*fSkin[17]-0.5303300858899105*coeff[1]*fSkin[11]-0.3061862178478971*coeff[0]*fSkin[11]-0.5303300858899105*coeff[7]*fSkin[10]-0.3061862178478971*coeff[6]*fSkin[10]-0.6123724356957944*coeff[5]*fSkin[9]-0.3535533905932737*coeff[3]*fSkin[9]-0.5303300858899105*fSkin[4]*coeff[5]-0.3061862178478971*coeff[3]*fSkin[4]; 
  boundSurf_incr[20] = -(0.6123724356957944*coeff[5]*fSkin[27])-0.3535533905932737*coeff[3]*fSkin[27]-0.5303300858899105*coeff[5]*fSkin[22]-0.3061862178478971*coeff[3]*fSkin[22]-0.6123724356957944*coeff[7]*fSkin[21]-0.3535533905932737*coeff[6]*fSkin[21]-0.6123724356957944*coeff[1]*fSkin[20]-0.3535533905932737*coeff[0]*fSkin[20]-0.5303300858899105*coeff[7]*fSkin[14]-0.3061862178478971*coeff[6]*fSkin[14]-0.5303300858899105*coeff[1]*fSkin[13]-0.3061862178478971*coeff[0]*fSkin[13]-0.6123724356957944*coeff[4]*fSkin[12]-0.3535533905932737*coeff[2]*fSkin[12]-0.5303300858899105*coeff[4]*fSkin[5]-0.3061862178478971*coeff[2]*fSkin[5]; 
  boundSurf_incr[21] = -(0.6123724356957944*coeff[4]*fSkin[27])-0.3535533905932737*coeff[2]*fSkin[27]-0.5303300858899105*coeff[4]*fSkin[22]-0.3061862178478971*coeff[2]*fSkin[22]-0.6123724356957944*coeff[1]*fSkin[21]-0.3535533905932737*coeff[0]*fSkin[21]-0.6123724356957944*coeff[7]*fSkin[20]-0.3535533905932737*coeff[6]*fSkin[20]-0.5303300858899105*coeff[1]*fSkin[14]-0.3061862178478971*coeff[0]*fSkin[14]-0.5303300858899105*coeff[7]*fSkin[13]-0.3061862178478971*coeff[6]*fSkin[13]-0.6123724356957944*coeff[5]*fSkin[12]-0.3535533905932737*coeff[3]*fSkin[12]-0.5303300858899105*coeff[5]*fSkin[5]-0.3061862178478971*coeff[3]*fSkin[5]; 
  boundSurf_incr[23] = -(0.6123724356957944*coeff[7]*fSkin[31])-0.3535533905932737*coeff[6]*fSkin[31]-0.5303300858899105*coeff[7]*fSkin[30]-0.3061862178478971*coeff[6]*fSkin[30]-0.6123724356957944*coeff[5]*fSkin[29]-0.3535533905932737*coeff[3]*fSkin[29]-0.6123724356957944*coeff[4]*fSkin[28]-0.3535533905932737*coeff[2]*fSkin[28]-0.5303300858899105*coeff[5]*fSkin[25]-0.3061862178478971*coeff[3]*fSkin[25]-0.5303300858899105*coeff[4]*fSkin[24]-0.3061862178478971*coeff[2]*fSkin[24]-0.6123724356957944*coeff[1]*fSkin[23]-0.3535533905932737*coeff[0]*fSkin[23]-0.5303300858899105*coeff[1]*fSkin[15]-0.3061862178478971*coeff[0]*fSkin[15]; 
  boundSurf_incr[26] = -(0.6123724356957944*coeff[1]*fSkin[26])-0.3535533905932737*coeff[0]*fSkin[26]-0.5303300858899105*coeff[1]*fSkin[19]-0.3061862178478971*coeff[0]*fSkin[19]-0.6123724356957944*coeff[4]*fSkin[18]-0.3535533905932737*coeff[2]*fSkin[18]-0.6123724356957944*coeff[5]*fSkin[17]-0.3535533905932737*coeff[3]*fSkin[17]-0.5303300858899105*coeff[4]*fSkin[11]-0.3061862178478971*coeff[2]*fSkin[11]-0.5303300858899105*coeff[5]*fSkin[10]-0.3061862178478971*coeff[3]*fSkin[10]-0.6123724356957944*coeff[7]*fSkin[9]-0.3535533905932737*coeff[6]*fSkin[9]-0.5303300858899105*fSkin[4]*coeff[7]-0.3061862178478971*fSkin[4]*coeff[6]; 
  boundSurf_incr[27] = -(0.6123724356957944*coeff[1]*fSkin[27])-0.3535533905932737*coeff[0]*fSkin[27]-0.5303300858899105*coeff[1]*fSkin[22]-0.3061862178478971*coeff[0]*fSkin[22]-0.6123724356957944*coeff[4]*fSkin[21]-0.3535533905932737*coeff[2]*fSkin[21]-0.6123724356957944*coeff[5]*fSkin[20]-0.3535533905932737*coeff[3]*fSkin[20]-0.5303300858899105*coeff[4]*fSkin[14]-0.3061862178478971*coeff[2]*fSkin[14]-0.5303300858899105*coeff[5]*fSkin[13]-0.3061862178478971*coeff[3]*fSkin[13]-0.6123724356957944*coeff[7]*fSkin[12]-0.3535533905932737*coeff[6]*fSkin[12]-0.5303300858899105*fSkin[5]*coeff[7]-0.3061862178478971*fSkin[5]*coeff[6]; 
  boundSurf_incr[28] = -(0.6123724356957944*coeff[5]*fSkin[31])-0.3535533905932737*coeff[3]*fSkin[31]-0.5303300858899105*coeff[5]*fSkin[30]-0.3061862178478971*coeff[3]*fSkin[30]-0.6123724356957944*coeff[7]*fSkin[29]-0.3535533905932737*coeff[6]*fSkin[29]-0.6123724356957944*coeff[1]*fSkin[28]-0.3535533905932737*coeff[0]*fSkin[28]-0.5303300858899105*coeff[7]*fSkin[25]-0.3061862178478971*coeff[6]*fSkin[25]-0.5303300858899105*coeff[1]*fSkin[24]-0.3061862178478971*coeff[0]*fSkin[24]-0.6123724356957944*coeff[4]*fSkin[23]-0.3535533905932737*coeff[2]*fSkin[23]-0.5303300858899105*coeff[4]*fSkin[15]-0.3061862178478971*coeff[2]*fSkin[15]; 
  boundSurf_incr[29] = -(0.6123724356957944*coeff[4]*fSkin[31])-0.3535533905932737*coeff[2]*fSkin[31]-0.5303300858899105*coeff[4]*fSkin[30]-0.3061862178478971*coeff[2]*fSkin[30]-0.6123724356957944*coeff[1]*fSkin[29]-0.3535533905932737*coeff[0]*fSkin[29]-0.6123724356957944*coeff[7]*fSkin[28]-0.3535533905932737*coeff[6]*fSkin[28]-0.5303300858899105*coeff[1]*fSkin[25]-0.3061862178478971*coeff[0]*fSkin[25]-0.5303300858899105*coeff[7]*fSkin[24]-0.3061862178478971*coeff[6]*fSkin[24]-0.6123724356957944*coeff[5]*fSkin[23]-0.3535533905932737*coeff[3]*fSkin[23]-0.5303300858899105*coeff[5]*fSkin[15]-0.3061862178478971*coeff[3]*fSkin[15]; 
  boundSurf_incr[31] = -(0.6123724356957944*coeff[1]*fSkin[31])-0.3535533905932737*coeff[0]*fSkin[31]-0.5303300858899105*coeff[1]*fSkin[30]-0.3061862178478971*coeff[0]*fSkin[30]-0.6123724356957944*coeff[4]*fSkin[29]-0.3535533905932737*coeff[2]*fSkin[29]-0.6123724356957944*coeff[5]*fSkin[28]-0.3535533905932737*coeff[3]*fSkin[28]-0.5303300858899105*coeff[4]*fSkin[25]-0.3061862178478971*coeff[2]*fSkin[25]-0.5303300858899105*coeff[5]*fSkin[24]-0.3061862178478971*coeff[3]*fSkin[24]-0.6123724356957944*coeff[7]*fSkin[23]-0.3535533905932737*coeff[6]*fSkin[23]-0.5303300858899105*coeff[7]*fSkin[15]-0.3061862178478971*coeff[6]*fSkin[15]; 
  boundSurf_incr[33] = -(0.6123724356957944*coeff[7]*fSkin[43])-0.3535533905932737*coeff[6]*fSkin[43]-0.5303300858899104*coeff[7]*fSkin[39]-0.3061862178478971*coeff[6]*fSkin[39]-0.6123724356957944*coeff[5]*fSkin[38]-0.3535533905932737*coeff[3]*fSkin[38]-0.6123724356957944*coeff[4]*fSkin[37]-0.3535533905932737*coeff[2]*fSkin[37]-0.5303300858899105*coeff[5]*fSkin[35]-0.3061862178478971*coeff[3]*fSkin[35]-0.5303300858899105*coeff[4]*fSkin[34]-0.3061862178478971*coeff[2]*fSkin[34]-0.6123724356957944*coeff[1]*fSkin[33]-0.3535533905932737*coeff[0]*fSkin[33]-0.5303300858899104*coeff[1]*fSkin[32]-0.3061862178478971*coeff[0]*fSkin[32]; 
  boundSurf_incr[37] = -(0.6123724356957944*coeff[5]*fSkin[43])-0.3535533905932737*coeff[3]*fSkin[43]-0.5303300858899105*coeff[5]*fSkin[39]-0.3061862178478971*coeff[3]*fSkin[39]-0.6123724356957944*coeff[7]*fSkin[38]-0.3535533905932737*coeff[6]*fSkin[38]-0.6123724356957944*coeff[1]*fSkin[37]-0.3535533905932737*coeff[0]*fSkin[37]-0.5303300858899104*coeff[7]*fSkin[35]-0.3061862178478971*coeff[6]*fSkin[35]-0.5303300858899104*coeff[1]*fSkin[34]-0.3061862178478971*coeff[0]*fSkin[34]-0.6123724356957944*coeff[4]*fSkin[33]-0.3535533905932737*coeff[2]*fSkin[33]-0.5303300858899105*coeff[4]*fSkin[32]-0.3061862178478971*coeff[2]*fSkin[32]; 
  boundSurf_incr[38] = -(0.6123724356957944*coeff[4]*fSkin[43])-0.3535533905932737*coeff[2]*fSkin[43]-0.5303300858899105*coeff[4]*fSkin[39]-0.3061862178478971*coeff[2]*fSkin[39]-0.6123724356957944*coeff[1]*fSkin[38]-0.3535533905932737*coeff[0]*fSkin[38]-0.6123724356957944*coeff[7]*fSkin[37]-0.3535533905932737*coeff[6]*fSkin[37]-0.5303300858899104*coeff[1]*fSkin[35]-0.3061862178478971*coeff[0]*fSkin[35]-0.5303300858899104*coeff[7]*fSkin[34]-0.3061862178478971*coeff[6]*fSkin[34]-0.6123724356957944*coeff[5]*fSkin[33]-0.3535533905932737*coeff[3]*fSkin[33]-0.5303300858899105*coeff[5]*fSkin[32]-0.3061862178478971*coeff[3]*fSkin[32]; 
  boundSurf_incr[40] = -(0.6123724356957944*coeff[7]*fSkin[47])-0.3535533905932737*coeff[6]*fSkin[47]-0.5303300858899104*coeff[7]*fSkin[46]-0.3061862178478971*coeff[6]*fSkin[46]-0.6123724356957944*coeff[5]*fSkin[45]-0.3535533905932737*coeff[3]*fSkin[45]-0.6123724356957944*coeff[4]*fSkin[44]-0.3535533905932737*coeff[2]*fSkin[44]-0.5303300858899105*coeff[5]*fSkin[42]-0.3061862178478971*coeff[3]*fSkin[42]-0.5303300858899105*coeff[4]*fSkin[41]-0.3061862178478971*coeff[2]*fSkin[41]-0.6123724356957944*coeff[1]*fSkin[40]-0.3535533905932737*coeff[0]*fSkin[40]-0.5303300858899104*coeff[1]*fSkin[36]-0.3061862178478971*coeff[0]*fSkin[36]; 
  boundSurf_incr[43] = -(0.6123724356957944*coeff[1]*fSkin[43])-0.3535533905932737*coeff[0]*fSkin[43]-0.5303300858899104*coeff[1]*fSkin[39]-0.3061862178478971*coeff[0]*fSkin[39]-0.6123724356957944*coeff[4]*fSkin[38]-0.3535533905932737*coeff[2]*fSkin[38]-0.6123724356957944*coeff[5]*fSkin[37]-0.3535533905932737*coeff[3]*fSkin[37]-0.5303300858899105*coeff[4]*fSkin[35]-0.3061862178478971*coeff[2]*fSkin[35]-0.5303300858899105*coeff[5]*fSkin[34]-0.3061862178478971*coeff[3]*fSkin[34]-0.6123724356957944*coeff[7]*fSkin[33]-0.3535533905932737*coeff[6]*fSkin[33]-0.5303300858899104*coeff[7]*fSkin[32]-0.3061862178478971*coeff[6]*fSkin[32]; 
  boundSurf_incr[44] = -(0.6123724356957944*coeff[5]*fSkin[47])-0.3535533905932737*coeff[3]*fSkin[47]-0.5303300858899105*coeff[5]*fSkin[46]-0.3061862178478971*coeff[3]*fSkin[46]-0.6123724356957944*coeff[7]*fSkin[45]-0.3535533905932737*coeff[6]*fSkin[45]-0.6123724356957944*coeff[1]*fSkin[44]-0.3535533905932737*coeff[0]*fSkin[44]-0.5303300858899104*coeff[7]*fSkin[42]-0.3061862178478971*coeff[6]*fSkin[42]-0.5303300858899104*coeff[1]*fSkin[41]-0.3061862178478971*coeff[0]*fSkin[41]-0.6123724356957944*coeff[4]*fSkin[40]-0.3535533905932737*coeff[2]*fSkin[40]-0.5303300858899105*coeff[4]*fSkin[36]-0.3061862178478971*coeff[2]*fSkin[36]; 
  boundSurf_incr[45] = -(0.6123724356957944*coeff[4]*fSkin[47])-0.3535533905932737*coeff[2]*fSkin[47]-0.5303300858899105*coeff[4]*fSkin[46]-0.3061862178478971*coeff[2]*fSkin[46]-0.6123724356957944*coeff[1]*fSkin[45]-0.3535533905932737*coeff[0]*fSkin[45]-0.6123724356957944*coeff[7]*fSkin[44]-0.3535533905932737*coeff[6]*fSkin[44]-0.5303300858899104*coeff[1]*fSkin[42]-0.3061862178478971*coeff[0]*fSkin[42]-0.5303300858899104*coeff[7]*fSkin[41]-0.3061862178478971*coeff[6]*fSkin[41]-0.6123724356957944*coeff[5]*fSkin[40]-0.3535533905932737*coeff[3]*fSkin[40]-0.5303300858899105*coeff[5]*fSkin[36]-0.3061862178478971*coeff[3]*fSkin[36]; 
  boundSurf_incr[47] = -(0.6123724356957944*coeff[1]*fSkin[47])-0.3535533905932737*coeff[0]*fSkin[47]-0.5303300858899104*coeff[1]*fSkin[46]-0.3061862178478971*coeff[0]*fSkin[46]-0.6123724356957944*coeff[4]*fSkin[45]-0.3535533905932737*coeff[2]*fSkin[45]-0.6123724356957944*coeff[5]*fSkin[44]-0.3535533905932737*coeff[3]*fSkin[44]-0.5303300858899105*coeff[4]*fSkin[42]-0.3061862178478971*coeff[2]*fSkin[42]-0.5303300858899105*coeff[5]*fSkin[41]-0.3061862178478971*coeff[3]*fSkin[41]-0.6123724356957944*coeff[7]*fSkin[40]-0.3535533905932737*coeff[6]*fSkin[40]-0.5303300858899104*coeff[7]*fSkin[36]-0.3061862178478971*coeff[6]*fSkin[36]; 

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
  out[24] += boundSurf_incr[24]*rdx2Sq; 
  out[25] += boundSurf_incr[25]*rdx2Sq; 
  out[26] += boundSurf_incr[26]*rdx2Sq; 
  out[27] += boundSurf_incr[27]*rdx2Sq; 
  out[28] += boundSurf_incr[28]*rdx2Sq; 
  out[29] += boundSurf_incr[29]*rdx2Sq; 
  out[30] += boundSurf_incr[30]*rdx2Sq; 
  out[31] += boundSurf_incr[31]*rdx2Sq; 
  out[32] += boundSurf_incr[32]*rdx2Sq; 
  out[33] += boundSurf_incr[33]*rdx2Sq; 
  out[34] += boundSurf_incr[34]*rdx2Sq; 
  out[35] += boundSurf_incr[35]*rdx2Sq; 
  out[36] += boundSurf_incr[36]*rdx2Sq; 
  out[37] += boundSurf_incr[37]*rdx2Sq; 
  out[38] += boundSurf_incr[38]*rdx2Sq; 
  out[39] += boundSurf_incr[39]*rdx2Sq; 
  out[40] += boundSurf_incr[40]*rdx2Sq; 
  out[41] += boundSurf_incr[41]*rdx2Sq; 
  out[42] += boundSurf_incr[42]*rdx2Sq; 
  out[43] += boundSurf_incr[43]*rdx2Sq; 
  out[44] += boundSurf_incr[44]*rdx2Sq; 
  out[45] += boundSurf_incr[45]*rdx2Sq; 
  out[46] += boundSurf_incr[46]*rdx2Sq; 
  out[47] += boundSurf_incr[47]*rdx2Sq; 

  return 0.;
}

