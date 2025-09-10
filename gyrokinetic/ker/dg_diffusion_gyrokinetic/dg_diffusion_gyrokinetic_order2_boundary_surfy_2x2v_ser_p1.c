#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],2.);

  double vol_incr[24] = {0.0}; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[1]*qSkin[2])-0.5412658773652741*coeff[1]*qEdge[2]-0.5625*qSkin[0]*coeff[1]+0.5625*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(0.5412658773652741*coeff[1]*qSkin[5])-0.5412658773652741*coeff[1]*qEdge[5]-0.5625*coeff[1]*qSkin[1]+0.5625*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*qSkin[2])-0.4375*coeff[1]*qEdge[2]-1.4072912811497125*qSkin[0]*coeff[1]+0.5412658773652739*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(0.5412658773652741*coeff[1]*qSkin[7])-0.5412658773652741*coeff[1]*qEdge[7]-0.5625*coeff[1]*qSkin[3]+0.5625*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = -(0.5412658773652741*coeff[1]*qSkin[9])-0.5412658773652741*coeff[1]*qEdge[9]-0.5625*coeff[1]*qSkin[4]+0.5625*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[1]*qSkin[5])-0.4375*coeff[1]*qEdge[5]-1.4072912811497125*coeff[1]*qSkin[1]+0.5412658773652739*coeff[1]*qEdge[1]; 
  edgeSurf_incr[6] = -(0.5412658773652741*coeff[1]*qSkin[11])-0.5412658773652741*coeff[1]*qEdge[11]-0.5625*coeff[1]*qSkin[6]+0.5625*coeff[1]*qEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[1]*qSkin[7])-0.4375*coeff[1]*qEdge[7]-1.4072912811497125*coeff[1]*qSkin[3]+0.5412658773652739*coeff[1]*qEdge[3]; 
  edgeSurf_incr[8] = -(0.5412658773652741*coeff[1]*qSkin[12])-0.5412658773652741*coeff[1]*qEdge[12]-0.5625*coeff[1]*qSkin[8]+0.5625*coeff[1]*qEdge[8]; 
  edgeSurf_incr[9] = -(1.4375*coeff[1]*qSkin[9])-0.4375*coeff[1]*qEdge[9]-1.4072912811497125*coeff[1]*qSkin[4]+0.5412658773652739*coeff[1]*qEdge[4]; 
  edgeSurf_incr[10] = -(0.5412658773652741*coeff[1]*qSkin[14])-0.5412658773652741*coeff[1]*qEdge[14]-0.5625*coeff[1]*qSkin[10]+0.5625*coeff[1]*qEdge[10]; 
  edgeSurf_incr[11] = -(1.4375*coeff[1]*qSkin[11])-0.4375*coeff[1]*qEdge[11]-1.4072912811497125*coeff[1]*qSkin[6]+0.5412658773652739*coeff[1]*qEdge[6]; 
  edgeSurf_incr[12] = -(1.4375*coeff[1]*qSkin[12])-0.4375*coeff[1]*qEdge[12]-1.4072912811497125*coeff[1]*qSkin[8]+0.5412658773652739*coeff[1]*qEdge[8]; 
  edgeSurf_incr[13] = -(0.5412658773652741*coeff[1]*qSkin[15])-0.5412658773652741*coeff[1]*qEdge[15]-0.5625*coeff[1]*qSkin[13]+0.5625*coeff[1]*qEdge[13]; 
  edgeSurf_incr[14] = -(1.4375*coeff[1]*qSkin[14])-0.4375*coeff[1]*qEdge[14]-1.4072912811497125*coeff[1]*qSkin[10]+0.5412658773652739*coeff[1]*qEdge[10]; 
  edgeSurf_incr[15] = -(1.4375*coeff[1]*qSkin[15])-0.4375*coeff[1]*qEdge[15]-1.4072912811497125*coeff[1]*qSkin[13]+0.5412658773652739*coeff[1]*qEdge[13]; 
  edgeSurf_incr[16] = -(0.5412658773652742*coeff[1]*qSkin[18])-0.5412658773652742*coeff[1]*qEdge[18]-0.5625*coeff[1]*qSkin[16]+0.5625*coeff[1]*qEdge[16]; 
  edgeSurf_incr[17] = -(0.5412658773652742*coeff[1]*qSkin[20])-0.5412658773652742*coeff[1]*qEdge[20]-0.5625*coeff[1]*qSkin[17]+0.5625*coeff[1]*qEdge[17]; 
  edgeSurf_incr[18] = -(1.4375*coeff[1]*qSkin[18])-0.4375*coeff[1]*qEdge[18]-1.4072912811497127*coeff[1]*qSkin[16]+0.5412658773652742*coeff[1]*qEdge[16]; 
  edgeSurf_incr[19] = -(0.5412658773652742*coeff[1]*qSkin[22])-0.5412658773652742*coeff[1]*qEdge[22]-0.5625*coeff[1]*qSkin[19]+0.5625*coeff[1]*qEdge[19]; 
  edgeSurf_incr[20] = -(1.4375*coeff[1]*qSkin[20])-0.4375*coeff[1]*qEdge[20]-1.4072912811497127*coeff[1]*qSkin[17]+0.5412658773652742*coeff[1]*qEdge[17]; 
  edgeSurf_incr[21] = -(0.5412658773652742*coeff[1]*qSkin[23])-0.5412658773652742*coeff[1]*qEdge[23]-0.5625*coeff[1]*qSkin[21]+0.5625*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = -(1.4375*coeff[1]*qSkin[22])-0.4375*coeff[1]*qEdge[22]-1.4072912811497127*coeff[1]*qSkin[19]+0.5412658773652742*coeff[1]*qEdge[19]; 
  edgeSurf_incr[23] = -(1.4375*coeff[1]*qSkin[23])-0.4375*coeff[1]*qEdge[23]-1.4072912811497127*coeff[1]*qSkin[21]+0.5412658773652742*coeff[1]*qEdge[21]; 

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

  edgeSurf_incr[0] = 0.5412658773652741*coeff[1]*qSkin[2]+0.5412658773652741*coeff[1]*qEdge[2]-0.5625*qSkin[0]*coeff[1]+0.5625*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 0.5412658773652741*coeff[1]*qSkin[5]+0.5412658773652741*coeff[1]*qEdge[5]-0.5625*coeff[1]*qSkin[1]+0.5625*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*qSkin[2])-0.4375*coeff[1]*qEdge[2]+1.4072912811497125*qSkin[0]*coeff[1]-0.5412658773652739*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[1]*qSkin[7]+0.5412658773652741*coeff[1]*qEdge[7]-0.5625*coeff[1]*qSkin[3]+0.5625*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = 0.5412658773652741*coeff[1]*qSkin[9]+0.5412658773652741*coeff[1]*qEdge[9]-0.5625*coeff[1]*qSkin[4]+0.5625*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[1]*qSkin[5])-0.4375*coeff[1]*qEdge[5]+1.4072912811497125*coeff[1]*qSkin[1]-0.5412658773652739*coeff[1]*qEdge[1]; 
  edgeSurf_incr[6] = 0.5412658773652741*coeff[1]*qSkin[11]+0.5412658773652741*coeff[1]*qEdge[11]-0.5625*coeff[1]*qSkin[6]+0.5625*coeff[1]*qEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[1]*qSkin[7])-0.4375*coeff[1]*qEdge[7]+1.4072912811497125*coeff[1]*qSkin[3]-0.5412658773652739*coeff[1]*qEdge[3]; 
  edgeSurf_incr[8] = 0.5412658773652741*coeff[1]*qSkin[12]+0.5412658773652741*coeff[1]*qEdge[12]-0.5625*coeff[1]*qSkin[8]+0.5625*coeff[1]*qEdge[8]; 
  edgeSurf_incr[9] = -(1.4375*coeff[1]*qSkin[9])-0.4375*coeff[1]*qEdge[9]+1.4072912811497125*coeff[1]*qSkin[4]-0.5412658773652739*coeff[1]*qEdge[4]; 
  edgeSurf_incr[10] = 0.5412658773652741*coeff[1]*qSkin[14]+0.5412658773652741*coeff[1]*qEdge[14]-0.5625*coeff[1]*qSkin[10]+0.5625*coeff[1]*qEdge[10]; 
  edgeSurf_incr[11] = -(1.4375*coeff[1]*qSkin[11])-0.4375*coeff[1]*qEdge[11]+1.4072912811497125*coeff[1]*qSkin[6]-0.5412658773652739*coeff[1]*qEdge[6]; 
  edgeSurf_incr[12] = -(1.4375*coeff[1]*qSkin[12])-0.4375*coeff[1]*qEdge[12]+1.4072912811497125*coeff[1]*qSkin[8]-0.5412658773652739*coeff[1]*qEdge[8]; 
  edgeSurf_incr[13] = 0.5412658773652741*coeff[1]*qSkin[15]+0.5412658773652741*coeff[1]*qEdge[15]-0.5625*coeff[1]*qSkin[13]+0.5625*coeff[1]*qEdge[13]; 
  edgeSurf_incr[14] = -(1.4375*coeff[1]*qSkin[14])-0.4375*coeff[1]*qEdge[14]+1.4072912811497125*coeff[1]*qSkin[10]-0.5412658773652739*coeff[1]*qEdge[10]; 
  edgeSurf_incr[15] = -(1.4375*coeff[1]*qSkin[15])-0.4375*coeff[1]*qEdge[15]+1.4072912811497125*coeff[1]*qSkin[13]-0.5412658773652739*coeff[1]*qEdge[13]; 
  edgeSurf_incr[16] = 0.5412658773652742*coeff[1]*qSkin[18]+0.5412658773652742*coeff[1]*qEdge[18]-0.5625*coeff[1]*qSkin[16]+0.5625*coeff[1]*qEdge[16]; 
  edgeSurf_incr[17] = 0.5412658773652742*coeff[1]*qSkin[20]+0.5412658773652742*coeff[1]*qEdge[20]-0.5625*coeff[1]*qSkin[17]+0.5625*coeff[1]*qEdge[17]; 
  edgeSurf_incr[18] = -(1.4375*coeff[1]*qSkin[18])-0.4375*coeff[1]*qEdge[18]+1.4072912811497127*coeff[1]*qSkin[16]-0.5412658773652742*coeff[1]*qEdge[16]; 
  edgeSurf_incr[19] = 0.5412658773652742*coeff[1]*qSkin[22]+0.5412658773652742*coeff[1]*qEdge[22]-0.5625*coeff[1]*qSkin[19]+0.5625*coeff[1]*qEdge[19]; 
  edgeSurf_incr[20] = -(1.4375*coeff[1]*qSkin[20])-0.4375*coeff[1]*qEdge[20]+1.4072912811497127*coeff[1]*qSkin[17]-0.5412658773652742*coeff[1]*qEdge[17]; 
  edgeSurf_incr[21] = 0.5412658773652742*coeff[1]*qSkin[23]+0.5412658773652742*coeff[1]*qEdge[23]-0.5625*coeff[1]*qSkin[21]+0.5625*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = -(1.4375*coeff[1]*qSkin[22])-0.4375*coeff[1]*qEdge[22]+1.4072912811497127*coeff[1]*qSkin[19]-0.5412658773652742*coeff[1]*qEdge[19]; 
  edgeSurf_incr[23] = -(1.4375*coeff[1]*qSkin[23])-0.4375*coeff[1]*qEdge[23]+1.4072912811497127*coeff[1]*qSkin[21]-0.5412658773652742*coeff[1]*qEdge[21]; 

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

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdx2Sq; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdx2Sq; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdx2Sq; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdx2Sq; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdx2Sq; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdx2Sq; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdx2Sq; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdx2Sq; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdx2Sq; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdx2Sq; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdx2Sq; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],2.);

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

  double fEdge[24];
  fEdge[0] = 0.5*(jacobgeo_inv[3]*qEdge[5]+jacobgeo_inv[2]*qEdge[2]+jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.5*(jacobgeo_inv[2]*qEdge[5]+qEdge[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1]); 
  fEdge[2] = 0.5*(jacobgeo_inv[1]*qEdge[5]+qEdge[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*qEdge[2]+qEdge[0]*jacobgeo_inv[2]); 
  fEdge[3] = 0.5*(jacobgeo_inv[3]*qEdge[11]+jacobgeo_inv[2]*qEdge[7]+jacobgeo_inv[1]*qEdge[6]+jacobgeo_inv[0]*qEdge[3]); 
  fEdge[4] = 0.5*(jacobgeo_inv[3]*qEdge[12]+jacobgeo_inv[2]*qEdge[9]+jacobgeo_inv[1]*qEdge[8]+jacobgeo_inv[0]*qEdge[4]); 
  fEdge[5] = 0.5*(jacobgeo_inv[0]*qEdge[5]+qEdge[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*qEdge[2]+qEdge[1]*jacobgeo_inv[2]); 
  fEdge[6] = 0.5*(jacobgeo_inv[2]*qEdge[11]+jacobgeo_inv[3]*qEdge[7]+jacobgeo_inv[0]*qEdge[6]+jacobgeo_inv[1]*qEdge[3]); 
  fEdge[7] = 0.5*(jacobgeo_inv[1]*qEdge[11]+jacobgeo_inv[0]*qEdge[7]+jacobgeo_inv[3]*qEdge[6]+jacobgeo_inv[2]*qEdge[3]); 
  fEdge[8] = 0.5*(jacobgeo_inv[2]*qEdge[12]+jacobgeo_inv[3]*qEdge[9]+jacobgeo_inv[0]*qEdge[8]+jacobgeo_inv[1]*qEdge[4]); 
  fEdge[9] = 0.5*(jacobgeo_inv[1]*qEdge[12]+jacobgeo_inv[0]*qEdge[9]+jacobgeo_inv[3]*qEdge[8]+jacobgeo_inv[2]*qEdge[4]); 
  fEdge[10] = 0.5*(jacobgeo_inv[3]*qEdge[15]+jacobgeo_inv[2]*qEdge[14]+jacobgeo_inv[1]*qEdge[13]+jacobgeo_inv[0]*qEdge[10]); 
  fEdge[11] = 0.5*(jacobgeo_inv[0]*qEdge[11]+jacobgeo_inv[1]*qEdge[7]+jacobgeo_inv[2]*qEdge[6]+jacobgeo_inv[3]*qEdge[3]); 
  fEdge[12] = 0.5*(jacobgeo_inv[0]*qEdge[12]+jacobgeo_inv[1]*qEdge[9]+jacobgeo_inv[2]*qEdge[8]+jacobgeo_inv[3]*qEdge[4]); 
  fEdge[13] = 0.5*(jacobgeo_inv[2]*qEdge[15]+jacobgeo_inv[3]*qEdge[14]+jacobgeo_inv[0]*qEdge[13]+jacobgeo_inv[1]*qEdge[10]); 
  fEdge[14] = 0.5*(jacobgeo_inv[1]*qEdge[15]+jacobgeo_inv[0]*qEdge[14]+jacobgeo_inv[3]*qEdge[13]+jacobgeo_inv[2]*qEdge[10]); 
  fEdge[15] = 0.5*(jacobgeo_inv[0]*qEdge[15]+jacobgeo_inv[1]*qEdge[14]+jacobgeo_inv[2]*qEdge[13]+jacobgeo_inv[3]*qEdge[10]); 
  fEdge[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qEdge[20]+15.000000000000002*(jacobgeo_inv[2]*qEdge[18]+jacobgeo_inv[1]*qEdge[17])+15.0*jacobgeo_inv[0]*qEdge[16]); 
  fEdge[17] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*qEdge[20]+15.0*(jacobgeo_inv[3]*qEdge[18]+jacobgeo_inv[0]*qEdge[17])+15.000000000000002*jacobgeo_inv[1]*qEdge[16]); 
  fEdge[18] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*qEdge[20]+15.0*(jacobgeo_inv[0]*qEdge[18]+jacobgeo_inv[3]*qEdge[17])+15.000000000000002*jacobgeo_inv[2]*qEdge[16]); 
  fEdge[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qEdge[23]+15.000000000000002*(jacobgeo_inv[2]*qEdge[22]+jacobgeo_inv[1]*qEdge[21])+15.0*jacobgeo_inv[0]*qEdge[19]); 
  fEdge[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qEdge[20]+15.000000000000002*(jacobgeo_inv[1]*qEdge[18]+jacobgeo_inv[2]*qEdge[17])+15.0*jacobgeo_inv[3]*qEdge[16]); 
  fEdge[21] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*qEdge[23]+15.0*(jacobgeo_inv[3]*qEdge[22]+jacobgeo_inv[0]*qEdge[21])+15.000000000000002*jacobgeo_inv[1]*qEdge[19]); 
  fEdge[22] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*qEdge[23]+15.0*(jacobgeo_inv[0]*qEdge[22]+jacobgeo_inv[3]*qEdge[21])+15.000000000000002*jacobgeo_inv[2]*qEdge[19]); 
  fEdge[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qEdge[23]+15.000000000000002*(jacobgeo_inv[1]*qEdge[22]+jacobgeo_inv[2]*qEdge[21])+15.0*jacobgeo_inv[3]*qEdge[19]); 

  double vol_incr[24] = {0.0}; 
  vol_incr[2] = 1.5*fSkin[1]*coeff[7]+1.5*fSkin[0]*coeff[6]; 
  vol_incr[5] = 1.5*fSkin[0]*coeff[7]+1.5*fSkin[1]*coeff[6]; 
  vol_incr[7] = 1.5*fSkin[6]*coeff[7]+1.5*fSkin[3]*coeff[6]; 
  vol_incr[9] = 1.5*coeff[7]*fSkin[8]+1.5*fSkin[4]*coeff[6]; 
  vol_incr[11] = 1.5*fSkin[3]*coeff[7]+1.5*coeff[6]*fSkin[6]; 
  vol_incr[12] = 1.5*coeff[6]*fSkin[8]+1.5*fSkin[4]*coeff[7]; 
  vol_incr[14] = 1.5*coeff[7]*fSkin[13]+1.5*coeff[6]*fSkin[10]; 
  vol_incr[15] = 1.5*coeff[6]*fSkin[13]+1.5*coeff[7]*fSkin[10]; 
  vol_incr[18] = 1.5*coeff[7]*fSkin[17]+1.5*coeff[6]*fSkin[16]; 
  vol_incr[20] = 1.5*coeff[6]*fSkin[17]+1.5*coeff[7]*fSkin[16]; 
  vol_incr[22] = 1.5*coeff[7]*fSkin[21]+1.5*coeff[6]*fSkin[19]; 
  vol_incr[23] = 1.5*coeff[6]*fSkin[21]+1.5*coeff[7]*fSkin[19]; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.27063293868263705*coeff[5]*fSkin[5])-0.27063293868263705*coeff[5]*fEdge[5]-0.28125*fSkin[1]*coeff[5]+0.28125*fEdge[1]*coeff[5]-0.27063293868263705*fSkin[2]*coeff[4]-0.27063293868263705*fEdge[2]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = -(0.27063293868263705*coeff[4]*fSkin[5])-0.27063293868263705*coeff[4]*fEdge[5]-0.27063293868263705*fSkin[2]*coeff[5]-0.27063293868263705*fEdge[2]*coeff[5]-0.28125*fSkin[0]*coeff[5]+0.28125*fEdge[0]*coeff[5]-0.28125*fSkin[1]*coeff[4]+0.28125*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = -(0.4330127018922193*fSkin[5]*coeff[7])+0.4330127018922193*fEdge[5]*coeff[7]-0.375*fSkin[1]*coeff[7]-0.375*fEdge[1]*coeff[7]-0.4330127018922193*fSkin[2]*coeff[6]+0.4330127018922193*fEdge[2]*coeff[6]-0.375*fSkin[0]*coeff[6]-0.375*fEdge[0]*coeff[6]-0.71875*coeff[5]*fSkin[5]-0.21875*coeff[5]*fEdge[5]-0.7036456405748562*fSkin[1]*coeff[5]+0.27063293868263694*fEdge[1]*coeff[5]-0.71875*fSkin[2]*coeff[4]-0.21875*fEdge[2]*coeff[4]-0.7036456405748562*fSkin[0]*coeff[4]+0.27063293868263694*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = -(0.27063293868263705*coeff[5]*fSkin[11])-0.27063293868263705*coeff[5]*fEdge[11]-0.27063293868263705*coeff[4]*fSkin[7]-0.27063293868263705*coeff[4]*fEdge[7]-0.28125*coeff[5]*fSkin[6]+0.28125*coeff[5]*fEdge[6]-0.28125*fSkin[3]*coeff[4]+0.28125*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = -(0.27063293868263705*coeff[5]*fSkin[12])-0.27063293868263705*coeff[5]*fEdge[12]-0.27063293868263705*coeff[4]*fSkin[9]-0.27063293868263705*coeff[4]*fEdge[9]-0.28125*coeff[5]*fSkin[8]+0.28125*coeff[5]*fEdge[8]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = -(0.4330127018922193*fSkin[2]*coeff[7])+0.4330127018922193*fEdge[2]*coeff[7]-0.375*fSkin[0]*coeff[7]-0.375*fEdge[0]*coeff[7]-0.4330127018922193*fSkin[5]*coeff[6]+0.4330127018922193*fEdge[5]*coeff[6]-0.375*fSkin[1]*coeff[6]-0.375*fEdge[1]*coeff[6]-0.71875*coeff[4]*fSkin[5]-0.21875*coeff[4]*fEdge[5]-0.71875*fSkin[2]*coeff[5]-0.21875*fEdge[2]*coeff[5]-0.7036456405748562*fSkin[0]*coeff[5]+0.27063293868263694*fEdge[0]*coeff[5]-0.7036456405748562*fSkin[1]*coeff[4]+0.27063293868263694*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = -(0.27063293868263705*coeff[4]*fSkin[11])-0.27063293868263705*coeff[4]*fEdge[11]-0.27063293868263705*coeff[5]*fSkin[7]-0.27063293868263705*coeff[5]*fEdge[7]-0.28125*coeff[4]*fSkin[6]+0.28125*coeff[4]*fEdge[6]-0.28125*fSkin[3]*coeff[5]+0.28125*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = -(0.4330127018922193*coeff[7]*fSkin[11])-0.71875*coeff[5]*fSkin[11]+0.4330127018922193*coeff[7]*fEdge[11]-0.21875*coeff[5]*fEdge[11]-0.4330127018922193*coeff[6]*fSkin[7]-0.71875*coeff[4]*fSkin[7]+0.4330127018922193*coeff[6]*fEdge[7]-0.21875*coeff[4]*fEdge[7]-0.375*fSkin[6]*coeff[7]-0.375*fEdge[6]*coeff[7]-0.7036456405748562*coeff[5]*fSkin[6]+0.27063293868263694*coeff[5]*fEdge[6]-0.375*fSkin[3]*coeff[6]-0.375*fEdge[3]*coeff[6]-0.7036456405748562*fSkin[3]*coeff[4]+0.27063293868263694*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = -(0.27063293868263705*coeff[4]*fSkin[12])-0.27063293868263705*coeff[4]*fEdge[12]-0.27063293868263705*coeff[5]*fSkin[9]-0.27063293868263705*coeff[5]*fEdge[9]-0.28125*coeff[4]*fSkin[8]+0.28125*coeff[4]*fEdge[8]-0.28125*fSkin[4]*coeff[5]+0.28125*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = -(0.4330127018922193*coeff[7]*fSkin[12])-0.71875*coeff[5]*fSkin[12]+0.4330127018922193*coeff[7]*fEdge[12]-0.21875*coeff[5]*fEdge[12]-0.4330127018922193*coeff[6]*fSkin[9]-0.71875*coeff[4]*fSkin[9]+0.4330127018922193*coeff[6]*fEdge[9]-0.21875*coeff[4]*fEdge[9]-0.375*coeff[7]*fSkin[8]-0.7036456405748562*coeff[5]*fSkin[8]-0.375*coeff[7]*fEdge[8]+0.27063293868263694*coeff[5]*fEdge[8]-0.375*fSkin[4]*coeff[6]-0.375*fEdge[4]*coeff[6]-0.7036456405748562*coeff[4]*fSkin[4]+0.27063293868263694*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = -(0.27063293868263705*coeff[5]*fSkin[15])-0.27063293868263705*coeff[5]*fEdge[15]-0.27063293868263705*coeff[4]*fSkin[14]-0.27063293868263705*coeff[4]*fEdge[14]-0.28125*coeff[5]*fSkin[13]+0.28125*coeff[5]*fEdge[13]-0.28125*coeff[4]*fSkin[10]+0.28125*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.4330127018922193*coeff[6]*fSkin[11])-0.71875*coeff[4]*fSkin[11]+0.4330127018922193*coeff[6]*fEdge[11]-0.21875*coeff[4]*fEdge[11]-0.4330127018922193*coeff[7]*fSkin[7]-0.71875*coeff[5]*fSkin[7]+0.4330127018922193*coeff[7]*fEdge[7]-0.21875*coeff[5]*fEdge[7]-0.375*fSkin[3]*coeff[7]-0.375*fEdge[3]*coeff[7]-0.375*coeff[6]*fSkin[6]-0.7036456405748562*coeff[4]*fSkin[6]-0.375*coeff[6]*fEdge[6]+0.27063293868263694*coeff[4]*fEdge[6]-0.7036456405748562*fSkin[3]*coeff[5]+0.27063293868263694*fEdge[3]*coeff[5]; 
  edgeSurf_incr[12] = -(0.4330127018922193*coeff[6]*fSkin[12])-0.71875*coeff[4]*fSkin[12]+0.4330127018922193*coeff[6]*fEdge[12]-0.21875*coeff[4]*fEdge[12]-0.4330127018922193*coeff[7]*fSkin[9]-0.71875*coeff[5]*fSkin[9]+0.4330127018922193*coeff[7]*fEdge[9]-0.21875*coeff[5]*fEdge[9]-0.375*coeff[6]*fSkin[8]-0.7036456405748562*coeff[4]*fSkin[8]-0.375*coeff[6]*fEdge[8]+0.27063293868263694*coeff[4]*fEdge[8]-0.375*fSkin[4]*coeff[7]-0.375*fEdge[4]*coeff[7]-0.7036456405748562*fSkin[4]*coeff[5]+0.27063293868263694*fEdge[4]*coeff[5]; 
  edgeSurf_incr[13] = -(0.27063293868263705*coeff[4]*fSkin[15])-0.27063293868263705*coeff[4]*fEdge[15]-0.27063293868263705*coeff[5]*fSkin[14]-0.27063293868263705*coeff[5]*fEdge[14]-0.28125*coeff[4]*fSkin[13]+0.28125*coeff[4]*fEdge[13]-0.28125*coeff[5]*fSkin[10]+0.28125*coeff[5]*fEdge[10]; 
  edgeSurf_incr[14] = -(0.4330127018922193*coeff[7]*fSkin[15])-0.71875*coeff[5]*fSkin[15]+0.4330127018922193*coeff[7]*fEdge[15]-0.21875*coeff[5]*fEdge[15]-0.4330127018922193*coeff[6]*fSkin[14]-0.71875*coeff[4]*fSkin[14]+0.4330127018922193*coeff[6]*fEdge[14]-0.21875*coeff[4]*fEdge[14]-0.375*coeff[7]*fSkin[13]-0.7036456405748562*coeff[5]*fSkin[13]-0.375*coeff[7]*fEdge[13]+0.27063293868263694*coeff[5]*fEdge[13]-0.375*coeff[6]*fSkin[10]-0.7036456405748562*coeff[4]*fSkin[10]-0.375*coeff[6]*fEdge[10]+0.27063293868263694*coeff[4]*fEdge[10]; 
  edgeSurf_incr[15] = -(0.4330127018922193*coeff[6]*fSkin[15])-0.71875*coeff[4]*fSkin[15]+0.4330127018922193*coeff[6]*fEdge[15]-0.21875*coeff[4]*fEdge[15]-0.4330127018922193*coeff[7]*fSkin[14]-0.71875*coeff[5]*fSkin[14]+0.4330127018922193*coeff[7]*fEdge[14]-0.21875*coeff[5]*fEdge[14]-0.375*coeff[6]*fSkin[13]-0.7036456405748562*coeff[4]*fSkin[13]-0.375*coeff[6]*fEdge[13]+0.27063293868263694*coeff[4]*fEdge[13]-0.375*coeff[7]*fSkin[10]-0.7036456405748562*coeff[5]*fSkin[10]-0.375*coeff[7]*fEdge[10]+0.27063293868263694*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = -(0.27063293868263705*coeff[5]*fSkin[20])-0.27063293868263705*coeff[5]*fEdge[20]-0.2706329386826371*coeff[4]*fSkin[18]-0.2706329386826371*coeff[4]*fEdge[18]-0.28124999999999994*coeff[5]*fSkin[17]+0.28124999999999994*coeff[5]*fEdge[17]-0.28125*coeff[4]*fSkin[16]+0.28125*coeff[4]*fEdge[16]; 
  edgeSurf_incr[17] = -(0.2706329386826371*coeff[4]*fSkin[20])-0.2706329386826371*coeff[4]*fEdge[20]-0.27063293868263705*coeff[5]*fSkin[18]-0.27063293868263705*coeff[5]*fEdge[18]-0.28125*coeff[4]*fSkin[17]+0.28125*coeff[4]*fEdge[17]-0.28124999999999994*coeff[5]*fSkin[16]+0.28124999999999994*coeff[5]*fEdge[16]; 
  edgeSurf_incr[18] = -(0.43301270189221935*coeff[7]*fSkin[20])-0.7187500000000001*coeff[5]*fSkin[20]+0.43301270189221935*coeff[7]*fEdge[20]-0.21875*coeff[5]*fEdge[20]-0.4330127018922193*coeff[6]*fSkin[18]-0.71875*coeff[4]*fSkin[18]+0.4330127018922193*coeff[6]*fEdge[18]-0.21875*coeff[4]*fEdge[18]-0.375*coeff[7]*fSkin[17]-0.7036456405748562*coeff[5]*fSkin[17]-0.375*coeff[7]*fEdge[17]+0.27063293868263694*coeff[5]*fEdge[17]-0.375*coeff[6]*fSkin[16]-0.7036456405748563*coeff[4]*fSkin[16]-0.375*coeff[6]*fEdge[16]+0.2706329386826371*coeff[4]*fEdge[16]; 
  edgeSurf_incr[19] = -(0.27063293868263705*coeff[5]*fSkin[23])-0.27063293868263705*coeff[5]*fEdge[23]-0.2706329386826371*coeff[4]*fSkin[22]-0.2706329386826371*coeff[4]*fEdge[22]-0.28124999999999994*coeff[5]*fSkin[21]+0.28124999999999994*coeff[5]*fEdge[21]-0.28125*coeff[4]*fSkin[19]+0.28125*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = -(0.4330127018922193*coeff[6]*fSkin[20])-0.71875*coeff[4]*fSkin[20]+0.4330127018922193*coeff[6]*fEdge[20]-0.21875*coeff[4]*fEdge[20]-0.43301270189221935*coeff[7]*fSkin[18]-0.7187500000000001*coeff[5]*fSkin[18]+0.43301270189221935*coeff[7]*fEdge[18]-0.21875*coeff[5]*fEdge[18]-0.375*coeff[6]*fSkin[17]-0.7036456405748563*coeff[4]*fSkin[17]-0.375*coeff[6]*fEdge[17]+0.2706329386826371*coeff[4]*fEdge[17]-0.375*coeff[7]*fSkin[16]-0.7036456405748562*coeff[5]*fSkin[16]-0.375*coeff[7]*fEdge[16]+0.27063293868263694*coeff[5]*fEdge[16]; 
  edgeSurf_incr[21] = -(0.2706329386826371*coeff[4]*fSkin[23])-0.2706329386826371*coeff[4]*fEdge[23]-0.27063293868263705*coeff[5]*fSkin[22]-0.27063293868263705*coeff[5]*fEdge[22]-0.28125*coeff[4]*fSkin[21]+0.28125*coeff[4]*fEdge[21]-0.28124999999999994*coeff[5]*fSkin[19]+0.28124999999999994*coeff[5]*fEdge[19]; 
  edgeSurf_incr[22] = -(0.43301270189221935*coeff[7]*fSkin[23])-0.7187500000000001*coeff[5]*fSkin[23]+0.43301270189221935*coeff[7]*fEdge[23]-0.21875*coeff[5]*fEdge[23]-0.4330127018922193*coeff[6]*fSkin[22]-0.71875*coeff[4]*fSkin[22]+0.4330127018922193*coeff[6]*fEdge[22]-0.21875*coeff[4]*fEdge[22]-0.375*coeff[7]*fSkin[21]-0.7036456405748562*coeff[5]*fSkin[21]-0.375*coeff[7]*fEdge[21]+0.27063293868263694*coeff[5]*fEdge[21]-0.375*coeff[6]*fSkin[19]-0.7036456405748563*coeff[4]*fSkin[19]-0.375*coeff[6]*fEdge[19]+0.2706329386826371*coeff[4]*fEdge[19]; 
  edgeSurf_incr[23] = -(0.4330127018922193*coeff[6]*fSkin[23])-0.71875*coeff[4]*fSkin[23]+0.4330127018922193*coeff[6]*fEdge[23]-0.21875*coeff[4]*fEdge[23]-0.43301270189221935*coeff[7]*fSkin[22]-0.7187500000000001*coeff[5]*fSkin[22]+0.43301270189221935*coeff[7]*fEdge[22]-0.21875*coeff[5]*fEdge[22]-0.375*coeff[6]*fSkin[21]-0.7036456405748563*coeff[4]*fSkin[21]-0.375*coeff[6]*fEdge[21]+0.2706329386826371*coeff[4]*fEdge[21]-0.375*coeff[7]*fSkin[19]-0.7036456405748562*coeff[5]*fSkin[19]-0.375*coeff[7]*fEdge[19]+0.27063293868263694*coeff[5]*fEdge[19]; 

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

  edgeSurf_incr[0] = 0.27063293868263705*coeff[5]*fSkin[5]+0.27063293868263705*coeff[5]*fEdge[5]-0.28125*fSkin[1]*coeff[5]+0.28125*fEdge[1]*coeff[5]+0.27063293868263705*fSkin[2]*coeff[4]+0.27063293868263705*fEdge[2]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = 0.27063293868263705*coeff[4]*fSkin[5]+0.27063293868263705*coeff[4]*fEdge[5]+0.27063293868263705*fSkin[2]*coeff[5]+0.27063293868263705*fEdge[2]*coeff[5]-0.28125*fSkin[0]*coeff[5]+0.28125*fEdge[0]*coeff[5]-0.28125*fSkin[1]*coeff[4]+0.28125*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = 0.4330127018922193*fSkin[5]*coeff[7]-0.4330127018922193*fEdge[5]*coeff[7]-0.375*fSkin[1]*coeff[7]-0.375*fEdge[1]*coeff[7]+0.4330127018922193*fSkin[2]*coeff[6]-0.4330127018922193*fEdge[2]*coeff[6]-0.375*fSkin[0]*coeff[6]-0.375*fEdge[0]*coeff[6]-0.71875*coeff[5]*fSkin[5]-0.21875*coeff[5]*fEdge[5]+0.7036456405748562*fSkin[1]*coeff[5]-0.27063293868263694*fEdge[1]*coeff[5]-0.71875*fSkin[2]*coeff[4]-0.21875*fEdge[2]*coeff[4]+0.7036456405748562*fSkin[0]*coeff[4]-0.27063293868263694*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = 0.27063293868263705*coeff[5]*fSkin[11]+0.27063293868263705*coeff[5]*fEdge[11]+0.27063293868263705*coeff[4]*fSkin[7]+0.27063293868263705*coeff[4]*fEdge[7]-0.28125*coeff[5]*fSkin[6]+0.28125*coeff[5]*fEdge[6]-0.28125*fSkin[3]*coeff[4]+0.28125*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = 0.27063293868263705*coeff[5]*fSkin[12]+0.27063293868263705*coeff[5]*fEdge[12]+0.27063293868263705*coeff[4]*fSkin[9]+0.27063293868263705*coeff[4]*fEdge[9]-0.28125*coeff[5]*fSkin[8]+0.28125*coeff[5]*fEdge[8]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = 0.4330127018922193*fSkin[2]*coeff[7]-0.4330127018922193*fEdge[2]*coeff[7]-0.375*fSkin[0]*coeff[7]-0.375*fEdge[0]*coeff[7]+0.4330127018922193*fSkin[5]*coeff[6]-0.4330127018922193*fEdge[5]*coeff[6]-0.375*fSkin[1]*coeff[6]-0.375*fEdge[1]*coeff[6]-0.71875*coeff[4]*fSkin[5]-0.21875*coeff[4]*fEdge[5]-0.71875*fSkin[2]*coeff[5]-0.21875*fEdge[2]*coeff[5]+0.7036456405748562*fSkin[0]*coeff[5]-0.27063293868263694*fEdge[0]*coeff[5]+0.7036456405748562*fSkin[1]*coeff[4]-0.27063293868263694*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = 0.27063293868263705*coeff[4]*fSkin[11]+0.27063293868263705*coeff[4]*fEdge[11]+0.27063293868263705*coeff[5]*fSkin[7]+0.27063293868263705*coeff[5]*fEdge[7]-0.28125*coeff[4]*fSkin[6]+0.28125*coeff[4]*fEdge[6]-0.28125*fSkin[3]*coeff[5]+0.28125*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = 0.4330127018922193*coeff[7]*fSkin[11]-0.71875*coeff[5]*fSkin[11]-0.4330127018922193*coeff[7]*fEdge[11]-0.21875*coeff[5]*fEdge[11]+0.4330127018922193*coeff[6]*fSkin[7]-0.71875*coeff[4]*fSkin[7]-0.4330127018922193*coeff[6]*fEdge[7]-0.21875*coeff[4]*fEdge[7]-0.375*fSkin[6]*coeff[7]-0.375*fEdge[6]*coeff[7]+0.7036456405748562*coeff[5]*fSkin[6]-0.27063293868263694*coeff[5]*fEdge[6]-0.375*fSkin[3]*coeff[6]-0.375*fEdge[3]*coeff[6]+0.7036456405748562*fSkin[3]*coeff[4]-0.27063293868263694*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = 0.27063293868263705*coeff[4]*fSkin[12]+0.27063293868263705*coeff[4]*fEdge[12]+0.27063293868263705*coeff[5]*fSkin[9]+0.27063293868263705*coeff[5]*fEdge[9]-0.28125*coeff[4]*fSkin[8]+0.28125*coeff[4]*fEdge[8]-0.28125*fSkin[4]*coeff[5]+0.28125*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = 0.4330127018922193*coeff[7]*fSkin[12]-0.71875*coeff[5]*fSkin[12]-0.4330127018922193*coeff[7]*fEdge[12]-0.21875*coeff[5]*fEdge[12]+0.4330127018922193*coeff[6]*fSkin[9]-0.71875*coeff[4]*fSkin[9]-0.4330127018922193*coeff[6]*fEdge[9]-0.21875*coeff[4]*fEdge[9]-0.375*coeff[7]*fSkin[8]+0.7036456405748562*coeff[5]*fSkin[8]-0.375*coeff[7]*fEdge[8]-0.27063293868263694*coeff[5]*fEdge[8]-0.375*fSkin[4]*coeff[6]-0.375*fEdge[4]*coeff[6]+0.7036456405748562*coeff[4]*fSkin[4]-0.27063293868263694*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = 0.27063293868263705*coeff[5]*fSkin[15]+0.27063293868263705*coeff[5]*fEdge[15]+0.27063293868263705*coeff[4]*fSkin[14]+0.27063293868263705*coeff[4]*fEdge[14]-0.28125*coeff[5]*fSkin[13]+0.28125*coeff[5]*fEdge[13]-0.28125*coeff[4]*fSkin[10]+0.28125*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = 0.4330127018922193*coeff[6]*fSkin[11]-0.71875*coeff[4]*fSkin[11]-0.4330127018922193*coeff[6]*fEdge[11]-0.21875*coeff[4]*fEdge[11]+0.4330127018922193*coeff[7]*fSkin[7]-0.71875*coeff[5]*fSkin[7]-0.4330127018922193*coeff[7]*fEdge[7]-0.21875*coeff[5]*fEdge[7]-0.375*fSkin[3]*coeff[7]-0.375*fEdge[3]*coeff[7]-0.375*coeff[6]*fSkin[6]+0.7036456405748562*coeff[4]*fSkin[6]-0.375*coeff[6]*fEdge[6]-0.27063293868263694*coeff[4]*fEdge[6]+0.7036456405748562*fSkin[3]*coeff[5]-0.27063293868263694*fEdge[3]*coeff[5]; 
  edgeSurf_incr[12] = 0.4330127018922193*coeff[6]*fSkin[12]-0.71875*coeff[4]*fSkin[12]-0.4330127018922193*coeff[6]*fEdge[12]-0.21875*coeff[4]*fEdge[12]+0.4330127018922193*coeff[7]*fSkin[9]-0.71875*coeff[5]*fSkin[9]-0.4330127018922193*coeff[7]*fEdge[9]-0.21875*coeff[5]*fEdge[9]-0.375*coeff[6]*fSkin[8]+0.7036456405748562*coeff[4]*fSkin[8]-0.375*coeff[6]*fEdge[8]-0.27063293868263694*coeff[4]*fEdge[8]-0.375*fSkin[4]*coeff[7]-0.375*fEdge[4]*coeff[7]+0.7036456405748562*fSkin[4]*coeff[5]-0.27063293868263694*fEdge[4]*coeff[5]; 
  edgeSurf_incr[13] = 0.27063293868263705*coeff[4]*fSkin[15]+0.27063293868263705*coeff[4]*fEdge[15]+0.27063293868263705*coeff[5]*fSkin[14]+0.27063293868263705*coeff[5]*fEdge[14]-0.28125*coeff[4]*fSkin[13]+0.28125*coeff[4]*fEdge[13]-0.28125*coeff[5]*fSkin[10]+0.28125*coeff[5]*fEdge[10]; 
  edgeSurf_incr[14] = 0.4330127018922193*coeff[7]*fSkin[15]-0.71875*coeff[5]*fSkin[15]-0.4330127018922193*coeff[7]*fEdge[15]-0.21875*coeff[5]*fEdge[15]+0.4330127018922193*coeff[6]*fSkin[14]-0.71875*coeff[4]*fSkin[14]-0.4330127018922193*coeff[6]*fEdge[14]-0.21875*coeff[4]*fEdge[14]-0.375*coeff[7]*fSkin[13]+0.7036456405748562*coeff[5]*fSkin[13]-0.375*coeff[7]*fEdge[13]-0.27063293868263694*coeff[5]*fEdge[13]-0.375*coeff[6]*fSkin[10]+0.7036456405748562*coeff[4]*fSkin[10]-0.375*coeff[6]*fEdge[10]-0.27063293868263694*coeff[4]*fEdge[10]; 
  edgeSurf_incr[15] = 0.4330127018922193*coeff[6]*fSkin[15]-0.71875*coeff[4]*fSkin[15]-0.4330127018922193*coeff[6]*fEdge[15]-0.21875*coeff[4]*fEdge[15]+0.4330127018922193*coeff[7]*fSkin[14]-0.71875*coeff[5]*fSkin[14]-0.4330127018922193*coeff[7]*fEdge[14]-0.21875*coeff[5]*fEdge[14]-0.375*coeff[6]*fSkin[13]+0.7036456405748562*coeff[4]*fSkin[13]-0.375*coeff[6]*fEdge[13]-0.27063293868263694*coeff[4]*fEdge[13]-0.375*coeff[7]*fSkin[10]+0.7036456405748562*coeff[5]*fSkin[10]-0.375*coeff[7]*fEdge[10]-0.27063293868263694*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = 0.27063293868263705*coeff[5]*fSkin[20]+0.27063293868263705*coeff[5]*fEdge[20]+0.2706329386826371*coeff[4]*fSkin[18]+0.2706329386826371*coeff[4]*fEdge[18]-0.28124999999999994*coeff[5]*fSkin[17]+0.28124999999999994*coeff[5]*fEdge[17]-0.28125*coeff[4]*fSkin[16]+0.28125*coeff[4]*fEdge[16]; 
  edgeSurf_incr[17] = 0.2706329386826371*coeff[4]*fSkin[20]+0.2706329386826371*coeff[4]*fEdge[20]+0.27063293868263705*coeff[5]*fSkin[18]+0.27063293868263705*coeff[5]*fEdge[18]-0.28125*coeff[4]*fSkin[17]+0.28125*coeff[4]*fEdge[17]-0.28124999999999994*coeff[5]*fSkin[16]+0.28124999999999994*coeff[5]*fEdge[16]; 
  edgeSurf_incr[18] = 0.43301270189221935*coeff[7]*fSkin[20]-0.7187500000000001*coeff[5]*fSkin[20]-0.43301270189221935*coeff[7]*fEdge[20]-0.21875*coeff[5]*fEdge[20]+0.4330127018922193*coeff[6]*fSkin[18]-0.71875*coeff[4]*fSkin[18]-0.4330127018922193*coeff[6]*fEdge[18]-0.21875*coeff[4]*fEdge[18]-0.375*coeff[7]*fSkin[17]+0.7036456405748562*coeff[5]*fSkin[17]-0.375*coeff[7]*fEdge[17]-0.27063293868263694*coeff[5]*fEdge[17]-0.375*coeff[6]*fSkin[16]+0.7036456405748563*coeff[4]*fSkin[16]-0.375*coeff[6]*fEdge[16]-0.2706329386826371*coeff[4]*fEdge[16]; 
  edgeSurf_incr[19] = 0.27063293868263705*coeff[5]*fSkin[23]+0.27063293868263705*coeff[5]*fEdge[23]+0.2706329386826371*coeff[4]*fSkin[22]+0.2706329386826371*coeff[4]*fEdge[22]-0.28124999999999994*coeff[5]*fSkin[21]+0.28124999999999994*coeff[5]*fEdge[21]-0.28125*coeff[4]*fSkin[19]+0.28125*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = 0.4330127018922193*coeff[6]*fSkin[20]-0.71875*coeff[4]*fSkin[20]-0.4330127018922193*coeff[6]*fEdge[20]-0.21875*coeff[4]*fEdge[20]+0.43301270189221935*coeff[7]*fSkin[18]-0.7187500000000001*coeff[5]*fSkin[18]-0.43301270189221935*coeff[7]*fEdge[18]-0.21875*coeff[5]*fEdge[18]-0.375*coeff[6]*fSkin[17]+0.7036456405748563*coeff[4]*fSkin[17]-0.375*coeff[6]*fEdge[17]-0.2706329386826371*coeff[4]*fEdge[17]-0.375*coeff[7]*fSkin[16]+0.7036456405748562*coeff[5]*fSkin[16]-0.375*coeff[7]*fEdge[16]-0.27063293868263694*coeff[5]*fEdge[16]; 
  edgeSurf_incr[21] = 0.2706329386826371*coeff[4]*fSkin[23]+0.2706329386826371*coeff[4]*fEdge[23]+0.27063293868263705*coeff[5]*fSkin[22]+0.27063293868263705*coeff[5]*fEdge[22]-0.28125*coeff[4]*fSkin[21]+0.28125*coeff[4]*fEdge[21]-0.28124999999999994*coeff[5]*fSkin[19]+0.28124999999999994*coeff[5]*fEdge[19]; 
  edgeSurf_incr[22] = 0.43301270189221935*coeff[7]*fSkin[23]-0.7187500000000001*coeff[5]*fSkin[23]-0.43301270189221935*coeff[7]*fEdge[23]-0.21875*coeff[5]*fEdge[23]+0.4330127018922193*coeff[6]*fSkin[22]-0.71875*coeff[4]*fSkin[22]-0.4330127018922193*coeff[6]*fEdge[22]-0.21875*coeff[4]*fEdge[22]-0.375*coeff[7]*fSkin[21]+0.7036456405748562*coeff[5]*fSkin[21]-0.375*coeff[7]*fEdge[21]-0.27063293868263694*coeff[5]*fEdge[21]-0.375*coeff[6]*fSkin[19]+0.7036456405748563*coeff[4]*fSkin[19]-0.375*coeff[6]*fEdge[19]-0.2706329386826371*coeff[4]*fEdge[19]; 
  edgeSurf_incr[23] = 0.4330127018922193*coeff[6]*fSkin[23]-0.71875*coeff[4]*fSkin[23]-0.4330127018922193*coeff[6]*fEdge[23]-0.21875*coeff[4]*fEdge[23]+0.43301270189221935*coeff[7]*fSkin[22]-0.7187500000000001*coeff[5]*fSkin[22]-0.43301270189221935*coeff[7]*fEdge[22]-0.21875*coeff[5]*fEdge[22]-0.375*coeff[6]*fSkin[21]+0.7036456405748563*coeff[4]*fSkin[21]-0.375*coeff[6]*fEdge[21]-0.2706329386826371*coeff[4]*fEdge[21]-0.375*coeff[7]*fSkin[19]+0.7036456405748562*coeff[5]*fSkin[19]-0.375*coeff[7]*fEdge[19]-0.27063293868263694*coeff[5]*fEdge[19]; 

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

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdx2Sq; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdx2Sq; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdx2Sq; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdx2Sq; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdx2Sq; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdx2Sq; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdx2Sq; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdx2Sq; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdx2Sq; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdx2Sq; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdx2Sq; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdx2Sq; 

  return 0.;
}

