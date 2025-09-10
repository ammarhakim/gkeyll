#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double vol_incr[24] = {0.0}; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.6237976320958223*coeff[0]*qSkin[1]+1.6237976320958223*coeff[0]*qEdge[1]+0.9375*coeff[0]*qSkin[0]-0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*qSkin[1]+2.0625*coeff[0]*qEdge[1]+1.6237976320958223*coeff[0]*qSkin[0]-1.6237976320958223*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 1.6237976320958223*coeff[0]*qSkin[5]+1.6237976320958223*coeff[0]*qEdge[5]+0.9375*coeff[0]*qSkin[2]-0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 1.6237976320958223*coeff[0]*qSkin[6]+1.6237976320958223*coeff[0]*qEdge[6]+0.9375*coeff[0]*qSkin[3]-0.9375*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = 1.6237976320958223*coeff[0]*qSkin[8]+1.6237976320958223*coeff[0]*qEdge[8]+0.9375*coeff[0]*qSkin[4]-0.9375*coeff[0]*qEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*qSkin[5]+2.0625*coeff[0]*qEdge[5]+1.6237976320958223*coeff[0]*qSkin[2]-1.6237976320958223*coeff[0]*qEdge[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[0]*qSkin[6]+2.0625*coeff[0]*qEdge[6]+1.6237976320958223*coeff[0]*qSkin[3]-1.6237976320958223*coeff[0]*qEdge[3]; 
  edgeSurf_incr[7] = 1.6237976320958223*coeff[0]*qSkin[11]+1.6237976320958223*coeff[0]*qEdge[11]+0.9375*coeff[0]*qSkin[7]-0.9375*coeff[0]*qEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[0]*qSkin[8]+2.0625*coeff[0]*qEdge[8]+1.6237976320958223*coeff[0]*qSkin[4]-1.6237976320958223*coeff[0]*qEdge[4]; 
  edgeSurf_incr[9] = 1.6237976320958223*coeff[0]*qSkin[12]+1.6237976320958223*coeff[0]*qEdge[12]+0.9375*coeff[0]*qSkin[9]-0.9375*coeff[0]*qEdge[9]; 
  edgeSurf_incr[10] = 1.6237976320958223*coeff[0]*qSkin[13]+1.6237976320958223*coeff[0]*qEdge[13]+0.9375*coeff[0]*qSkin[10]-0.9375*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*qSkin[11]+2.0625*coeff[0]*qEdge[11]+1.6237976320958223*coeff[0]*qSkin[7]-1.6237976320958223*coeff[0]*qEdge[7]; 
  edgeSurf_incr[12] = 3.5625*coeff[0]*qSkin[12]+2.0625*coeff[0]*qEdge[12]+1.6237976320958223*coeff[0]*qSkin[9]-1.6237976320958223*coeff[0]*qEdge[9]; 
  edgeSurf_incr[13] = 3.5625*coeff[0]*qSkin[13]+2.0625*coeff[0]*qEdge[13]+1.6237976320958223*coeff[0]*qSkin[10]-1.6237976320958223*coeff[0]*qEdge[10]; 
  edgeSurf_incr[14] = 1.6237976320958223*coeff[0]*qSkin[15]+1.6237976320958223*coeff[0]*qEdge[15]+0.9375*coeff[0]*qSkin[14]-0.9375*coeff[0]*qEdge[14]; 
  edgeSurf_incr[15] = 3.5625*coeff[0]*qSkin[15]+2.0625*coeff[0]*qEdge[15]+1.6237976320958223*coeff[0]*qSkin[14]-1.6237976320958223*coeff[0]*qEdge[14]; 
  edgeSurf_incr[16] = 1.6237976320958225*coeff[0]*qSkin[17]+1.6237976320958225*coeff[0]*qEdge[17]+0.9375*coeff[0]*qSkin[16]-0.9375*coeff[0]*qEdge[16]; 
  edgeSurf_incr[17] = 3.5625*coeff[0]*qSkin[17]+2.0625*coeff[0]*qEdge[17]+1.6237976320958225*coeff[0]*qSkin[16]-1.6237976320958225*coeff[0]*qEdge[16]; 
  edgeSurf_incr[18] = 1.6237976320958225*coeff[0]*qSkin[20]+1.6237976320958225*coeff[0]*qEdge[20]+0.9375*coeff[0]*qSkin[18]-0.9375*coeff[0]*qEdge[18]; 
  edgeSurf_incr[19] = 1.6237976320958225*coeff[0]*qSkin[21]+1.6237976320958225*coeff[0]*qEdge[21]+0.9375*coeff[0]*qSkin[19]-0.9375*coeff[0]*qEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[0]*qSkin[20]+2.0625*coeff[0]*qEdge[20]+1.6237976320958225*coeff[0]*qSkin[18]-1.6237976320958225*coeff[0]*qEdge[18]; 
  edgeSurf_incr[21] = 3.5625*coeff[0]*qSkin[21]+2.0625*coeff[0]*qEdge[21]+1.6237976320958225*coeff[0]*qSkin[19]-1.6237976320958225*coeff[0]*qEdge[19]; 
  edgeSurf_incr[22] = 1.6237976320958225*coeff[0]*qSkin[23]+1.6237976320958225*coeff[0]*qEdge[23]+0.9375*coeff[0]*qSkin[22]-0.9375*coeff[0]*qEdge[22]; 
  edgeSurf_incr[23] = 3.5625*coeff[0]*qSkin[23]+2.0625*coeff[0]*qEdge[23]+1.6237976320958225*coeff[0]*qSkin[22]-1.6237976320958225*coeff[0]*qEdge[22]; 

  boundSurf_incr[1] = 1.5*coeff[0]*qSkin[1]; 
  boundSurf_incr[5] = 1.5*coeff[0]*qSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[0]*qSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[0]*qSkin[8]; 
  boundSurf_incr[11] = 1.5*coeff[0]*qSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[0]*qSkin[12]; 
  boundSurf_incr[13] = 1.5*coeff[0]*qSkin[13]; 
  boundSurf_incr[15] = 1.5*coeff[0]*qSkin[15]; 
  boundSurf_incr[17] = 1.5*coeff[0]*qSkin[17]; 
  boundSurf_incr[20] = 1.5*coeff[0]*qSkin[20]; 
  boundSurf_incr[21] = 1.5*coeff[0]*qSkin[21]; 
  boundSurf_incr[23] = 1.5*coeff[0]*qSkin[23]; 

  } else { 

  edgeSurf_incr[0] = -(1.6237976320958223*coeff[0]*qSkin[1])-1.6237976320958223*coeff[0]*qEdge[1]+0.9375*coeff[0]*qSkin[0]-0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*qSkin[1]+2.0625*coeff[0]*qEdge[1]-1.6237976320958223*coeff[0]*qSkin[0]+1.6237976320958223*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = -(1.6237976320958223*coeff[0]*qSkin[5])-1.6237976320958223*coeff[0]*qEdge[5]+0.9375*coeff[0]*qSkin[2]-0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = -(1.6237976320958223*coeff[0]*qSkin[6])-1.6237976320958223*coeff[0]*qEdge[6]+0.9375*coeff[0]*qSkin[3]-0.9375*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = -(1.6237976320958223*coeff[0]*qSkin[8])-1.6237976320958223*coeff[0]*qEdge[8]+0.9375*coeff[0]*qSkin[4]-0.9375*coeff[0]*qEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*qSkin[5]+2.0625*coeff[0]*qEdge[5]-1.6237976320958223*coeff[0]*qSkin[2]+1.6237976320958223*coeff[0]*qEdge[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[0]*qSkin[6]+2.0625*coeff[0]*qEdge[6]-1.6237976320958223*coeff[0]*qSkin[3]+1.6237976320958223*coeff[0]*qEdge[3]; 
  edgeSurf_incr[7] = -(1.6237976320958223*coeff[0]*qSkin[11])-1.6237976320958223*coeff[0]*qEdge[11]+0.9375*coeff[0]*qSkin[7]-0.9375*coeff[0]*qEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[0]*qSkin[8]+2.0625*coeff[0]*qEdge[8]-1.6237976320958223*coeff[0]*qSkin[4]+1.6237976320958223*coeff[0]*qEdge[4]; 
  edgeSurf_incr[9] = -(1.6237976320958223*coeff[0]*qSkin[12])-1.6237976320958223*coeff[0]*qEdge[12]+0.9375*coeff[0]*qSkin[9]-0.9375*coeff[0]*qEdge[9]; 
  edgeSurf_incr[10] = -(1.6237976320958223*coeff[0]*qSkin[13])-1.6237976320958223*coeff[0]*qEdge[13]+0.9375*coeff[0]*qSkin[10]-0.9375*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*qSkin[11]+2.0625*coeff[0]*qEdge[11]-1.6237976320958223*coeff[0]*qSkin[7]+1.6237976320958223*coeff[0]*qEdge[7]; 
  edgeSurf_incr[12] = 3.5625*coeff[0]*qSkin[12]+2.0625*coeff[0]*qEdge[12]-1.6237976320958223*coeff[0]*qSkin[9]+1.6237976320958223*coeff[0]*qEdge[9]; 
  edgeSurf_incr[13] = 3.5625*coeff[0]*qSkin[13]+2.0625*coeff[0]*qEdge[13]-1.6237976320958223*coeff[0]*qSkin[10]+1.6237976320958223*coeff[0]*qEdge[10]; 
  edgeSurf_incr[14] = -(1.6237976320958223*coeff[0]*qSkin[15])-1.6237976320958223*coeff[0]*qEdge[15]+0.9375*coeff[0]*qSkin[14]-0.9375*coeff[0]*qEdge[14]; 
  edgeSurf_incr[15] = 3.5625*coeff[0]*qSkin[15]+2.0625*coeff[0]*qEdge[15]-1.6237976320958223*coeff[0]*qSkin[14]+1.6237976320958223*coeff[0]*qEdge[14]; 
  edgeSurf_incr[16] = -(1.6237976320958225*coeff[0]*qSkin[17])-1.6237976320958225*coeff[0]*qEdge[17]+0.9375*coeff[0]*qSkin[16]-0.9375*coeff[0]*qEdge[16]; 
  edgeSurf_incr[17] = 3.5625*coeff[0]*qSkin[17]+2.0625*coeff[0]*qEdge[17]-1.6237976320958225*coeff[0]*qSkin[16]+1.6237976320958225*coeff[0]*qEdge[16]; 
  edgeSurf_incr[18] = -(1.6237976320958225*coeff[0]*qSkin[20])-1.6237976320958225*coeff[0]*qEdge[20]+0.9375*coeff[0]*qSkin[18]-0.9375*coeff[0]*qEdge[18]; 
  edgeSurf_incr[19] = -(1.6237976320958225*coeff[0]*qSkin[21])-1.6237976320958225*coeff[0]*qEdge[21]+0.9375*coeff[0]*qSkin[19]-0.9375*coeff[0]*qEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[0]*qSkin[20]+2.0625*coeff[0]*qEdge[20]-1.6237976320958225*coeff[0]*qSkin[18]+1.6237976320958225*coeff[0]*qEdge[18]; 
  edgeSurf_incr[21] = 3.5625*coeff[0]*qSkin[21]+2.0625*coeff[0]*qEdge[21]-1.6237976320958225*coeff[0]*qSkin[19]+1.6237976320958225*coeff[0]*qEdge[19]; 
  edgeSurf_incr[22] = -(1.6237976320958225*coeff[0]*qSkin[23])-1.6237976320958225*coeff[0]*qEdge[23]+0.9375*coeff[0]*qSkin[22]-0.9375*coeff[0]*qEdge[22]; 
  edgeSurf_incr[23] = 3.5625*coeff[0]*qSkin[23]+2.0625*coeff[0]*qEdge[23]-1.6237976320958225*coeff[0]*qSkin[22]+1.6237976320958225*coeff[0]*qEdge[22]; 

  boundSurf_incr[1] = 1.5*coeff[0]*qSkin[1]; 
  boundSurf_incr[5] = 1.5*coeff[0]*qSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[0]*qSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[0]*qSkin[8]; 
  boundSurf_incr[11] = 1.5*coeff[0]*qSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[0]*qSkin[12]; 
  boundSurf_incr[13] = 1.5*coeff[0]*qSkin[13]; 
  boundSurf_incr[15] = 1.5*coeff[0]*qSkin[15]; 
  boundSurf_incr[17] = 1.5*coeff[0]*qSkin[17]; 
  boundSurf_incr[20] = 1.5*coeff[0]*qSkin[20]; 
  boundSurf_incr[21] = 1.5*coeff[0]*qSkin[21]; 
  boundSurf_incr[23] = 1.5*coeff[0]*qSkin[23]; 

  }

  out[0] += -(1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq); 
  out[1] += -(1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq); 
  out[2] += -(1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq); 
  out[3] += -(1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq); 
  out[4] += -(1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq); 
  out[5] += -(1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq); 
  out[6] += -(1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq); 
  out[7] += -(1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq); 
  out[8] += -(1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq); 
  out[9] += -(1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq); 
  out[10] += -(1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq); 
  out[11] += -(1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq); 
  out[12] += -(1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdx2Sq); 
  out[13] += -(1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdx2Sq); 
  out[14] += -(1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdx2Sq); 
  out[15] += -(1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdx2Sq); 
  out[16] += -(1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdx2Sq); 
  out[17] += -(1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdx2Sq); 
  out[18] += -(1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdx2Sq); 
  out[19] += -(1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdx2Sq); 
  out[20] += -(1.0*(vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdx2Sq); 
  out[21] += -(1.0*(vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdx2Sq); 
  out[22] += -(1.0*(vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdx2Sq); 
  out[23] += -(1.0*(vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdx2Sq); 

  return 0.;
}

