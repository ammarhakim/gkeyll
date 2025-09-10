#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfy_3x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],4.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.6237976320958223*coeff[1]*qSkin[2]+1.6237976320958223*coeff[1]*qEdge[2]+0.9375*qSkin[0]*coeff[1]-0.9375*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 1.6237976320958223*coeff[1]*qSkin[6]+1.6237976320958223*coeff[1]*qEdge[6]+0.9375*coeff[1]*qSkin[1]-0.9375*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = 3.5625*coeff[1]*qSkin[2]+2.0625*coeff[1]*qEdge[2]+1.6237976320958223*qSkin[0]*coeff[1]-1.6237976320958223*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 1.6237976320958223*coeff[1]*qSkin[8]+1.6237976320958223*coeff[1]*qEdge[8]+0.9375*coeff[1]*qSkin[3]-0.9375*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = 1.6237976320958223*coeff[1]*qSkin[10]+1.6237976320958223*coeff[1]*qEdge[10]+0.9375*coeff[1]*qSkin[4]-0.9375*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = 1.6237976320958223*coeff[1]*qSkin[13]+1.6237976320958223*coeff[1]*qEdge[13]+0.9375*coeff[1]*qSkin[5]-0.9375*coeff[1]*qEdge[5]; 
  edgeSurf_incr[6] = 3.5625*coeff[1]*qSkin[6]+2.0625*coeff[1]*qEdge[6]+1.6237976320958223*coeff[1]*qSkin[1]-1.6237976320958223*coeff[1]*qEdge[1]; 
  edgeSurf_incr[7] = 1.6237976320958223*coeff[1]*qSkin[16]+1.6237976320958223*coeff[1]*qEdge[16]+0.9375*coeff[1]*qSkin[7]-0.9375*coeff[1]*qEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[1]*qSkin[8]+2.0625*coeff[1]*qEdge[8]+1.6237976320958223*coeff[1]*qSkin[3]-1.6237976320958223*coeff[1]*qEdge[3]; 
  edgeSurf_incr[9] = 1.6237976320958223*coeff[1]*qSkin[17]+1.6237976320958223*coeff[1]*qEdge[17]+0.9375*coeff[1]*qSkin[9]-0.9375*coeff[1]*qEdge[9]; 
  edgeSurf_incr[10] = 3.5625*coeff[1]*qSkin[10]+2.0625*coeff[1]*qEdge[10]+1.6237976320958223*coeff[1]*qSkin[4]-1.6237976320958223*coeff[1]*qEdge[4]; 
  edgeSurf_incr[11] = 1.6237976320958223*coeff[1]*qSkin[19]+1.6237976320958223*coeff[1]*qEdge[19]+0.9375*coeff[1]*qSkin[11]-0.9375*coeff[1]*qEdge[11]; 
  edgeSurf_incr[12] = 1.6237976320958223*coeff[1]*qSkin[20]+1.6237976320958223*coeff[1]*qEdge[20]+0.9375*coeff[1]*qSkin[12]-0.9375*coeff[1]*qEdge[12]; 
  edgeSurf_incr[13] = 3.5625*coeff[1]*qSkin[13]+2.0625*coeff[1]*qEdge[13]+1.6237976320958223*coeff[1]*qSkin[5]-1.6237976320958223*coeff[1]*qEdge[5]; 
  edgeSurf_incr[14] = 1.6237976320958223*coeff[1]*qSkin[22]+1.6237976320958223*coeff[1]*qEdge[22]+0.9375*coeff[1]*qSkin[14]-0.9375*coeff[1]*qEdge[14]; 
  edgeSurf_incr[15] = 1.6237976320958223*coeff[1]*qSkin[24]+1.6237976320958223*coeff[1]*qEdge[24]+0.9375*coeff[1]*qSkin[15]-0.9375*coeff[1]*qEdge[15]; 
  edgeSurf_incr[16] = 3.5625*coeff[1]*qSkin[16]+2.0625*coeff[1]*qEdge[16]+1.6237976320958223*coeff[1]*qSkin[7]-1.6237976320958223*coeff[1]*qEdge[7]; 
  edgeSurf_incr[17] = 3.5625*coeff[1]*qSkin[17]+2.0625*coeff[1]*qEdge[17]+1.6237976320958223*coeff[1]*qSkin[9]-1.6237976320958223*coeff[1]*qEdge[9]; 
  edgeSurf_incr[18] = 1.6237976320958223*coeff[1]*qSkin[26]+1.6237976320958223*coeff[1]*qEdge[26]+0.9375*coeff[1]*qSkin[18]-0.9375*coeff[1]*qEdge[18]; 
  edgeSurf_incr[19] = 3.5625*coeff[1]*qSkin[19]+2.0625*coeff[1]*qEdge[19]+1.6237976320958223*coeff[1]*qSkin[11]-1.6237976320958223*coeff[1]*qEdge[11]; 
  edgeSurf_incr[20] = 3.5625*coeff[1]*qSkin[20]+2.0625*coeff[1]*qEdge[20]+1.6237976320958223*coeff[1]*qSkin[12]-1.6237976320958223*coeff[1]*qEdge[12]; 
  edgeSurf_incr[21] = 1.6237976320958223*coeff[1]*qSkin[27]+1.6237976320958223*coeff[1]*qEdge[27]+0.9375*coeff[1]*qSkin[21]-0.9375*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = 3.5625*coeff[1]*qSkin[22]+2.0625*coeff[1]*qEdge[22]+1.6237976320958223*coeff[1]*qSkin[14]-1.6237976320958223*coeff[1]*qEdge[14]; 
  edgeSurf_incr[23] = 1.6237976320958223*coeff[1]*qSkin[28]+1.6237976320958223*coeff[1]*qEdge[28]+0.9375*coeff[1]*qSkin[23]-0.9375*coeff[1]*qEdge[23]; 
  edgeSurf_incr[24] = 3.5625*coeff[1]*qSkin[24]+2.0625*coeff[1]*qEdge[24]+1.6237976320958223*coeff[1]*qSkin[15]-1.6237976320958223*coeff[1]*qEdge[15]; 
  edgeSurf_incr[25] = 1.6237976320958223*coeff[1]*qSkin[30]+1.6237976320958223*coeff[1]*qEdge[30]+0.9375*coeff[1]*qSkin[25]-0.9375*coeff[1]*qEdge[25]; 
  edgeSurf_incr[26] = 3.5625*coeff[1]*qSkin[26]+2.0625*coeff[1]*qEdge[26]+1.6237976320958223*coeff[1]*qSkin[18]-1.6237976320958223*coeff[1]*qEdge[18]; 
  edgeSurf_incr[27] = 3.5625*coeff[1]*qSkin[27]+2.0625*coeff[1]*qEdge[27]+1.6237976320958223*coeff[1]*qSkin[21]-1.6237976320958223*coeff[1]*qEdge[21]; 
  edgeSurf_incr[28] = 3.5625*coeff[1]*qSkin[28]+2.0625*coeff[1]*qEdge[28]+1.6237976320958223*coeff[1]*qSkin[23]-1.6237976320958223*coeff[1]*qEdge[23]; 
  edgeSurf_incr[29] = 1.6237976320958223*coeff[1]*qSkin[31]+1.6237976320958223*coeff[1]*qEdge[31]+0.9375*coeff[1]*qSkin[29]-0.9375*coeff[1]*qEdge[29]; 
  edgeSurf_incr[30] = 3.5625*coeff[1]*qSkin[30]+2.0625*coeff[1]*qEdge[30]+1.6237976320958223*coeff[1]*qSkin[25]-1.6237976320958223*coeff[1]*qEdge[25]; 
  edgeSurf_incr[31] = 3.5625*coeff[1]*qSkin[31]+2.0625*coeff[1]*qEdge[31]+1.6237976320958223*coeff[1]*qSkin[29]-1.6237976320958223*coeff[1]*qEdge[29]; 
  edgeSurf_incr[32] = 1.6237976320958225*coeff[1]*qSkin[34]+1.6237976320958225*coeff[1]*qEdge[34]+0.9375*coeff[1]*qSkin[32]-0.9375*coeff[1]*qEdge[32]; 
  edgeSurf_incr[33] = 1.6237976320958225*coeff[1]*qSkin[37]+1.6237976320958225*coeff[1]*qEdge[37]+0.9375*coeff[1]*qSkin[33]-0.9375*coeff[1]*qEdge[33]; 
  edgeSurf_incr[34] = 3.5625*coeff[1]*qSkin[34]+2.0625*coeff[1]*qEdge[34]+1.6237976320958225*coeff[1]*qSkin[32]-1.6237976320958225*coeff[1]*qEdge[32]; 
  edgeSurf_incr[35] = 1.6237976320958225*coeff[1]*qSkin[39]+1.6237976320958225*coeff[1]*qEdge[39]+0.9375*coeff[1]*qSkin[35]-0.9375*coeff[1]*qEdge[35]; 
  edgeSurf_incr[36] = 1.6237976320958225*coeff[1]*qSkin[41]+1.6237976320958225*coeff[1]*qEdge[41]+0.9375*coeff[1]*qSkin[36]-0.9375*coeff[1]*qEdge[36]; 
  edgeSurf_incr[37] = 3.5625*coeff[1]*qSkin[37]+2.0625*coeff[1]*qEdge[37]+1.6237976320958225*coeff[1]*qSkin[33]-1.6237976320958225*coeff[1]*qEdge[33]; 
  edgeSurf_incr[38] = 1.6237976320958225*coeff[1]*qSkin[43]+1.6237976320958225*coeff[1]*qEdge[43]+0.9375*coeff[1]*qSkin[38]-0.9375*coeff[1]*qEdge[38]; 
  edgeSurf_incr[39] = 3.5625*coeff[1]*qSkin[39]+2.0625*coeff[1]*qEdge[39]+1.6237976320958225*coeff[1]*qSkin[35]-1.6237976320958225*coeff[1]*qEdge[35]; 
  edgeSurf_incr[40] = 1.6237976320958225*coeff[1]*qSkin[44]+1.6237976320958225*coeff[1]*qEdge[44]+0.9375*coeff[1]*qSkin[40]-0.9375*coeff[1]*qEdge[40]; 
  edgeSurf_incr[41] = 3.5625*coeff[1]*qSkin[41]+2.0625*coeff[1]*qEdge[41]+1.6237976320958225*coeff[1]*qSkin[36]-1.6237976320958225*coeff[1]*qEdge[36]; 
  edgeSurf_incr[42] = 1.6237976320958225*coeff[1]*qSkin[46]+1.6237976320958225*coeff[1]*qEdge[46]+0.9375*coeff[1]*qSkin[42]-0.9375*coeff[1]*qEdge[42]; 
  edgeSurf_incr[43] = 3.5625*coeff[1]*qSkin[43]+2.0625*coeff[1]*qEdge[43]+1.6237976320958225*coeff[1]*qSkin[38]-1.6237976320958225*coeff[1]*qEdge[38]; 
  edgeSurf_incr[44] = 3.5625*coeff[1]*qSkin[44]+2.0625*coeff[1]*qEdge[44]+1.6237976320958225*coeff[1]*qSkin[40]-1.6237976320958225*coeff[1]*qEdge[40]; 
  edgeSurf_incr[45] = 1.6237976320958225*coeff[1]*qSkin[47]+1.6237976320958225*coeff[1]*qEdge[47]+0.9375*coeff[1]*qSkin[45]-0.9375*coeff[1]*qEdge[45]; 
  edgeSurf_incr[46] = 3.5625*coeff[1]*qSkin[46]+2.0625*coeff[1]*qEdge[46]+1.6237976320958225*coeff[1]*qSkin[42]-1.6237976320958225*coeff[1]*qEdge[42]; 
  edgeSurf_incr[47] = 3.5625*coeff[1]*qSkin[47]+2.0625*coeff[1]*qEdge[47]+1.6237976320958225*coeff[1]*qSkin[45]-1.6237976320958225*coeff[1]*qEdge[45]; 

  boundSurf_incr[2] = 1.5*coeff[1]*qSkin[2]; 
  boundSurf_incr[6] = 1.5*coeff[1]*qSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[1]*qSkin[8]; 
  boundSurf_incr[10] = 1.5*coeff[1]*qSkin[10]; 
  boundSurf_incr[13] = 1.5*coeff[1]*qSkin[13]; 
  boundSurf_incr[16] = 1.5*coeff[1]*qSkin[16]; 
  boundSurf_incr[17] = 1.5*coeff[1]*qSkin[17]; 
  boundSurf_incr[19] = 1.5*coeff[1]*qSkin[19]; 
  boundSurf_incr[20] = 1.5*coeff[1]*qSkin[20]; 
  boundSurf_incr[22] = 1.5*coeff[1]*qSkin[22]; 
  boundSurf_incr[24] = 1.5*coeff[1]*qSkin[24]; 
  boundSurf_incr[26] = 1.5*coeff[1]*qSkin[26]; 
  boundSurf_incr[27] = 1.5*coeff[1]*qSkin[27]; 
  boundSurf_incr[28] = 1.5*coeff[1]*qSkin[28]; 
  boundSurf_incr[30] = 1.5*coeff[1]*qSkin[30]; 
  boundSurf_incr[31] = 1.5*coeff[1]*qSkin[31]; 
  boundSurf_incr[34] = 1.5*coeff[1]*qSkin[34]; 
  boundSurf_incr[37] = 1.5*coeff[1]*qSkin[37]; 
  boundSurf_incr[39] = 1.5*coeff[1]*qSkin[39]; 
  boundSurf_incr[41] = 1.5*coeff[1]*qSkin[41]; 
  boundSurf_incr[43] = 1.5*coeff[1]*qSkin[43]; 
  boundSurf_incr[44] = 1.5*coeff[1]*qSkin[44]; 
  boundSurf_incr[46] = 1.5*coeff[1]*qSkin[46]; 
  boundSurf_incr[47] = 1.5*coeff[1]*qSkin[47]; 

  } else { 

  edgeSurf_incr[0] = -(1.6237976320958223*coeff[1]*qSkin[2])-1.6237976320958223*coeff[1]*qEdge[2]+0.9375*qSkin[0]*coeff[1]-0.9375*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(1.6237976320958223*coeff[1]*qSkin[6])-1.6237976320958223*coeff[1]*qEdge[6]+0.9375*coeff[1]*qSkin[1]-0.9375*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = 3.5625*coeff[1]*qSkin[2]+2.0625*coeff[1]*qEdge[2]-1.6237976320958223*qSkin[0]*coeff[1]+1.6237976320958223*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(1.6237976320958223*coeff[1]*qSkin[8])-1.6237976320958223*coeff[1]*qEdge[8]+0.9375*coeff[1]*qSkin[3]-0.9375*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = -(1.6237976320958223*coeff[1]*qSkin[10])-1.6237976320958223*coeff[1]*qEdge[10]+0.9375*coeff[1]*qSkin[4]-0.9375*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = -(1.6237976320958223*coeff[1]*qSkin[13])-1.6237976320958223*coeff[1]*qEdge[13]+0.9375*coeff[1]*qSkin[5]-0.9375*coeff[1]*qEdge[5]; 
  edgeSurf_incr[6] = 3.5625*coeff[1]*qSkin[6]+2.0625*coeff[1]*qEdge[6]-1.6237976320958223*coeff[1]*qSkin[1]+1.6237976320958223*coeff[1]*qEdge[1]; 
  edgeSurf_incr[7] = -(1.6237976320958223*coeff[1]*qSkin[16])-1.6237976320958223*coeff[1]*qEdge[16]+0.9375*coeff[1]*qSkin[7]-0.9375*coeff[1]*qEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[1]*qSkin[8]+2.0625*coeff[1]*qEdge[8]-1.6237976320958223*coeff[1]*qSkin[3]+1.6237976320958223*coeff[1]*qEdge[3]; 
  edgeSurf_incr[9] = -(1.6237976320958223*coeff[1]*qSkin[17])-1.6237976320958223*coeff[1]*qEdge[17]+0.9375*coeff[1]*qSkin[9]-0.9375*coeff[1]*qEdge[9]; 
  edgeSurf_incr[10] = 3.5625*coeff[1]*qSkin[10]+2.0625*coeff[1]*qEdge[10]-1.6237976320958223*coeff[1]*qSkin[4]+1.6237976320958223*coeff[1]*qEdge[4]; 
  edgeSurf_incr[11] = -(1.6237976320958223*coeff[1]*qSkin[19])-1.6237976320958223*coeff[1]*qEdge[19]+0.9375*coeff[1]*qSkin[11]-0.9375*coeff[1]*qEdge[11]; 
  edgeSurf_incr[12] = -(1.6237976320958223*coeff[1]*qSkin[20])-1.6237976320958223*coeff[1]*qEdge[20]+0.9375*coeff[1]*qSkin[12]-0.9375*coeff[1]*qEdge[12]; 
  edgeSurf_incr[13] = 3.5625*coeff[1]*qSkin[13]+2.0625*coeff[1]*qEdge[13]-1.6237976320958223*coeff[1]*qSkin[5]+1.6237976320958223*coeff[1]*qEdge[5]; 
  edgeSurf_incr[14] = -(1.6237976320958223*coeff[1]*qSkin[22])-1.6237976320958223*coeff[1]*qEdge[22]+0.9375*coeff[1]*qSkin[14]-0.9375*coeff[1]*qEdge[14]; 
  edgeSurf_incr[15] = -(1.6237976320958223*coeff[1]*qSkin[24])-1.6237976320958223*coeff[1]*qEdge[24]+0.9375*coeff[1]*qSkin[15]-0.9375*coeff[1]*qEdge[15]; 
  edgeSurf_incr[16] = 3.5625*coeff[1]*qSkin[16]+2.0625*coeff[1]*qEdge[16]-1.6237976320958223*coeff[1]*qSkin[7]+1.6237976320958223*coeff[1]*qEdge[7]; 
  edgeSurf_incr[17] = 3.5625*coeff[1]*qSkin[17]+2.0625*coeff[1]*qEdge[17]-1.6237976320958223*coeff[1]*qSkin[9]+1.6237976320958223*coeff[1]*qEdge[9]; 
  edgeSurf_incr[18] = -(1.6237976320958223*coeff[1]*qSkin[26])-1.6237976320958223*coeff[1]*qEdge[26]+0.9375*coeff[1]*qSkin[18]-0.9375*coeff[1]*qEdge[18]; 
  edgeSurf_incr[19] = 3.5625*coeff[1]*qSkin[19]+2.0625*coeff[1]*qEdge[19]-1.6237976320958223*coeff[1]*qSkin[11]+1.6237976320958223*coeff[1]*qEdge[11]; 
  edgeSurf_incr[20] = 3.5625*coeff[1]*qSkin[20]+2.0625*coeff[1]*qEdge[20]-1.6237976320958223*coeff[1]*qSkin[12]+1.6237976320958223*coeff[1]*qEdge[12]; 
  edgeSurf_incr[21] = -(1.6237976320958223*coeff[1]*qSkin[27])-1.6237976320958223*coeff[1]*qEdge[27]+0.9375*coeff[1]*qSkin[21]-0.9375*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = 3.5625*coeff[1]*qSkin[22]+2.0625*coeff[1]*qEdge[22]-1.6237976320958223*coeff[1]*qSkin[14]+1.6237976320958223*coeff[1]*qEdge[14]; 
  edgeSurf_incr[23] = -(1.6237976320958223*coeff[1]*qSkin[28])-1.6237976320958223*coeff[1]*qEdge[28]+0.9375*coeff[1]*qSkin[23]-0.9375*coeff[1]*qEdge[23]; 
  edgeSurf_incr[24] = 3.5625*coeff[1]*qSkin[24]+2.0625*coeff[1]*qEdge[24]-1.6237976320958223*coeff[1]*qSkin[15]+1.6237976320958223*coeff[1]*qEdge[15]; 
  edgeSurf_incr[25] = -(1.6237976320958223*coeff[1]*qSkin[30])-1.6237976320958223*coeff[1]*qEdge[30]+0.9375*coeff[1]*qSkin[25]-0.9375*coeff[1]*qEdge[25]; 
  edgeSurf_incr[26] = 3.5625*coeff[1]*qSkin[26]+2.0625*coeff[1]*qEdge[26]-1.6237976320958223*coeff[1]*qSkin[18]+1.6237976320958223*coeff[1]*qEdge[18]; 
  edgeSurf_incr[27] = 3.5625*coeff[1]*qSkin[27]+2.0625*coeff[1]*qEdge[27]-1.6237976320958223*coeff[1]*qSkin[21]+1.6237976320958223*coeff[1]*qEdge[21]; 
  edgeSurf_incr[28] = 3.5625*coeff[1]*qSkin[28]+2.0625*coeff[1]*qEdge[28]-1.6237976320958223*coeff[1]*qSkin[23]+1.6237976320958223*coeff[1]*qEdge[23]; 
  edgeSurf_incr[29] = -(1.6237976320958223*coeff[1]*qSkin[31])-1.6237976320958223*coeff[1]*qEdge[31]+0.9375*coeff[1]*qSkin[29]-0.9375*coeff[1]*qEdge[29]; 
  edgeSurf_incr[30] = 3.5625*coeff[1]*qSkin[30]+2.0625*coeff[1]*qEdge[30]-1.6237976320958223*coeff[1]*qSkin[25]+1.6237976320958223*coeff[1]*qEdge[25]; 
  edgeSurf_incr[31] = 3.5625*coeff[1]*qSkin[31]+2.0625*coeff[1]*qEdge[31]-1.6237976320958223*coeff[1]*qSkin[29]+1.6237976320958223*coeff[1]*qEdge[29]; 
  edgeSurf_incr[32] = -(1.6237976320958225*coeff[1]*qSkin[34])-1.6237976320958225*coeff[1]*qEdge[34]+0.9375*coeff[1]*qSkin[32]-0.9375*coeff[1]*qEdge[32]; 
  edgeSurf_incr[33] = -(1.6237976320958225*coeff[1]*qSkin[37])-1.6237976320958225*coeff[1]*qEdge[37]+0.9375*coeff[1]*qSkin[33]-0.9375*coeff[1]*qEdge[33]; 
  edgeSurf_incr[34] = 3.5625*coeff[1]*qSkin[34]+2.0625*coeff[1]*qEdge[34]-1.6237976320958225*coeff[1]*qSkin[32]+1.6237976320958225*coeff[1]*qEdge[32]; 
  edgeSurf_incr[35] = -(1.6237976320958225*coeff[1]*qSkin[39])-1.6237976320958225*coeff[1]*qEdge[39]+0.9375*coeff[1]*qSkin[35]-0.9375*coeff[1]*qEdge[35]; 
  edgeSurf_incr[36] = -(1.6237976320958225*coeff[1]*qSkin[41])-1.6237976320958225*coeff[1]*qEdge[41]+0.9375*coeff[1]*qSkin[36]-0.9375*coeff[1]*qEdge[36]; 
  edgeSurf_incr[37] = 3.5625*coeff[1]*qSkin[37]+2.0625*coeff[1]*qEdge[37]-1.6237976320958225*coeff[1]*qSkin[33]+1.6237976320958225*coeff[1]*qEdge[33]; 
  edgeSurf_incr[38] = -(1.6237976320958225*coeff[1]*qSkin[43])-1.6237976320958225*coeff[1]*qEdge[43]+0.9375*coeff[1]*qSkin[38]-0.9375*coeff[1]*qEdge[38]; 
  edgeSurf_incr[39] = 3.5625*coeff[1]*qSkin[39]+2.0625*coeff[1]*qEdge[39]-1.6237976320958225*coeff[1]*qSkin[35]+1.6237976320958225*coeff[1]*qEdge[35]; 
  edgeSurf_incr[40] = -(1.6237976320958225*coeff[1]*qSkin[44])-1.6237976320958225*coeff[1]*qEdge[44]+0.9375*coeff[1]*qSkin[40]-0.9375*coeff[1]*qEdge[40]; 
  edgeSurf_incr[41] = 3.5625*coeff[1]*qSkin[41]+2.0625*coeff[1]*qEdge[41]-1.6237976320958225*coeff[1]*qSkin[36]+1.6237976320958225*coeff[1]*qEdge[36]; 
  edgeSurf_incr[42] = -(1.6237976320958225*coeff[1]*qSkin[46])-1.6237976320958225*coeff[1]*qEdge[46]+0.9375*coeff[1]*qSkin[42]-0.9375*coeff[1]*qEdge[42]; 
  edgeSurf_incr[43] = 3.5625*coeff[1]*qSkin[43]+2.0625*coeff[1]*qEdge[43]-1.6237976320958225*coeff[1]*qSkin[38]+1.6237976320958225*coeff[1]*qEdge[38]; 
  edgeSurf_incr[44] = 3.5625*coeff[1]*qSkin[44]+2.0625*coeff[1]*qEdge[44]-1.6237976320958225*coeff[1]*qSkin[40]+1.6237976320958225*coeff[1]*qEdge[40]; 
  edgeSurf_incr[45] = -(1.6237976320958225*coeff[1]*qSkin[47])-1.6237976320958225*coeff[1]*qEdge[47]+0.9375*coeff[1]*qSkin[45]-0.9375*coeff[1]*qEdge[45]; 
  edgeSurf_incr[46] = 3.5625*coeff[1]*qSkin[46]+2.0625*coeff[1]*qEdge[46]-1.6237976320958225*coeff[1]*qSkin[42]+1.6237976320958225*coeff[1]*qEdge[42]; 
  edgeSurf_incr[47] = 3.5625*coeff[1]*qSkin[47]+2.0625*coeff[1]*qEdge[47]-1.6237976320958225*coeff[1]*qSkin[45]+1.6237976320958225*coeff[1]*qEdge[45]; 

  boundSurf_incr[2] = 1.5*coeff[1]*qSkin[2]; 
  boundSurf_incr[6] = 1.5*coeff[1]*qSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[1]*qSkin[8]; 
  boundSurf_incr[10] = 1.5*coeff[1]*qSkin[10]; 
  boundSurf_incr[13] = 1.5*coeff[1]*qSkin[13]; 
  boundSurf_incr[16] = 1.5*coeff[1]*qSkin[16]; 
  boundSurf_incr[17] = 1.5*coeff[1]*qSkin[17]; 
  boundSurf_incr[19] = 1.5*coeff[1]*qSkin[19]; 
  boundSurf_incr[20] = 1.5*coeff[1]*qSkin[20]; 
  boundSurf_incr[22] = 1.5*coeff[1]*qSkin[22]; 
  boundSurf_incr[24] = 1.5*coeff[1]*qSkin[24]; 
  boundSurf_incr[26] = 1.5*coeff[1]*qSkin[26]; 
  boundSurf_incr[27] = 1.5*coeff[1]*qSkin[27]; 
  boundSurf_incr[28] = 1.5*coeff[1]*qSkin[28]; 
  boundSurf_incr[30] = 1.5*coeff[1]*qSkin[30]; 
  boundSurf_incr[31] = 1.5*coeff[1]*qSkin[31]; 
  boundSurf_incr[34] = 1.5*coeff[1]*qSkin[34]; 
  boundSurf_incr[37] = 1.5*coeff[1]*qSkin[37]; 
  boundSurf_incr[39] = 1.5*coeff[1]*qSkin[39]; 
  boundSurf_incr[41] = 1.5*coeff[1]*qSkin[41]; 
  boundSurf_incr[43] = 1.5*coeff[1]*qSkin[43]; 
  boundSurf_incr[44] = 1.5*coeff[1]*qSkin[44]; 
  boundSurf_incr[46] = 1.5*coeff[1]*qSkin[46]; 
  boundSurf_incr[47] = 1.5*coeff[1]*qSkin[47]; 

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
  out[24] += -(1.0*(vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*rdx2Sq); 
  out[25] += -(1.0*(vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*rdx2Sq); 
  out[26] += -(1.0*(vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*rdx2Sq); 
  out[27] += -(1.0*(vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*rdx2Sq); 
  out[28] += -(1.0*(vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*rdx2Sq); 
  out[29] += -(1.0*(vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*rdx2Sq); 
  out[30] += -(1.0*(vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*rdx2Sq); 
  out[31] += -(1.0*(vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*rdx2Sq); 
  out[32] += -(1.0*(vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*rdx2Sq); 
  out[33] += -(1.0*(vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*rdx2Sq); 
  out[34] += -(1.0*(vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*rdx2Sq); 
  out[35] += -(1.0*(vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*rdx2Sq); 
  out[36] += -(1.0*(vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*rdx2Sq); 
  out[37] += -(1.0*(vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*rdx2Sq); 
  out[38] += -(1.0*(vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*rdx2Sq); 
  out[39] += -(1.0*(vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*rdx2Sq); 
  out[40] += -(1.0*(vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*rdx2Sq); 
  out[41] += -(1.0*(vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*rdx2Sq); 
  out[42] += -(1.0*(vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*rdx2Sq); 
  out[43] += -(1.0*(vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*rdx2Sq); 
  out[44] += -(1.0*(vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*rdx2Sq); 
  out[45] += -(1.0*(vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*rdx2Sq); 
  out[46] += -(1.0*(vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*rdx2Sq); 
  out[47] += -(1.0*(vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*rdx2Sq); 

  return 0.;
}

