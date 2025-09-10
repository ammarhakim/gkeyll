#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfy_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
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

  edgeSurf_incr[0] = 6.708203932499369*coeff[1]*qSkin[12]-6.708203932499369*coeff[1]*qEdge[12]+8.11898816047911*coeff[1]*qSkin[2]+8.11898816047911*coeff[1]*qEdge[2]+4.6875*qSkin[0]*coeff[1]-4.6875*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[1]*qSkin[20]-6.7082039324993685*coeff[1]*qEdge[20]+8.11898816047911*coeff[1]*qSkin[5]+8.11898816047911*coeff[1]*qEdge[5]+4.6875*coeff[1]*qSkin[1]-4.6875*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = 14.160595359570864*coeff[1]*qSkin[12]-9.077304717673634*coeff[1]*qEdge[12]+15.46875*coeff[1]*qSkin[2]+12.65625*coeff[1]*qEdge[2]+8.11898816047911*qSkin[0]*coeff[1]-8.11898816047911*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[1]*qSkin[22]-6.7082039324993685*coeff[1]*qEdge[22]+8.11898816047911*coeff[1]*qSkin[7]+8.11898816047911*coeff[1]*qEdge[7]+4.6875*coeff[1]*qSkin[3]-4.6875*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = 6.7082039324993685*coeff[1]*qSkin[26]-6.7082039324993685*coeff[1]*qEdge[26]+8.11898816047911*coeff[1]*qSkin[9]+8.11898816047911*coeff[1]*qEdge[9]+4.6875*coeff[1]*qSkin[4]-4.6875*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = 14.160595359570868*coeff[1]*qSkin[20]-9.077304717673634*coeff[1]*qEdge[20]+15.46875*coeff[1]*qSkin[5]+12.65625*coeff[1]*qEdge[5]+8.11898816047911*coeff[1]*qSkin[1]-8.11898816047911*coeff[1]*qEdge[1]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[1]*qSkin[33]-6.708203932499369*coeff[1]*qEdge[33]+8.11898816047911*coeff[1]*qSkin[15]+8.11898816047911*coeff[1]*qEdge[15]+4.6875*coeff[1]*qSkin[6]-4.6875*coeff[1]*qEdge[6]; 
  edgeSurf_incr[7] = 14.160595359570868*coeff[1]*qSkin[22]-9.077304717673634*coeff[1]*qEdge[22]+15.46875*coeff[1]*qSkin[7]+12.65625*coeff[1]*qEdge[7]+8.11898816047911*coeff[1]*qSkin[3]-8.11898816047911*coeff[1]*qEdge[3]; 
  edgeSurf_incr[8] = 6.708203932499369*coeff[1]*qSkin[36]-6.708203932499369*coeff[1]*qEdge[36]+8.11898816047911*coeff[1]*qSkin[16]+8.11898816047911*coeff[1]*qEdge[16]+4.6875*coeff[1]*qSkin[8]-4.6875*coeff[1]*qEdge[8]; 
  edgeSurf_incr[9] = 14.160595359570868*coeff[1]*qSkin[26]-9.077304717673634*coeff[1]*qEdge[26]+15.46875*coeff[1]*qSkin[9]+12.65625*coeff[1]*qEdge[9]+8.11898816047911*coeff[1]*qSkin[4]-8.11898816047911*coeff[1]*qEdge[4]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[1]*qSkin[38]-6.708203932499369*coeff[1]*qEdge[38]+8.11898816047911*coeff[1]*qSkin[18]+8.11898816047911*coeff[1]*qEdge[18]+4.6875*coeff[1]*qSkin[10]-4.6875*coeff[1]*qEdge[10]; 
  edgeSurf_incr[11] = 8.118988160479114*coeff[1]*qSkin[19]+8.118988160479114*coeff[1]*qEdge[19]+4.6875*coeff[1]*qSkin[11]-4.6875*coeff[1]*qEdge[11]; 
  edgeSurf_incr[12] = 20.34375*coeff[1]*qSkin[12]-0.65625*coeff[1]*qEdge[12]+15.61296411439865*coeff[1]*qSkin[2]+4.720198453190292*coeff[1]*qEdge[2]+4.192627457812107*qSkin[0]*coeff[1]-4.192627457812107*qEdge[0]*coeff[1]; 
  edgeSurf_incr[13] = 8.118988160479114*coeff[1]*qSkin[24]+8.118988160479114*coeff[1]*qEdge[24]+4.6875*coeff[1]*qSkin[13]-4.6875*coeff[1]*qEdge[13]; 
  edgeSurf_incr[14] = 8.118988160479114*coeff[1]*qSkin[29]+8.118988160479114*coeff[1]*qEdge[29]+4.6875*coeff[1]*qSkin[14]-4.6875*coeff[1]*qEdge[14]; 
  edgeSurf_incr[15] = 14.160595359570864*coeff[1]*qSkin[33]-9.077304717673634*coeff[1]*qEdge[33]+15.46875*coeff[1]*qSkin[15]+12.65625*coeff[1]*qEdge[15]+8.11898816047911*coeff[1]*qSkin[6]-8.11898816047911*coeff[1]*qEdge[6]; 
  edgeSurf_incr[16] = 14.160595359570864*coeff[1]*qSkin[36]-9.077304717673634*coeff[1]*qEdge[36]+15.46875*coeff[1]*qSkin[16]+12.65625*coeff[1]*qEdge[16]+8.11898816047911*coeff[1]*qSkin[8]-8.11898816047911*coeff[1]*qEdge[8]; 
  edgeSurf_incr[17] = 6.7082039324993685*coeff[1]*qSkin[45]-6.7082039324993685*coeff[1]*qEdge[45]+8.11898816047911*coeff[1]*qSkin[31]+8.11898816047911*coeff[1]*qEdge[31]+4.6875*coeff[1]*qSkin[17]-4.6875*coeff[1]*qEdge[17]; 
  edgeSurf_incr[18] = 14.160595359570864*coeff[1]*qSkin[38]-9.077304717673634*coeff[1]*qEdge[38]+15.46875*coeff[1]*qSkin[18]+12.65625*coeff[1]*qEdge[18]+8.11898816047911*coeff[1]*qSkin[10]-8.11898816047911*coeff[1]*qEdge[10]; 
  edgeSurf_incr[19] = 15.46875*coeff[1]*qSkin[19]+12.65625*coeff[1]*qEdge[19]+8.118988160479114*coeff[1]*qSkin[11]-8.118988160479114*coeff[1]*qEdge[11]; 
  edgeSurf_incr[20] = 20.34375*coeff[1]*qSkin[20]-0.65625*coeff[1]*qEdge[20]+15.612964114398654*coeff[1]*qSkin[5]+4.72019845319029*coeff[1]*qEdge[5]+4.192627457812107*coeff[1]*qSkin[1]-4.192627457812107*coeff[1]*qEdge[1]; 
  edgeSurf_incr[21] = 8.118988160479114*coeff[1]*qSkin[32]+8.118988160479114*coeff[1]*qEdge[32]+4.6875*coeff[1]*qSkin[21]-4.6875*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = 20.34375*coeff[1]*qSkin[22]-0.65625*coeff[1]*qEdge[22]+15.612964114398654*coeff[1]*qSkin[7]+4.72019845319029*coeff[1]*qEdge[7]+4.192627457812107*coeff[1]*qSkin[3]-4.192627457812107*coeff[1]*qEdge[3]; 
  edgeSurf_incr[23] = 8.118988160479114*coeff[1]*qSkin[34]+8.118988160479114*coeff[1]*qEdge[34]+4.6875*coeff[1]*qSkin[23]-4.6875*coeff[1]*qEdge[23]; 
  edgeSurf_incr[24] = 15.46875*coeff[1]*qSkin[24]+12.65625*coeff[1]*qEdge[24]+8.118988160479114*coeff[1]*qSkin[13]-8.118988160479114*coeff[1]*qEdge[13]; 
  edgeSurf_incr[25] = 8.118988160479114*coeff[1]*qSkin[35]+8.118988160479114*coeff[1]*qEdge[35]+4.6875*coeff[1]*qSkin[25]-4.6875*coeff[1]*qEdge[25]; 
  edgeSurf_incr[26] = 20.34375*coeff[1]*qSkin[26]-0.65625*coeff[1]*qEdge[26]+15.612964114398654*coeff[1]*qSkin[9]+4.72019845319029*coeff[1]*qEdge[9]+4.192627457812107*coeff[1]*qSkin[4]-4.192627457812107*coeff[1]*qEdge[4]; 
  edgeSurf_incr[27] = 8.118988160479114*coeff[1]*qSkin[40]+8.118988160479114*coeff[1]*qEdge[40]+4.6875*coeff[1]*qSkin[27]-4.6875*coeff[1]*qEdge[27]; 
  edgeSurf_incr[28] = 8.118988160479114*coeff[1]*qSkin[41]+8.118988160479114*coeff[1]*qEdge[41]+4.6875*coeff[1]*qSkin[28]-4.6875*coeff[1]*qEdge[28]; 
  edgeSurf_incr[29] = 15.46875*coeff[1]*qSkin[29]+12.65625*coeff[1]*qEdge[29]+8.118988160479114*coeff[1]*qSkin[14]-8.118988160479114*coeff[1]*qEdge[14]; 
  edgeSurf_incr[30] = 8.118988160479114*coeff[1]*qSkin[43]+8.118988160479114*coeff[1]*qEdge[43]+4.6875*coeff[1]*qSkin[30]-4.6875*coeff[1]*qEdge[30]; 
  edgeSurf_incr[31] = 14.160595359570868*coeff[1]*qSkin[45]-9.077304717673634*coeff[1]*qEdge[45]+15.46875*coeff[1]*qSkin[31]+12.65625*coeff[1]*qEdge[31]+8.11898816047911*coeff[1]*qSkin[17]-8.11898816047911*coeff[1]*qEdge[17]; 
  edgeSurf_incr[32] = 15.46875*coeff[1]*qSkin[32]+12.65625*coeff[1]*qEdge[32]+8.118988160479114*coeff[1]*qSkin[21]-8.118988160479114*coeff[1]*qEdge[21]; 
  edgeSurf_incr[33] = 20.34375*coeff[1]*qSkin[33]-0.65625*coeff[1]*qEdge[33]+15.61296411439865*coeff[1]*qSkin[15]+4.720198453190292*coeff[1]*qEdge[15]+4.192627457812107*coeff[1]*qSkin[6]-4.192627457812107*coeff[1]*qEdge[6]; 
  edgeSurf_incr[34] = 15.46875*coeff[1]*qSkin[34]+12.65625*coeff[1]*qEdge[34]+8.118988160479114*coeff[1]*qSkin[23]-8.118988160479114*coeff[1]*qEdge[23]; 
  edgeSurf_incr[35] = 15.46875*coeff[1]*qSkin[35]+12.65625*coeff[1]*qEdge[35]+8.118988160479114*coeff[1]*qSkin[25]-8.118988160479114*coeff[1]*qEdge[25]; 
  edgeSurf_incr[36] = 20.34375*coeff[1]*qSkin[36]-0.65625*coeff[1]*qEdge[36]+15.61296411439865*coeff[1]*qSkin[16]+4.720198453190292*coeff[1]*qEdge[16]+4.192627457812107*coeff[1]*qSkin[8]-4.192627457812107*coeff[1]*qEdge[8]; 
  edgeSurf_incr[37] = 8.118988160479114*coeff[1]*qSkin[44]+8.118988160479114*coeff[1]*qEdge[44]+4.6875*coeff[1]*qSkin[37]-4.6875*coeff[1]*qEdge[37]; 
  edgeSurf_incr[38] = 20.34375*coeff[1]*qSkin[38]-0.65625*coeff[1]*qEdge[38]+15.61296411439865*coeff[1]*qSkin[18]+4.720198453190292*coeff[1]*qEdge[18]+4.192627457812107*coeff[1]*qSkin[10]-4.192627457812107*coeff[1]*qEdge[10]; 
  edgeSurf_incr[39] = 8.118988160479114*coeff[1]*qSkin[46]+8.118988160479114*coeff[1]*qEdge[46]+4.6875*coeff[1]*qSkin[39]-4.6875*coeff[1]*qEdge[39]; 
  edgeSurf_incr[40] = 15.46875*coeff[1]*qSkin[40]+12.65625*coeff[1]*qEdge[40]+8.118988160479114*coeff[1]*qSkin[27]-8.118988160479114*coeff[1]*qEdge[27]; 
  edgeSurf_incr[41] = 15.46875*coeff[1]*qSkin[41]+12.65625*coeff[1]*qEdge[41]+8.118988160479114*coeff[1]*qSkin[28]-8.118988160479114*coeff[1]*qEdge[28]; 
  edgeSurf_incr[42] = 8.118988160479114*coeff[1]*qSkin[47]+8.118988160479114*coeff[1]*qEdge[47]+4.6875*coeff[1]*qSkin[42]-4.6875*coeff[1]*qEdge[42]; 
  edgeSurf_incr[43] = 15.46875*coeff[1]*qSkin[43]+12.65625*coeff[1]*qEdge[43]+8.118988160479114*coeff[1]*qSkin[30]-8.118988160479114*coeff[1]*qEdge[30]; 
  edgeSurf_incr[44] = 15.46875*coeff[1]*qSkin[44]+12.65625*coeff[1]*qEdge[44]+8.118988160479114*coeff[1]*qSkin[37]-8.118988160479114*coeff[1]*qEdge[37]; 
  edgeSurf_incr[45] = 20.34375*coeff[1]*qSkin[45]-0.65625*coeff[1]*qEdge[45]+15.612964114398654*coeff[1]*qSkin[31]+4.72019845319029*coeff[1]*qEdge[31]+4.192627457812107*coeff[1]*qSkin[17]-4.192627457812107*coeff[1]*qEdge[17]; 
  edgeSurf_incr[46] = 15.46875*coeff[1]*qSkin[46]+12.65625*coeff[1]*qEdge[46]+8.118988160479114*coeff[1]*qSkin[39]-8.118988160479114*coeff[1]*qEdge[39]; 
  edgeSurf_incr[47] = 15.46875*coeff[1]*qSkin[47]+12.65625*coeff[1]*qEdge[47]+8.118988160479114*coeff[1]*qSkin[42]-8.118988160479114*coeff[1]*qEdge[42]; 

  boundSurf_incr[2] = 2.8125*coeff[1]*qSkin[2]-5.083290641897234*coeff[1]*qSkin[12]; 
  boundSurf_incr[5] = 2.8125*coeff[1]*qSkin[5]-5.083290641897235*coeff[1]*qSkin[20]; 
  boundSurf_incr[7] = 2.8125*coeff[1]*qSkin[7]-5.083290641897235*coeff[1]*qSkin[22]; 
  boundSurf_incr[9] = 2.8125*coeff[1]*qSkin[9]-5.083290641897235*coeff[1]*qSkin[26]; 
  boundSurf_incr[12] = 19.6875*coeff[1]*qSkin[12]-10.892765661208358*coeff[1]*qSkin[2]; 
  boundSurf_incr[15] = 2.8125*coeff[1]*qSkin[15]-5.083290641897234*coeff[1]*qSkin[33]; 
  boundSurf_incr[16] = 2.8125*coeff[1]*qSkin[16]-5.083290641897234*coeff[1]*qSkin[36]; 
  boundSurf_incr[18] = 2.8125*coeff[1]*qSkin[18]-5.083290641897234*coeff[1]*qSkin[38]; 
  boundSurf_incr[19] = 2.8125*coeff[1]*qSkin[19]; 
  boundSurf_incr[20] = 19.6875*coeff[1]*qSkin[20]-10.892765661208362*coeff[1]*qSkin[5]; 
  boundSurf_incr[22] = 19.6875*coeff[1]*qSkin[22]-10.892765661208362*coeff[1]*qSkin[7]; 
  boundSurf_incr[24] = 2.8125*coeff[1]*qSkin[24]; 
  boundSurf_incr[26] = 19.6875*coeff[1]*qSkin[26]-10.892765661208362*coeff[1]*qSkin[9]; 
  boundSurf_incr[29] = 2.8125*coeff[1]*qSkin[29]; 
  boundSurf_incr[31] = 2.8125*coeff[1]*qSkin[31]-5.083290641897235*coeff[1]*qSkin[45]; 
  boundSurf_incr[32] = 2.8125*coeff[1]*qSkin[32]; 
  boundSurf_incr[33] = 19.6875*coeff[1]*qSkin[33]-10.892765661208358*coeff[1]*qSkin[15]; 
  boundSurf_incr[34] = 2.8125*coeff[1]*qSkin[34]; 
  boundSurf_incr[35] = 2.8125*coeff[1]*qSkin[35]; 
  boundSurf_incr[36] = 19.6875*coeff[1]*qSkin[36]-10.892765661208358*coeff[1]*qSkin[16]; 
  boundSurf_incr[38] = 19.6875*coeff[1]*qSkin[38]-10.892765661208358*coeff[1]*qSkin[18]; 
  boundSurf_incr[40] = 2.8125*coeff[1]*qSkin[40]; 
  boundSurf_incr[41] = 2.8125*coeff[1]*qSkin[41]; 
  boundSurf_incr[43] = 2.8125*coeff[1]*qSkin[43]; 
  boundSurf_incr[44] = 2.8125*coeff[1]*qSkin[44]; 
  boundSurf_incr[45] = 19.6875*coeff[1]*qSkin[45]-10.892765661208362*coeff[1]*qSkin[31]; 
  boundSurf_incr[46] = 2.8125*coeff[1]*qSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[1]*qSkin[47]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[1]*qSkin[12]-6.708203932499369*coeff[1]*qEdge[12]-8.11898816047911*coeff[1]*qSkin[2]-8.11898816047911*coeff[1]*qEdge[2]+4.6875*qSkin[0]*coeff[1]-4.6875*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[1]*qSkin[20]-6.7082039324993685*coeff[1]*qEdge[20]-8.11898816047911*coeff[1]*qSkin[5]-8.11898816047911*coeff[1]*qEdge[5]+4.6875*coeff[1]*qSkin[1]-4.6875*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = -(14.160595359570864*coeff[1]*qSkin[12])+9.077304717673634*coeff[1]*qEdge[12]+15.46875*coeff[1]*qSkin[2]+12.65625*coeff[1]*qEdge[2]-8.11898816047911*qSkin[0]*coeff[1]+8.11898816047911*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[1]*qSkin[22]-6.7082039324993685*coeff[1]*qEdge[22]-8.11898816047911*coeff[1]*qSkin[7]-8.11898816047911*coeff[1]*qEdge[7]+4.6875*coeff[1]*qSkin[3]-4.6875*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = 6.7082039324993685*coeff[1]*qSkin[26]-6.7082039324993685*coeff[1]*qEdge[26]-8.11898816047911*coeff[1]*qSkin[9]-8.11898816047911*coeff[1]*qEdge[9]+4.6875*coeff[1]*qSkin[4]-4.6875*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = -(14.160595359570868*coeff[1]*qSkin[20])+9.077304717673634*coeff[1]*qEdge[20]+15.46875*coeff[1]*qSkin[5]+12.65625*coeff[1]*qEdge[5]-8.11898816047911*coeff[1]*qSkin[1]+8.11898816047911*coeff[1]*qEdge[1]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[1]*qSkin[33]-6.708203932499369*coeff[1]*qEdge[33]-8.11898816047911*coeff[1]*qSkin[15]-8.11898816047911*coeff[1]*qEdge[15]+4.6875*coeff[1]*qSkin[6]-4.6875*coeff[1]*qEdge[6]; 
  edgeSurf_incr[7] = -(14.160595359570868*coeff[1]*qSkin[22])+9.077304717673634*coeff[1]*qEdge[22]+15.46875*coeff[1]*qSkin[7]+12.65625*coeff[1]*qEdge[7]-8.11898816047911*coeff[1]*qSkin[3]+8.11898816047911*coeff[1]*qEdge[3]; 
  edgeSurf_incr[8] = 6.708203932499369*coeff[1]*qSkin[36]-6.708203932499369*coeff[1]*qEdge[36]-8.11898816047911*coeff[1]*qSkin[16]-8.11898816047911*coeff[1]*qEdge[16]+4.6875*coeff[1]*qSkin[8]-4.6875*coeff[1]*qEdge[8]; 
  edgeSurf_incr[9] = -(14.160595359570868*coeff[1]*qSkin[26])+9.077304717673634*coeff[1]*qEdge[26]+15.46875*coeff[1]*qSkin[9]+12.65625*coeff[1]*qEdge[9]-8.11898816047911*coeff[1]*qSkin[4]+8.11898816047911*coeff[1]*qEdge[4]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[1]*qSkin[38]-6.708203932499369*coeff[1]*qEdge[38]-8.11898816047911*coeff[1]*qSkin[18]-8.11898816047911*coeff[1]*qEdge[18]+4.6875*coeff[1]*qSkin[10]-4.6875*coeff[1]*qEdge[10]; 
  edgeSurf_incr[11] = -(8.118988160479114*coeff[1]*qSkin[19])-8.118988160479114*coeff[1]*qEdge[19]+4.6875*coeff[1]*qSkin[11]-4.6875*coeff[1]*qEdge[11]; 
  edgeSurf_incr[12] = 20.34375*coeff[1]*qSkin[12]-0.65625*coeff[1]*qEdge[12]-15.61296411439865*coeff[1]*qSkin[2]-4.720198453190292*coeff[1]*qEdge[2]+4.192627457812107*qSkin[0]*coeff[1]-4.192627457812107*qEdge[0]*coeff[1]; 
  edgeSurf_incr[13] = -(8.118988160479114*coeff[1]*qSkin[24])-8.118988160479114*coeff[1]*qEdge[24]+4.6875*coeff[1]*qSkin[13]-4.6875*coeff[1]*qEdge[13]; 
  edgeSurf_incr[14] = -(8.118988160479114*coeff[1]*qSkin[29])-8.118988160479114*coeff[1]*qEdge[29]+4.6875*coeff[1]*qSkin[14]-4.6875*coeff[1]*qEdge[14]; 
  edgeSurf_incr[15] = -(14.160595359570864*coeff[1]*qSkin[33])+9.077304717673634*coeff[1]*qEdge[33]+15.46875*coeff[1]*qSkin[15]+12.65625*coeff[1]*qEdge[15]-8.11898816047911*coeff[1]*qSkin[6]+8.11898816047911*coeff[1]*qEdge[6]; 
  edgeSurf_incr[16] = -(14.160595359570864*coeff[1]*qSkin[36])+9.077304717673634*coeff[1]*qEdge[36]+15.46875*coeff[1]*qSkin[16]+12.65625*coeff[1]*qEdge[16]-8.11898816047911*coeff[1]*qSkin[8]+8.11898816047911*coeff[1]*qEdge[8]; 
  edgeSurf_incr[17] = 6.7082039324993685*coeff[1]*qSkin[45]-6.7082039324993685*coeff[1]*qEdge[45]-8.11898816047911*coeff[1]*qSkin[31]-8.11898816047911*coeff[1]*qEdge[31]+4.6875*coeff[1]*qSkin[17]-4.6875*coeff[1]*qEdge[17]; 
  edgeSurf_incr[18] = -(14.160595359570864*coeff[1]*qSkin[38])+9.077304717673634*coeff[1]*qEdge[38]+15.46875*coeff[1]*qSkin[18]+12.65625*coeff[1]*qEdge[18]-8.11898816047911*coeff[1]*qSkin[10]+8.11898816047911*coeff[1]*qEdge[10]; 
  edgeSurf_incr[19] = 15.46875*coeff[1]*qSkin[19]+12.65625*coeff[1]*qEdge[19]-8.118988160479114*coeff[1]*qSkin[11]+8.118988160479114*coeff[1]*qEdge[11]; 
  edgeSurf_incr[20] = 20.34375*coeff[1]*qSkin[20]-0.65625*coeff[1]*qEdge[20]-15.612964114398654*coeff[1]*qSkin[5]-4.72019845319029*coeff[1]*qEdge[5]+4.192627457812107*coeff[1]*qSkin[1]-4.192627457812107*coeff[1]*qEdge[1]; 
  edgeSurf_incr[21] = -(8.118988160479114*coeff[1]*qSkin[32])-8.118988160479114*coeff[1]*qEdge[32]+4.6875*coeff[1]*qSkin[21]-4.6875*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = 20.34375*coeff[1]*qSkin[22]-0.65625*coeff[1]*qEdge[22]-15.612964114398654*coeff[1]*qSkin[7]-4.72019845319029*coeff[1]*qEdge[7]+4.192627457812107*coeff[1]*qSkin[3]-4.192627457812107*coeff[1]*qEdge[3]; 
  edgeSurf_incr[23] = -(8.118988160479114*coeff[1]*qSkin[34])-8.118988160479114*coeff[1]*qEdge[34]+4.6875*coeff[1]*qSkin[23]-4.6875*coeff[1]*qEdge[23]; 
  edgeSurf_incr[24] = 15.46875*coeff[1]*qSkin[24]+12.65625*coeff[1]*qEdge[24]-8.118988160479114*coeff[1]*qSkin[13]+8.118988160479114*coeff[1]*qEdge[13]; 
  edgeSurf_incr[25] = -(8.118988160479114*coeff[1]*qSkin[35])-8.118988160479114*coeff[1]*qEdge[35]+4.6875*coeff[1]*qSkin[25]-4.6875*coeff[1]*qEdge[25]; 
  edgeSurf_incr[26] = 20.34375*coeff[1]*qSkin[26]-0.65625*coeff[1]*qEdge[26]-15.612964114398654*coeff[1]*qSkin[9]-4.72019845319029*coeff[1]*qEdge[9]+4.192627457812107*coeff[1]*qSkin[4]-4.192627457812107*coeff[1]*qEdge[4]; 
  edgeSurf_incr[27] = -(8.118988160479114*coeff[1]*qSkin[40])-8.118988160479114*coeff[1]*qEdge[40]+4.6875*coeff[1]*qSkin[27]-4.6875*coeff[1]*qEdge[27]; 
  edgeSurf_incr[28] = -(8.118988160479114*coeff[1]*qSkin[41])-8.118988160479114*coeff[1]*qEdge[41]+4.6875*coeff[1]*qSkin[28]-4.6875*coeff[1]*qEdge[28]; 
  edgeSurf_incr[29] = 15.46875*coeff[1]*qSkin[29]+12.65625*coeff[1]*qEdge[29]-8.118988160479114*coeff[1]*qSkin[14]+8.118988160479114*coeff[1]*qEdge[14]; 
  edgeSurf_incr[30] = -(8.118988160479114*coeff[1]*qSkin[43])-8.118988160479114*coeff[1]*qEdge[43]+4.6875*coeff[1]*qSkin[30]-4.6875*coeff[1]*qEdge[30]; 
  edgeSurf_incr[31] = -(14.160595359570868*coeff[1]*qSkin[45])+9.077304717673634*coeff[1]*qEdge[45]+15.46875*coeff[1]*qSkin[31]+12.65625*coeff[1]*qEdge[31]-8.11898816047911*coeff[1]*qSkin[17]+8.11898816047911*coeff[1]*qEdge[17]; 
  edgeSurf_incr[32] = 15.46875*coeff[1]*qSkin[32]+12.65625*coeff[1]*qEdge[32]-8.118988160479114*coeff[1]*qSkin[21]+8.118988160479114*coeff[1]*qEdge[21]; 
  edgeSurf_incr[33] = 20.34375*coeff[1]*qSkin[33]-0.65625*coeff[1]*qEdge[33]-15.61296411439865*coeff[1]*qSkin[15]-4.720198453190292*coeff[1]*qEdge[15]+4.192627457812107*coeff[1]*qSkin[6]-4.192627457812107*coeff[1]*qEdge[6]; 
  edgeSurf_incr[34] = 15.46875*coeff[1]*qSkin[34]+12.65625*coeff[1]*qEdge[34]-8.118988160479114*coeff[1]*qSkin[23]+8.118988160479114*coeff[1]*qEdge[23]; 
  edgeSurf_incr[35] = 15.46875*coeff[1]*qSkin[35]+12.65625*coeff[1]*qEdge[35]-8.118988160479114*coeff[1]*qSkin[25]+8.118988160479114*coeff[1]*qEdge[25]; 
  edgeSurf_incr[36] = 20.34375*coeff[1]*qSkin[36]-0.65625*coeff[1]*qEdge[36]-15.61296411439865*coeff[1]*qSkin[16]-4.720198453190292*coeff[1]*qEdge[16]+4.192627457812107*coeff[1]*qSkin[8]-4.192627457812107*coeff[1]*qEdge[8]; 
  edgeSurf_incr[37] = -(8.118988160479114*coeff[1]*qSkin[44])-8.118988160479114*coeff[1]*qEdge[44]+4.6875*coeff[1]*qSkin[37]-4.6875*coeff[1]*qEdge[37]; 
  edgeSurf_incr[38] = 20.34375*coeff[1]*qSkin[38]-0.65625*coeff[1]*qEdge[38]-15.61296411439865*coeff[1]*qSkin[18]-4.720198453190292*coeff[1]*qEdge[18]+4.192627457812107*coeff[1]*qSkin[10]-4.192627457812107*coeff[1]*qEdge[10]; 
  edgeSurf_incr[39] = -(8.118988160479114*coeff[1]*qSkin[46])-8.118988160479114*coeff[1]*qEdge[46]+4.6875*coeff[1]*qSkin[39]-4.6875*coeff[1]*qEdge[39]; 
  edgeSurf_incr[40] = 15.46875*coeff[1]*qSkin[40]+12.65625*coeff[1]*qEdge[40]-8.118988160479114*coeff[1]*qSkin[27]+8.118988160479114*coeff[1]*qEdge[27]; 
  edgeSurf_incr[41] = 15.46875*coeff[1]*qSkin[41]+12.65625*coeff[1]*qEdge[41]-8.118988160479114*coeff[1]*qSkin[28]+8.118988160479114*coeff[1]*qEdge[28]; 
  edgeSurf_incr[42] = -(8.118988160479114*coeff[1]*qSkin[47])-8.118988160479114*coeff[1]*qEdge[47]+4.6875*coeff[1]*qSkin[42]-4.6875*coeff[1]*qEdge[42]; 
  edgeSurf_incr[43] = 15.46875*coeff[1]*qSkin[43]+12.65625*coeff[1]*qEdge[43]-8.118988160479114*coeff[1]*qSkin[30]+8.118988160479114*coeff[1]*qEdge[30]; 
  edgeSurf_incr[44] = 15.46875*coeff[1]*qSkin[44]+12.65625*coeff[1]*qEdge[44]-8.118988160479114*coeff[1]*qSkin[37]+8.118988160479114*coeff[1]*qEdge[37]; 
  edgeSurf_incr[45] = 20.34375*coeff[1]*qSkin[45]-0.65625*coeff[1]*qEdge[45]-15.612964114398654*coeff[1]*qSkin[31]-4.72019845319029*coeff[1]*qEdge[31]+4.192627457812107*coeff[1]*qSkin[17]-4.192627457812107*coeff[1]*qEdge[17]; 
  edgeSurf_incr[46] = 15.46875*coeff[1]*qSkin[46]+12.65625*coeff[1]*qEdge[46]-8.118988160479114*coeff[1]*qSkin[39]+8.118988160479114*coeff[1]*qEdge[39]; 
  edgeSurf_incr[47] = 15.46875*coeff[1]*qSkin[47]+12.65625*coeff[1]*qEdge[47]-8.118988160479114*coeff[1]*qSkin[42]+8.118988160479114*coeff[1]*qEdge[42]; 

  boundSurf_incr[2] = 5.083290641897234*coeff[1]*qSkin[12]+2.8125*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[1]*qSkin[20]+2.8125*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = 5.083290641897235*coeff[1]*qSkin[22]+2.8125*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = 5.083290641897235*coeff[1]*qSkin[26]+2.8125*coeff[1]*qSkin[9]; 
  boundSurf_incr[12] = 19.6875*coeff[1]*qSkin[12]+10.892765661208358*coeff[1]*qSkin[2]; 
  boundSurf_incr[15] = 5.083290641897234*coeff[1]*qSkin[33]+2.8125*coeff[1]*qSkin[15]; 
  boundSurf_incr[16] = 5.083290641897234*coeff[1]*qSkin[36]+2.8125*coeff[1]*qSkin[16]; 
  boundSurf_incr[18] = 5.083290641897234*coeff[1]*qSkin[38]+2.8125*coeff[1]*qSkin[18]; 
  boundSurf_incr[19] = 2.8125*coeff[1]*qSkin[19]; 
  boundSurf_incr[20] = 19.6875*coeff[1]*qSkin[20]+10.892765661208362*coeff[1]*qSkin[5]; 
  boundSurf_incr[22] = 19.6875*coeff[1]*qSkin[22]+10.892765661208362*coeff[1]*qSkin[7]; 
  boundSurf_incr[24] = 2.8125*coeff[1]*qSkin[24]; 
  boundSurf_incr[26] = 19.6875*coeff[1]*qSkin[26]+10.892765661208362*coeff[1]*qSkin[9]; 
  boundSurf_incr[29] = 2.8125*coeff[1]*qSkin[29]; 
  boundSurf_incr[31] = 5.083290641897235*coeff[1]*qSkin[45]+2.8125*coeff[1]*qSkin[31]; 
  boundSurf_incr[32] = 2.8125*coeff[1]*qSkin[32]; 
  boundSurf_incr[33] = 19.6875*coeff[1]*qSkin[33]+10.892765661208358*coeff[1]*qSkin[15]; 
  boundSurf_incr[34] = 2.8125*coeff[1]*qSkin[34]; 
  boundSurf_incr[35] = 2.8125*coeff[1]*qSkin[35]; 
  boundSurf_incr[36] = 19.6875*coeff[1]*qSkin[36]+10.892765661208358*coeff[1]*qSkin[16]; 
  boundSurf_incr[38] = 19.6875*coeff[1]*qSkin[38]+10.892765661208358*coeff[1]*qSkin[18]; 
  boundSurf_incr[40] = 2.8125*coeff[1]*qSkin[40]; 
  boundSurf_incr[41] = 2.8125*coeff[1]*qSkin[41]; 
  boundSurf_incr[43] = 2.8125*coeff[1]*qSkin[43]; 
  boundSurf_incr[44] = 2.8125*coeff[1]*qSkin[44]; 
  boundSurf_incr[45] = 19.6875*coeff[1]*qSkin[45]+10.892765661208362*coeff[1]*qSkin[31]; 
  boundSurf_incr[46] = 2.8125*coeff[1]*qSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[1]*qSkin[47]; 

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

