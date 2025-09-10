#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*qSkin[11]-6.708203932499369*coeff[0]*qEdge[11]+8.11898816047911*coeff[0]*qSkin[1]+8.11898816047911*coeff[0]*qEdge[1]+4.6875*coeff[0]*qSkin[0]-4.6875*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 14.160595359570864*coeff[0]*qSkin[11]-9.077304717673634*coeff[0]*qEdge[11]+15.46875*coeff[0]*qSkin[1]+12.65625*coeff[0]*qEdge[1]+8.11898816047911*coeff[0]*qSkin[0]-8.11898816047911*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*qSkin[19]-6.7082039324993685*coeff[0]*qEdge[19]+8.11898816047911*coeff[0]*qSkin[5]+8.11898816047911*coeff[0]*qEdge[5]+4.6875*coeff[0]*qSkin[2]-4.6875*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*qSkin[21]-6.7082039324993685*coeff[0]*qEdge[21]+8.11898816047911*coeff[0]*qSkin[6]+8.11898816047911*coeff[0]*qEdge[6]+4.6875*coeff[0]*qSkin[3]-4.6875*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = 6.7082039324993685*coeff[0]*qSkin[25]-6.7082039324993685*coeff[0]*qEdge[25]+8.11898816047911*coeff[0]*qSkin[8]+8.11898816047911*coeff[0]*qEdge[8]+4.6875*coeff[0]*qSkin[4]-4.6875*coeff[0]*qEdge[4]; 
  edgeSurf_incr[5] = 14.160595359570868*coeff[0]*qSkin[19]-9.077304717673634*coeff[0]*qEdge[19]+15.46875*coeff[0]*qSkin[5]+12.65625*coeff[0]*qEdge[5]+8.11898816047911*coeff[0]*qSkin[2]-8.11898816047911*coeff[0]*qEdge[2]; 
  edgeSurf_incr[6] = 14.160595359570868*coeff[0]*qSkin[21]-9.077304717673634*coeff[0]*qEdge[21]+15.46875*coeff[0]*qSkin[6]+12.65625*coeff[0]*qEdge[6]+8.11898816047911*coeff[0]*qSkin[3]-8.11898816047911*coeff[0]*qEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[0]*qSkin[32]-6.708203932499369*coeff[0]*qEdge[32]+8.11898816047911*coeff[0]*qSkin[15]+8.11898816047911*coeff[0]*qEdge[15]+4.6875*coeff[0]*qSkin[7]-4.6875*coeff[0]*qEdge[7]; 
  edgeSurf_incr[8] = 14.160595359570868*coeff[0]*qSkin[25]-9.077304717673634*coeff[0]*qEdge[25]+15.46875*coeff[0]*qSkin[8]+12.65625*coeff[0]*qEdge[8]+8.11898816047911*coeff[0]*qSkin[4]-8.11898816047911*coeff[0]*qEdge[4]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*qSkin[35]-6.708203932499369*coeff[0]*qEdge[35]+8.11898816047911*coeff[0]*qSkin[16]+8.11898816047911*coeff[0]*qEdge[16]+4.6875*coeff[0]*qSkin[9]-4.6875*coeff[0]*qEdge[9]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[0]*qSkin[37]-6.708203932499369*coeff[0]*qEdge[37]+8.11898816047911*coeff[0]*qSkin[17]+8.11898816047911*coeff[0]*qEdge[17]+4.6875*coeff[0]*qSkin[10]-4.6875*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*qSkin[11]-0.65625*coeff[0]*qEdge[11]+15.61296411439865*coeff[0]*qSkin[1]+4.720198453190292*coeff[0]*qEdge[1]+4.192627457812107*coeff[0]*qSkin[0]-4.192627457812107*coeff[0]*qEdge[0]; 
  edgeSurf_incr[12] = 8.118988160479114*coeff[0]*qSkin[20]+8.118988160479114*coeff[0]*qEdge[20]+4.6875*coeff[0]*qSkin[12]-4.6875*coeff[0]*qEdge[12]; 
  edgeSurf_incr[13] = 8.118988160479114*coeff[0]*qSkin[23]+8.118988160479114*coeff[0]*qEdge[23]+4.6875*coeff[0]*qSkin[13]-4.6875*coeff[0]*qEdge[13]; 
  edgeSurf_incr[14] = 8.118988160479114*coeff[0]*qSkin[28]+8.118988160479114*coeff[0]*qEdge[28]+4.6875*coeff[0]*qSkin[14]-4.6875*coeff[0]*qEdge[14]; 
  edgeSurf_incr[15] = 14.160595359570864*coeff[0]*qSkin[32]-9.077304717673634*coeff[0]*qEdge[32]+15.46875*coeff[0]*qSkin[15]+12.65625*coeff[0]*qEdge[15]+8.11898816047911*coeff[0]*qSkin[7]-8.11898816047911*coeff[0]*qEdge[7]; 
  edgeSurf_incr[16] = 14.160595359570864*coeff[0]*qSkin[35]-9.077304717673634*coeff[0]*qEdge[35]+15.46875*coeff[0]*qSkin[16]+12.65625*coeff[0]*qEdge[16]+8.11898816047911*coeff[0]*qSkin[9]-8.11898816047911*coeff[0]*qEdge[9]; 
  edgeSurf_incr[17] = 14.160595359570864*coeff[0]*qSkin[37]-9.077304717673634*coeff[0]*qEdge[37]+15.46875*coeff[0]*qSkin[17]+12.65625*coeff[0]*qEdge[17]+8.11898816047911*coeff[0]*qSkin[10]-8.11898816047911*coeff[0]*qEdge[10]; 
  edgeSurf_incr[18] = 6.7082039324993685*coeff[0]*qSkin[44]-6.7082039324993685*coeff[0]*qEdge[44]+8.11898816047911*coeff[0]*qSkin[31]+8.11898816047911*coeff[0]*qEdge[31]+4.6875*coeff[0]*qSkin[18]-4.6875*coeff[0]*qEdge[18]; 
  edgeSurf_incr[19] = 20.34375*coeff[0]*qSkin[19]-0.65625*coeff[0]*qEdge[19]+15.612964114398654*coeff[0]*qSkin[5]+4.72019845319029*coeff[0]*qEdge[5]+4.192627457812107*coeff[0]*qSkin[2]-4.192627457812107*coeff[0]*qEdge[2]; 
  edgeSurf_incr[20] = 15.46875*coeff[0]*qSkin[20]+12.65625*coeff[0]*qEdge[20]+8.118988160479114*coeff[0]*qSkin[12]-8.118988160479114*coeff[0]*qEdge[12]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*qSkin[21]-0.65625*coeff[0]*qEdge[21]+15.612964114398654*coeff[0]*qSkin[6]+4.72019845319029*coeff[0]*qEdge[6]+4.192627457812107*coeff[0]*qSkin[3]-4.192627457812107*coeff[0]*qEdge[3]; 
  edgeSurf_incr[22] = 8.118988160479114*coeff[0]*qSkin[33]+8.118988160479114*coeff[0]*qEdge[33]+4.6875*coeff[0]*qSkin[22]-4.6875*coeff[0]*qEdge[22]; 
  edgeSurf_incr[23] = 15.46875*coeff[0]*qSkin[23]+12.65625*coeff[0]*qEdge[23]+8.118988160479114*coeff[0]*qSkin[13]-8.118988160479114*coeff[0]*qEdge[13]; 
  edgeSurf_incr[24] = 8.118988160479114*coeff[0]*qSkin[34]+8.118988160479114*coeff[0]*qEdge[34]+4.6875*coeff[0]*qSkin[24]-4.6875*coeff[0]*qEdge[24]; 
  edgeSurf_incr[25] = 20.34375*coeff[0]*qSkin[25]-0.65625*coeff[0]*qEdge[25]+15.612964114398654*coeff[0]*qSkin[8]+4.72019845319029*coeff[0]*qEdge[8]+4.192627457812107*coeff[0]*qSkin[4]-4.192627457812107*coeff[0]*qEdge[4]; 
  edgeSurf_incr[26] = 8.118988160479114*coeff[0]*qSkin[36]+8.118988160479114*coeff[0]*qEdge[36]+4.6875*coeff[0]*qSkin[26]-4.6875*coeff[0]*qEdge[26]; 
  edgeSurf_incr[27] = 8.118988160479114*coeff[0]*qSkin[39]+8.118988160479114*coeff[0]*qEdge[39]+4.6875*coeff[0]*qSkin[27]-4.6875*coeff[0]*qEdge[27]; 
  edgeSurf_incr[28] = 15.46875*coeff[0]*qSkin[28]+12.65625*coeff[0]*qEdge[28]+8.118988160479114*coeff[0]*qSkin[14]-8.118988160479114*coeff[0]*qEdge[14]; 
  edgeSurf_incr[29] = 8.118988160479114*coeff[0]*qSkin[41]+8.118988160479114*coeff[0]*qEdge[41]+4.6875*coeff[0]*qSkin[29]-4.6875*coeff[0]*qEdge[29]; 
  edgeSurf_incr[30] = 8.118988160479114*coeff[0]*qSkin[42]+8.118988160479114*coeff[0]*qEdge[42]+4.6875*coeff[0]*qSkin[30]-4.6875*coeff[0]*qEdge[30]; 
  edgeSurf_incr[31] = 14.160595359570868*coeff[0]*qSkin[44]-9.077304717673634*coeff[0]*qEdge[44]+15.46875*coeff[0]*qSkin[31]+12.65625*coeff[0]*qEdge[31]+8.11898816047911*coeff[0]*qSkin[18]-8.11898816047911*coeff[0]*qEdge[18]; 
  edgeSurf_incr[32] = 20.34375*coeff[0]*qSkin[32]-0.65625*coeff[0]*qEdge[32]+15.61296411439865*coeff[0]*qSkin[15]+4.720198453190292*coeff[0]*qEdge[15]+4.192627457812107*coeff[0]*qSkin[7]-4.192627457812107*coeff[0]*qEdge[7]; 
  edgeSurf_incr[33] = 15.46875*coeff[0]*qSkin[33]+12.65625*coeff[0]*qEdge[33]+8.118988160479114*coeff[0]*qSkin[22]-8.118988160479114*coeff[0]*qEdge[22]; 
  edgeSurf_incr[34] = 15.46875*coeff[0]*qSkin[34]+12.65625*coeff[0]*qEdge[34]+8.118988160479114*coeff[0]*qSkin[24]-8.118988160479114*coeff[0]*qEdge[24]; 
  edgeSurf_incr[35] = 20.34375*coeff[0]*qSkin[35]-0.65625*coeff[0]*qEdge[35]+15.61296411439865*coeff[0]*qSkin[16]+4.720198453190292*coeff[0]*qEdge[16]+4.192627457812107*coeff[0]*qSkin[9]-4.192627457812107*coeff[0]*qEdge[9]; 
  edgeSurf_incr[36] = 15.46875*coeff[0]*qSkin[36]+12.65625*coeff[0]*qEdge[36]+8.118988160479114*coeff[0]*qSkin[26]-8.118988160479114*coeff[0]*qEdge[26]; 
  edgeSurf_incr[37] = 20.34375*coeff[0]*qSkin[37]-0.65625*coeff[0]*qEdge[37]+15.61296411439865*coeff[0]*qSkin[17]+4.720198453190292*coeff[0]*qEdge[17]+4.192627457812107*coeff[0]*qSkin[10]-4.192627457812107*coeff[0]*qEdge[10]; 
  edgeSurf_incr[38] = 8.118988160479114*coeff[0]*qSkin[45]+8.118988160479114*coeff[0]*qEdge[45]+4.6875*coeff[0]*qSkin[38]-4.6875*coeff[0]*qEdge[38]; 
  edgeSurf_incr[39] = 15.46875*coeff[0]*qSkin[39]+12.65625*coeff[0]*qEdge[39]+8.118988160479114*coeff[0]*qSkin[27]-8.118988160479114*coeff[0]*qEdge[27]; 
  edgeSurf_incr[40] = 8.118988160479114*coeff[0]*qSkin[46]+8.118988160479114*coeff[0]*qEdge[46]+4.6875*coeff[0]*qSkin[40]-4.6875*coeff[0]*qEdge[40]; 
  edgeSurf_incr[41] = 15.46875*coeff[0]*qSkin[41]+12.65625*coeff[0]*qEdge[41]+8.118988160479114*coeff[0]*qSkin[29]-8.118988160479114*coeff[0]*qEdge[29]; 
  edgeSurf_incr[42] = 15.46875*coeff[0]*qSkin[42]+12.65625*coeff[0]*qEdge[42]+8.118988160479114*coeff[0]*qSkin[30]-8.118988160479114*coeff[0]*qEdge[30]; 
  edgeSurf_incr[43] = 8.118988160479114*coeff[0]*qSkin[47]+8.118988160479114*coeff[0]*qEdge[47]+4.6875*coeff[0]*qSkin[43]-4.6875*coeff[0]*qEdge[43]; 
  edgeSurf_incr[44] = 20.34375*coeff[0]*qSkin[44]-0.65625*coeff[0]*qEdge[44]+15.612964114398654*coeff[0]*qSkin[31]+4.72019845319029*coeff[0]*qEdge[31]+4.192627457812107*coeff[0]*qSkin[18]-4.192627457812107*coeff[0]*qEdge[18]; 
  edgeSurf_incr[45] = 15.46875*coeff[0]*qSkin[45]+12.65625*coeff[0]*qEdge[45]+8.118988160479114*coeff[0]*qSkin[38]-8.118988160479114*coeff[0]*qEdge[38]; 
  edgeSurf_incr[46] = 15.46875*coeff[0]*qSkin[46]+12.65625*coeff[0]*qEdge[46]+8.118988160479114*coeff[0]*qSkin[40]-8.118988160479114*coeff[0]*qEdge[40]; 
  edgeSurf_incr[47] = 15.46875*coeff[0]*qSkin[47]+12.65625*coeff[0]*qEdge[47]+8.118988160479114*coeff[0]*qSkin[43]-8.118988160479114*coeff[0]*qEdge[43]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*qSkin[1]-5.083290641897234*coeff[0]*qSkin[11]; 
  boundSurf_incr[5] = 2.8125*coeff[0]*qSkin[5]-5.083290641897235*coeff[0]*qSkin[19]; 
  boundSurf_incr[6] = 2.8125*coeff[0]*qSkin[6]-5.083290641897235*coeff[0]*qSkin[21]; 
  boundSurf_incr[8] = 2.8125*coeff[0]*qSkin[8]-5.083290641897235*coeff[0]*qSkin[25]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*qSkin[11]-10.892765661208358*coeff[0]*qSkin[1]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*qSkin[15]-5.083290641897234*coeff[0]*qSkin[32]; 
  boundSurf_incr[16] = 2.8125*coeff[0]*qSkin[16]-5.083290641897234*coeff[0]*qSkin[35]; 
  boundSurf_incr[17] = 2.8125*coeff[0]*qSkin[17]-5.083290641897234*coeff[0]*qSkin[37]; 
  boundSurf_incr[19] = 19.6875*coeff[0]*qSkin[19]-10.892765661208362*coeff[0]*qSkin[5]; 
  boundSurf_incr[20] = 2.8125*coeff[0]*qSkin[20]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*qSkin[21]-10.892765661208362*coeff[0]*qSkin[6]; 
  boundSurf_incr[23] = 2.8125*coeff[0]*qSkin[23]; 
  boundSurf_incr[25] = 19.6875*coeff[0]*qSkin[25]-10.892765661208362*coeff[0]*qSkin[8]; 
  boundSurf_incr[28] = 2.8125*coeff[0]*qSkin[28]; 
  boundSurf_incr[31] = 2.8125*coeff[0]*qSkin[31]-5.083290641897235*coeff[0]*qSkin[44]; 
  boundSurf_incr[32] = 19.6875*coeff[0]*qSkin[32]-10.892765661208358*coeff[0]*qSkin[15]; 
  boundSurf_incr[33] = 2.8125*coeff[0]*qSkin[33]; 
  boundSurf_incr[34] = 2.8125*coeff[0]*qSkin[34]; 
  boundSurf_incr[35] = 19.6875*coeff[0]*qSkin[35]-10.892765661208358*coeff[0]*qSkin[16]; 
  boundSurf_incr[36] = 2.8125*coeff[0]*qSkin[36]; 
  boundSurf_incr[37] = 19.6875*coeff[0]*qSkin[37]-10.892765661208358*coeff[0]*qSkin[17]; 
  boundSurf_incr[39] = 2.8125*coeff[0]*qSkin[39]; 
  boundSurf_incr[41] = 2.8125*coeff[0]*qSkin[41]; 
  boundSurf_incr[42] = 2.8125*coeff[0]*qSkin[42]; 
  boundSurf_incr[44] = 19.6875*coeff[0]*qSkin[44]-10.892765661208362*coeff[0]*qSkin[31]; 
  boundSurf_incr[45] = 2.8125*coeff[0]*qSkin[45]; 
  boundSurf_incr[46] = 2.8125*coeff[0]*qSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[0]*qSkin[47]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*qSkin[11]-6.708203932499369*coeff[0]*qEdge[11]-8.11898816047911*coeff[0]*qSkin[1]-8.11898816047911*coeff[0]*qEdge[1]+4.6875*coeff[0]*qSkin[0]-4.6875*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(14.160595359570864*coeff[0]*qSkin[11])+9.077304717673634*coeff[0]*qEdge[11]+15.46875*coeff[0]*qSkin[1]+12.65625*coeff[0]*qEdge[1]-8.11898816047911*coeff[0]*qSkin[0]+8.11898816047911*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*qSkin[19]-6.7082039324993685*coeff[0]*qEdge[19]-8.11898816047911*coeff[0]*qSkin[5]-8.11898816047911*coeff[0]*qEdge[5]+4.6875*coeff[0]*qSkin[2]-4.6875*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*qSkin[21]-6.7082039324993685*coeff[0]*qEdge[21]-8.11898816047911*coeff[0]*qSkin[6]-8.11898816047911*coeff[0]*qEdge[6]+4.6875*coeff[0]*qSkin[3]-4.6875*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = 6.7082039324993685*coeff[0]*qSkin[25]-6.7082039324993685*coeff[0]*qEdge[25]-8.11898816047911*coeff[0]*qSkin[8]-8.11898816047911*coeff[0]*qEdge[8]+4.6875*coeff[0]*qSkin[4]-4.6875*coeff[0]*qEdge[4]; 
  edgeSurf_incr[5] = -(14.160595359570868*coeff[0]*qSkin[19])+9.077304717673634*coeff[0]*qEdge[19]+15.46875*coeff[0]*qSkin[5]+12.65625*coeff[0]*qEdge[5]-8.11898816047911*coeff[0]*qSkin[2]+8.11898816047911*coeff[0]*qEdge[2]; 
  edgeSurf_incr[6] = -(14.160595359570868*coeff[0]*qSkin[21])+9.077304717673634*coeff[0]*qEdge[21]+15.46875*coeff[0]*qSkin[6]+12.65625*coeff[0]*qEdge[6]-8.11898816047911*coeff[0]*qSkin[3]+8.11898816047911*coeff[0]*qEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[0]*qSkin[32]-6.708203932499369*coeff[0]*qEdge[32]-8.11898816047911*coeff[0]*qSkin[15]-8.11898816047911*coeff[0]*qEdge[15]+4.6875*coeff[0]*qSkin[7]-4.6875*coeff[0]*qEdge[7]; 
  edgeSurf_incr[8] = -(14.160595359570868*coeff[0]*qSkin[25])+9.077304717673634*coeff[0]*qEdge[25]+15.46875*coeff[0]*qSkin[8]+12.65625*coeff[0]*qEdge[8]-8.11898816047911*coeff[0]*qSkin[4]+8.11898816047911*coeff[0]*qEdge[4]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*qSkin[35]-6.708203932499369*coeff[0]*qEdge[35]-8.11898816047911*coeff[0]*qSkin[16]-8.11898816047911*coeff[0]*qEdge[16]+4.6875*coeff[0]*qSkin[9]-4.6875*coeff[0]*qEdge[9]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[0]*qSkin[37]-6.708203932499369*coeff[0]*qEdge[37]-8.11898816047911*coeff[0]*qSkin[17]-8.11898816047911*coeff[0]*qEdge[17]+4.6875*coeff[0]*qSkin[10]-4.6875*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*qSkin[11]-0.65625*coeff[0]*qEdge[11]-15.61296411439865*coeff[0]*qSkin[1]-4.720198453190292*coeff[0]*qEdge[1]+4.192627457812107*coeff[0]*qSkin[0]-4.192627457812107*coeff[0]*qEdge[0]; 
  edgeSurf_incr[12] = -(8.118988160479114*coeff[0]*qSkin[20])-8.118988160479114*coeff[0]*qEdge[20]+4.6875*coeff[0]*qSkin[12]-4.6875*coeff[0]*qEdge[12]; 
  edgeSurf_incr[13] = -(8.118988160479114*coeff[0]*qSkin[23])-8.118988160479114*coeff[0]*qEdge[23]+4.6875*coeff[0]*qSkin[13]-4.6875*coeff[0]*qEdge[13]; 
  edgeSurf_incr[14] = -(8.118988160479114*coeff[0]*qSkin[28])-8.118988160479114*coeff[0]*qEdge[28]+4.6875*coeff[0]*qSkin[14]-4.6875*coeff[0]*qEdge[14]; 
  edgeSurf_incr[15] = -(14.160595359570864*coeff[0]*qSkin[32])+9.077304717673634*coeff[0]*qEdge[32]+15.46875*coeff[0]*qSkin[15]+12.65625*coeff[0]*qEdge[15]-8.11898816047911*coeff[0]*qSkin[7]+8.11898816047911*coeff[0]*qEdge[7]; 
  edgeSurf_incr[16] = -(14.160595359570864*coeff[0]*qSkin[35])+9.077304717673634*coeff[0]*qEdge[35]+15.46875*coeff[0]*qSkin[16]+12.65625*coeff[0]*qEdge[16]-8.11898816047911*coeff[0]*qSkin[9]+8.11898816047911*coeff[0]*qEdge[9]; 
  edgeSurf_incr[17] = -(14.160595359570864*coeff[0]*qSkin[37])+9.077304717673634*coeff[0]*qEdge[37]+15.46875*coeff[0]*qSkin[17]+12.65625*coeff[0]*qEdge[17]-8.11898816047911*coeff[0]*qSkin[10]+8.11898816047911*coeff[0]*qEdge[10]; 
  edgeSurf_incr[18] = 6.7082039324993685*coeff[0]*qSkin[44]-6.7082039324993685*coeff[0]*qEdge[44]-8.11898816047911*coeff[0]*qSkin[31]-8.11898816047911*coeff[0]*qEdge[31]+4.6875*coeff[0]*qSkin[18]-4.6875*coeff[0]*qEdge[18]; 
  edgeSurf_incr[19] = 20.34375*coeff[0]*qSkin[19]-0.65625*coeff[0]*qEdge[19]-15.612964114398654*coeff[0]*qSkin[5]-4.72019845319029*coeff[0]*qEdge[5]+4.192627457812107*coeff[0]*qSkin[2]-4.192627457812107*coeff[0]*qEdge[2]; 
  edgeSurf_incr[20] = 15.46875*coeff[0]*qSkin[20]+12.65625*coeff[0]*qEdge[20]-8.118988160479114*coeff[0]*qSkin[12]+8.118988160479114*coeff[0]*qEdge[12]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*qSkin[21]-0.65625*coeff[0]*qEdge[21]-15.612964114398654*coeff[0]*qSkin[6]-4.72019845319029*coeff[0]*qEdge[6]+4.192627457812107*coeff[0]*qSkin[3]-4.192627457812107*coeff[0]*qEdge[3]; 
  edgeSurf_incr[22] = -(8.118988160479114*coeff[0]*qSkin[33])-8.118988160479114*coeff[0]*qEdge[33]+4.6875*coeff[0]*qSkin[22]-4.6875*coeff[0]*qEdge[22]; 
  edgeSurf_incr[23] = 15.46875*coeff[0]*qSkin[23]+12.65625*coeff[0]*qEdge[23]-8.118988160479114*coeff[0]*qSkin[13]+8.118988160479114*coeff[0]*qEdge[13]; 
  edgeSurf_incr[24] = -(8.118988160479114*coeff[0]*qSkin[34])-8.118988160479114*coeff[0]*qEdge[34]+4.6875*coeff[0]*qSkin[24]-4.6875*coeff[0]*qEdge[24]; 
  edgeSurf_incr[25] = 20.34375*coeff[0]*qSkin[25]-0.65625*coeff[0]*qEdge[25]-15.612964114398654*coeff[0]*qSkin[8]-4.72019845319029*coeff[0]*qEdge[8]+4.192627457812107*coeff[0]*qSkin[4]-4.192627457812107*coeff[0]*qEdge[4]; 
  edgeSurf_incr[26] = -(8.118988160479114*coeff[0]*qSkin[36])-8.118988160479114*coeff[0]*qEdge[36]+4.6875*coeff[0]*qSkin[26]-4.6875*coeff[0]*qEdge[26]; 
  edgeSurf_incr[27] = -(8.118988160479114*coeff[0]*qSkin[39])-8.118988160479114*coeff[0]*qEdge[39]+4.6875*coeff[0]*qSkin[27]-4.6875*coeff[0]*qEdge[27]; 
  edgeSurf_incr[28] = 15.46875*coeff[0]*qSkin[28]+12.65625*coeff[0]*qEdge[28]-8.118988160479114*coeff[0]*qSkin[14]+8.118988160479114*coeff[0]*qEdge[14]; 
  edgeSurf_incr[29] = -(8.118988160479114*coeff[0]*qSkin[41])-8.118988160479114*coeff[0]*qEdge[41]+4.6875*coeff[0]*qSkin[29]-4.6875*coeff[0]*qEdge[29]; 
  edgeSurf_incr[30] = -(8.118988160479114*coeff[0]*qSkin[42])-8.118988160479114*coeff[0]*qEdge[42]+4.6875*coeff[0]*qSkin[30]-4.6875*coeff[0]*qEdge[30]; 
  edgeSurf_incr[31] = -(14.160595359570868*coeff[0]*qSkin[44])+9.077304717673634*coeff[0]*qEdge[44]+15.46875*coeff[0]*qSkin[31]+12.65625*coeff[0]*qEdge[31]-8.11898816047911*coeff[0]*qSkin[18]+8.11898816047911*coeff[0]*qEdge[18]; 
  edgeSurf_incr[32] = 20.34375*coeff[0]*qSkin[32]-0.65625*coeff[0]*qEdge[32]-15.61296411439865*coeff[0]*qSkin[15]-4.720198453190292*coeff[0]*qEdge[15]+4.192627457812107*coeff[0]*qSkin[7]-4.192627457812107*coeff[0]*qEdge[7]; 
  edgeSurf_incr[33] = 15.46875*coeff[0]*qSkin[33]+12.65625*coeff[0]*qEdge[33]-8.118988160479114*coeff[0]*qSkin[22]+8.118988160479114*coeff[0]*qEdge[22]; 
  edgeSurf_incr[34] = 15.46875*coeff[0]*qSkin[34]+12.65625*coeff[0]*qEdge[34]-8.118988160479114*coeff[0]*qSkin[24]+8.118988160479114*coeff[0]*qEdge[24]; 
  edgeSurf_incr[35] = 20.34375*coeff[0]*qSkin[35]-0.65625*coeff[0]*qEdge[35]-15.61296411439865*coeff[0]*qSkin[16]-4.720198453190292*coeff[0]*qEdge[16]+4.192627457812107*coeff[0]*qSkin[9]-4.192627457812107*coeff[0]*qEdge[9]; 
  edgeSurf_incr[36] = 15.46875*coeff[0]*qSkin[36]+12.65625*coeff[0]*qEdge[36]-8.118988160479114*coeff[0]*qSkin[26]+8.118988160479114*coeff[0]*qEdge[26]; 
  edgeSurf_incr[37] = 20.34375*coeff[0]*qSkin[37]-0.65625*coeff[0]*qEdge[37]-15.61296411439865*coeff[0]*qSkin[17]-4.720198453190292*coeff[0]*qEdge[17]+4.192627457812107*coeff[0]*qSkin[10]-4.192627457812107*coeff[0]*qEdge[10]; 
  edgeSurf_incr[38] = -(8.118988160479114*coeff[0]*qSkin[45])-8.118988160479114*coeff[0]*qEdge[45]+4.6875*coeff[0]*qSkin[38]-4.6875*coeff[0]*qEdge[38]; 
  edgeSurf_incr[39] = 15.46875*coeff[0]*qSkin[39]+12.65625*coeff[0]*qEdge[39]-8.118988160479114*coeff[0]*qSkin[27]+8.118988160479114*coeff[0]*qEdge[27]; 
  edgeSurf_incr[40] = -(8.118988160479114*coeff[0]*qSkin[46])-8.118988160479114*coeff[0]*qEdge[46]+4.6875*coeff[0]*qSkin[40]-4.6875*coeff[0]*qEdge[40]; 
  edgeSurf_incr[41] = 15.46875*coeff[0]*qSkin[41]+12.65625*coeff[0]*qEdge[41]-8.118988160479114*coeff[0]*qSkin[29]+8.118988160479114*coeff[0]*qEdge[29]; 
  edgeSurf_incr[42] = 15.46875*coeff[0]*qSkin[42]+12.65625*coeff[0]*qEdge[42]-8.118988160479114*coeff[0]*qSkin[30]+8.118988160479114*coeff[0]*qEdge[30]; 
  edgeSurf_incr[43] = -(8.118988160479114*coeff[0]*qSkin[47])-8.118988160479114*coeff[0]*qEdge[47]+4.6875*coeff[0]*qSkin[43]-4.6875*coeff[0]*qEdge[43]; 
  edgeSurf_incr[44] = 20.34375*coeff[0]*qSkin[44]-0.65625*coeff[0]*qEdge[44]-15.612964114398654*coeff[0]*qSkin[31]-4.72019845319029*coeff[0]*qEdge[31]+4.192627457812107*coeff[0]*qSkin[18]-4.192627457812107*coeff[0]*qEdge[18]; 
  edgeSurf_incr[45] = 15.46875*coeff[0]*qSkin[45]+12.65625*coeff[0]*qEdge[45]-8.118988160479114*coeff[0]*qSkin[38]+8.118988160479114*coeff[0]*qEdge[38]; 
  edgeSurf_incr[46] = 15.46875*coeff[0]*qSkin[46]+12.65625*coeff[0]*qEdge[46]-8.118988160479114*coeff[0]*qSkin[40]+8.118988160479114*coeff[0]*qEdge[40]; 
  edgeSurf_incr[47] = 15.46875*coeff[0]*qSkin[47]+12.65625*coeff[0]*qEdge[47]-8.118988160479114*coeff[0]*qSkin[43]+8.118988160479114*coeff[0]*qEdge[43]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*qSkin[11]+2.8125*coeff[0]*qSkin[1]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[0]*qSkin[19]+2.8125*coeff[0]*qSkin[5]; 
  boundSurf_incr[6] = 5.083290641897235*coeff[0]*qSkin[21]+2.8125*coeff[0]*qSkin[6]; 
  boundSurf_incr[8] = 5.083290641897235*coeff[0]*qSkin[25]+2.8125*coeff[0]*qSkin[8]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*qSkin[11]+10.892765661208358*coeff[0]*qSkin[1]; 
  boundSurf_incr[15] = 5.083290641897234*coeff[0]*qSkin[32]+2.8125*coeff[0]*qSkin[15]; 
  boundSurf_incr[16] = 5.083290641897234*coeff[0]*qSkin[35]+2.8125*coeff[0]*qSkin[16]; 
  boundSurf_incr[17] = 5.083290641897234*coeff[0]*qSkin[37]+2.8125*coeff[0]*qSkin[17]; 
  boundSurf_incr[19] = 19.6875*coeff[0]*qSkin[19]+10.892765661208362*coeff[0]*qSkin[5]; 
  boundSurf_incr[20] = 2.8125*coeff[0]*qSkin[20]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*qSkin[21]+10.892765661208362*coeff[0]*qSkin[6]; 
  boundSurf_incr[23] = 2.8125*coeff[0]*qSkin[23]; 
  boundSurf_incr[25] = 19.6875*coeff[0]*qSkin[25]+10.892765661208362*coeff[0]*qSkin[8]; 
  boundSurf_incr[28] = 2.8125*coeff[0]*qSkin[28]; 
  boundSurf_incr[31] = 5.083290641897235*coeff[0]*qSkin[44]+2.8125*coeff[0]*qSkin[31]; 
  boundSurf_incr[32] = 19.6875*coeff[0]*qSkin[32]+10.892765661208358*coeff[0]*qSkin[15]; 
  boundSurf_incr[33] = 2.8125*coeff[0]*qSkin[33]; 
  boundSurf_incr[34] = 2.8125*coeff[0]*qSkin[34]; 
  boundSurf_incr[35] = 19.6875*coeff[0]*qSkin[35]+10.892765661208358*coeff[0]*qSkin[16]; 
  boundSurf_incr[36] = 2.8125*coeff[0]*qSkin[36]; 
  boundSurf_incr[37] = 19.6875*coeff[0]*qSkin[37]+10.892765661208358*coeff[0]*qSkin[17]; 
  boundSurf_incr[39] = 2.8125*coeff[0]*qSkin[39]; 
  boundSurf_incr[41] = 2.8125*coeff[0]*qSkin[41]; 
  boundSurf_incr[42] = 2.8125*coeff[0]*qSkin[42]; 
  boundSurf_incr[44] = 19.6875*coeff[0]*qSkin[44]+10.892765661208362*coeff[0]*qSkin[31]; 
  boundSurf_incr[45] = 2.8125*coeff[0]*qSkin[45]; 
  boundSurf_incr[46] = 2.8125*coeff[0]*qSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[0]*qSkin[47]; 

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

