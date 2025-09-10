#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_boundary_surfx_1x3v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.160595359570864*coeff[0]*fSkin[11]-9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*fSkin[19]-6.7082039324993685*coeff[0]*fEdge[19]+8.11898816047911*coeff[0]*fSkin[5]+8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*fSkin[21]-6.7082039324993685*coeff[0]*fEdge[21]+8.11898816047911*coeff[0]*fSkin[6]+8.11898816047911*coeff[0]*fEdge[6]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 6.7082039324993685*coeff[0]*fSkin[25]-6.7082039324993685*coeff[0]*fEdge[25]+8.11898816047911*coeff[0]*fSkin[8]+8.11898816047911*coeff[0]*fEdge[8]+4.6875*coeff[0]*fSkin[4]-4.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 14.160595359570868*coeff[0]*fSkin[19]-9.077304717673634*coeff[0]*fEdge[19]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 14.160595359570868*coeff[0]*fSkin[21]-9.077304717673634*coeff[0]*fEdge[21]+15.46875*coeff[0]*fSkin[6]+12.65625*coeff[0]*fEdge[6]+8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[0]*fSkin[32]-6.708203932499369*coeff[0]*fEdge[32]+8.11898816047911*coeff[0]*fSkin[15]+8.11898816047911*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[7]-4.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 14.160595359570868*coeff[0]*fSkin[25]-9.077304717673634*coeff[0]*fEdge[25]+15.46875*coeff[0]*fSkin[8]+12.65625*coeff[0]*fEdge[8]+8.11898816047911*coeff[0]*fSkin[4]-8.11898816047911*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*fSkin[35]-6.708203932499369*coeff[0]*fEdge[35]+8.11898816047911*coeff[0]*fSkin[16]+8.11898816047911*coeff[0]*fEdge[16]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[0]*fSkin[37]-6.708203932499369*coeff[0]*fEdge[37]+8.11898816047911*coeff[0]*fSkin[17]+8.11898816047911*coeff[0]*fEdge[17]+4.6875*coeff[0]*fSkin[10]-4.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = 8.118988160479114*coeff[0]*fSkin[20]+8.118988160479114*coeff[0]*fEdge[20]+4.6875*coeff[0]*fSkin[12]-4.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 8.118988160479114*coeff[0]*fSkin[23]+8.118988160479114*coeff[0]*fEdge[23]+4.6875*coeff[0]*fSkin[13]-4.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = 8.118988160479114*coeff[0]*fSkin[28]+8.118988160479114*coeff[0]*fEdge[28]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 14.160595359570864*coeff[0]*fSkin[32]-9.077304717673634*coeff[0]*fEdge[32]+15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]+8.11898816047911*coeff[0]*fSkin[7]-8.11898816047911*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = 14.160595359570864*coeff[0]*fSkin[35]-9.077304717673634*coeff[0]*fEdge[35]+15.46875*coeff[0]*fSkin[16]+12.65625*coeff[0]*fEdge[16]+8.11898816047911*coeff[0]*fSkin[9]-8.11898816047911*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = 14.160595359570864*coeff[0]*fSkin[37]-9.077304717673634*coeff[0]*fEdge[37]+15.46875*coeff[0]*fSkin[17]+12.65625*coeff[0]*fEdge[17]+8.11898816047911*coeff[0]*fSkin[10]-8.11898816047911*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = 6.7082039324993685*coeff[0]*fSkin[44]-6.7082039324993685*coeff[0]*fEdge[44]+8.11898816047911*coeff[0]*fSkin[31]+8.11898816047911*coeff[0]*fEdge[31]+4.6875*coeff[0]*fSkin[18]-4.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 20.34375*coeff[0]*fSkin[19]-0.65625*coeff[0]*fEdge[19]+15.612964114398654*coeff[0]*fSkin[5]+4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = 15.46875*coeff[0]*fSkin[20]+12.65625*coeff[0]*fEdge[20]+8.118988160479114*coeff[0]*fSkin[12]-8.118988160479114*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*fSkin[21]-0.65625*coeff[0]*fEdge[21]+15.612964114398654*coeff[0]*fSkin[6]+4.72019845319029*coeff[0]*fEdge[6]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = 8.118988160479114*coeff[0]*fSkin[33]+8.118988160479114*coeff[0]*fEdge[33]+4.6875*coeff[0]*fSkin[22]-4.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 15.46875*coeff[0]*fSkin[23]+12.65625*coeff[0]*fEdge[23]+8.118988160479114*coeff[0]*fSkin[13]-8.118988160479114*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = 8.118988160479114*coeff[0]*fSkin[34]+8.118988160479114*coeff[0]*fEdge[34]+4.6875*coeff[0]*fSkin[24]-4.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = 20.34375*coeff[0]*fSkin[25]-0.65625*coeff[0]*fEdge[25]+15.612964114398654*coeff[0]*fSkin[8]+4.72019845319029*coeff[0]*fEdge[8]+4.192627457812107*coeff[0]*fSkin[4]-4.192627457812107*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = 8.118988160479114*coeff[0]*fSkin[36]+8.118988160479114*coeff[0]*fEdge[36]+4.6875*coeff[0]*fSkin[26]-4.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 8.118988160479114*coeff[0]*fSkin[39]+8.118988160479114*coeff[0]*fEdge[39]+4.6875*coeff[0]*fSkin[27]-4.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = 15.46875*coeff[0]*fSkin[28]+12.65625*coeff[0]*fEdge[28]+8.118988160479114*coeff[0]*fSkin[14]-8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = 8.118988160479114*coeff[0]*fSkin[41]+8.118988160479114*coeff[0]*fEdge[41]+4.6875*coeff[0]*fSkin[29]-4.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = 8.118988160479114*coeff[0]*fSkin[42]+8.118988160479114*coeff[0]*fEdge[42]+4.6875*coeff[0]*fSkin[30]-4.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 14.160595359570868*coeff[0]*fSkin[44]-9.077304717673634*coeff[0]*fEdge[44]+15.46875*coeff[0]*fSkin[31]+12.65625*coeff[0]*fEdge[31]+8.11898816047911*coeff[0]*fSkin[18]-8.11898816047911*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = 20.34375*coeff[0]*fSkin[32]-0.65625*coeff[0]*fEdge[32]+15.61296411439865*coeff[0]*fSkin[15]+4.720198453190292*coeff[0]*fEdge[15]+4.192627457812107*coeff[0]*fSkin[7]-4.192627457812107*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = 15.46875*coeff[0]*fSkin[33]+12.65625*coeff[0]*fEdge[33]+8.118988160479114*coeff[0]*fSkin[22]-8.118988160479114*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = 15.46875*coeff[0]*fSkin[34]+12.65625*coeff[0]*fEdge[34]+8.118988160479114*coeff[0]*fSkin[24]-8.118988160479114*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = 20.34375*coeff[0]*fSkin[35]-0.65625*coeff[0]*fEdge[35]+15.61296411439865*coeff[0]*fSkin[16]+4.720198453190292*coeff[0]*fEdge[16]+4.192627457812107*coeff[0]*fSkin[9]-4.192627457812107*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = 15.46875*coeff[0]*fSkin[36]+12.65625*coeff[0]*fEdge[36]+8.118988160479114*coeff[0]*fSkin[26]-8.118988160479114*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = 20.34375*coeff[0]*fSkin[37]-0.65625*coeff[0]*fEdge[37]+15.61296411439865*coeff[0]*fSkin[17]+4.720198453190292*coeff[0]*fEdge[17]+4.192627457812107*coeff[0]*fSkin[10]-4.192627457812107*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = 8.118988160479114*coeff[0]*fSkin[45]+8.118988160479114*coeff[0]*fEdge[45]+4.6875*coeff[0]*fSkin[38]-4.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 15.46875*coeff[0]*fSkin[39]+12.65625*coeff[0]*fEdge[39]+8.118988160479114*coeff[0]*fSkin[27]-8.118988160479114*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = 8.118988160479114*coeff[0]*fSkin[46]+8.118988160479114*coeff[0]*fEdge[46]+4.6875*coeff[0]*fSkin[40]-4.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = 15.46875*coeff[0]*fSkin[41]+12.65625*coeff[0]*fEdge[41]+8.118988160479114*coeff[0]*fSkin[29]-8.118988160479114*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = 15.46875*coeff[0]*fSkin[42]+12.65625*coeff[0]*fEdge[42]+8.118988160479114*coeff[0]*fSkin[30]-8.118988160479114*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = 8.118988160479114*coeff[0]*fSkin[47]+8.118988160479114*coeff[0]*fEdge[47]+4.6875*coeff[0]*fSkin[43]-4.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = 20.34375*coeff[0]*fSkin[44]-0.65625*coeff[0]*fEdge[44]+15.612964114398654*coeff[0]*fSkin[31]+4.72019845319029*coeff[0]*fEdge[31]+4.192627457812107*coeff[0]*fSkin[18]-4.192627457812107*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = 15.46875*coeff[0]*fSkin[45]+12.65625*coeff[0]*fEdge[45]+8.118988160479114*coeff[0]*fSkin[38]-8.118988160479114*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 15.46875*coeff[0]*fSkin[46]+12.65625*coeff[0]*fEdge[46]+8.118988160479114*coeff[0]*fSkin[40]-8.118988160479114*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = 15.46875*coeff[0]*fSkin[47]+12.65625*coeff[0]*fEdge[47]+8.118988160479114*coeff[0]*fSkin[43]-8.118988160479114*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[11]; 
  boundSurf_incr[5] = 2.8125*coeff[0]*fSkin[5]-5.083290641897235*coeff[0]*fSkin[19]; 
  boundSurf_incr[6] = 2.8125*coeff[0]*fSkin[6]-5.083290641897235*coeff[0]*fSkin[21]; 
  boundSurf_incr[8] = 2.8125*coeff[0]*fSkin[8]-5.083290641897235*coeff[0]*fSkin[25]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]-10.892765661208358*coeff[0]*fSkin[1]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*fSkin[15]-5.083290641897234*coeff[0]*fSkin[32]; 
  boundSurf_incr[16] = 2.8125*coeff[0]*fSkin[16]-5.083290641897234*coeff[0]*fSkin[35]; 
  boundSurf_incr[17] = 2.8125*coeff[0]*fSkin[17]-5.083290641897234*coeff[0]*fSkin[37]; 
  boundSurf_incr[19] = 19.6875*coeff[0]*fSkin[19]-10.892765661208362*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = 2.8125*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*fSkin[21]-10.892765661208362*coeff[0]*fSkin[6]; 
  boundSurf_incr[23] = 2.8125*coeff[0]*fSkin[23]; 
  boundSurf_incr[25] = 19.6875*coeff[0]*fSkin[25]-10.892765661208362*coeff[0]*fSkin[8]; 
  boundSurf_incr[28] = 2.8125*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = 2.8125*coeff[0]*fSkin[31]-5.083290641897235*coeff[0]*fSkin[44]; 
  boundSurf_incr[32] = 19.6875*coeff[0]*fSkin[32]-10.892765661208358*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = 2.8125*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = 2.8125*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 19.6875*coeff[0]*fSkin[35]-10.892765661208358*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = 2.8125*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = 19.6875*coeff[0]*fSkin[37]-10.892765661208358*coeff[0]*fSkin[17]; 
  boundSurf_incr[39] = 2.8125*coeff[0]*fSkin[39]; 
  boundSurf_incr[41] = 2.8125*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = 2.8125*coeff[0]*fSkin[42]; 
  boundSurf_incr[44] = 19.6875*coeff[0]*fSkin[44]-10.892765661208362*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = 2.8125*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = 2.8125*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[0]*fSkin[47]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(14.160595359570864*coeff[0]*fSkin[11])+9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*fSkin[19]-6.7082039324993685*coeff[0]*fEdge[19]-8.11898816047911*coeff[0]*fSkin[5]-8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*fSkin[21]-6.7082039324993685*coeff[0]*fEdge[21]-8.11898816047911*coeff[0]*fSkin[6]-8.11898816047911*coeff[0]*fEdge[6]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 6.7082039324993685*coeff[0]*fSkin[25]-6.7082039324993685*coeff[0]*fEdge[25]-8.11898816047911*coeff[0]*fSkin[8]-8.11898816047911*coeff[0]*fEdge[8]+4.6875*coeff[0]*fSkin[4]-4.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = -(14.160595359570868*coeff[0]*fSkin[19])+9.077304717673634*coeff[0]*fEdge[19]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = -(14.160595359570868*coeff[0]*fSkin[21])+9.077304717673634*coeff[0]*fEdge[21]+15.46875*coeff[0]*fSkin[6]+12.65625*coeff[0]*fEdge[6]-8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[0]*fSkin[32]-6.708203932499369*coeff[0]*fEdge[32]-8.11898816047911*coeff[0]*fSkin[15]-8.11898816047911*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[7]-4.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = -(14.160595359570868*coeff[0]*fSkin[25])+9.077304717673634*coeff[0]*fEdge[25]+15.46875*coeff[0]*fSkin[8]+12.65625*coeff[0]*fEdge[8]-8.11898816047911*coeff[0]*fSkin[4]+8.11898816047911*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*fSkin[35]-6.708203932499369*coeff[0]*fEdge[35]-8.11898816047911*coeff[0]*fSkin[16]-8.11898816047911*coeff[0]*fEdge[16]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[0]*fSkin[37]-6.708203932499369*coeff[0]*fEdge[37]-8.11898816047911*coeff[0]*fSkin[17]-8.11898816047911*coeff[0]*fEdge[17]+4.6875*coeff[0]*fSkin[10]-4.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = -(8.118988160479114*coeff[0]*fSkin[20])-8.118988160479114*coeff[0]*fEdge[20]+4.6875*coeff[0]*fSkin[12]-4.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = -(8.118988160479114*coeff[0]*fSkin[23])-8.118988160479114*coeff[0]*fEdge[23]+4.6875*coeff[0]*fSkin[13]-4.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = -(8.118988160479114*coeff[0]*fSkin[28])-8.118988160479114*coeff[0]*fEdge[28]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = -(14.160595359570864*coeff[0]*fSkin[32])+9.077304717673634*coeff[0]*fEdge[32]+15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]-8.11898816047911*coeff[0]*fSkin[7]+8.11898816047911*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = -(14.160595359570864*coeff[0]*fSkin[35])+9.077304717673634*coeff[0]*fEdge[35]+15.46875*coeff[0]*fSkin[16]+12.65625*coeff[0]*fEdge[16]-8.11898816047911*coeff[0]*fSkin[9]+8.11898816047911*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = -(14.160595359570864*coeff[0]*fSkin[37])+9.077304717673634*coeff[0]*fEdge[37]+15.46875*coeff[0]*fSkin[17]+12.65625*coeff[0]*fEdge[17]-8.11898816047911*coeff[0]*fSkin[10]+8.11898816047911*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = 6.7082039324993685*coeff[0]*fSkin[44]-6.7082039324993685*coeff[0]*fEdge[44]-8.11898816047911*coeff[0]*fSkin[31]-8.11898816047911*coeff[0]*fEdge[31]+4.6875*coeff[0]*fSkin[18]-4.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 20.34375*coeff[0]*fSkin[19]-0.65625*coeff[0]*fEdge[19]-15.612964114398654*coeff[0]*fSkin[5]-4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = 15.46875*coeff[0]*fSkin[20]+12.65625*coeff[0]*fEdge[20]-8.118988160479114*coeff[0]*fSkin[12]+8.118988160479114*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*fSkin[21]-0.65625*coeff[0]*fEdge[21]-15.612964114398654*coeff[0]*fSkin[6]-4.72019845319029*coeff[0]*fEdge[6]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = -(8.118988160479114*coeff[0]*fSkin[33])-8.118988160479114*coeff[0]*fEdge[33]+4.6875*coeff[0]*fSkin[22]-4.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 15.46875*coeff[0]*fSkin[23]+12.65625*coeff[0]*fEdge[23]-8.118988160479114*coeff[0]*fSkin[13]+8.118988160479114*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = -(8.118988160479114*coeff[0]*fSkin[34])-8.118988160479114*coeff[0]*fEdge[34]+4.6875*coeff[0]*fSkin[24]-4.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = 20.34375*coeff[0]*fSkin[25]-0.65625*coeff[0]*fEdge[25]-15.612964114398654*coeff[0]*fSkin[8]-4.72019845319029*coeff[0]*fEdge[8]+4.192627457812107*coeff[0]*fSkin[4]-4.192627457812107*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = -(8.118988160479114*coeff[0]*fSkin[36])-8.118988160479114*coeff[0]*fEdge[36]+4.6875*coeff[0]*fSkin[26]-4.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = -(8.118988160479114*coeff[0]*fSkin[39])-8.118988160479114*coeff[0]*fEdge[39]+4.6875*coeff[0]*fSkin[27]-4.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = 15.46875*coeff[0]*fSkin[28]+12.65625*coeff[0]*fEdge[28]-8.118988160479114*coeff[0]*fSkin[14]+8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = -(8.118988160479114*coeff[0]*fSkin[41])-8.118988160479114*coeff[0]*fEdge[41]+4.6875*coeff[0]*fSkin[29]-4.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = -(8.118988160479114*coeff[0]*fSkin[42])-8.118988160479114*coeff[0]*fEdge[42]+4.6875*coeff[0]*fSkin[30]-4.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = -(14.160595359570868*coeff[0]*fSkin[44])+9.077304717673634*coeff[0]*fEdge[44]+15.46875*coeff[0]*fSkin[31]+12.65625*coeff[0]*fEdge[31]-8.11898816047911*coeff[0]*fSkin[18]+8.11898816047911*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = 20.34375*coeff[0]*fSkin[32]-0.65625*coeff[0]*fEdge[32]-15.61296411439865*coeff[0]*fSkin[15]-4.720198453190292*coeff[0]*fEdge[15]+4.192627457812107*coeff[0]*fSkin[7]-4.192627457812107*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = 15.46875*coeff[0]*fSkin[33]+12.65625*coeff[0]*fEdge[33]-8.118988160479114*coeff[0]*fSkin[22]+8.118988160479114*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = 15.46875*coeff[0]*fSkin[34]+12.65625*coeff[0]*fEdge[34]-8.118988160479114*coeff[0]*fSkin[24]+8.118988160479114*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = 20.34375*coeff[0]*fSkin[35]-0.65625*coeff[0]*fEdge[35]-15.61296411439865*coeff[0]*fSkin[16]-4.720198453190292*coeff[0]*fEdge[16]+4.192627457812107*coeff[0]*fSkin[9]-4.192627457812107*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = 15.46875*coeff[0]*fSkin[36]+12.65625*coeff[0]*fEdge[36]-8.118988160479114*coeff[0]*fSkin[26]+8.118988160479114*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = 20.34375*coeff[0]*fSkin[37]-0.65625*coeff[0]*fEdge[37]-15.61296411439865*coeff[0]*fSkin[17]-4.720198453190292*coeff[0]*fEdge[17]+4.192627457812107*coeff[0]*fSkin[10]-4.192627457812107*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = -(8.118988160479114*coeff[0]*fSkin[45])-8.118988160479114*coeff[0]*fEdge[45]+4.6875*coeff[0]*fSkin[38]-4.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 15.46875*coeff[0]*fSkin[39]+12.65625*coeff[0]*fEdge[39]-8.118988160479114*coeff[0]*fSkin[27]+8.118988160479114*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = -(8.118988160479114*coeff[0]*fSkin[46])-8.118988160479114*coeff[0]*fEdge[46]+4.6875*coeff[0]*fSkin[40]-4.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = 15.46875*coeff[0]*fSkin[41]+12.65625*coeff[0]*fEdge[41]-8.118988160479114*coeff[0]*fSkin[29]+8.118988160479114*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = 15.46875*coeff[0]*fSkin[42]+12.65625*coeff[0]*fEdge[42]-8.118988160479114*coeff[0]*fSkin[30]+8.118988160479114*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = -(8.118988160479114*coeff[0]*fSkin[47])-8.118988160479114*coeff[0]*fEdge[47]+4.6875*coeff[0]*fSkin[43]-4.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = 20.34375*coeff[0]*fSkin[44]-0.65625*coeff[0]*fEdge[44]-15.612964114398654*coeff[0]*fSkin[31]-4.72019845319029*coeff[0]*fEdge[31]+4.192627457812107*coeff[0]*fSkin[18]-4.192627457812107*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = 15.46875*coeff[0]*fSkin[45]+12.65625*coeff[0]*fEdge[45]-8.118988160479114*coeff[0]*fSkin[38]+8.118988160479114*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 15.46875*coeff[0]*fSkin[46]+12.65625*coeff[0]*fEdge[46]-8.118988160479114*coeff[0]*fSkin[40]+8.118988160479114*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = 15.46875*coeff[0]*fSkin[47]+12.65625*coeff[0]*fEdge[47]-8.118988160479114*coeff[0]*fSkin[43]+8.118988160479114*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[11]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[0]*fSkin[19]+2.8125*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 5.083290641897235*coeff[0]*fSkin[21]+2.8125*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 5.083290641897235*coeff[0]*fSkin[25]+2.8125*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]+10.892765661208358*coeff[0]*fSkin[1]; 
  boundSurf_incr[15] = 5.083290641897234*coeff[0]*fSkin[32]+2.8125*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 5.083290641897234*coeff[0]*fSkin[35]+2.8125*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = 5.083290641897234*coeff[0]*fSkin[37]+2.8125*coeff[0]*fSkin[17]; 
  boundSurf_incr[19] = 19.6875*coeff[0]*fSkin[19]+10.892765661208362*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = 2.8125*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*fSkin[21]+10.892765661208362*coeff[0]*fSkin[6]; 
  boundSurf_incr[23] = 2.8125*coeff[0]*fSkin[23]; 
  boundSurf_incr[25] = 19.6875*coeff[0]*fSkin[25]+10.892765661208362*coeff[0]*fSkin[8]; 
  boundSurf_incr[28] = 2.8125*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = 5.083290641897235*coeff[0]*fSkin[44]+2.8125*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 19.6875*coeff[0]*fSkin[32]+10.892765661208358*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = 2.8125*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = 2.8125*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 19.6875*coeff[0]*fSkin[35]+10.892765661208358*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = 2.8125*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = 19.6875*coeff[0]*fSkin[37]+10.892765661208358*coeff[0]*fSkin[17]; 
  boundSurf_incr[39] = 2.8125*coeff[0]*fSkin[39]; 
  boundSurf_incr[41] = 2.8125*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = 2.8125*coeff[0]*fSkin[42]; 
  boundSurf_incr[44] = 19.6875*coeff[0]*fSkin[44]+10.892765661208362*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = 2.8125*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = 2.8125*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[0]*fSkin[47]; 

  }

  out[0] += -(1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac); 
  out[1] += -(1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac); 
  out[2] += -(1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac); 
  out[3] += -(1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac); 
  out[4] += -(1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac); 
  out[5] += -(1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac); 
  out[6] += -(1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac); 
  out[7] += -(1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac); 
  out[8] += -(1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac); 
  out[9] += -(1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac); 
  out[10] += -(1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac); 
  out[11] += -(1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac); 
  out[12] += -(1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac); 
  out[13] += -(1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac); 
  out[14] += -(1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac); 
  out[15] += -(1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac); 
  out[16] += -(1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac); 
  out[17] += -(1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac); 
  out[18] += -(1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac); 
  out[19] += -(1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac); 
  out[20] += -(1.0*(vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac); 
  out[21] += -(1.0*(vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac); 
  out[22] += -(1.0*(vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac); 
  out[23] += -(1.0*(vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac); 
  out[24] += -(1.0*(vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac); 
  out[25] += -(1.0*(vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac); 
  out[26] += -(1.0*(vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac); 
  out[27] += -(1.0*(vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac); 
  out[28] += -(1.0*(vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac); 
  out[29] += -(1.0*(vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac); 
  out[30] += -(1.0*(vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac); 
  out[31] += -(1.0*(vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac); 
  out[32] += -(1.0*(vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac); 
  out[33] += -(1.0*(vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac); 
  out[34] += -(1.0*(vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac); 
  out[35] += -(1.0*(vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac); 
  out[36] += -(1.0*(vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac); 
  out[37] += -(1.0*(vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac); 
  out[38] += -(1.0*(vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac); 
  out[39] += -(1.0*(vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac); 
  out[40] += -(1.0*(vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac); 
  out[41] += -(1.0*(vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac); 
  out[42] += -(1.0*(vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac); 
  out[43] += -(1.0*(vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac); 
  out[44] += -(1.0*(vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac); 
  out[45] += -(1.0*(vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac); 
  out[46] += -(1.0*(vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac); 
  out[47] += -(1.0*(vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac); 

  return 0.;
}

