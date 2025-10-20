#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  double vol_incr[27] = {0.0}; 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[1]*fSkin[8]-6.708203932499369*coeff[1]*fEdge[8]+8.11898816047911*coeff[1]*fSkin[2]+8.11898816047911*coeff[1]*fEdge[2]+4.6875*fSkin[0]*coeff[1]-4.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[1]*fSkin[12]-6.7082039324993685*coeff[1]*fEdge[12]+8.11898816047911*coeff[1]*fSkin[4]+8.11898816047911*coeff[1]*fEdge[4]+4.6875*coeff[1]*fSkin[1]-4.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 14.160595359570864*coeff[1]*fSkin[8]-9.077304717673634*coeff[1]*fEdge[8]+15.46875*coeff[1]*fSkin[2]+12.65625*coeff[1]*fEdge[2]+8.11898816047911*fSkin[0]*coeff[1]-8.11898816047911*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[1]*fSkin[14]-6.7082039324993685*coeff[1]*fEdge[14]+8.11898816047911*coeff[1]*fSkin[6]+8.11898816047911*coeff[1]*fEdge[6]+4.6875*coeff[1]*fSkin[3]-4.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 14.160595359570868*coeff[1]*fSkin[12]-9.077304717673634*coeff[1]*fEdge[12]+15.46875*coeff[1]*fSkin[4]+12.65625*coeff[1]*fEdge[4]+8.11898816047911*coeff[1]*fSkin[1]-8.11898816047911*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = 6.708203932499369*coeff[1]*fSkin[18]-6.708203932499369*coeff[1]*fEdge[18]+8.11898816047911*coeff[1]*fSkin[10]+8.11898816047911*coeff[1]*fEdge[10]+4.6875*coeff[1]*fSkin[5]-4.6875*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = 14.160595359570868*coeff[1]*fSkin[14]-9.077304717673634*coeff[1]*fEdge[14]+15.46875*coeff[1]*fSkin[6]+12.65625*coeff[1]*fEdge[6]+8.11898816047911*coeff[1]*fSkin[3]-8.11898816047911*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[1]*fSkin[20]-6.708203932499369*coeff[1]*fEdge[20]+8.118988160479114*coeff[1]*fSkin[11]+8.118988160479114*coeff[1]*fEdge[11]+4.6875*coeff[1]*fSkin[7]-4.6875*coeff[1]*fEdge[7]; 
  edgeSurf_incr[8] = 20.34375*coeff[1]*fSkin[8]-0.65625*coeff[1]*fEdge[8]+15.61296411439865*coeff[1]*fSkin[2]+4.720198453190292*coeff[1]*fEdge[2]+4.192627457812107*fSkin[0]*coeff[1]-4.192627457812107*fEdge[0]*coeff[1]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[1]*fSkin[22]-6.708203932499369*coeff[1]*fEdge[22]+8.118988160479114*coeff[1]*fSkin[16]+8.118988160479114*coeff[1]*fEdge[16]+4.6875*coeff[1]*fSkin[9]-4.6875*coeff[1]*fEdge[9]; 
  edgeSurf_incr[10] = 14.160595359570864*coeff[1]*fSkin[18]-9.077304717673634*coeff[1]*fEdge[18]+15.46875*coeff[1]*fSkin[10]+12.65625*coeff[1]*fEdge[10]+8.11898816047911*coeff[1]*fSkin[5]-8.11898816047911*coeff[1]*fEdge[5]; 
  edgeSurf_incr[11] = 14.160595359570868*coeff[1]*fSkin[20]-9.077304717673634*coeff[1]*fEdge[20]+15.46875*coeff[1]*fSkin[11]+12.65625*coeff[1]*fEdge[11]+8.118988160479114*coeff[1]*fSkin[7]-8.118988160479114*coeff[1]*fEdge[7]; 
  edgeSurf_incr[12] = 20.34375*coeff[1]*fSkin[12]-0.65625*coeff[1]*fEdge[12]+15.612964114398654*coeff[1]*fSkin[4]+4.72019845319029*coeff[1]*fEdge[4]+4.192627457812107*coeff[1]*fSkin[1]-4.192627457812107*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = 6.7082039324993685*coeff[1]*fSkin[23]-6.7082039324993685*coeff[1]*fEdge[23]+8.118988160479114*coeff[1]*fSkin[17]+8.118988160479114*coeff[1]*fEdge[17]+4.6875*coeff[1]*fSkin[13]-4.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = 20.34375*coeff[1]*fSkin[14]-0.65625*coeff[1]*fEdge[14]+15.612964114398654*coeff[1]*fSkin[6]+4.72019845319029*coeff[1]*fEdge[6]+4.192627457812107*coeff[1]*fSkin[3]-4.192627457812107*coeff[1]*fEdge[3]; 
  edgeSurf_incr[15] = 6.7082039324993685*coeff[1]*fSkin[25]-6.7082039324993685*coeff[1]*fEdge[25]+8.118988160479114*coeff[1]*fSkin[19]+8.118988160479114*coeff[1]*fEdge[19]+4.6875*coeff[1]*fSkin[15]-4.6875*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = 14.160595359570868*coeff[1]*fSkin[22]-9.077304717673634*coeff[1]*fEdge[22]+15.46875*coeff[1]*fSkin[16]+12.65625*coeff[1]*fEdge[16]+8.118988160479114*coeff[1]*fSkin[9]-8.118988160479114*coeff[1]*fEdge[9]; 
  edgeSurf_incr[17] = 14.160595359570864*coeff[1]*fSkin[23]-9.077304717673634*coeff[1]*fEdge[23]+15.46875*coeff[1]*fSkin[17]+12.65625*coeff[1]*fEdge[17]+8.118988160479114*coeff[1]*fSkin[13]-8.118988160479114*coeff[1]*fEdge[13]; 
  edgeSurf_incr[18] = 20.34375*coeff[1]*fSkin[18]-0.65625*coeff[1]*fEdge[18]+15.61296411439865*coeff[1]*fSkin[10]+4.720198453190292*coeff[1]*fEdge[10]+4.192627457812107*coeff[1]*fSkin[5]-4.192627457812107*coeff[1]*fEdge[5]; 
  edgeSurf_incr[19] = 14.160595359570864*coeff[1]*fSkin[25]-9.077304717673634*coeff[1]*fEdge[25]+15.46875*coeff[1]*fSkin[19]+12.65625*coeff[1]*fEdge[19]+8.118988160479114*coeff[1]*fSkin[15]-8.118988160479114*coeff[1]*fEdge[15]; 
  edgeSurf_incr[20] = 20.34375*coeff[1]*fSkin[20]-0.65625*coeff[1]*fEdge[20]+15.612964114398654*coeff[1]*fSkin[11]+4.72019845319029*coeff[1]*fEdge[11]+4.192627457812107*coeff[1]*fSkin[7]-4.192627457812107*coeff[1]*fEdge[7]; 
  edgeSurf_incr[21] = 6.708203932499369*coeff[1]*fSkin[26]-6.708203932499369*coeff[1]*fEdge[26]+8.11898816047911*coeff[1]*fSkin[24]+8.11898816047911*coeff[1]*fEdge[24]+4.6875*coeff[1]*fSkin[21]-4.6875*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = 20.34375*coeff[1]*fSkin[22]-0.65625*coeff[1]*fEdge[22]+15.612964114398654*coeff[1]*fSkin[16]+4.72019845319029*coeff[1]*fEdge[16]+4.192627457812107*coeff[1]*fSkin[9]-4.192627457812107*coeff[1]*fEdge[9]; 
  edgeSurf_incr[23] = 20.34375*coeff[1]*fSkin[23]-0.65625*coeff[1]*fEdge[23]+15.61296411439865*coeff[1]*fSkin[17]+4.720198453190292*coeff[1]*fEdge[17]+4.192627457812107*coeff[1]*fSkin[13]-4.192627457812107*coeff[1]*fEdge[13]; 
  edgeSurf_incr[24] = 14.160595359570864*coeff[1]*fSkin[26]-9.077304717673634*coeff[1]*fEdge[26]+15.46875*coeff[1]*fSkin[24]+12.65625*coeff[1]*fEdge[24]+8.11898816047911*coeff[1]*fSkin[21]-8.11898816047911*coeff[1]*fEdge[21]; 
  edgeSurf_incr[25] = 20.34375*coeff[1]*fSkin[25]-0.65625*coeff[1]*fEdge[25]+15.61296411439865*coeff[1]*fSkin[19]+4.720198453190292*coeff[1]*fEdge[19]+4.192627457812107*coeff[1]*fSkin[15]-4.192627457812107*coeff[1]*fEdge[15]; 
  edgeSurf_incr[26] = 20.34375*coeff[1]*fSkin[26]-0.65625*coeff[1]*fEdge[26]+15.61296411439865*coeff[1]*fSkin[24]+4.720198453190292*coeff[1]*fEdge[24]+4.192627457812107*coeff[1]*fSkin[21]-4.192627457812107*coeff[1]*fEdge[21]; 

  boundSurf_incr[2] = 2.8125*coeff[1]*fSkin[2]-5.083290641897234*coeff[1]*fSkin[8]; 
  boundSurf_incr[4] = 2.8125*coeff[1]*fSkin[4]-5.083290641897235*coeff[1]*fSkin[12]; 
  boundSurf_incr[6] = 2.8125*coeff[1]*fSkin[6]-5.083290641897235*coeff[1]*fSkin[14]; 
  boundSurf_incr[8] = 19.6875*coeff[1]*fSkin[8]-10.892765661208358*coeff[1]*fSkin[2]; 
  boundSurf_incr[10] = 2.8125*coeff[1]*fSkin[10]-5.083290641897234*coeff[1]*fSkin[18]; 
  boundSurf_incr[11] = 2.8125*coeff[1]*fSkin[11]-5.083290641897235*coeff[1]*fSkin[20]; 
  boundSurf_incr[12] = 19.6875*coeff[1]*fSkin[12]-10.892765661208362*coeff[1]*fSkin[4]; 
  boundSurf_incr[14] = 19.6875*coeff[1]*fSkin[14]-10.892765661208362*coeff[1]*fSkin[6]; 
  boundSurf_incr[16] = 2.8125*coeff[1]*fSkin[16]-5.083290641897235*coeff[1]*fSkin[22]; 
  boundSurf_incr[17] = 2.8125*coeff[1]*fSkin[17]-5.083290641897234*coeff[1]*fSkin[23]; 
  boundSurf_incr[18] = 19.6875*coeff[1]*fSkin[18]-10.892765661208358*coeff[1]*fSkin[10]; 
  boundSurf_incr[19] = 2.8125*coeff[1]*fSkin[19]-5.083290641897234*coeff[1]*fSkin[25]; 
  boundSurf_incr[20] = 19.6875*coeff[1]*fSkin[20]-10.892765661208362*coeff[1]*fSkin[11]; 
  boundSurf_incr[22] = 19.6875*coeff[1]*fSkin[22]-10.892765661208362*coeff[1]*fSkin[16]; 
  boundSurf_incr[23] = 19.6875*coeff[1]*fSkin[23]-10.892765661208358*coeff[1]*fSkin[17]; 
  boundSurf_incr[24] = 2.8125*coeff[1]*fSkin[24]-5.083290641897234*coeff[1]*fSkin[26]; 
  boundSurf_incr[25] = 19.6875*coeff[1]*fSkin[25]-10.892765661208358*coeff[1]*fSkin[19]; 
  boundSurf_incr[26] = 19.6875*coeff[1]*fSkin[26]-10.892765661208358*coeff[1]*fSkin[24]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[1]*fSkin[8]-6.708203932499369*coeff[1]*fEdge[8]-8.11898816047911*coeff[1]*fSkin[2]-8.11898816047911*coeff[1]*fEdge[2]+4.6875*fSkin[0]*coeff[1]-4.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[1]*fSkin[12]-6.7082039324993685*coeff[1]*fEdge[12]-8.11898816047911*coeff[1]*fSkin[4]-8.11898816047911*coeff[1]*fEdge[4]+4.6875*coeff[1]*fSkin[1]-4.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(14.160595359570864*coeff[1]*fSkin[8])+9.077304717673634*coeff[1]*fEdge[8]+15.46875*coeff[1]*fSkin[2]+12.65625*coeff[1]*fEdge[2]-8.11898816047911*fSkin[0]*coeff[1]+8.11898816047911*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[1]*fSkin[14]-6.7082039324993685*coeff[1]*fEdge[14]-8.11898816047911*coeff[1]*fSkin[6]-8.11898816047911*coeff[1]*fEdge[6]+4.6875*coeff[1]*fSkin[3]-4.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(14.160595359570868*coeff[1]*fSkin[12])+9.077304717673634*coeff[1]*fEdge[12]+15.46875*coeff[1]*fSkin[4]+12.65625*coeff[1]*fEdge[4]-8.11898816047911*coeff[1]*fSkin[1]+8.11898816047911*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = 6.708203932499369*coeff[1]*fSkin[18]-6.708203932499369*coeff[1]*fEdge[18]-8.11898816047911*coeff[1]*fSkin[10]-8.11898816047911*coeff[1]*fEdge[10]+4.6875*coeff[1]*fSkin[5]-4.6875*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = -(14.160595359570868*coeff[1]*fSkin[14])+9.077304717673634*coeff[1]*fEdge[14]+15.46875*coeff[1]*fSkin[6]+12.65625*coeff[1]*fEdge[6]-8.11898816047911*coeff[1]*fSkin[3]+8.11898816047911*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[1]*fSkin[20]-6.708203932499369*coeff[1]*fEdge[20]-8.118988160479114*coeff[1]*fSkin[11]-8.118988160479114*coeff[1]*fEdge[11]+4.6875*coeff[1]*fSkin[7]-4.6875*coeff[1]*fEdge[7]; 
  edgeSurf_incr[8] = 20.34375*coeff[1]*fSkin[8]-0.65625*coeff[1]*fEdge[8]-15.61296411439865*coeff[1]*fSkin[2]-4.720198453190292*coeff[1]*fEdge[2]+4.192627457812107*fSkin[0]*coeff[1]-4.192627457812107*fEdge[0]*coeff[1]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[1]*fSkin[22]-6.708203932499369*coeff[1]*fEdge[22]-8.118988160479114*coeff[1]*fSkin[16]-8.118988160479114*coeff[1]*fEdge[16]+4.6875*coeff[1]*fSkin[9]-4.6875*coeff[1]*fEdge[9]; 
  edgeSurf_incr[10] = -(14.160595359570864*coeff[1]*fSkin[18])+9.077304717673634*coeff[1]*fEdge[18]+15.46875*coeff[1]*fSkin[10]+12.65625*coeff[1]*fEdge[10]-8.11898816047911*coeff[1]*fSkin[5]+8.11898816047911*coeff[1]*fEdge[5]; 
  edgeSurf_incr[11] = -(14.160595359570868*coeff[1]*fSkin[20])+9.077304717673634*coeff[1]*fEdge[20]+15.46875*coeff[1]*fSkin[11]+12.65625*coeff[1]*fEdge[11]-8.118988160479114*coeff[1]*fSkin[7]+8.118988160479114*coeff[1]*fEdge[7]; 
  edgeSurf_incr[12] = 20.34375*coeff[1]*fSkin[12]-0.65625*coeff[1]*fEdge[12]-15.612964114398654*coeff[1]*fSkin[4]-4.72019845319029*coeff[1]*fEdge[4]+4.192627457812107*coeff[1]*fSkin[1]-4.192627457812107*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = 6.7082039324993685*coeff[1]*fSkin[23]-6.7082039324993685*coeff[1]*fEdge[23]-8.118988160479114*coeff[1]*fSkin[17]-8.118988160479114*coeff[1]*fEdge[17]+4.6875*coeff[1]*fSkin[13]-4.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = 20.34375*coeff[1]*fSkin[14]-0.65625*coeff[1]*fEdge[14]-15.612964114398654*coeff[1]*fSkin[6]-4.72019845319029*coeff[1]*fEdge[6]+4.192627457812107*coeff[1]*fSkin[3]-4.192627457812107*coeff[1]*fEdge[3]; 
  edgeSurf_incr[15] = 6.7082039324993685*coeff[1]*fSkin[25]-6.7082039324993685*coeff[1]*fEdge[25]-8.118988160479114*coeff[1]*fSkin[19]-8.118988160479114*coeff[1]*fEdge[19]+4.6875*coeff[1]*fSkin[15]-4.6875*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = -(14.160595359570868*coeff[1]*fSkin[22])+9.077304717673634*coeff[1]*fEdge[22]+15.46875*coeff[1]*fSkin[16]+12.65625*coeff[1]*fEdge[16]-8.118988160479114*coeff[1]*fSkin[9]+8.118988160479114*coeff[1]*fEdge[9]; 
  edgeSurf_incr[17] = -(14.160595359570864*coeff[1]*fSkin[23])+9.077304717673634*coeff[1]*fEdge[23]+15.46875*coeff[1]*fSkin[17]+12.65625*coeff[1]*fEdge[17]-8.118988160479114*coeff[1]*fSkin[13]+8.118988160479114*coeff[1]*fEdge[13]; 
  edgeSurf_incr[18] = 20.34375*coeff[1]*fSkin[18]-0.65625*coeff[1]*fEdge[18]-15.61296411439865*coeff[1]*fSkin[10]-4.720198453190292*coeff[1]*fEdge[10]+4.192627457812107*coeff[1]*fSkin[5]-4.192627457812107*coeff[1]*fEdge[5]; 
  edgeSurf_incr[19] = -(14.160595359570864*coeff[1]*fSkin[25])+9.077304717673634*coeff[1]*fEdge[25]+15.46875*coeff[1]*fSkin[19]+12.65625*coeff[1]*fEdge[19]-8.118988160479114*coeff[1]*fSkin[15]+8.118988160479114*coeff[1]*fEdge[15]; 
  edgeSurf_incr[20] = 20.34375*coeff[1]*fSkin[20]-0.65625*coeff[1]*fEdge[20]-15.612964114398654*coeff[1]*fSkin[11]-4.72019845319029*coeff[1]*fEdge[11]+4.192627457812107*coeff[1]*fSkin[7]-4.192627457812107*coeff[1]*fEdge[7]; 
  edgeSurf_incr[21] = 6.708203932499369*coeff[1]*fSkin[26]-6.708203932499369*coeff[1]*fEdge[26]-8.11898816047911*coeff[1]*fSkin[24]-8.11898816047911*coeff[1]*fEdge[24]+4.6875*coeff[1]*fSkin[21]-4.6875*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = 20.34375*coeff[1]*fSkin[22]-0.65625*coeff[1]*fEdge[22]-15.612964114398654*coeff[1]*fSkin[16]-4.72019845319029*coeff[1]*fEdge[16]+4.192627457812107*coeff[1]*fSkin[9]-4.192627457812107*coeff[1]*fEdge[9]; 
  edgeSurf_incr[23] = 20.34375*coeff[1]*fSkin[23]-0.65625*coeff[1]*fEdge[23]-15.61296411439865*coeff[1]*fSkin[17]-4.720198453190292*coeff[1]*fEdge[17]+4.192627457812107*coeff[1]*fSkin[13]-4.192627457812107*coeff[1]*fEdge[13]; 
  edgeSurf_incr[24] = -(14.160595359570864*coeff[1]*fSkin[26])+9.077304717673634*coeff[1]*fEdge[26]+15.46875*coeff[1]*fSkin[24]+12.65625*coeff[1]*fEdge[24]-8.11898816047911*coeff[1]*fSkin[21]+8.11898816047911*coeff[1]*fEdge[21]; 
  edgeSurf_incr[25] = 20.34375*coeff[1]*fSkin[25]-0.65625*coeff[1]*fEdge[25]-15.61296411439865*coeff[1]*fSkin[19]-4.720198453190292*coeff[1]*fEdge[19]+4.192627457812107*coeff[1]*fSkin[15]-4.192627457812107*coeff[1]*fEdge[15]; 
  edgeSurf_incr[26] = 20.34375*coeff[1]*fSkin[26]-0.65625*coeff[1]*fEdge[26]-15.61296411439865*coeff[1]*fSkin[24]-4.720198453190292*coeff[1]*fEdge[24]+4.192627457812107*coeff[1]*fSkin[21]-4.192627457812107*coeff[1]*fEdge[21]; 

  boundSurf_incr[2] = 5.083290641897234*coeff[1]*fSkin[8]+2.8125*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = 5.083290641897235*coeff[1]*fSkin[12]+2.8125*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = 5.083290641897235*coeff[1]*fSkin[14]+2.8125*coeff[1]*fSkin[6]; 
  boundSurf_incr[8] = 19.6875*coeff[1]*fSkin[8]+10.892765661208358*coeff[1]*fSkin[2]; 
  boundSurf_incr[10] = 5.083290641897234*coeff[1]*fSkin[18]+2.8125*coeff[1]*fSkin[10]; 
  boundSurf_incr[11] = 5.083290641897235*coeff[1]*fSkin[20]+2.8125*coeff[1]*fSkin[11]; 
  boundSurf_incr[12] = 19.6875*coeff[1]*fSkin[12]+10.892765661208362*coeff[1]*fSkin[4]; 
  boundSurf_incr[14] = 19.6875*coeff[1]*fSkin[14]+10.892765661208362*coeff[1]*fSkin[6]; 
  boundSurf_incr[16] = 5.083290641897235*coeff[1]*fSkin[22]+2.8125*coeff[1]*fSkin[16]; 
  boundSurf_incr[17] = 5.083290641897234*coeff[1]*fSkin[23]+2.8125*coeff[1]*fSkin[17]; 
  boundSurf_incr[18] = 19.6875*coeff[1]*fSkin[18]+10.892765661208358*coeff[1]*fSkin[10]; 
  boundSurf_incr[19] = 5.083290641897234*coeff[1]*fSkin[25]+2.8125*coeff[1]*fSkin[19]; 
  boundSurf_incr[20] = 19.6875*coeff[1]*fSkin[20]+10.892765661208362*coeff[1]*fSkin[11]; 
  boundSurf_incr[22] = 19.6875*coeff[1]*fSkin[22]+10.892765661208362*coeff[1]*fSkin[16]; 
  boundSurf_incr[23] = 19.6875*coeff[1]*fSkin[23]+10.892765661208358*coeff[1]*fSkin[17]; 
  boundSurf_incr[24] = 5.083290641897234*coeff[1]*fSkin[26]+2.8125*coeff[1]*fSkin[24]; 
  boundSurf_incr[25] = 19.6875*coeff[1]*fSkin[25]+10.892765661208358*coeff[1]*fSkin[19]; 
  boundSurf_incr[26] = 19.6875*coeff[1]*fSkin[26]+10.892765661208358*coeff[1]*fSkin[24]; 

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

  return 0.;
}

