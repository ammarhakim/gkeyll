#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfz_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],4.);

  double vol_incr[27] = {0.0}; 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[2]*fSkin[9]-6.708203932499369*coeff[2]*fEdge[9]+8.11898816047911*coeff[2]*fSkin[3]+8.11898816047911*coeff[2]*fEdge[3]+4.6875*fSkin[0]*coeff[2]-4.6875*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[2]*fSkin[15]-6.7082039324993685*coeff[2]*fEdge[15]+8.11898816047911*coeff[2]*fSkin[5]+8.11898816047911*coeff[2]*fEdge[5]+4.6875*fSkin[1]*coeff[2]-4.6875*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[2]*fSkin[16]-6.7082039324993685*coeff[2]*fEdge[16]+8.11898816047911*coeff[2]*fSkin[6]+8.11898816047911*coeff[2]*fEdge[6]+4.6875*coeff[2]*fSkin[2]-4.6875*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = 14.160595359570864*coeff[2]*fSkin[9]-9.077304717673634*coeff[2]*fEdge[9]+15.46875*coeff[2]*fSkin[3]+12.65625*coeff[2]*fEdge[3]+8.11898816047911*fSkin[0]*coeff[2]-8.11898816047911*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = 6.708203932499369*coeff[2]*fSkin[19]-6.708203932499369*coeff[2]*fEdge[19]+8.11898816047911*coeff[2]*fSkin[10]+8.11898816047911*coeff[2]*fEdge[10]+4.6875*coeff[2]*fSkin[4]-4.6875*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = 14.160595359570868*coeff[2]*fSkin[15]-9.077304717673634*coeff[2]*fEdge[15]+15.46875*coeff[2]*fSkin[5]+12.65625*coeff[2]*fEdge[5]+8.11898816047911*fSkin[1]*coeff[2]-8.11898816047911*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = 14.160595359570868*coeff[2]*fSkin[16]-9.077304717673634*coeff[2]*fEdge[16]+15.46875*coeff[2]*fSkin[6]+12.65625*coeff[2]*fEdge[6]+8.11898816047911*coeff[2]*fSkin[2]-8.11898816047911*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[2]*fSkin[21]-6.708203932499369*coeff[2]*fEdge[21]+8.118988160479114*coeff[2]*fSkin[13]+8.118988160479114*coeff[2]*fEdge[13]+4.6875*coeff[2]*fSkin[7]-4.6875*coeff[2]*fEdge[7]; 
  edgeSurf_incr[8] = 6.708203932499369*coeff[2]*fSkin[22]-6.708203932499369*coeff[2]*fEdge[22]+8.118988160479114*coeff[2]*fSkin[14]+8.118988160479114*coeff[2]*fEdge[14]+4.6875*coeff[2]*fSkin[8]-4.6875*coeff[2]*fEdge[8]; 
  edgeSurf_incr[9] = 20.34375*coeff[2]*fSkin[9]-0.65625*coeff[2]*fEdge[9]+15.61296411439865*coeff[2]*fSkin[3]+4.720198453190292*coeff[2]*fEdge[3]+4.192627457812107*fSkin[0]*coeff[2]-4.192627457812107*fEdge[0]*coeff[2]; 
  edgeSurf_incr[10] = 14.160595359570864*coeff[2]*fSkin[19]-9.077304717673634*coeff[2]*fEdge[19]+15.46875*coeff[2]*fSkin[10]+12.65625*coeff[2]*fEdge[10]+8.11898816047911*coeff[2]*fSkin[4]-8.11898816047911*coeff[2]*fEdge[4]; 
  edgeSurf_incr[11] = 6.7082039324993685*coeff[2]*fSkin[24]-6.7082039324993685*coeff[2]*fEdge[24]+8.118988160479114*coeff[2]*fSkin[17]+8.118988160479114*coeff[2]*fEdge[17]+4.6875*coeff[2]*fSkin[11]-4.6875*coeff[2]*fEdge[11]; 
  edgeSurf_incr[12] = 6.7082039324993685*coeff[2]*fSkin[25]-6.7082039324993685*coeff[2]*fEdge[25]+8.118988160479114*coeff[2]*fSkin[18]+8.118988160479114*coeff[2]*fEdge[18]+4.6875*coeff[2]*fSkin[12]-4.6875*coeff[2]*fEdge[12]; 
  edgeSurf_incr[13] = 14.160595359570868*coeff[2]*fSkin[21]-9.077304717673634*coeff[2]*fEdge[21]+15.46875*coeff[2]*fSkin[13]+12.65625*coeff[2]*fEdge[13]+8.118988160479114*coeff[2]*fSkin[7]-8.118988160479114*coeff[2]*fEdge[7]; 
  edgeSurf_incr[14] = 14.160595359570868*coeff[2]*fSkin[22]-9.077304717673634*coeff[2]*fEdge[22]+15.46875*coeff[2]*fSkin[14]+12.65625*coeff[2]*fEdge[14]+8.118988160479114*coeff[2]*fSkin[8]-8.118988160479114*coeff[2]*fEdge[8]; 
  edgeSurf_incr[15] = 20.34375*coeff[2]*fSkin[15]-0.65625*coeff[2]*fEdge[15]+15.612964114398654*coeff[2]*fSkin[5]+4.72019845319029*coeff[2]*fEdge[5]+4.192627457812107*fSkin[1]*coeff[2]-4.192627457812107*fEdge[1]*coeff[2]; 
  edgeSurf_incr[16] = 20.34375*coeff[2]*fSkin[16]-0.65625*coeff[2]*fEdge[16]+15.612964114398654*coeff[2]*fSkin[6]+4.72019845319029*coeff[2]*fEdge[6]+4.192627457812107*coeff[2]*fSkin[2]-4.192627457812107*coeff[2]*fEdge[2]; 
  edgeSurf_incr[17] = 14.160595359570864*coeff[2]*fSkin[24]-9.077304717673634*coeff[2]*fEdge[24]+15.46875*coeff[2]*fSkin[17]+12.65625*coeff[2]*fEdge[17]+8.118988160479114*coeff[2]*fSkin[11]-8.118988160479114*coeff[2]*fEdge[11]; 
  edgeSurf_incr[18] = 14.160595359570864*coeff[2]*fSkin[25]-9.077304717673634*coeff[2]*fEdge[25]+15.46875*coeff[2]*fSkin[18]+12.65625*coeff[2]*fEdge[18]+8.118988160479114*coeff[2]*fSkin[12]-8.118988160479114*coeff[2]*fEdge[12]; 
  edgeSurf_incr[19] = 20.34375*coeff[2]*fSkin[19]-0.65625*coeff[2]*fEdge[19]+15.61296411439865*coeff[2]*fSkin[10]+4.720198453190292*coeff[2]*fEdge[10]+4.192627457812107*coeff[2]*fSkin[4]-4.192627457812107*coeff[2]*fEdge[4]; 
  edgeSurf_incr[20] = 6.708203932499369*coeff[2]*fSkin[26]-6.708203932499369*coeff[2]*fEdge[26]+8.11898816047911*coeff[2]*fSkin[23]+8.11898816047911*coeff[2]*fEdge[23]+4.6875*coeff[2]*fSkin[20]-4.6875*coeff[2]*fEdge[20]; 
  edgeSurf_incr[21] = 20.34375*coeff[2]*fSkin[21]-0.65625*coeff[2]*fEdge[21]+15.612964114398654*coeff[2]*fSkin[13]+4.72019845319029*coeff[2]*fEdge[13]+4.192627457812107*coeff[2]*fSkin[7]-4.192627457812107*coeff[2]*fEdge[7]; 
  edgeSurf_incr[22] = 20.34375*coeff[2]*fSkin[22]-0.65625*coeff[2]*fEdge[22]+15.612964114398654*coeff[2]*fSkin[14]+4.72019845319029*coeff[2]*fEdge[14]+4.192627457812107*coeff[2]*fSkin[8]-4.192627457812107*coeff[2]*fEdge[8]; 
  edgeSurf_incr[23] = 14.160595359570864*coeff[2]*fSkin[26]-9.077304717673634*coeff[2]*fEdge[26]+15.46875*coeff[2]*fSkin[23]+12.65625*coeff[2]*fEdge[23]+8.11898816047911*coeff[2]*fSkin[20]-8.11898816047911*coeff[2]*fEdge[20]; 
  edgeSurf_incr[24] = 20.34375*coeff[2]*fSkin[24]-0.65625*coeff[2]*fEdge[24]+15.61296411439865*coeff[2]*fSkin[17]+4.720198453190292*coeff[2]*fEdge[17]+4.192627457812107*coeff[2]*fSkin[11]-4.192627457812107*coeff[2]*fEdge[11]; 
  edgeSurf_incr[25] = 20.34375*coeff[2]*fSkin[25]-0.65625*coeff[2]*fEdge[25]+15.61296411439865*coeff[2]*fSkin[18]+4.720198453190292*coeff[2]*fEdge[18]+4.192627457812107*coeff[2]*fSkin[12]-4.192627457812107*coeff[2]*fEdge[12]; 
  edgeSurf_incr[26] = 20.34375*coeff[2]*fSkin[26]-0.65625*coeff[2]*fEdge[26]+15.61296411439865*coeff[2]*fSkin[23]+4.720198453190292*coeff[2]*fEdge[23]+4.192627457812107*coeff[2]*fSkin[20]-4.192627457812107*coeff[2]*fEdge[20]; 

  boundSurf_incr[3] = 2.8125*coeff[2]*fSkin[3]-5.083290641897234*coeff[2]*fSkin[9]; 
  boundSurf_incr[5] = 2.8125*coeff[2]*fSkin[5]-5.083290641897235*coeff[2]*fSkin[15]; 
  boundSurf_incr[6] = 2.8125*coeff[2]*fSkin[6]-5.083290641897235*coeff[2]*fSkin[16]; 
  boundSurf_incr[9] = 19.6875*coeff[2]*fSkin[9]-10.892765661208358*coeff[2]*fSkin[3]; 
  boundSurf_incr[10] = 2.8125*coeff[2]*fSkin[10]-5.083290641897234*coeff[2]*fSkin[19]; 
  boundSurf_incr[13] = 2.8125*coeff[2]*fSkin[13]-5.083290641897235*coeff[2]*fSkin[21]; 
  boundSurf_incr[14] = 2.8125*coeff[2]*fSkin[14]-5.083290641897235*coeff[2]*fSkin[22]; 
  boundSurf_incr[15] = 19.6875*coeff[2]*fSkin[15]-10.892765661208362*coeff[2]*fSkin[5]; 
  boundSurf_incr[16] = 19.6875*coeff[2]*fSkin[16]-10.892765661208362*coeff[2]*fSkin[6]; 
  boundSurf_incr[17] = 2.8125*coeff[2]*fSkin[17]-5.083290641897234*coeff[2]*fSkin[24]; 
  boundSurf_incr[18] = 2.8125*coeff[2]*fSkin[18]-5.083290641897234*coeff[2]*fSkin[25]; 
  boundSurf_incr[19] = 19.6875*coeff[2]*fSkin[19]-10.892765661208358*coeff[2]*fSkin[10]; 
  boundSurf_incr[21] = 19.6875*coeff[2]*fSkin[21]-10.892765661208362*coeff[2]*fSkin[13]; 
  boundSurf_incr[22] = 19.6875*coeff[2]*fSkin[22]-10.892765661208362*coeff[2]*fSkin[14]; 
  boundSurf_incr[23] = 2.8125*coeff[2]*fSkin[23]-5.083290641897234*coeff[2]*fSkin[26]; 
  boundSurf_incr[24] = 19.6875*coeff[2]*fSkin[24]-10.892765661208358*coeff[2]*fSkin[17]; 
  boundSurf_incr[25] = 19.6875*coeff[2]*fSkin[25]-10.892765661208358*coeff[2]*fSkin[18]; 
  boundSurf_incr[26] = 19.6875*coeff[2]*fSkin[26]-10.892765661208358*coeff[2]*fSkin[23]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[2]*fSkin[9]-6.708203932499369*coeff[2]*fEdge[9]-8.11898816047911*coeff[2]*fSkin[3]-8.11898816047911*coeff[2]*fEdge[3]+4.6875*fSkin[0]*coeff[2]-4.6875*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[2]*fSkin[15]-6.7082039324993685*coeff[2]*fEdge[15]-8.11898816047911*coeff[2]*fSkin[5]-8.11898816047911*coeff[2]*fEdge[5]+4.6875*fSkin[1]*coeff[2]-4.6875*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[2]*fSkin[16]-6.7082039324993685*coeff[2]*fEdge[16]-8.11898816047911*coeff[2]*fSkin[6]-8.11898816047911*coeff[2]*fEdge[6]+4.6875*coeff[2]*fSkin[2]-4.6875*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = -(14.160595359570864*coeff[2]*fSkin[9])+9.077304717673634*coeff[2]*fEdge[9]+15.46875*coeff[2]*fSkin[3]+12.65625*coeff[2]*fEdge[3]-8.11898816047911*fSkin[0]*coeff[2]+8.11898816047911*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = 6.708203932499369*coeff[2]*fSkin[19]-6.708203932499369*coeff[2]*fEdge[19]-8.11898816047911*coeff[2]*fSkin[10]-8.11898816047911*coeff[2]*fEdge[10]+4.6875*coeff[2]*fSkin[4]-4.6875*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = -(14.160595359570868*coeff[2]*fSkin[15])+9.077304717673634*coeff[2]*fEdge[15]+15.46875*coeff[2]*fSkin[5]+12.65625*coeff[2]*fEdge[5]-8.11898816047911*fSkin[1]*coeff[2]+8.11898816047911*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = -(14.160595359570868*coeff[2]*fSkin[16])+9.077304717673634*coeff[2]*fEdge[16]+15.46875*coeff[2]*fSkin[6]+12.65625*coeff[2]*fEdge[6]-8.11898816047911*coeff[2]*fSkin[2]+8.11898816047911*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[2]*fSkin[21]-6.708203932499369*coeff[2]*fEdge[21]-8.118988160479114*coeff[2]*fSkin[13]-8.118988160479114*coeff[2]*fEdge[13]+4.6875*coeff[2]*fSkin[7]-4.6875*coeff[2]*fEdge[7]; 
  edgeSurf_incr[8] = 6.708203932499369*coeff[2]*fSkin[22]-6.708203932499369*coeff[2]*fEdge[22]-8.118988160479114*coeff[2]*fSkin[14]-8.118988160479114*coeff[2]*fEdge[14]+4.6875*coeff[2]*fSkin[8]-4.6875*coeff[2]*fEdge[8]; 
  edgeSurf_incr[9] = 20.34375*coeff[2]*fSkin[9]-0.65625*coeff[2]*fEdge[9]-15.61296411439865*coeff[2]*fSkin[3]-4.720198453190292*coeff[2]*fEdge[3]+4.192627457812107*fSkin[0]*coeff[2]-4.192627457812107*fEdge[0]*coeff[2]; 
  edgeSurf_incr[10] = -(14.160595359570864*coeff[2]*fSkin[19])+9.077304717673634*coeff[2]*fEdge[19]+15.46875*coeff[2]*fSkin[10]+12.65625*coeff[2]*fEdge[10]-8.11898816047911*coeff[2]*fSkin[4]+8.11898816047911*coeff[2]*fEdge[4]; 
  edgeSurf_incr[11] = 6.7082039324993685*coeff[2]*fSkin[24]-6.7082039324993685*coeff[2]*fEdge[24]-8.118988160479114*coeff[2]*fSkin[17]-8.118988160479114*coeff[2]*fEdge[17]+4.6875*coeff[2]*fSkin[11]-4.6875*coeff[2]*fEdge[11]; 
  edgeSurf_incr[12] = 6.7082039324993685*coeff[2]*fSkin[25]-6.7082039324993685*coeff[2]*fEdge[25]-8.118988160479114*coeff[2]*fSkin[18]-8.118988160479114*coeff[2]*fEdge[18]+4.6875*coeff[2]*fSkin[12]-4.6875*coeff[2]*fEdge[12]; 
  edgeSurf_incr[13] = -(14.160595359570868*coeff[2]*fSkin[21])+9.077304717673634*coeff[2]*fEdge[21]+15.46875*coeff[2]*fSkin[13]+12.65625*coeff[2]*fEdge[13]-8.118988160479114*coeff[2]*fSkin[7]+8.118988160479114*coeff[2]*fEdge[7]; 
  edgeSurf_incr[14] = -(14.160595359570868*coeff[2]*fSkin[22])+9.077304717673634*coeff[2]*fEdge[22]+15.46875*coeff[2]*fSkin[14]+12.65625*coeff[2]*fEdge[14]-8.118988160479114*coeff[2]*fSkin[8]+8.118988160479114*coeff[2]*fEdge[8]; 
  edgeSurf_incr[15] = 20.34375*coeff[2]*fSkin[15]-0.65625*coeff[2]*fEdge[15]-15.612964114398654*coeff[2]*fSkin[5]-4.72019845319029*coeff[2]*fEdge[5]+4.192627457812107*fSkin[1]*coeff[2]-4.192627457812107*fEdge[1]*coeff[2]; 
  edgeSurf_incr[16] = 20.34375*coeff[2]*fSkin[16]-0.65625*coeff[2]*fEdge[16]-15.612964114398654*coeff[2]*fSkin[6]-4.72019845319029*coeff[2]*fEdge[6]+4.192627457812107*coeff[2]*fSkin[2]-4.192627457812107*coeff[2]*fEdge[2]; 
  edgeSurf_incr[17] = -(14.160595359570864*coeff[2]*fSkin[24])+9.077304717673634*coeff[2]*fEdge[24]+15.46875*coeff[2]*fSkin[17]+12.65625*coeff[2]*fEdge[17]-8.118988160479114*coeff[2]*fSkin[11]+8.118988160479114*coeff[2]*fEdge[11]; 
  edgeSurf_incr[18] = -(14.160595359570864*coeff[2]*fSkin[25])+9.077304717673634*coeff[2]*fEdge[25]+15.46875*coeff[2]*fSkin[18]+12.65625*coeff[2]*fEdge[18]-8.118988160479114*coeff[2]*fSkin[12]+8.118988160479114*coeff[2]*fEdge[12]; 
  edgeSurf_incr[19] = 20.34375*coeff[2]*fSkin[19]-0.65625*coeff[2]*fEdge[19]-15.61296411439865*coeff[2]*fSkin[10]-4.720198453190292*coeff[2]*fEdge[10]+4.192627457812107*coeff[2]*fSkin[4]-4.192627457812107*coeff[2]*fEdge[4]; 
  edgeSurf_incr[20] = 6.708203932499369*coeff[2]*fSkin[26]-6.708203932499369*coeff[2]*fEdge[26]-8.11898816047911*coeff[2]*fSkin[23]-8.11898816047911*coeff[2]*fEdge[23]+4.6875*coeff[2]*fSkin[20]-4.6875*coeff[2]*fEdge[20]; 
  edgeSurf_incr[21] = 20.34375*coeff[2]*fSkin[21]-0.65625*coeff[2]*fEdge[21]-15.612964114398654*coeff[2]*fSkin[13]-4.72019845319029*coeff[2]*fEdge[13]+4.192627457812107*coeff[2]*fSkin[7]-4.192627457812107*coeff[2]*fEdge[7]; 
  edgeSurf_incr[22] = 20.34375*coeff[2]*fSkin[22]-0.65625*coeff[2]*fEdge[22]-15.612964114398654*coeff[2]*fSkin[14]-4.72019845319029*coeff[2]*fEdge[14]+4.192627457812107*coeff[2]*fSkin[8]-4.192627457812107*coeff[2]*fEdge[8]; 
  edgeSurf_incr[23] = -(14.160595359570864*coeff[2]*fSkin[26])+9.077304717673634*coeff[2]*fEdge[26]+15.46875*coeff[2]*fSkin[23]+12.65625*coeff[2]*fEdge[23]-8.11898816047911*coeff[2]*fSkin[20]+8.11898816047911*coeff[2]*fEdge[20]; 
  edgeSurf_incr[24] = 20.34375*coeff[2]*fSkin[24]-0.65625*coeff[2]*fEdge[24]-15.61296411439865*coeff[2]*fSkin[17]-4.720198453190292*coeff[2]*fEdge[17]+4.192627457812107*coeff[2]*fSkin[11]-4.192627457812107*coeff[2]*fEdge[11]; 
  edgeSurf_incr[25] = 20.34375*coeff[2]*fSkin[25]-0.65625*coeff[2]*fEdge[25]-15.61296411439865*coeff[2]*fSkin[18]-4.720198453190292*coeff[2]*fEdge[18]+4.192627457812107*coeff[2]*fSkin[12]-4.192627457812107*coeff[2]*fEdge[12]; 
  edgeSurf_incr[26] = 20.34375*coeff[2]*fSkin[26]-0.65625*coeff[2]*fEdge[26]-15.61296411439865*coeff[2]*fSkin[23]-4.720198453190292*coeff[2]*fEdge[23]+4.192627457812107*coeff[2]*fSkin[20]-4.192627457812107*coeff[2]*fEdge[20]; 

  boundSurf_incr[3] = 5.083290641897234*coeff[2]*fSkin[9]+2.8125*coeff[2]*fSkin[3]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[2]*fSkin[15]+2.8125*coeff[2]*fSkin[5]; 
  boundSurf_incr[6] = 5.083290641897235*coeff[2]*fSkin[16]+2.8125*coeff[2]*fSkin[6]; 
  boundSurf_incr[9] = 19.6875*coeff[2]*fSkin[9]+10.892765661208358*coeff[2]*fSkin[3]; 
  boundSurf_incr[10] = 5.083290641897234*coeff[2]*fSkin[19]+2.8125*coeff[2]*fSkin[10]; 
  boundSurf_incr[13] = 5.083290641897235*coeff[2]*fSkin[21]+2.8125*coeff[2]*fSkin[13]; 
  boundSurf_incr[14] = 5.083290641897235*coeff[2]*fSkin[22]+2.8125*coeff[2]*fSkin[14]; 
  boundSurf_incr[15] = 19.6875*coeff[2]*fSkin[15]+10.892765661208362*coeff[2]*fSkin[5]; 
  boundSurf_incr[16] = 19.6875*coeff[2]*fSkin[16]+10.892765661208362*coeff[2]*fSkin[6]; 
  boundSurf_incr[17] = 5.083290641897234*coeff[2]*fSkin[24]+2.8125*coeff[2]*fSkin[17]; 
  boundSurf_incr[18] = 5.083290641897234*coeff[2]*fSkin[25]+2.8125*coeff[2]*fSkin[18]; 
  boundSurf_incr[19] = 19.6875*coeff[2]*fSkin[19]+10.892765661208358*coeff[2]*fSkin[10]; 
  boundSurf_incr[21] = 19.6875*coeff[2]*fSkin[21]+10.892765661208362*coeff[2]*fSkin[13]; 
  boundSurf_incr[22] = 19.6875*coeff[2]*fSkin[22]+10.892765661208362*coeff[2]*fSkin[14]; 
  boundSurf_incr[23] = 5.083290641897234*coeff[2]*fSkin[26]+2.8125*coeff[2]*fSkin[23]; 
  boundSurf_incr[24] = 19.6875*coeff[2]*fSkin[24]+10.892765661208358*coeff[2]*fSkin[17]; 
  boundSurf_incr[25] = 19.6875*coeff[2]*fSkin[25]+10.892765661208358*coeff[2]*fSkin[18]; 
  boundSurf_incr[26] = 19.6875*coeff[2]*fSkin[26]+10.892765661208358*coeff[2]*fSkin[23]; 

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

