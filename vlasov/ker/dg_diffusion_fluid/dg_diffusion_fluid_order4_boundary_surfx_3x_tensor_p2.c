#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[27] = {0.0}; 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[7]-6.708203932499369*coeff[0]*fEdge[7]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.160595359570864*coeff[0]*fSkin[7]-9.077304717673634*coeff[0]*fEdge[7]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*fSkin[11]-6.7082039324993685*coeff[0]*fEdge[11]+8.11898816047911*coeff[0]*fSkin[4]+8.11898816047911*coeff[0]*fEdge[4]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*fSkin[13]-6.7082039324993685*coeff[0]*fEdge[13]+8.11898816047911*coeff[0]*fSkin[5]+8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 14.160595359570868*coeff[0]*fSkin[11]-9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[4]+12.65625*coeff[0]*fEdge[4]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 14.160595359570868*coeff[0]*fSkin[13]-9.077304717673634*coeff[0]*fEdge[13]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]+8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*fSkin[17]-6.708203932499369*coeff[0]*fEdge[17]+8.11898816047911*coeff[0]*fSkin[10]+8.11898816047911*coeff[0]*fEdge[10]+4.6875*coeff[0]*fSkin[6]-4.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*fSkin[7]-0.65625*coeff[0]*fEdge[7]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 6.708203932499369*coeff[0]*fSkin[20]-6.708203932499369*coeff[0]*fEdge[20]+8.118988160479114*coeff[0]*fSkin[12]+8.118988160479114*coeff[0]*fEdge[12]+4.6875*coeff[0]*fSkin[8]-4.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*fSkin[21]-6.708203932499369*coeff[0]*fEdge[21]+8.118988160479114*coeff[0]*fSkin[15]+8.118988160479114*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 14.160595359570864*coeff[0]*fSkin[17]-9.077304717673634*coeff[0]*fEdge[17]+15.46875*coeff[0]*fSkin[10]+12.65625*coeff[0]*fEdge[10]+8.11898816047911*coeff[0]*fSkin[6]-8.11898816047911*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]+15.612964114398654*coeff[0]*fSkin[4]+4.72019845319029*coeff[0]*fEdge[4]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 14.160595359570868*coeff[0]*fSkin[20]-9.077304717673634*coeff[0]*fEdge[20]+15.46875*coeff[0]*fSkin[12]+12.65625*coeff[0]*fEdge[12]+8.118988160479114*coeff[0]*fSkin[8]-8.118988160479114*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*fSkin[13]-0.65625*coeff[0]*fEdge[13]+15.612964114398654*coeff[0]*fSkin[5]+4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 6.7082039324993685*coeff[0]*fSkin[23]-6.7082039324993685*coeff[0]*fEdge[23]+8.118988160479114*coeff[0]*fSkin[18]+8.118988160479114*coeff[0]*fEdge[18]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 14.160595359570868*coeff[0]*fSkin[21]-9.077304717673634*coeff[0]*fEdge[21]+15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]+8.118988160479114*coeff[0]*fSkin[9]-8.118988160479114*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 6.7082039324993685*coeff[0]*fSkin[24]-6.7082039324993685*coeff[0]*fEdge[24]+8.118988160479114*coeff[0]*fSkin[19]+8.118988160479114*coeff[0]*fEdge[19]+4.6875*coeff[0]*fSkin[16]-4.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*fSkin[17]-0.65625*coeff[0]*fEdge[17]+15.61296411439865*coeff[0]*fSkin[10]+4.720198453190292*coeff[0]*fEdge[10]+4.192627457812107*coeff[0]*fSkin[6]-4.192627457812107*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 14.160595359570864*coeff[0]*fSkin[23]-9.077304717673634*coeff[0]*fEdge[23]+15.46875*coeff[0]*fSkin[18]+12.65625*coeff[0]*fEdge[18]+8.118988160479114*coeff[0]*fSkin[14]-8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 14.160595359570864*coeff[0]*fSkin[24]-9.077304717673634*coeff[0]*fEdge[24]+15.46875*coeff[0]*fSkin[19]+12.65625*coeff[0]*fEdge[19]+8.118988160479114*coeff[0]*fSkin[16]-8.118988160479114*coeff[0]*fEdge[16]; 
  edgeSurf_incr[20] = 20.34375*coeff[0]*fSkin[20]-0.65625*coeff[0]*fEdge[20]+15.612964114398654*coeff[0]*fSkin[12]+4.72019845319029*coeff[0]*fEdge[12]+4.192627457812107*coeff[0]*fSkin[8]-4.192627457812107*coeff[0]*fEdge[8]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*fSkin[21]-0.65625*coeff[0]*fEdge[21]+15.612964114398654*coeff[0]*fSkin[15]+4.72019845319029*coeff[0]*fEdge[15]+4.192627457812107*coeff[0]*fSkin[9]-4.192627457812107*coeff[0]*fEdge[9]; 
  edgeSurf_incr[22] = 6.708203932499369*coeff[0]*fSkin[26]-6.708203932499369*coeff[0]*fEdge[26]+8.11898816047911*coeff[0]*fSkin[25]+8.11898816047911*coeff[0]*fEdge[25]+4.6875*coeff[0]*fSkin[22]-4.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 20.34375*coeff[0]*fSkin[23]-0.65625*coeff[0]*fEdge[23]+15.61296411439865*coeff[0]*fSkin[18]+4.720198453190292*coeff[0]*fEdge[18]+4.192627457812107*coeff[0]*fSkin[14]-4.192627457812107*coeff[0]*fEdge[14]; 
  edgeSurf_incr[24] = 20.34375*coeff[0]*fSkin[24]-0.65625*coeff[0]*fEdge[24]+15.61296411439865*coeff[0]*fSkin[19]+4.720198453190292*coeff[0]*fEdge[19]+4.192627457812107*coeff[0]*fSkin[16]-4.192627457812107*coeff[0]*fEdge[16]; 
  edgeSurf_incr[25] = 14.160595359570864*coeff[0]*fSkin[26]-9.077304717673634*coeff[0]*fEdge[26]+15.46875*coeff[0]*fSkin[25]+12.65625*coeff[0]*fEdge[25]+8.11898816047911*coeff[0]*fSkin[22]-8.11898816047911*coeff[0]*fEdge[22]; 
  edgeSurf_incr[26] = 20.34375*coeff[0]*fSkin[26]-0.65625*coeff[0]*fEdge[26]+15.61296411439865*coeff[0]*fSkin[25]+4.720198453190292*coeff[0]*fEdge[25]+4.192627457812107*coeff[0]*fSkin[22]-4.192627457812107*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[7]; 
  boundSurf_incr[4] = 2.8125*coeff[0]*fSkin[4]-5.083290641897235*coeff[0]*fSkin[11]; 
  boundSurf_incr[5] = 2.8125*coeff[0]*fSkin[5]-5.083290641897235*coeff[0]*fSkin[13]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*fSkin[7]-10.892765661208358*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = 2.8125*coeff[0]*fSkin[10]-5.083290641897234*coeff[0]*fSkin[17]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]-10.892765661208362*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 2.8125*coeff[0]*fSkin[12]-5.083290641897235*coeff[0]*fSkin[20]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*fSkin[13]-10.892765661208362*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*fSkin[15]-5.083290641897235*coeff[0]*fSkin[21]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*fSkin[17]-10.892765661208358*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 2.8125*coeff[0]*fSkin[18]-5.083290641897234*coeff[0]*fSkin[23]; 
  boundSurf_incr[19] = 2.8125*coeff[0]*fSkin[19]-5.083290641897234*coeff[0]*fSkin[24]; 
  boundSurf_incr[20] = 19.6875*coeff[0]*fSkin[20]-10.892765661208362*coeff[0]*fSkin[12]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*fSkin[21]-10.892765661208362*coeff[0]*fSkin[15]; 
  boundSurf_incr[23] = 19.6875*coeff[0]*fSkin[23]-10.892765661208358*coeff[0]*fSkin[18]; 
  boundSurf_incr[24] = 19.6875*coeff[0]*fSkin[24]-10.892765661208358*coeff[0]*fSkin[19]; 
  boundSurf_incr[25] = 2.8125*coeff[0]*fSkin[25]-5.083290641897234*coeff[0]*fSkin[26]; 
  boundSurf_incr[26] = 19.6875*coeff[0]*fSkin[26]-10.892765661208358*coeff[0]*fSkin[25]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[7]-6.708203932499369*coeff[0]*fEdge[7]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(14.160595359570864*coeff[0]*fSkin[7])+9.077304717673634*coeff[0]*fEdge[7]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*fSkin[11]-6.7082039324993685*coeff[0]*fEdge[11]-8.11898816047911*coeff[0]*fSkin[4]-8.11898816047911*coeff[0]*fEdge[4]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*fSkin[13]-6.7082039324993685*coeff[0]*fEdge[13]-8.11898816047911*coeff[0]*fSkin[5]-8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(14.160595359570868*coeff[0]*fSkin[11])+9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[4]+12.65625*coeff[0]*fEdge[4]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = -(14.160595359570868*coeff[0]*fSkin[13])+9.077304717673634*coeff[0]*fEdge[13]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]-8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*fSkin[17]-6.708203932499369*coeff[0]*fEdge[17]-8.11898816047911*coeff[0]*fSkin[10]-8.11898816047911*coeff[0]*fEdge[10]+4.6875*coeff[0]*fSkin[6]-4.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*fSkin[7]-0.65625*coeff[0]*fEdge[7]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 6.708203932499369*coeff[0]*fSkin[20]-6.708203932499369*coeff[0]*fEdge[20]-8.118988160479114*coeff[0]*fSkin[12]-8.118988160479114*coeff[0]*fEdge[12]+4.6875*coeff[0]*fSkin[8]-4.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*fSkin[21]-6.708203932499369*coeff[0]*fEdge[21]-8.118988160479114*coeff[0]*fSkin[15]-8.118988160479114*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(14.160595359570864*coeff[0]*fSkin[17])+9.077304717673634*coeff[0]*fEdge[17]+15.46875*coeff[0]*fSkin[10]+12.65625*coeff[0]*fEdge[10]-8.11898816047911*coeff[0]*fSkin[6]+8.11898816047911*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]-15.612964114398654*coeff[0]*fSkin[4]-4.72019845319029*coeff[0]*fEdge[4]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = -(14.160595359570868*coeff[0]*fSkin[20])+9.077304717673634*coeff[0]*fEdge[20]+15.46875*coeff[0]*fSkin[12]+12.65625*coeff[0]*fEdge[12]-8.118988160479114*coeff[0]*fSkin[8]+8.118988160479114*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*fSkin[13]-0.65625*coeff[0]*fEdge[13]-15.612964114398654*coeff[0]*fSkin[5]-4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 6.7082039324993685*coeff[0]*fSkin[23]-6.7082039324993685*coeff[0]*fEdge[23]-8.118988160479114*coeff[0]*fSkin[18]-8.118988160479114*coeff[0]*fEdge[18]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = -(14.160595359570868*coeff[0]*fSkin[21])+9.077304717673634*coeff[0]*fEdge[21]+15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]-8.118988160479114*coeff[0]*fSkin[9]+8.118988160479114*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 6.7082039324993685*coeff[0]*fSkin[24]-6.7082039324993685*coeff[0]*fEdge[24]-8.118988160479114*coeff[0]*fSkin[19]-8.118988160479114*coeff[0]*fEdge[19]+4.6875*coeff[0]*fSkin[16]-4.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*fSkin[17]-0.65625*coeff[0]*fEdge[17]-15.61296411439865*coeff[0]*fSkin[10]-4.720198453190292*coeff[0]*fEdge[10]+4.192627457812107*coeff[0]*fSkin[6]-4.192627457812107*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = -(14.160595359570864*coeff[0]*fSkin[23])+9.077304717673634*coeff[0]*fEdge[23]+15.46875*coeff[0]*fSkin[18]+12.65625*coeff[0]*fEdge[18]-8.118988160479114*coeff[0]*fSkin[14]+8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = -(14.160595359570864*coeff[0]*fSkin[24])+9.077304717673634*coeff[0]*fEdge[24]+15.46875*coeff[0]*fSkin[19]+12.65625*coeff[0]*fEdge[19]-8.118988160479114*coeff[0]*fSkin[16]+8.118988160479114*coeff[0]*fEdge[16]; 
  edgeSurf_incr[20] = 20.34375*coeff[0]*fSkin[20]-0.65625*coeff[0]*fEdge[20]-15.612964114398654*coeff[0]*fSkin[12]-4.72019845319029*coeff[0]*fEdge[12]+4.192627457812107*coeff[0]*fSkin[8]-4.192627457812107*coeff[0]*fEdge[8]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*fSkin[21]-0.65625*coeff[0]*fEdge[21]-15.612964114398654*coeff[0]*fSkin[15]-4.72019845319029*coeff[0]*fEdge[15]+4.192627457812107*coeff[0]*fSkin[9]-4.192627457812107*coeff[0]*fEdge[9]; 
  edgeSurf_incr[22] = 6.708203932499369*coeff[0]*fSkin[26]-6.708203932499369*coeff[0]*fEdge[26]-8.11898816047911*coeff[0]*fSkin[25]-8.11898816047911*coeff[0]*fEdge[25]+4.6875*coeff[0]*fSkin[22]-4.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 20.34375*coeff[0]*fSkin[23]-0.65625*coeff[0]*fEdge[23]-15.61296411439865*coeff[0]*fSkin[18]-4.720198453190292*coeff[0]*fEdge[18]+4.192627457812107*coeff[0]*fSkin[14]-4.192627457812107*coeff[0]*fEdge[14]; 
  edgeSurf_incr[24] = 20.34375*coeff[0]*fSkin[24]-0.65625*coeff[0]*fEdge[24]-15.61296411439865*coeff[0]*fSkin[19]-4.720198453190292*coeff[0]*fEdge[19]+4.192627457812107*coeff[0]*fSkin[16]-4.192627457812107*coeff[0]*fEdge[16]; 
  edgeSurf_incr[25] = -(14.160595359570864*coeff[0]*fSkin[26])+9.077304717673634*coeff[0]*fEdge[26]+15.46875*coeff[0]*fSkin[25]+12.65625*coeff[0]*fEdge[25]-8.11898816047911*coeff[0]*fSkin[22]+8.11898816047911*coeff[0]*fEdge[22]; 
  edgeSurf_incr[26] = 20.34375*coeff[0]*fSkin[26]-0.65625*coeff[0]*fEdge[26]-15.61296411439865*coeff[0]*fSkin[25]-4.720198453190292*coeff[0]*fEdge[25]+4.192627457812107*coeff[0]*fSkin[22]-4.192627457812107*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[7]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 5.083290641897235*coeff[0]*fSkin[11]+2.8125*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[0]*fSkin[13]+2.8125*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*fSkin[7]+10.892765661208358*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = 5.083290641897234*coeff[0]*fSkin[17]+2.8125*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]+10.892765661208362*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 5.083290641897235*coeff[0]*fSkin[20]+2.8125*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*fSkin[13]+10.892765661208362*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = 5.083290641897235*coeff[0]*fSkin[21]+2.8125*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*fSkin[17]+10.892765661208358*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 5.083290641897234*coeff[0]*fSkin[23]+2.8125*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 5.083290641897234*coeff[0]*fSkin[24]+2.8125*coeff[0]*fSkin[19]; 
  boundSurf_incr[20] = 19.6875*coeff[0]*fSkin[20]+10.892765661208362*coeff[0]*fSkin[12]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*fSkin[21]+10.892765661208362*coeff[0]*fSkin[15]; 
  boundSurf_incr[23] = 19.6875*coeff[0]*fSkin[23]+10.892765661208358*coeff[0]*fSkin[18]; 
  boundSurf_incr[24] = 19.6875*coeff[0]*fSkin[24]+10.892765661208358*coeff[0]*fSkin[19]; 
  boundSurf_incr[25] = 5.083290641897234*coeff[0]*fSkin[26]+2.8125*coeff[0]*fSkin[25]; 
  boundSurf_incr[26] = 19.6875*coeff[0]*fSkin[26]+10.892765661208358*coeff[0]*fSkin[25]; 

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

