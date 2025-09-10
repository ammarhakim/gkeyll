#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double vol_incr[20] = {0.0}; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*qSkin[7]-6.708203932499369*coeff[0]*qEdge[7]+8.11898816047911*coeff[0]*qSkin[1]+8.11898816047911*coeff[0]*qEdge[1]+4.6875*coeff[0]*qSkin[0]-4.6875*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 14.160595359570864*coeff[0]*qSkin[7]-9.077304717673634*coeff[0]*qEdge[7]+15.46875*coeff[0]*qSkin[1]+12.65625*coeff[0]*qEdge[1]+8.11898816047911*coeff[0]*qSkin[0]-8.11898816047911*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*qSkin[11]-6.7082039324993685*coeff[0]*qEdge[11]+8.11898816047911*coeff[0]*qSkin[4]+8.11898816047911*coeff[0]*qEdge[4]+4.6875*coeff[0]*qSkin[2]-4.6875*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*qSkin[13]-6.7082039324993685*coeff[0]*qEdge[13]+8.11898816047911*coeff[0]*qSkin[5]+8.11898816047911*coeff[0]*qEdge[5]+4.6875*coeff[0]*qSkin[3]-4.6875*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = 14.160595359570868*coeff[0]*qSkin[11]-9.077304717673634*coeff[0]*qEdge[11]+15.46875*coeff[0]*qSkin[4]+12.65625*coeff[0]*qEdge[4]+8.11898816047911*coeff[0]*qSkin[2]-8.11898816047911*coeff[0]*qEdge[2]; 
  edgeSurf_incr[5] = 14.160595359570868*coeff[0]*qSkin[13]-9.077304717673634*coeff[0]*qEdge[13]+15.46875*coeff[0]*qSkin[5]+12.65625*coeff[0]*qEdge[5]+8.11898816047911*coeff[0]*qSkin[3]-8.11898816047911*coeff[0]*qEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*qSkin[17]-6.708203932499369*coeff[0]*qEdge[17]+8.11898816047911*coeff[0]*qSkin[10]+8.11898816047911*coeff[0]*qEdge[10]+4.6875*coeff[0]*qSkin[6]-4.6875*coeff[0]*qEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*qSkin[7]-0.65625*coeff[0]*qEdge[7]+15.61296411439865*coeff[0]*qSkin[1]+4.720198453190292*coeff[0]*qEdge[1]+4.192627457812107*coeff[0]*qSkin[0]-4.192627457812107*coeff[0]*qEdge[0]; 
  edgeSurf_incr[8] = 8.118988160479114*coeff[0]*qSkin[12]+8.118988160479114*coeff[0]*qEdge[12]+4.6875*coeff[0]*qSkin[8]-4.6875*coeff[0]*qEdge[8]; 
  edgeSurf_incr[9] = 8.118988160479114*coeff[0]*qSkin[15]+8.118988160479114*coeff[0]*qEdge[15]+4.6875*coeff[0]*qSkin[9]-4.6875*coeff[0]*qEdge[9]; 
  edgeSurf_incr[10] = 14.160595359570864*coeff[0]*qSkin[17]-9.077304717673634*coeff[0]*qEdge[17]+15.46875*coeff[0]*qSkin[10]+12.65625*coeff[0]*qEdge[10]+8.11898816047911*coeff[0]*qSkin[6]-8.11898816047911*coeff[0]*qEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*qSkin[11]-0.65625*coeff[0]*qEdge[11]+15.612964114398654*coeff[0]*qSkin[4]+4.72019845319029*coeff[0]*qEdge[4]+4.192627457812107*coeff[0]*qSkin[2]-4.192627457812107*coeff[0]*qEdge[2]; 
  edgeSurf_incr[12] = 15.46875*coeff[0]*qSkin[12]+12.65625*coeff[0]*qEdge[12]+8.118988160479114*coeff[0]*qSkin[8]-8.118988160479114*coeff[0]*qEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*qSkin[13]-0.65625*coeff[0]*qEdge[13]+15.612964114398654*coeff[0]*qSkin[5]+4.72019845319029*coeff[0]*qEdge[5]+4.192627457812107*coeff[0]*qSkin[3]-4.192627457812107*coeff[0]*qEdge[3]; 
  edgeSurf_incr[14] = 8.118988160479114*coeff[0]*qSkin[18]+8.118988160479114*coeff[0]*qEdge[18]+4.6875*coeff[0]*qSkin[14]-4.6875*coeff[0]*qEdge[14]; 
  edgeSurf_incr[15] = 15.46875*coeff[0]*qSkin[15]+12.65625*coeff[0]*qEdge[15]+8.118988160479114*coeff[0]*qSkin[9]-8.118988160479114*coeff[0]*qEdge[9]; 
  edgeSurf_incr[16] = 8.118988160479114*coeff[0]*qSkin[19]+8.118988160479114*coeff[0]*qEdge[19]+4.6875*coeff[0]*qSkin[16]-4.6875*coeff[0]*qEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*qSkin[17]-0.65625*coeff[0]*qEdge[17]+15.61296411439865*coeff[0]*qSkin[10]+4.720198453190292*coeff[0]*qEdge[10]+4.192627457812107*coeff[0]*qSkin[6]-4.192627457812107*coeff[0]*qEdge[6]; 
  edgeSurf_incr[18] = 15.46875*coeff[0]*qSkin[18]+12.65625*coeff[0]*qEdge[18]+8.118988160479114*coeff[0]*qSkin[14]-8.118988160479114*coeff[0]*qEdge[14]; 
  edgeSurf_incr[19] = 15.46875*coeff[0]*qSkin[19]+12.65625*coeff[0]*qEdge[19]+8.118988160479114*coeff[0]*qSkin[16]-8.118988160479114*coeff[0]*qEdge[16]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*qSkin[1]-5.083290641897234*coeff[0]*qSkin[7]; 
  boundSurf_incr[4] = 2.8125*coeff[0]*qSkin[4]-5.083290641897235*coeff[0]*qSkin[11]; 
  boundSurf_incr[5] = 2.8125*coeff[0]*qSkin[5]-5.083290641897235*coeff[0]*qSkin[13]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*qSkin[7]-10.892765661208358*coeff[0]*qSkin[1]; 
  boundSurf_incr[10] = 2.8125*coeff[0]*qSkin[10]-5.083290641897234*coeff[0]*qSkin[17]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*qSkin[11]-10.892765661208362*coeff[0]*qSkin[4]; 
  boundSurf_incr[12] = 2.8125*coeff[0]*qSkin[12]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*qSkin[13]-10.892765661208362*coeff[0]*qSkin[5]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*qSkin[15]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*qSkin[17]-10.892765661208358*coeff[0]*qSkin[10]; 
  boundSurf_incr[18] = 2.8125*coeff[0]*qSkin[18]; 
  boundSurf_incr[19] = 2.8125*coeff[0]*qSkin[19]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*qSkin[7]-6.708203932499369*coeff[0]*qEdge[7]-8.11898816047911*coeff[0]*qSkin[1]-8.11898816047911*coeff[0]*qEdge[1]+4.6875*coeff[0]*qSkin[0]-4.6875*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(14.160595359570864*coeff[0]*qSkin[7])+9.077304717673634*coeff[0]*qEdge[7]+15.46875*coeff[0]*qSkin[1]+12.65625*coeff[0]*qEdge[1]-8.11898816047911*coeff[0]*qSkin[0]+8.11898816047911*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*qSkin[11]-6.7082039324993685*coeff[0]*qEdge[11]-8.11898816047911*coeff[0]*qSkin[4]-8.11898816047911*coeff[0]*qEdge[4]+4.6875*coeff[0]*qSkin[2]-4.6875*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 6.7082039324993685*coeff[0]*qSkin[13]-6.7082039324993685*coeff[0]*qEdge[13]-8.11898816047911*coeff[0]*qSkin[5]-8.11898816047911*coeff[0]*qEdge[5]+4.6875*coeff[0]*qSkin[3]-4.6875*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = -(14.160595359570868*coeff[0]*qSkin[11])+9.077304717673634*coeff[0]*qEdge[11]+15.46875*coeff[0]*qSkin[4]+12.65625*coeff[0]*qEdge[4]-8.11898816047911*coeff[0]*qSkin[2]+8.11898816047911*coeff[0]*qEdge[2]; 
  edgeSurf_incr[5] = -(14.160595359570868*coeff[0]*qSkin[13])+9.077304717673634*coeff[0]*qEdge[13]+15.46875*coeff[0]*qSkin[5]+12.65625*coeff[0]*qEdge[5]-8.11898816047911*coeff[0]*qSkin[3]+8.11898816047911*coeff[0]*qEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*qSkin[17]-6.708203932499369*coeff[0]*qEdge[17]-8.11898816047911*coeff[0]*qSkin[10]-8.11898816047911*coeff[0]*qEdge[10]+4.6875*coeff[0]*qSkin[6]-4.6875*coeff[0]*qEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*qSkin[7]-0.65625*coeff[0]*qEdge[7]-15.61296411439865*coeff[0]*qSkin[1]-4.720198453190292*coeff[0]*qEdge[1]+4.192627457812107*coeff[0]*qSkin[0]-4.192627457812107*coeff[0]*qEdge[0]; 
  edgeSurf_incr[8] = -(8.118988160479114*coeff[0]*qSkin[12])-8.118988160479114*coeff[0]*qEdge[12]+4.6875*coeff[0]*qSkin[8]-4.6875*coeff[0]*qEdge[8]; 
  edgeSurf_incr[9] = -(8.118988160479114*coeff[0]*qSkin[15])-8.118988160479114*coeff[0]*qEdge[15]+4.6875*coeff[0]*qSkin[9]-4.6875*coeff[0]*qEdge[9]; 
  edgeSurf_incr[10] = -(14.160595359570864*coeff[0]*qSkin[17])+9.077304717673634*coeff[0]*qEdge[17]+15.46875*coeff[0]*qSkin[10]+12.65625*coeff[0]*qEdge[10]-8.11898816047911*coeff[0]*qSkin[6]+8.11898816047911*coeff[0]*qEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*qSkin[11]-0.65625*coeff[0]*qEdge[11]-15.612964114398654*coeff[0]*qSkin[4]-4.72019845319029*coeff[0]*qEdge[4]+4.192627457812107*coeff[0]*qSkin[2]-4.192627457812107*coeff[0]*qEdge[2]; 
  edgeSurf_incr[12] = 15.46875*coeff[0]*qSkin[12]+12.65625*coeff[0]*qEdge[12]-8.118988160479114*coeff[0]*qSkin[8]+8.118988160479114*coeff[0]*qEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*qSkin[13]-0.65625*coeff[0]*qEdge[13]-15.612964114398654*coeff[0]*qSkin[5]-4.72019845319029*coeff[0]*qEdge[5]+4.192627457812107*coeff[0]*qSkin[3]-4.192627457812107*coeff[0]*qEdge[3]; 
  edgeSurf_incr[14] = -(8.118988160479114*coeff[0]*qSkin[18])-8.118988160479114*coeff[0]*qEdge[18]+4.6875*coeff[0]*qSkin[14]-4.6875*coeff[0]*qEdge[14]; 
  edgeSurf_incr[15] = 15.46875*coeff[0]*qSkin[15]+12.65625*coeff[0]*qEdge[15]-8.118988160479114*coeff[0]*qSkin[9]+8.118988160479114*coeff[0]*qEdge[9]; 
  edgeSurf_incr[16] = -(8.118988160479114*coeff[0]*qSkin[19])-8.118988160479114*coeff[0]*qEdge[19]+4.6875*coeff[0]*qSkin[16]-4.6875*coeff[0]*qEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*qSkin[17]-0.65625*coeff[0]*qEdge[17]-15.61296411439865*coeff[0]*qSkin[10]-4.720198453190292*coeff[0]*qEdge[10]+4.192627457812107*coeff[0]*qSkin[6]-4.192627457812107*coeff[0]*qEdge[6]; 
  edgeSurf_incr[18] = 15.46875*coeff[0]*qSkin[18]+12.65625*coeff[0]*qEdge[18]-8.118988160479114*coeff[0]*qSkin[14]+8.118988160479114*coeff[0]*qEdge[14]; 
  edgeSurf_incr[19] = 15.46875*coeff[0]*qSkin[19]+12.65625*coeff[0]*qEdge[19]-8.118988160479114*coeff[0]*qSkin[16]+8.118988160479114*coeff[0]*qEdge[16]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*qSkin[7]+2.8125*coeff[0]*qSkin[1]; 
  boundSurf_incr[4] = 5.083290641897235*coeff[0]*qSkin[11]+2.8125*coeff[0]*qSkin[4]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[0]*qSkin[13]+2.8125*coeff[0]*qSkin[5]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*qSkin[7]+10.892765661208358*coeff[0]*qSkin[1]; 
  boundSurf_incr[10] = 5.083290641897234*coeff[0]*qSkin[17]+2.8125*coeff[0]*qSkin[10]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*qSkin[11]+10.892765661208362*coeff[0]*qSkin[4]; 
  boundSurf_incr[12] = 2.8125*coeff[0]*qSkin[12]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*qSkin[13]+10.892765661208362*coeff[0]*qSkin[5]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*qSkin[15]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*qSkin[17]+10.892765661208358*coeff[0]*qSkin[10]; 
  boundSurf_incr[18] = 2.8125*coeff[0]*qSkin[18]; 
  boundSurf_incr[19] = 2.8125*coeff[0]*qSkin[19]; 

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

  return 0.;
}

