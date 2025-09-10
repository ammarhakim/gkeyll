#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*qSkin[4]-6.708203932499369*coeff[0]*qEdge[4]+8.11898816047911*coeff[0]*qSkin[1]+8.11898816047911*coeff[0]*qEdge[1]+4.6875*coeff[0]*qSkin[0]-4.6875*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 14.160595359570864*coeff[0]*qSkin[4]-9.077304717673634*coeff[0]*qEdge[4]+15.46875*coeff[0]*qSkin[1]+12.65625*coeff[0]*qEdge[1]+8.11898816047911*coeff[0]*qSkin[0]-8.11898816047911*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*qSkin[6]-6.7082039324993685*coeff[0]*qEdge[6]+8.11898816047911*coeff[0]*qSkin[3]+8.11898816047911*coeff[0]*qEdge[3]+4.6875*coeff[0]*qSkin[2]-4.6875*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 14.160595359570868*coeff[0]*qSkin[6]-9.077304717673634*coeff[0]*qEdge[6]+15.46875*coeff[0]*qSkin[3]+12.65625*coeff[0]*qEdge[3]+8.11898816047911*coeff[0]*qSkin[2]-8.11898816047911*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*qSkin[4]-0.65625*coeff[0]*qEdge[4]+15.61296411439865*coeff[0]*qSkin[1]+4.720198453190292*coeff[0]*qEdge[1]+4.192627457812107*coeff[0]*qSkin[0]-4.192627457812107*coeff[0]*qEdge[0]; 
  edgeSurf_incr[5] = 8.118988160479114*coeff[0]*qSkin[7]+8.118988160479114*coeff[0]*qEdge[7]+4.6875*coeff[0]*qSkin[5]-4.6875*coeff[0]*qEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*qSkin[6]-0.65625*coeff[0]*qEdge[6]+15.612964114398654*coeff[0]*qSkin[3]+4.72019845319029*coeff[0]*qEdge[3]+4.192627457812107*coeff[0]*qSkin[2]-4.192627457812107*coeff[0]*qEdge[2]; 
  edgeSurf_incr[7] = 15.46875*coeff[0]*qSkin[7]+12.65625*coeff[0]*qEdge[7]+8.118988160479114*coeff[0]*qSkin[5]-8.118988160479114*coeff[0]*qEdge[5]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*qSkin[1]-5.083290641897234*coeff[0]*qSkin[4]; 
  boundSurf_incr[3] = 2.8125*coeff[0]*qSkin[3]-5.083290641897235*coeff[0]*qSkin[6]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*qSkin[4]-10.892765661208358*coeff[0]*qSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*qSkin[6]-10.892765661208362*coeff[0]*qSkin[3]; 
  boundSurf_incr[7] = 2.8125*coeff[0]*qSkin[7]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*qSkin[4]-6.708203932499369*coeff[0]*qEdge[4]-8.11898816047911*coeff[0]*qSkin[1]-8.11898816047911*coeff[0]*qEdge[1]+4.6875*coeff[0]*qSkin[0]-4.6875*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(14.160595359570864*coeff[0]*qSkin[4])+9.077304717673634*coeff[0]*qEdge[4]+15.46875*coeff[0]*qSkin[1]+12.65625*coeff[0]*qEdge[1]-8.11898816047911*coeff[0]*qSkin[0]+8.11898816047911*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*qSkin[6]-6.7082039324993685*coeff[0]*qEdge[6]-8.11898816047911*coeff[0]*qSkin[3]-8.11898816047911*coeff[0]*qEdge[3]+4.6875*coeff[0]*qSkin[2]-4.6875*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = -(14.160595359570868*coeff[0]*qSkin[6])+9.077304717673634*coeff[0]*qEdge[6]+15.46875*coeff[0]*qSkin[3]+12.65625*coeff[0]*qEdge[3]-8.11898816047911*coeff[0]*qSkin[2]+8.11898816047911*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*qSkin[4]-0.65625*coeff[0]*qEdge[4]-15.61296411439865*coeff[0]*qSkin[1]-4.720198453190292*coeff[0]*qEdge[1]+4.192627457812107*coeff[0]*qSkin[0]-4.192627457812107*coeff[0]*qEdge[0]; 
  edgeSurf_incr[5] = -(8.118988160479114*coeff[0]*qSkin[7])-8.118988160479114*coeff[0]*qEdge[7]+4.6875*coeff[0]*qSkin[5]-4.6875*coeff[0]*qEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*qSkin[6]-0.65625*coeff[0]*qEdge[6]-15.612964114398654*coeff[0]*qSkin[3]-4.72019845319029*coeff[0]*qEdge[3]+4.192627457812107*coeff[0]*qSkin[2]-4.192627457812107*coeff[0]*qEdge[2]; 
  edgeSurf_incr[7] = 15.46875*coeff[0]*qSkin[7]+12.65625*coeff[0]*qEdge[7]-8.118988160479114*coeff[0]*qSkin[5]+8.118988160479114*coeff[0]*qEdge[5]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*qSkin[4]+2.8125*coeff[0]*qSkin[1]; 
  boundSurf_incr[3] = 5.083290641897235*coeff[0]*qSkin[6]+2.8125*coeff[0]*qSkin[3]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*qSkin[4]+10.892765661208358*coeff[0]*qSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*qSkin[6]+10.892765661208362*coeff[0]*qSkin[3]; 
  boundSurf_incr[7] = 2.8125*coeff[0]*qSkin[7]; 

  }

  out[0] += -(1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq); 
  out[1] += -(1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq); 
  out[2] += -(1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq); 
  out[3] += -(1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq); 
  out[4] += -(1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq); 
  out[5] += -(1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq); 
  out[6] += -(1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq); 
  out[7] += -(1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq); 

  return 0.;
}

