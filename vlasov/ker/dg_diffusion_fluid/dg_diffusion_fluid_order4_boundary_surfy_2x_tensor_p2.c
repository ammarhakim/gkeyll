#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfy_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  double vol_incr[9] = {0.0}; 

  double edgeSurf_incr[9] = {0.0}; 
  double boundSurf_incr[9] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[1]*fSkin[5]-6.708203932499369*coeff[1]*fEdge[5]+8.11898816047911*coeff[1]*fSkin[2]+8.11898816047911*coeff[1]*fEdge[2]+4.6875*fSkin[0]*coeff[1]-4.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[1]*fSkin[7]-6.7082039324993685*coeff[1]*fEdge[7]+8.11898816047911*coeff[1]*fSkin[3]+8.11898816047911*coeff[1]*fEdge[3]+4.6875*coeff[1]*fSkin[1]-4.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 14.160595359570864*coeff[1]*fSkin[5]-9.077304717673634*coeff[1]*fEdge[5]+15.46875*coeff[1]*fSkin[2]+12.65625*coeff[1]*fEdge[2]+8.11898816047911*fSkin[0]*coeff[1]-8.11898816047911*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 14.160595359570868*coeff[1]*fSkin[7]-9.077304717673634*coeff[1]*fEdge[7]+15.46875*coeff[1]*fSkin[3]+12.65625*coeff[1]*fEdge[3]+8.11898816047911*coeff[1]*fSkin[1]-8.11898816047911*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = 6.708203932499369*coeff[1]*fSkin[8]-6.708203932499369*coeff[1]*fEdge[8]+8.118988160479114*coeff[1]*fSkin[6]+8.118988160479114*coeff[1]*fEdge[6]+4.6875*coeff[1]*fSkin[4]-4.6875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 20.34375*coeff[1]*fSkin[5]-0.65625*coeff[1]*fEdge[5]+15.61296411439865*coeff[1]*fSkin[2]+4.720198453190292*coeff[1]*fEdge[2]+4.192627457812107*fSkin[0]*coeff[1]-4.192627457812107*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = 14.160595359570868*coeff[1]*fSkin[8]-9.077304717673634*coeff[1]*fEdge[8]+15.46875*coeff[1]*fSkin[6]+12.65625*coeff[1]*fEdge[6]+8.118988160479114*coeff[1]*fSkin[4]-8.118988160479114*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = 20.34375*coeff[1]*fSkin[7]-0.65625*coeff[1]*fEdge[7]+15.612964114398654*coeff[1]*fSkin[3]+4.72019845319029*coeff[1]*fEdge[3]+4.192627457812107*coeff[1]*fSkin[1]-4.192627457812107*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = 20.34375*coeff[1]*fSkin[8]-0.65625*coeff[1]*fEdge[8]+15.612964114398654*coeff[1]*fSkin[6]+4.72019845319029*coeff[1]*fEdge[6]+4.192627457812107*coeff[1]*fSkin[4]-4.192627457812107*coeff[1]*fEdge[4]; 

  boundSurf_incr[2] = 2.8125*coeff[1]*fSkin[2]-5.083290641897234*coeff[1]*fSkin[5]; 
  boundSurf_incr[3] = 2.8125*coeff[1]*fSkin[3]-5.083290641897235*coeff[1]*fSkin[7]; 
  boundSurf_incr[5] = 19.6875*coeff[1]*fSkin[5]-10.892765661208358*coeff[1]*fSkin[2]; 
  boundSurf_incr[6] = 2.8125*coeff[1]*fSkin[6]-5.083290641897235*coeff[1]*fSkin[8]; 
  boundSurf_incr[7] = 19.6875*coeff[1]*fSkin[7]-10.892765661208362*coeff[1]*fSkin[3]; 
  boundSurf_incr[8] = 19.6875*coeff[1]*fSkin[8]-10.892765661208362*coeff[1]*fSkin[6]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[1]*fSkin[5]-6.708203932499369*coeff[1]*fEdge[5]-8.11898816047911*coeff[1]*fSkin[2]-8.11898816047911*coeff[1]*fEdge[2]+4.6875*fSkin[0]*coeff[1]-4.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 6.7082039324993685*coeff[1]*fSkin[7]-6.7082039324993685*coeff[1]*fEdge[7]-8.11898816047911*coeff[1]*fSkin[3]-8.11898816047911*coeff[1]*fEdge[3]+4.6875*coeff[1]*fSkin[1]-4.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(14.160595359570864*coeff[1]*fSkin[5])+9.077304717673634*coeff[1]*fEdge[5]+15.46875*coeff[1]*fSkin[2]+12.65625*coeff[1]*fEdge[2]-8.11898816047911*fSkin[0]*coeff[1]+8.11898816047911*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(14.160595359570868*coeff[1]*fSkin[7])+9.077304717673634*coeff[1]*fEdge[7]+15.46875*coeff[1]*fSkin[3]+12.65625*coeff[1]*fEdge[3]-8.11898816047911*coeff[1]*fSkin[1]+8.11898816047911*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = 6.708203932499369*coeff[1]*fSkin[8]-6.708203932499369*coeff[1]*fEdge[8]-8.118988160479114*coeff[1]*fSkin[6]-8.118988160479114*coeff[1]*fEdge[6]+4.6875*coeff[1]*fSkin[4]-4.6875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 20.34375*coeff[1]*fSkin[5]-0.65625*coeff[1]*fEdge[5]-15.61296411439865*coeff[1]*fSkin[2]-4.720198453190292*coeff[1]*fEdge[2]+4.192627457812107*fSkin[0]*coeff[1]-4.192627457812107*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = -(14.160595359570868*coeff[1]*fSkin[8])+9.077304717673634*coeff[1]*fEdge[8]+15.46875*coeff[1]*fSkin[6]+12.65625*coeff[1]*fEdge[6]-8.118988160479114*coeff[1]*fSkin[4]+8.118988160479114*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = 20.34375*coeff[1]*fSkin[7]-0.65625*coeff[1]*fEdge[7]-15.612964114398654*coeff[1]*fSkin[3]-4.72019845319029*coeff[1]*fEdge[3]+4.192627457812107*coeff[1]*fSkin[1]-4.192627457812107*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = 20.34375*coeff[1]*fSkin[8]-0.65625*coeff[1]*fEdge[8]-15.612964114398654*coeff[1]*fSkin[6]-4.72019845319029*coeff[1]*fEdge[6]+4.192627457812107*coeff[1]*fSkin[4]-4.192627457812107*coeff[1]*fEdge[4]; 

  boundSurf_incr[2] = 5.083290641897234*coeff[1]*fSkin[5]+2.8125*coeff[1]*fSkin[2]; 
  boundSurf_incr[3] = 5.083290641897235*coeff[1]*fSkin[7]+2.8125*coeff[1]*fSkin[3]; 
  boundSurf_incr[5] = 19.6875*coeff[1]*fSkin[5]+10.892765661208358*coeff[1]*fSkin[2]; 
  boundSurf_incr[6] = 5.083290641897235*coeff[1]*fSkin[8]+2.8125*coeff[1]*fSkin[6]; 
  boundSurf_incr[7] = 19.6875*coeff[1]*fSkin[7]+10.892765661208362*coeff[1]*fSkin[3]; 
  boundSurf_incr[8] = 19.6875*coeff[1]*fSkin[8]+10.892765661208362*coeff[1]*fSkin[6]; 

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

  return 0.;
}

