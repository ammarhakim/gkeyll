#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[9] = {0.0}; 

  double edgeSurf_incr[9] = {0.0}; 
  double boundSurf_incr[9] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[4]-6.708203932499369*coeff[0]*fEdge[4]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.160595359570864*coeff[0]*fSkin[4]-9.077304717673634*coeff[0]*fEdge[4]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*fSkin[6]-6.7082039324993685*coeff[0]*fEdge[6]+8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 14.160595359570868*coeff[0]*fSkin[6]-9.077304717673634*coeff[0]*fEdge[6]+15.46875*coeff[0]*fSkin[3]+12.65625*coeff[0]*fEdge[3]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*fSkin[4]-0.65625*coeff[0]*fEdge[4]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 6.708203932499369*coeff[0]*fSkin[8]-6.708203932499369*coeff[0]*fEdge[8]+8.118988160479114*coeff[0]*fSkin[7]+8.118988160479114*coeff[0]*fEdge[7]+4.6875*coeff[0]*fSkin[5]-4.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*fSkin[6]-0.65625*coeff[0]*fEdge[6]+15.612964114398654*coeff[0]*fSkin[3]+4.72019845319029*coeff[0]*fEdge[3]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 14.160595359570868*coeff[0]*fSkin[8]-9.077304717673634*coeff[0]*fEdge[8]+15.46875*coeff[0]*fSkin[7]+12.65625*coeff[0]*fEdge[7]+8.118988160479114*coeff[0]*fSkin[5]-8.118988160479114*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = 20.34375*coeff[0]*fSkin[8]-0.65625*coeff[0]*fEdge[8]+15.612964114398654*coeff[0]*fSkin[7]+4.72019845319029*coeff[0]*fEdge[7]+4.192627457812107*coeff[0]*fSkin[5]-4.192627457812107*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[4]; 
  boundSurf_incr[3] = 2.8125*coeff[0]*fSkin[3]-5.083290641897235*coeff[0]*fSkin[6]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*fSkin[4]-10.892765661208358*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*fSkin[6]-10.892765661208362*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 2.8125*coeff[0]*fSkin[7]-5.083290641897235*coeff[0]*fSkin[8]; 
  boundSurf_incr[8] = 19.6875*coeff[0]*fSkin[8]-10.892765661208362*coeff[0]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[4]-6.708203932499369*coeff[0]*fEdge[4]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(14.160595359570864*coeff[0]*fSkin[4])+9.077304717673634*coeff[0]*fEdge[4]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.7082039324993685*coeff[0]*fSkin[6]-6.7082039324993685*coeff[0]*fEdge[6]-8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(14.160595359570868*coeff[0]*fSkin[6])+9.077304717673634*coeff[0]*fEdge[6]+15.46875*coeff[0]*fSkin[3]+12.65625*coeff[0]*fEdge[3]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*fSkin[4]-0.65625*coeff[0]*fEdge[4]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 6.708203932499369*coeff[0]*fSkin[8]-6.708203932499369*coeff[0]*fEdge[8]-8.118988160479114*coeff[0]*fSkin[7]-8.118988160479114*coeff[0]*fEdge[7]+4.6875*coeff[0]*fSkin[5]-4.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*fSkin[6]-0.65625*coeff[0]*fEdge[6]-15.612964114398654*coeff[0]*fSkin[3]-4.72019845319029*coeff[0]*fEdge[3]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = -(14.160595359570868*coeff[0]*fSkin[8])+9.077304717673634*coeff[0]*fEdge[8]+15.46875*coeff[0]*fSkin[7]+12.65625*coeff[0]*fEdge[7]-8.118988160479114*coeff[0]*fSkin[5]+8.118988160479114*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = 20.34375*coeff[0]*fSkin[8]-0.65625*coeff[0]*fEdge[8]-15.612964114398654*coeff[0]*fSkin[7]-4.72019845319029*coeff[0]*fEdge[7]+4.192627457812107*coeff[0]*fSkin[5]-4.192627457812107*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[4]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 5.083290641897235*coeff[0]*fSkin[6]+2.8125*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*fSkin[4]+10.892765661208358*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*fSkin[6]+10.892765661208362*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 5.083290641897235*coeff[0]*fSkin[8]+2.8125*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = 19.6875*coeff[0]*fSkin[8]+10.892765661208362*coeff[0]*fSkin[7]; 

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

