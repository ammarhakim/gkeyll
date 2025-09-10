#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_2x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[9] = {0.0}; 

  double edgeSurf_incr[9] = {0.0}; 
  double boundSurf_incr[9] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[0]*fSkin[4])+35.21807064562169*coeff[0]*fEdge[4]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(70.53065765632411*coeff[0]*fSkin[4])+51.46831774920949*coeff[0]*fEdge[4]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(35.21807064562168*coeff[0]*fSkin[6])+35.21807064562168*coeff[0]*fEdge[6]-34.09975027401226*coeff[0]*fSkin[3]-34.09975027401226*coeff[0]*fEdge[3]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(70.53065765632414*coeff[0]*fSkin[6])+51.468317749209504*coeff[0]*fEdge[6]-61.5234375*coeff[0]*fSkin[3]-56.6015625*coeff[0]*fEdge[3]-34.09975027401226*coeff[0]*fSkin[2]+34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(70.6640625*coeff[0]*fSkin[4])-3.1640625*coeff[0]*fEdge[4]-31.316701275974005*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = -(35.21807064562169*coeff[0]*fSkin[8])+35.21807064562169*coeff[0]*fEdge[8]-34.09975027401227*coeff[0]*fSkin[7]-34.09975027401227*coeff[0]*fEdge[7]-19.6875*coeff[0]*fSkin[5]+19.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = -(70.6640625*coeff[0]*fSkin[6])-3.1640625*coeff[0]*fEdge[6]-31.316701275974033*coeff[0]*fSkin[3]-12.2543613688594*coeff[0]*fEdge[3]-12.577882373436315*coeff[0]*fSkin[2]+12.577882373436315*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = -(70.53065765632414*coeff[0]*fSkin[8])+51.468317749209504*coeff[0]*fEdge[8]-61.5234375*coeff[0]*fSkin[7]-56.6015625*coeff[0]*fEdge[7]-34.09975027401227*coeff[0]*fSkin[5]+34.09975027401227*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = -(70.6640625*coeff[0]*fSkin[8])-3.1640625*coeff[0]*fEdge[8]-31.316701275974033*coeff[0]*fSkin[7]-12.2543613688594*coeff[0]*fEdge[7]-12.57788237343632*coeff[0]*fSkin[5]+12.57788237343632*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 19.062339907114627*coeff[0]*fSkin[4]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 19.062339907114634*coeff[0]*fSkin[6]-4.921875*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = 19.062339907114627*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[4]; 
  boundSurf_incr[6] = 19.062339907114634*coeff[0]*fSkin[3]-73.828125*coeff[0]*fSkin[6]; 
  boundSurf_incr[7] = 19.062339907114634*coeff[0]*fSkin[8]-4.921875*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = 19.062339907114634*coeff[0]*fSkin[7]-73.828125*coeff[0]*fSkin[8]; 

  } else { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[0]*fSkin[4])+35.21807064562169*coeff[0]*fEdge[4]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[4]-51.46831774920949*coeff[0]*fEdge[4]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(35.21807064562168*coeff[0]*fSkin[6])+35.21807064562168*coeff[0]*fEdge[6]+34.09975027401226*coeff[0]*fSkin[3]+34.09975027401226*coeff[0]*fEdge[3]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 70.53065765632414*coeff[0]*fSkin[6]-51.468317749209504*coeff[0]*fEdge[6]-61.5234375*coeff[0]*fSkin[3]-56.6015625*coeff[0]*fEdge[3]+34.09975027401226*coeff[0]*fSkin[2]-34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(70.6640625*coeff[0]*fSkin[4])-3.1640625*coeff[0]*fEdge[4]+31.316701275974005*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = -(35.21807064562169*coeff[0]*fSkin[8])+35.21807064562169*coeff[0]*fEdge[8]+34.09975027401227*coeff[0]*fSkin[7]+34.09975027401227*coeff[0]*fEdge[7]-19.6875*coeff[0]*fSkin[5]+19.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = -(70.6640625*coeff[0]*fSkin[6])-3.1640625*coeff[0]*fEdge[6]+31.316701275974033*coeff[0]*fSkin[3]+12.2543613688594*coeff[0]*fEdge[3]-12.577882373436315*coeff[0]*fSkin[2]+12.577882373436315*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 70.53065765632414*coeff[0]*fSkin[8]-51.468317749209504*coeff[0]*fEdge[8]-61.5234375*coeff[0]*fSkin[7]-56.6015625*coeff[0]*fEdge[7]+34.09975027401227*coeff[0]*fSkin[5]-34.09975027401227*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = -(70.6640625*coeff[0]*fSkin[8])-3.1640625*coeff[0]*fEdge[8]+31.316701275974033*coeff[0]*fSkin[7]+12.2543613688594*coeff[0]*fEdge[7]-12.57788237343632*coeff[0]*fSkin[5]+12.57788237343632*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = -(19.062339907114627*coeff[0]*fSkin[4])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = -(19.062339907114634*coeff[0]*fSkin[6])-4.921875*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = -(73.828125*coeff[0]*fSkin[4])-19.062339907114627*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = -(73.828125*coeff[0]*fSkin[6])-19.062339907114634*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(19.062339907114634*coeff[0]*fSkin[8])-4.921875*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = -(73.828125*coeff[0]*fSkin[8])-19.062339907114634*coeff[0]*fSkin[7]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 

  return 0.;
}

