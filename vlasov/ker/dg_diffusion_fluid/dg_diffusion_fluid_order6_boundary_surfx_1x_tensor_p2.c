#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[3] = {0.0}; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[0]*fSkin[2])+35.21807064562169*coeff[0]*fEdge[2]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(70.53065765632411*coeff[0]*fSkin[2])+51.46831774920949*coeff[0]*fEdge[2]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(70.6640625*coeff[0]*fSkin[2])-3.1640625*coeff[0]*fEdge[2]-31.316701275974005*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 19.062339907114627*coeff[0]*fSkin[2]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 19.062339907114627*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[2]; 

  } else { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[0]*fSkin[2])+35.21807064562169*coeff[0]*fEdge[2]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[2]-51.46831774920949*coeff[0]*fEdge[2]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(70.6640625*coeff[0]*fSkin[2])-3.1640625*coeff[0]*fEdge[2]+31.316701275974005*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = -(19.062339907114627*coeff[0]*fSkin[2])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = -(73.828125*coeff[0]*fSkin[2])-19.062339907114627*coeff[0]*fSkin[1]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  return 0.;
}

