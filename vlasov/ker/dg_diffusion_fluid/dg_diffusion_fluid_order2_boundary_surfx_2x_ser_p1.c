#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[4] = {0.0}; 

  double edgeSurf_incr[4] = {0.0}; 
  double boundSurf_incr[4] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[0]*fSkin[1])-0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]-1.4072912811497125*coeff[0]*fSkin[0]+0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.5412658773652741*coeff[0]*fSkin[3])-0.5412658773652741*coeff[0]*fEdge[3]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(1.4375*coeff[0]*fSkin[3])-0.4375*coeff[0]*fEdge[3]-1.4072912811497125*coeff[0]*fSkin[2]+0.5412658773652739*coeff[0]*fEdge[2]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*fSkin[0]-1.0*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 0.8660254037844386*coeff[0]*fSkin[2]-1.0*coeff[0]*fSkin[3]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[0]*fSkin[1]+0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]+1.4072912811497125*coeff[0]*fSkin[0]-0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[0]*fSkin[3]+0.5412658773652741*coeff[0]*fEdge[3]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(1.4375*coeff[0]*fSkin[3])-0.4375*coeff[0]*fEdge[3]+1.4072912811497125*coeff[0]*fSkin[2]-0.5412658773652739*coeff[0]*fEdge[2]; 

  boundSurf_incr[1] = -(1.0*coeff[0]*fSkin[1])-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(1.0*coeff[0]*fSkin[3])-0.8660254037844386*coeff[0]*fSkin[2]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[4] = {0.0}; 
  vol_incr[1] = 1.5*fSkin[2]*coeff[3]+1.5*fSkin[0]*coeff[1]; 
  vol_incr[3] = 1.5*fSkin[0]*coeff[3]+1.5*coeff[1]*fSkin[2]; 

  double edgeSurf_incr[4] = {0.0}; 
  double boundSurf_incr[4] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.27063293868263705*coeff[2]*fSkin[3])-0.27063293868263705*coeff[2]*fEdge[3]-0.28125*coeff[2]*fSkin[2]+0.28125*coeff[2]*fEdge[2]-0.27063293868263705*coeff[0]*fSkin[1]-0.27063293868263705*coeff[0]*fEdge[1]-0.28125*coeff[0]*fSkin[0]+0.28125*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(0.4330127018922193*coeff[3]*fSkin[3])-0.71875*coeff[2]*fSkin[3]+0.4330127018922193*coeff[3]*fEdge[3]-0.21875*coeff[2]*fEdge[3]-0.375*fSkin[2]*coeff[3]-0.375*fEdge[2]*coeff[3]-0.7036456405748562*coeff[2]*fSkin[2]+0.27063293868263694*coeff[2]*fEdge[2]-0.4330127018922193*coeff[1]*fSkin[1]-0.71875*coeff[0]*fSkin[1]+0.4330127018922193*coeff[1]*fEdge[1]-0.21875*coeff[0]*fEdge[1]-0.375*fSkin[0]*coeff[1]-0.375*fEdge[0]*coeff[1]-0.7036456405748562*coeff[0]*fSkin[0]+0.27063293868263694*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.27063293868263705*coeff[0]*fSkin[3])-0.27063293868263705*coeff[0]*fEdge[3]-0.28125*coeff[0]*fSkin[2]+0.28125*coeff[0]*fEdge[2]-0.27063293868263705*fSkin[1]*coeff[2]-0.27063293868263705*fEdge[1]*coeff[2]-0.28125*fSkin[0]*coeff[2]+0.28125*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = -(0.4330127018922193*coeff[1]*fSkin[3])-0.71875*coeff[0]*fSkin[3]+0.4330127018922193*coeff[1]*fEdge[3]-0.21875*coeff[0]*fEdge[3]-0.4330127018922193*fSkin[1]*coeff[3]+0.4330127018922193*fEdge[1]*coeff[3]-0.375*fSkin[0]*coeff[3]-0.375*fEdge[0]*coeff[3]-0.375*coeff[1]*fSkin[2]-0.7036456405748562*coeff[0]*fSkin[2]-0.375*coeff[1]*fEdge[2]+0.27063293868263694*coeff[0]*fEdge[2]-0.71875*fSkin[1]*coeff[2]-0.21875*fEdge[1]*coeff[2]-0.7036456405748562*fSkin[0]*coeff[2]+0.27063293868263694*fEdge[0]*coeff[2]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[3]*fSkin[3]-0.5*coeff[2]*fSkin[3]-0.75*fSkin[2]*coeff[3]+0.4330127018922193*coeff[2]*fSkin[2]+0.8660254037844386*coeff[1]*fSkin[1]-0.5*coeff[0]*fSkin[1]-0.75*fSkin[0]*coeff[1]+0.4330127018922193*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = 0.8660254037844386*coeff[1]*fSkin[3]-0.5*coeff[0]*fSkin[3]+0.8660254037844386*fSkin[1]*coeff[3]-0.75*fSkin[0]*coeff[3]-0.75*coeff[1]*fSkin[2]+0.4330127018922193*coeff[0]*fSkin[2]-0.5*fSkin[1]*coeff[2]+0.4330127018922193*fSkin[0]*coeff[2]; 

  } else { 

  edgeSurf_incr[0] = 0.27063293868263705*coeff[2]*fSkin[3]+0.27063293868263705*coeff[2]*fEdge[3]-0.28125*coeff[2]*fSkin[2]+0.28125*coeff[2]*fEdge[2]+0.27063293868263705*coeff[0]*fSkin[1]+0.27063293868263705*coeff[0]*fEdge[1]-0.28125*coeff[0]*fSkin[0]+0.28125*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.4330127018922193*coeff[3]*fSkin[3]-0.71875*coeff[2]*fSkin[3]-0.4330127018922193*coeff[3]*fEdge[3]-0.21875*coeff[2]*fEdge[3]-0.375*fSkin[2]*coeff[3]-0.375*fEdge[2]*coeff[3]+0.7036456405748562*coeff[2]*fSkin[2]-0.27063293868263694*coeff[2]*fEdge[2]+0.4330127018922193*coeff[1]*fSkin[1]-0.71875*coeff[0]*fSkin[1]-0.4330127018922193*coeff[1]*fEdge[1]-0.21875*coeff[0]*fEdge[1]-0.375*fSkin[0]*coeff[1]-0.375*fEdge[0]*coeff[1]+0.7036456405748562*coeff[0]*fSkin[0]-0.27063293868263694*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.27063293868263705*coeff[0]*fSkin[3]+0.27063293868263705*coeff[0]*fEdge[3]-0.28125*coeff[0]*fSkin[2]+0.28125*coeff[0]*fEdge[2]+0.27063293868263705*fSkin[1]*coeff[2]+0.27063293868263705*fEdge[1]*coeff[2]-0.28125*fSkin[0]*coeff[2]+0.28125*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = 0.4330127018922193*coeff[1]*fSkin[3]-0.71875*coeff[0]*fSkin[3]-0.4330127018922193*coeff[1]*fEdge[3]-0.21875*coeff[0]*fEdge[3]+0.4330127018922193*fSkin[1]*coeff[3]-0.4330127018922193*fEdge[1]*coeff[3]-0.375*fSkin[0]*coeff[3]-0.375*fEdge[0]*coeff[3]-0.375*coeff[1]*fSkin[2]+0.7036456405748562*coeff[0]*fSkin[2]-0.375*coeff[1]*fEdge[2]-0.27063293868263694*coeff[0]*fEdge[2]-0.71875*fSkin[1]*coeff[2]-0.21875*fEdge[1]*coeff[2]+0.7036456405748562*fSkin[0]*coeff[2]-0.27063293868263694*fEdge[0]*coeff[2]; 

  boundSurf_incr[1] = -(0.8660254037844386*coeff[3]*fSkin[3])-0.5*coeff[2]*fSkin[3]-0.75*fSkin[2]*coeff[3]-0.4330127018922193*coeff[2]*fSkin[2]-0.8660254037844386*coeff[1]*fSkin[1]-0.5*coeff[0]*fSkin[1]-0.75*fSkin[0]*coeff[1]-0.4330127018922193*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(0.8660254037844386*coeff[1]*fSkin[3])-0.5*coeff[0]*fSkin[3]-0.8660254037844386*fSkin[1]*coeff[3]-0.75*fSkin[0]*coeff[3]-0.75*coeff[1]*fSkin[2]-0.4330127018922193*coeff[0]*fSkin[2]-0.5*fSkin[1]*coeff[2]-0.4330127018922193*fSkin[0]*coeff[2]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 

  return 0.;
}

