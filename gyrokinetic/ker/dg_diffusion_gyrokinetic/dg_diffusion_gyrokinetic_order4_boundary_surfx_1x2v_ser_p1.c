#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double vol_incr[12] = {0.0}; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.6237976320958223*coeff[0]*qSkin[1]+1.6237976320958223*coeff[0]*qEdge[1]+0.9375*coeff[0]*qSkin[0]-0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*qSkin[1]+2.0625*coeff[0]*qEdge[1]+1.6237976320958223*coeff[0]*qSkin[0]-1.6237976320958223*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 1.6237976320958223*coeff[0]*qSkin[4]+1.6237976320958223*coeff[0]*qEdge[4]+0.9375*coeff[0]*qSkin[2]-0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 1.6237976320958223*coeff[0]*qSkin[5]+1.6237976320958223*coeff[0]*qEdge[5]+0.9375*coeff[0]*qSkin[3]-0.9375*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[0]*qSkin[4]+2.0625*coeff[0]*qEdge[4]+1.6237976320958223*coeff[0]*qSkin[2]-1.6237976320958223*coeff[0]*qEdge[2]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*qSkin[5]+2.0625*coeff[0]*qEdge[5]+1.6237976320958223*coeff[0]*qSkin[3]-1.6237976320958223*coeff[0]*qEdge[3]; 
  edgeSurf_incr[6] = 1.6237976320958223*coeff[0]*qSkin[7]+1.6237976320958223*coeff[0]*qEdge[7]+0.9375*coeff[0]*qSkin[6]-0.9375*coeff[0]*qEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[0]*qSkin[7]+2.0625*coeff[0]*qEdge[7]+1.6237976320958223*coeff[0]*qSkin[6]-1.6237976320958223*coeff[0]*qEdge[6]; 
  edgeSurf_incr[8] = 1.6237976320958225*coeff[0]*qSkin[9]+1.6237976320958225*coeff[0]*qEdge[9]+0.9375*coeff[0]*qSkin[8]-0.9375*coeff[0]*qEdge[8]; 
  edgeSurf_incr[9] = 3.5625*coeff[0]*qSkin[9]+2.0625*coeff[0]*qEdge[9]+1.6237976320958225*coeff[0]*qSkin[8]-1.6237976320958225*coeff[0]*qEdge[8]; 
  edgeSurf_incr[10] = 1.6237976320958225*coeff[0]*qSkin[11]+1.6237976320958225*coeff[0]*qEdge[11]+0.9375*coeff[0]*qSkin[10]-0.9375*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*qSkin[11]+2.0625*coeff[0]*qEdge[11]+1.6237976320958225*coeff[0]*qSkin[10]-1.6237976320958225*coeff[0]*qEdge[10]; 

  boundSurf_incr[1] = 1.5*coeff[0]*qSkin[1]; 
  boundSurf_incr[4] = 1.5*coeff[0]*qSkin[4]; 
  boundSurf_incr[5] = 1.5*coeff[0]*qSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[0]*qSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[0]*qSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[0]*qSkin[11]; 

  } else { 

  edgeSurf_incr[0] = -(1.6237976320958223*coeff[0]*qSkin[1])-1.6237976320958223*coeff[0]*qEdge[1]+0.9375*coeff[0]*qSkin[0]-0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*qSkin[1]+2.0625*coeff[0]*qEdge[1]-1.6237976320958223*coeff[0]*qSkin[0]+1.6237976320958223*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = -(1.6237976320958223*coeff[0]*qSkin[4])-1.6237976320958223*coeff[0]*qEdge[4]+0.9375*coeff[0]*qSkin[2]-0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = -(1.6237976320958223*coeff[0]*qSkin[5])-1.6237976320958223*coeff[0]*qEdge[5]+0.9375*coeff[0]*qSkin[3]-0.9375*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[0]*qSkin[4]+2.0625*coeff[0]*qEdge[4]-1.6237976320958223*coeff[0]*qSkin[2]+1.6237976320958223*coeff[0]*qEdge[2]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*qSkin[5]+2.0625*coeff[0]*qEdge[5]-1.6237976320958223*coeff[0]*qSkin[3]+1.6237976320958223*coeff[0]*qEdge[3]; 
  edgeSurf_incr[6] = -(1.6237976320958223*coeff[0]*qSkin[7])-1.6237976320958223*coeff[0]*qEdge[7]+0.9375*coeff[0]*qSkin[6]-0.9375*coeff[0]*qEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[0]*qSkin[7]+2.0625*coeff[0]*qEdge[7]-1.6237976320958223*coeff[0]*qSkin[6]+1.6237976320958223*coeff[0]*qEdge[6]; 
  edgeSurf_incr[8] = -(1.6237976320958225*coeff[0]*qSkin[9])-1.6237976320958225*coeff[0]*qEdge[9]+0.9375*coeff[0]*qSkin[8]-0.9375*coeff[0]*qEdge[8]; 
  edgeSurf_incr[9] = 3.5625*coeff[0]*qSkin[9]+2.0625*coeff[0]*qEdge[9]-1.6237976320958225*coeff[0]*qSkin[8]+1.6237976320958225*coeff[0]*qEdge[8]; 
  edgeSurf_incr[10] = -(1.6237976320958225*coeff[0]*qSkin[11])-1.6237976320958225*coeff[0]*qEdge[11]+0.9375*coeff[0]*qSkin[10]-0.9375*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*qSkin[11]+2.0625*coeff[0]*qEdge[11]-1.6237976320958225*coeff[0]*qSkin[10]+1.6237976320958225*coeff[0]*qEdge[10]; 

  boundSurf_incr[1] = 1.5*coeff[0]*qSkin[1]; 
  boundSurf_incr[4] = 1.5*coeff[0]*qSkin[4]; 
  boundSurf_incr[5] = 1.5*coeff[0]*qSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[0]*qSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[0]*qSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[0]*qSkin[11]; 

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

  return 0.;
}

