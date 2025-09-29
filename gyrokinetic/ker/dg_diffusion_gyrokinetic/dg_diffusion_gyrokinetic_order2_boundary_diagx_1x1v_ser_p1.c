#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double boundSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*qSkin[0]-1.0*coeff[0]*qSkin[1]; 
  boundSurf_incr[3] = 0.8660254037844386*coeff[0]*qSkin[2]-1.0*coeff[0]*qSkin[3]; 
  boundSurf_incr[5] = 0.8660254037844387*coeff[0]*qSkin[4]-1.0*coeff[0]*qSkin[5]; 

  } else { 

  boundSurf_incr[1] = -(1.0*coeff[0]*qSkin[1])-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = -(1.0*coeff[0]*qSkin[3])-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[5] = -(1.0*coeff[0]*qSkin[5])-0.8660254037844387*coeff[0]*qSkin[4]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fGhost[6];
  fGhost[0] = 0.7071067811865476*(jacobgeo_inv[1]*qGhost[1]+jacobgeo_inv[0]*qGhost[0]); 
  fGhost[1] = 0.7071067811865476*(jacobgeo_inv[0]*qGhost[1]+qGhost[0]*jacobgeo_inv[1]); 
  fGhost[2] = 0.7071067811865476*(jacobgeo_inv[1]*qGhost[3]+jacobgeo_inv[0]*qGhost[2]); 
  fGhost[3] = 0.7071067811865476*(jacobgeo_inv[0]*qGhost[3]+jacobgeo_inv[1]*qGhost[2]); 
  fGhost[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qGhost[5]+21.213203435596427*jacobgeo_inv[0]*qGhost[4]); 
  fGhost[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_inv[0]*qGhost[5]+21.21320343559643*jacobgeo_inv[1]*qGhost[4]); 

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.7071067811865476*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_inv[1]*qSkin[3]+jacobgeo_inv[0]*qSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_inv[0]*qSkin[3]+jacobgeo_inv[1]*qSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qSkin[5]+21.213203435596427*jacobgeo_inv[0]*qSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_inv[0]*qSkin[5]+21.21320343559643*jacobgeo_inv[1]*qSkin[4]); 

  double boundSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 1.224744871391589*coeff[1]*fSkin[1]-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = 1.224744871391589*coeff[1]*fSkin[3]-0.7071067811865475*coeff[0]*fSkin[3]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 1.224744871391589*coeff[1]*fSkin[5]-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[4]+0.6123724356957944*coeff[0]*fSkin[4]; 

  } else { 

  boundSurf_incr[1] = -(1.224744871391589*coeff[1]*fSkin[1])-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(1.224744871391589*coeff[1]*fSkin[3])-0.7071067811865475*coeff[0]*fSkin[3]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = -(1.224744871391589*coeff[1]*fSkin[5])-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[4]-0.6123724356957944*coeff[0]*fSkin[4]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 

  return 0.;
}

