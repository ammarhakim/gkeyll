#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*qSkin[0]-1.0*coeff[0]*qSkin[1]; 
  boundSurf_incr[4] = 0.8660254037844386*coeff[0]*qSkin[2]-1.0*coeff[0]*qSkin[4]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[0]*qSkin[3]-1.0*coeff[0]*qSkin[5]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[0]*qSkin[6]-1.0*coeff[0]*qSkin[7]; 
  boundSurf_incr[9] = 0.8660254037844387*coeff[0]*qSkin[8]-1.0*coeff[0]*qSkin[9]; 
  boundSurf_incr[11] = 0.8660254037844387*coeff[0]*qSkin[10]-1.0*coeff[0]*qSkin[11]; 

  } else { 

  boundSurf_incr[1] = -(1.0*coeff[0]*qSkin[1])-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[4] = -(1.0*coeff[0]*qSkin[4])-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[5] = -(1.0*coeff[0]*qSkin[5])-0.8660254037844386*coeff[0]*qSkin[3]; 
  boundSurf_incr[7] = -(1.0*coeff[0]*qSkin[7])-0.8660254037844386*coeff[0]*qSkin[6]; 
  boundSurf_incr[9] = -(1.0*coeff[0]*qSkin[9])-0.8660254037844387*coeff[0]*qSkin[8]; 
  boundSurf_incr[11] = -(1.0*coeff[0]*qSkin[11])-0.8660254037844387*coeff[0]*qSkin[10]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 
  out[8] += boundSurf_incr[8]*rdx2Sq; 
  out[9] += boundSurf_incr[9]*rdx2Sq; 
  out[10] += boundSurf_incr[10]*rdx2Sq; 
  out[11] += boundSurf_incr[11]*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fGhost[12];
  fGhost[0] = 0.7071067811865475*(jacobgeo_inv[1]*qGhost[1]+jacobgeo_inv[0]*qGhost[0]); 
  fGhost[1] = 0.7071067811865475*(jacobgeo_inv[0]*qGhost[1]+qGhost[0]*jacobgeo_inv[1]); 
  fGhost[2] = 0.7071067811865475*(jacobgeo_inv[1]*qGhost[4]+jacobgeo_inv[0]*qGhost[2]); 
  fGhost[3] = 0.7071067811865475*(jacobgeo_inv[1]*qGhost[5]+jacobgeo_inv[0]*qGhost[3]); 
  fGhost[4] = 0.7071067811865475*(jacobgeo_inv[0]*qGhost[4]+jacobgeo_inv[1]*qGhost[2]); 
  fGhost[5] = 0.7071067811865475*(jacobgeo_inv[0]*qGhost[5]+jacobgeo_inv[1]*qGhost[3]); 
  fGhost[6] = 0.7071067811865475*(jacobgeo_inv[1]*qGhost[7]+jacobgeo_inv[0]*qGhost[6]); 
  fGhost[7] = 0.7071067811865475*(jacobgeo_inv[0]*qGhost[7]+jacobgeo_inv[1]*qGhost[6]); 
  fGhost[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qGhost[9]+15.0*jacobgeo_inv[0]*qGhost[8]); 
  fGhost[9] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qGhost[9]+15.000000000000002*jacobgeo_inv[1]*qGhost[8]); 
  fGhost[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qGhost[11]+15.0*jacobgeo_inv[0]*qGhost[10]); 
  fGhost[11] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qGhost[11]+15.000000000000002*jacobgeo_inv[1]*qGhost[10]); 

  double fSkin[12];
  fSkin[0] = 0.7071067811865475*(jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.7071067811865475*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1]); 
  fSkin[2] = 0.7071067811865475*(jacobgeo_inv[1]*qSkin[4]+jacobgeo_inv[0]*qSkin[2]); 
  fSkin[3] = 0.7071067811865475*(jacobgeo_inv[1]*qSkin[5]+jacobgeo_inv[0]*qSkin[3]); 
  fSkin[4] = 0.7071067811865475*(jacobgeo_inv[0]*qSkin[4]+jacobgeo_inv[1]*qSkin[2]); 
  fSkin[5] = 0.7071067811865475*(jacobgeo_inv[0]*qSkin[5]+jacobgeo_inv[1]*qSkin[3]); 
  fSkin[6] = 0.7071067811865475*(jacobgeo_inv[1]*qSkin[7]+jacobgeo_inv[0]*qSkin[6]); 
  fSkin[7] = 0.7071067811865475*(jacobgeo_inv[0]*qSkin[7]+jacobgeo_inv[1]*qSkin[6]); 
  fSkin[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qSkin[9]+15.0*jacobgeo_inv[0]*qSkin[8]); 
  fSkin[9] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qSkin[9]+15.000000000000002*jacobgeo_inv[1]*qSkin[8]); 
  fSkin[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qSkin[11]+15.0*jacobgeo_inv[0]*qSkin[10]); 
  fSkin[11] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qSkin[11]+15.000000000000002*jacobgeo_inv[1]*qSkin[10]); 

  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 1.224744871391589*coeff[1]*fSkin[1]-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 1.224744871391589*coeff[1]*fSkin[4]-0.7071067811865475*coeff[0]*fSkin[4]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 1.224744871391589*coeff[1]*fSkin[5]-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[3]+0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 1.224744871391589*coeff[1]*fSkin[7]-0.7071067811865475*coeff[0]*fSkin[7]-1.060660171779821*coeff[1]*fSkin[6]+0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[9] = 1.224744871391589*coeff[1]*fSkin[9]-0.7071067811865475*coeff[0]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[8]+0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 1.224744871391589*coeff[1]*fSkin[11]-0.7071067811865475*coeff[0]*fSkin[11]-1.060660171779821*coeff[1]*fSkin[10]+0.6123724356957944*coeff[0]*fSkin[10]; 

  } else { 

  boundSurf_incr[1] = -(1.224744871391589*coeff[1]*fSkin[1])-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = -(1.224744871391589*coeff[1]*fSkin[4])-0.7071067811865475*coeff[0]*fSkin[4]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = -(1.224744871391589*coeff[1]*fSkin[5])-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[3]-0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(1.224744871391589*coeff[1]*fSkin[7])-0.7071067811865475*coeff[0]*fSkin[7]-1.060660171779821*coeff[1]*fSkin[6]-0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[9] = -(1.224744871391589*coeff[1]*fSkin[9])-0.7071067811865475*coeff[0]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[8]-0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = -(1.224744871391589*coeff[1]*fSkin[11])-0.7071067811865475*coeff[0]*fSkin[11]-1.060660171779821*coeff[1]*fSkin[10]-0.6123724356957944*coeff[0]*fSkin[10]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 
  out[8] += boundSurf_incr[8]*rdx2Sq; 
  out[9] += boundSurf_incr[9]*rdx2Sq; 
  out[10] += boundSurf_incr[10]*rdx2Sq; 
  out[11] += boundSurf_incr[11]*rdx2Sq; 

  return 0.;
}

