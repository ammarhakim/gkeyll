#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double vol_incr[12] = {0.0}; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[0]*qSkin[1])-0.5412658773652741*coeff[0]*qEdge[1]-0.5625*coeff[0]*qSkin[0]+0.5625*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*qSkin[1])-0.4375*coeff[0]*qEdge[1]-1.4072912811497125*coeff[0]*qSkin[0]+0.5412658773652739*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = -(0.5412658773652741*coeff[0]*qSkin[4])-0.5412658773652741*coeff[0]*qEdge[4]-0.5625*coeff[0]*qSkin[2]+0.5625*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = -(0.5412658773652741*coeff[0]*qSkin[5])-0.5412658773652741*coeff[0]*qEdge[5]-0.5625*coeff[0]*qSkin[3]+0.5625*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = -(1.4375*coeff[0]*qSkin[4])-0.4375*coeff[0]*qEdge[4]-1.4072912811497125*coeff[0]*qSkin[2]+0.5412658773652739*coeff[0]*qEdge[2]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*qSkin[5])-0.4375*coeff[0]*qEdge[5]-1.4072912811497125*coeff[0]*qSkin[3]+0.5412658773652739*coeff[0]*qEdge[3]; 
  edgeSurf_incr[6] = -(0.5412658773652741*coeff[0]*qSkin[7])-0.5412658773652741*coeff[0]*qEdge[7]-0.5625*coeff[0]*qSkin[6]+0.5625*coeff[0]*qEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[0]*qSkin[7])-0.4375*coeff[0]*qEdge[7]-1.4072912811497125*coeff[0]*qSkin[6]+0.5412658773652739*coeff[0]*qEdge[6]; 
  edgeSurf_incr[8] = -(0.5412658773652742*coeff[0]*qSkin[9])-0.5412658773652742*coeff[0]*qEdge[9]-0.5625*coeff[0]*qSkin[8]+0.5625*coeff[0]*qEdge[8]; 
  edgeSurf_incr[9] = -(1.4375*coeff[0]*qSkin[9])-0.4375*coeff[0]*qEdge[9]-1.4072912811497127*coeff[0]*qSkin[8]+0.5412658773652742*coeff[0]*qEdge[8]; 
  edgeSurf_incr[10] = -(0.5412658773652742*coeff[0]*qSkin[11])-0.5412658773652742*coeff[0]*qEdge[11]-0.5625*coeff[0]*qSkin[10]+0.5625*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = -(1.4375*coeff[0]*qSkin[11])-0.4375*coeff[0]*qEdge[11]-1.4072912811497127*coeff[0]*qSkin[10]+0.5412658773652742*coeff[0]*qEdge[10]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*qSkin[0]-1.0*coeff[0]*qSkin[1]; 
  boundSurf_incr[4] = 0.8660254037844386*coeff[0]*qSkin[2]-1.0*coeff[0]*qSkin[4]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[0]*qSkin[3]-1.0*coeff[0]*qSkin[5]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[0]*qSkin[6]-1.0*coeff[0]*qSkin[7]; 
  boundSurf_incr[9] = 0.8660254037844387*coeff[0]*qSkin[8]-1.0*coeff[0]*qSkin[9]; 
  boundSurf_incr[11] = 0.8660254037844387*coeff[0]*qSkin[10]-1.0*coeff[0]*qSkin[11]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[0]*qSkin[1]+0.5412658773652741*coeff[0]*qEdge[1]-0.5625*coeff[0]*qSkin[0]+0.5625*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*qSkin[1])-0.4375*coeff[0]*qEdge[1]+1.4072912811497125*coeff[0]*qSkin[0]-0.5412658773652739*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[0]*qSkin[4]+0.5412658773652741*coeff[0]*qEdge[4]-0.5625*coeff[0]*qSkin[2]+0.5625*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[0]*qSkin[5]+0.5412658773652741*coeff[0]*qEdge[5]-0.5625*coeff[0]*qSkin[3]+0.5625*coeff[0]*qEdge[3]; 
  edgeSurf_incr[4] = -(1.4375*coeff[0]*qSkin[4])-0.4375*coeff[0]*qEdge[4]+1.4072912811497125*coeff[0]*qSkin[2]-0.5412658773652739*coeff[0]*qEdge[2]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*qSkin[5])-0.4375*coeff[0]*qEdge[5]+1.4072912811497125*coeff[0]*qSkin[3]-0.5412658773652739*coeff[0]*qEdge[3]; 
  edgeSurf_incr[6] = 0.5412658773652741*coeff[0]*qSkin[7]+0.5412658773652741*coeff[0]*qEdge[7]-0.5625*coeff[0]*qSkin[6]+0.5625*coeff[0]*qEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[0]*qSkin[7])-0.4375*coeff[0]*qEdge[7]+1.4072912811497125*coeff[0]*qSkin[6]-0.5412658773652739*coeff[0]*qEdge[6]; 
  edgeSurf_incr[8] = 0.5412658773652742*coeff[0]*qSkin[9]+0.5412658773652742*coeff[0]*qEdge[9]-0.5625*coeff[0]*qSkin[8]+0.5625*coeff[0]*qEdge[8]; 
  edgeSurf_incr[9] = -(1.4375*coeff[0]*qSkin[9])-0.4375*coeff[0]*qEdge[9]+1.4072912811497127*coeff[0]*qSkin[8]-0.5412658773652742*coeff[0]*qEdge[8]; 
  edgeSurf_incr[10] = 0.5412658773652742*coeff[0]*qSkin[11]+0.5412658773652742*coeff[0]*qEdge[11]-0.5625*coeff[0]*qSkin[10]+0.5625*coeff[0]*qEdge[10]; 
  edgeSurf_incr[11] = -(1.4375*coeff[0]*qSkin[11])-0.4375*coeff[0]*qEdge[11]+1.4072912811497127*coeff[0]*qSkin[10]-0.5412658773652742*coeff[0]*qEdge[10]; 

  boundSurf_incr[1] = -(1.0*coeff[0]*qSkin[1])-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[4] = -(1.0*coeff[0]*qSkin[4])-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[5] = -(1.0*coeff[0]*qSkin[5])-0.8660254037844386*coeff[0]*qSkin[3]; 
  boundSurf_incr[7] = -(1.0*coeff[0]*qSkin[7])-0.8660254037844386*coeff[0]*qSkin[6]; 
  boundSurf_incr[9] = -(1.0*coeff[0]*qSkin[9])-0.8660254037844387*coeff[0]*qSkin[8]; 
  boundSurf_incr[11] = -(1.0*coeff[0]*qSkin[11])-0.8660254037844387*coeff[0]*qSkin[10]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

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

  double fEdge[12];
  fEdge[0] = 0.7071067811865475*(jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.7071067811865475*(jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1]); 
  fEdge[2] = 0.7071067811865475*(jacobgeo_inv[1]*qEdge[4]+jacobgeo_inv[0]*qEdge[2]); 
  fEdge[3] = 0.7071067811865475*(jacobgeo_inv[1]*qEdge[5]+jacobgeo_inv[0]*qEdge[3]); 
  fEdge[4] = 0.7071067811865475*(jacobgeo_inv[0]*qEdge[4]+jacobgeo_inv[1]*qEdge[2]); 
  fEdge[5] = 0.7071067811865475*(jacobgeo_inv[0]*qEdge[5]+jacobgeo_inv[1]*qEdge[3]); 
  fEdge[6] = 0.7071067811865475*(jacobgeo_inv[1]*qEdge[7]+jacobgeo_inv[0]*qEdge[6]); 
  fEdge[7] = 0.7071067811865475*(jacobgeo_inv[0]*qEdge[7]+jacobgeo_inv[1]*qEdge[6]); 
  fEdge[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qEdge[9]+15.0*jacobgeo_inv[0]*qEdge[8]); 
  fEdge[9] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qEdge[9]+15.000000000000002*jacobgeo_inv[1]*qEdge[8]); 
  fEdge[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qEdge[11]+15.0*jacobgeo_inv[0]*qEdge[10]); 
  fEdge[11] = 0.04714045207910316*(15.0*jacobgeo_inv[0]*qEdge[11]+15.000000000000002*jacobgeo_inv[1]*qEdge[10]); 

  double vol_incr[12] = {0.0}; 
  vol_incr[1] = 2.1213203435596424*fSkin[0]*coeff[1]; 
  vol_incr[4] = 2.1213203435596424*coeff[1]*fSkin[2]; 
  vol_incr[5] = 2.1213203435596424*coeff[1]*fSkin[3]; 
  vol_incr[7] = 2.1213203435596424*coeff[1]*fSkin[6]; 
  vol_incr[9] = 2.1213203435596424*coeff[1]*fSkin[8]; 
  vol_incr[11] = 2.1213203435596424*coeff[1]*fSkin[10]; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.38273277230987135*coeff[0]*fSkin[1])-0.38273277230987135*coeff[0]*fEdge[1]-0.39774756441743275*coeff[0]*fSkin[0]+0.39774756441743275*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(0.6123724356957944*coeff[1]*fSkin[1])-1.0164659979556616*coeff[0]*fSkin[1]+0.6123724356957944*coeff[1]*fEdge[1]-0.3093592167691142*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]-0.9951052080056654*coeff[0]*fSkin[0]+0.3827327723098711*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.38273277230987135*coeff[0]*fSkin[4])-0.38273277230987135*coeff[0]*fEdge[4]-0.39774756441743275*coeff[0]*fSkin[2]+0.39774756441743275*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.38273277230987135*coeff[0]*fSkin[5])-0.38273277230987135*coeff[0]*fEdge[5]-0.39774756441743275*coeff[0]*fSkin[3]+0.39774756441743275*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(0.6123724356957944*coeff[1]*fSkin[4])-1.0164659979556616*coeff[0]*fSkin[4]+0.6123724356957944*coeff[1]*fEdge[4]-0.3093592167691142*coeff[0]*fEdge[4]-0.5303300858899105*coeff[1]*fSkin[2]-0.9951052080056654*coeff[0]*fSkin[2]-0.5303300858899105*coeff[1]*fEdge[2]+0.3827327723098711*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = -(0.6123724356957944*coeff[1]*fSkin[5])-1.0164659979556616*coeff[0]*fSkin[5]+0.6123724356957944*coeff[1]*fEdge[5]-0.3093592167691142*coeff[0]*fEdge[5]-0.5303300858899105*coeff[1]*fSkin[3]-0.9951052080056654*coeff[0]*fSkin[3]-0.5303300858899105*coeff[1]*fEdge[3]+0.3827327723098711*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.38273277230987135*coeff[0]*fSkin[7])-0.38273277230987135*coeff[0]*fEdge[7]-0.39774756441743275*coeff[0]*fSkin[6]+0.39774756441743275*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = -(0.6123724356957944*coeff[1]*fSkin[7])-1.0164659979556616*coeff[0]*fSkin[7]+0.6123724356957944*coeff[1]*fEdge[7]-0.3093592167691142*coeff[0]*fEdge[7]-0.5303300858899105*coeff[1]*fSkin[6]-0.9951052080056654*coeff[0]*fSkin[6]-0.5303300858899105*coeff[1]*fEdge[6]+0.3827327723098711*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = -(0.3827327723098714*coeff[0]*fSkin[9])-0.3827327723098714*coeff[0]*fEdge[9]-0.39774756441743275*coeff[0]*fSkin[8]+0.39774756441743275*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = -(0.6123724356957944*coeff[1]*fSkin[9])-1.0164659979556616*coeff[0]*fSkin[9]+0.6123724356957944*coeff[1]*fEdge[9]-0.3093592167691142*coeff[0]*fEdge[9]-0.5303300858899104*coeff[1]*fSkin[8]-0.9951052080056656*coeff[0]*fSkin[8]-0.5303300858899104*coeff[1]*fEdge[8]+0.38273277230987135*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = -(0.3827327723098714*coeff[0]*fSkin[11])-0.3827327723098714*coeff[0]*fEdge[11]-0.39774756441743275*coeff[0]*fSkin[10]+0.39774756441743275*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.6123724356957944*coeff[1]*fSkin[11])-1.0164659979556616*coeff[0]*fSkin[11]+0.6123724356957944*coeff[1]*fEdge[11]-0.3093592167691142*coeff[0]*fEdge[11]-0.5303300858899104*coeff[1]*fSkin[10]-0.9951052080056656*coeff[0]*fSkin[10]-0.5303300858899104*coeff[1]*fEdge[10]+0.38273277230987135*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = 1.224744871391589*coeff[1]*fSkin[1]-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 1.224744871391589*coeff[1]*fSkin[4]-0.7071067811865475*coeff[0]*fSkin[4]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 1.224744871391589*coeff[1]*fSkin[5]-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[3]+0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 1.224744871391589*coeff[1]*fSkin[7]-0.7071067811865475*coeff[0]*fSkin[7]-1.060660171779821*coeff[1]*fSkin[6]+0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[9] = 1.224744871391589*coeff[1]*fSkin[9]-0.7071067811865475*coeff[0]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[8]+0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 1.224744871391589*coeff[1]*fSkin[11]-0.7071067811865475*coeff[0]*fSkin[11]-1.060660171779821*coeff[1]*fSkin[10]+0.6123724356957944*coeff[0]*fSkin[10]; 

  } else { 

  edgeSurf_incr[0] = 0.38273277230987135*coeff[0]*fSkin[1]+0.38273277230987135*coeff[0]*fEdge[1]-0.39774756441743275*coeff[0]*fSkin[0]+0.39774756441743275*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.6123724356957944*coeff[1]*fSkin[1]-1.0164659979556616*coeff[0]*fSkin[1]-0.6123724356957944*coeff[1]*fEdge[1]-0.3093592167691142*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]+0.9951052080056654*coeff[0]*fSkin[0]-0.3827327723098711*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.38273277230987135*coeff[0]*fSkin[4]+0.38273277230987135*coeff[0]*fEdge[4]-0.39774756441743275*coeff[0]*fSkin[2]+0.39774756441743275*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.38273277230987135*coeff[0]*fSkin[5]+0.38273277230987135*coeff[0]*fEdge[5]-0.39774756441743275*coeff[0]*fSkin[3]+0.39774756441743275*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.6123724356957944*coeff[1]*fSkin[4]-1.0164659979556616*coeff[0]*fSkin[4]-0.6123724356957944*coeff[1]*fEdge[4]-0.3093592167691142*coeff[0]*fEdge[4]-0.5303300858899105*coeff[1]*fSkin[2]+0.9951052080056654*coeff[0]*fSkin[2]-0.5303300858899105*coeff[1]*fEdge[2]-0.3827327723098711*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 0.6123724356957944*coeff[1]*fSkin[5]-1.0164659979556616*coeff[0]*fSkin[5]-0.6123724356957944*coeff[1]*fEdge[5]-0.3093592167691142*coeff[0]*fEdge[5]-0.5303300858899105*coeff[1]*fSkin[3]+0.9951052080056654*coeff[0]*fSkin[3]-0.5303300858899105*coeff[1]*fEdge[3]-0.3827327723098711*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 0.38273277230987135*coeff[0]*fSkin[7]+0.38273277230987135*coeff[0]*fEdge[7]-0.39774756441743275*coeff[0]*fSkin[6]+0.39774756441743275*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 0.6123724356957944*coeff[1]*fSkin[7]-1.0164659979556616*coeff[0]*fSkin[7]-0.6123724356957944*coeff[1]*fEdge[7]-0.3093592167691142*coeff[0]*fEdge[7]-0.5303300858899105*coeff[1]*fSkin[6]+0.9951052080056654*coeff[0]*fSkin[6]-0.5303300858899105*coeff[1]*fEdge[6]-0.3827327723098711*coeff[0]*fEdge[6]; 
  edgeSurf_incr[8] = 0.3827327723098714*coeff[0]*fSkin[9]+0.3827327723098714*coeff[0]*fEdge[9]-0.39774756441743275*coeff[0]*fSkin[8]+0.39774756441743275*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 0.6123724356957944*coeff[1]*fSkin[9]-1.0164659979556616*coeff[0]*fSkin[9]-0.6123724356957944*coeff[1]*fEdge[9]-0.3093592167691142*coeff[0]*fEdge[9]-0.5303300858899104*coeff[1]*fSkin[8]+0.9951052080056656*coeff[0]*fSkin[8]-0.5303300858899104*coeff[1]*fEdge[8]-0.38273277230987135*coeff[0]*fEdge[8]; 
  edgeSurf_incr[10] = 0.3827327723098714*coeff[0]*fSkin[11]+0.3827327723098714*coeff[0]*fEdge[11]-0.39774756441743275*coeff[0]*fSkin[10]+0.39774756441743275*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.6123724356957944*coeff[1]*fSkin[11]-1.0164659979556616*coeff[0]*fSkin[11]-0.6123724356957944*coeff[1]*fEdge[11]-0.3093592167691142*coeff[0]*fEdge[11]-0.5303300858899104*coeff[1]*fSkin[10]+0.9951052080056656*coeff[0]*fSkin[10]-0.5303300858899104*coeff[1]*fEdge[10]-0.38273277230987135*coeff[0]*fEdge[10]; 

  boundSurf_incr[1] = -(1.224744871391589*coeff[1]*fSkin[1])-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = -(1.224744871391589*coeff[1]*fSkin[4])-0.7071067811865475*coeff[0]*fSkin[4]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = -(1.224744871391589*coeff[1]*fSkin[5])-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[3]-0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(1.224744871391589*coeff[1]*fSkin[7])-0.7071067811865475*coeff[0]*fSkin[7]-1.060660171779821*coeff[1]*fSkin[6]-0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[9] = -(1.224744871391589*coeff[1]*fSkin[9])-0.7071067811865475*coeff[0]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[8]-0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = -(1.224744871391589*coeff[1]*fSkin[11])-0.7071067811865475*coeff[0]*fSkin[11]-1.060660171779821*coeff[1]*fSkin[10]-0.6123724356957944*coeff[0]*fSkin[10]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq; 

  return 0.;
}

