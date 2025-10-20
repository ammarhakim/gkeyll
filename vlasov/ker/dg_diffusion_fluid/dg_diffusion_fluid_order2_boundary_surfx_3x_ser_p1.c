#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[0]*fSkin[1])-0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]-1.4072912811497125*coeff[0]*fSkin[0]+0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.5412658773652741*coeff[0]*fSkin[4])-0.5412658773652741*coeff[0]*fEdge[4]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.5412658773652741*coeff[0]*fSkin[5])-0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(1.4375*coeff[0]*fSkin[4])-0.4375*coeff[0]*fEdge[4]-1.4072912811497125*coeff[0]*fSkin[2]+0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]-1.4072912811497125*coeff[0]*fSkin[3]+0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.5412658773652741*coeff[0]*fSkin[7])-0.5412658773652741*coeff[0]*fEdge[7]-0.5625*coeff[0]*fSkin[6]+0.5625*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[0]*fSkin[7])-0.4375*coeff[0]*fEdge[7]-1.4072912811497125*coeff[0]*fSkin[6]+0.5412658773652739*coeff[0]*fEdge[6]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*fSkin[0]-1.0*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 0.8660254037844386*coeff[0]*fSkin[2]-1.0*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[0]*fSkin[3]-1.0*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[0]*fSkin[6]-1.0*coeff[0]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[0]*fSkin[1]+0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]+1.4072912811497125*coeff[0]*fSkin[0]-0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[0]*fSkin[4]+0.5412658773652741*coeff[0]*fEdge[4]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[0]*fSkin[5]+0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(1.4375*coeff[0]*fSkin[4])-0.4375*coeff[0]*fEdge[4]+1.4072912811497125*coeff[0]*fSkin[2]-0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]+1.4072912811497125*coeff[0]*fSkin[3]-0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 0.5412658773652741*coeff[0]*fSkin[7]+0.5412658773652741*coeff[0]*fEdge[7]-0.5625*coeff[0]*fSkin[6]+0.5625*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[0]*fSkin[7])-0.4375*coeff[0]*fEdge[7]+1.4072912811497125*coeff[0]*fSkin[6]-0.5412658773652739*coeff[0]*fEdge[6]; 

  boundSurf_incr[1] = -(1.0*coeff[0]*fSkin[1])-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = -(1.0*coeff[0]*fSkin[4])-0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = -(1.0*coeff[0]*fSkin[5])-0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(1.0*coeff[0]*fSkin[7])-0.8660254037844386*coeff[0]*fSkin[6]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_3x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[8] = {0.0}; 
  vol_incr[1] = 1.060660171779821*fSkin[6]*coeff[7]+1.060660171779821*fSkin[3]*coeff[5]+1.060660171779821*fSkin[2]*coeff[4]+1.060660171779821*fSkin[0]*coeff[1]; 
  vol_incr[4] = 1.060660171779821*fSkin[3]*coeff[7]+1.060660171779821*coeff[5]*fSkin[6]+1.060660171779821*fSkin[0]*coeff[4]+1.060660171779821*coeff[1]*fSkin[2]; 
  vol_incr[5] = 1.060660171779821*fSkin[2]*coeff[7]+1.060660171779821*coeff[4]*fSkin[6]+1.060660171779821*fSkin[0]*coeff[5]+1.060660171779821*coeff[1]*fSkin[3]; 
  vol_incr[7] = 1.060660171779821*fSkin[0]*coeff[7]+1.060660171779821*coeff[1]*fSkin[6]+1.060660171779821*fSkin[2]*coeff[5]+1.060660171779821*fSkin[3]*coeff[4]; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.19136638615493565*coeff[6]*fSkin[7])-0.19136638615493565*coeff[6]*fEdge[7]-0.19887378220871635*coeff[6]*fSkin[6]+0.19887378220871635*coeff[6]*fEdge[6]-0.19136638615493565*coeff[3]*fSkin[5]-0.19136638615493565*coeff[3]*fEdge[5]-0.19136638615493565*coeff[2]*fSkin[4]-0.19136638615493565*coeff[2]*fEdge[4]-0.19887378220871635*coeff[3]*fSkin[3]+0.19887378220871635*coeff[3]*fEdge[3]-0.19887378220871635*coeff[2]*fSkin[2]+0.19887378220871635*coeff[2]*fEdge[2]-0.19136638615493565*coeff[0]*fSkin[1]-0.19136638615493565*coeff[0]*fEdge[1]-0.19887378220871635*coeff[0]*fSkin[0]+0.19887378220871635*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(0.3061862178478971*coeff[7]*fSkin[7])-0.5082329989778307*coeff[6]*fSkin[7]+0.3061862178478971*coeff[7]*fEdge[7]-0.15467960838455713*coeff[6]*fEdge[7]-0.26516504294495524*fSkin[6]*coeff[7]-0.26516504294495524*fEdge[6]*coeff[7]-0.49755260400283263*coeff[6]*fSkin[6]+0.1913663861549355*coeff[6]*fEdge[6]-0.3061862178478971*coeff[5]*fSkin[5]-0.5082329989778307*coeff[3]*fSkin[5]+0.3061862178478971*coeff[5]*fEdge[5]-0.15467960838455713*coeff[3]*fEdge[5]-0.26516504294495524*fSkin[3]*coeff[5]-0.26516504294495524*fEdge[3]*coeff[5]-0.3061862178478971*coeff[4]*fSkin[4]-0.5082329989778307*coeff[2]*fSkin[4]+0.3061862178478971*coeff[4]*fEdge[4]-0.15467960838455713*coeff[2]*fEdge[4]-0.26516504294495524*fSkin[2]*coeff[4]-0.26516504294495524*fEdge[2]*coeff[4]-0.49755260400283263*coeff[3]*fSkin[3]+0.1913663861549355*coeff[3]*fEdge[3]-0.49755260400283263*coeff[2]*fSkin[2]+0.1913663861549355*coeff[2]*fEdge[2]-0.3061862178478971*coeff[1]*fSkin[1]-0.5082329989778307*coeff[0]*fSkin[1]+0.3061862178478971*coeff[1]*fEdge[1]-0.15467960838455713*coeff[0]*fEdge[1]-0.26516504294495524*fSkin[0]*coeff[1]-0.26516504294495524*fEdge[0]*coeff[1]-0.49755260400283263*coeff[0]*fSkin[0]+0.1913663861549355*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.19136638615493565*coeff[3]*fSkin[7])-0.19136638615493565*coeff[3]*fEdge[7]-0.19887378220871635*coeff[3]*fSkin[6]+0.19887378220871635*coeff[3]*fEdge[6]-0.19136638615493565*fSkin[5]*coeff[6]-0.19136638615493565*fEdge[5]*coeff[6]-0.19887378220871635*fSkin[3]*coeff[6]+0.19887378220871635*fEdge[3]*coeff[6]-0.19136638615493565*coeff[0]*fSkin[4]-0.19136638615493565*coeff[0]*fEdge[4]-0.19887378220871635*coeff[0]*fSkin[2]+0.19887378220871635*coeff[0]*fEdge[2]-0.19136638615493565*fSkin[1]*coeff[2]-0.19136638615493565*fEdge[1]*coeff[2]-0.19887378220871635*fSkin[0]*coeff[2]+0.19887378220871635*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = -(0.19136638615493565*coeff[2]*fSkin[7])-0.19136638615493565*coeff[2]*fEdge[7]-0.19887378220871635*coeff[2]*fSkin[6]+0.19887378220871635*coeff[2]*fEdge[6]-0.19136638615493565*fSkin[4]*coeff[6]-0.19136638615493565*fEdge[4]*coeff[6]-0.19887378220871635*fSkin[2]*coeff[6]+0.19887378220871635*fEdge[2]*coeff[6]-0.19136638615493565*coeff[0]*fSkin[5]-0.19136638615493565*coeff[0]*fEdge[5]-0.19887378220871635*coeff[0]*fSkin[3]+0.19887378220871635*coeff[0]*fEdge[3]-0.19136638615493565*fSkin[1]*coeff[3]-0.19136638615493565*fEdge[1]*coeff[3]-0.19887378220871635*fSkin[0]*coeff[3]+0.19887378220871635*fEdge[0]*coeff[3]; 
  edgeSurf_incr[4] = -(0.3061862178478971*coeff[5]*fSkin[7])-0.5082329989778307*coeff[3]*fSkin[7]+0.3061862178478971*coeff[5]*fEdge[7]-0.15467960838455713*coeff[3]*fEdge[7]-0.3061862178478971*fSkin[5]*coeff[7]+0.3061862178478971*fEdge[5]*coeff[7]-0.26516504294495524*fSkin[3]*coeff[7]-0.26516504294495524*fEdge[3]*coeff[7]-0.26516504294495524*coeff[5]*fSkin[6]-0.49755260400283263*coeff[3]*fSkin[6]-0.26516504294495524*coeff[5]*fEdge[6]+0.1913663861549355*coeff[3]*fEdge[6]-0.5082329989778307*fSkin[5]*coeff[6]-0.15467960838455713*fEdge[5]*coeff[6]-0.49755260400283263*fSkin[3]*coeff[6]+0.1913663861549355*fEdge[3]*coeff[6]-0.3061862178478971*coeff[1]*fSkin[4]-0.5082329989778307*coeff[0]*fSkin[4]+0.3061862178478971*coeff[1]*fEdge[4]-0.15467960838455713*coeff[0]*fEdge[4]-0.3061862178478971*fSkin[1]*coeff[4]+0.3061862178478971*fEdge[1]*coeff[4]-0.26516504294495524*fSkin[0]*coeff[4]-0.26516504294495524*fEdge[0]*coeff[4]-0.26516504294495524*coeff[1]*fSkin[2]-0.49755260400283263*coeff[0]*fSkin[2]-0.26516504294495524*coeff[1]*fEdge[2]+0.1913663861549355*coeff[0]*fEdge[2]-0.5082329989778307*fSkin[1]*coeff[2]-0.15467960838455713*fEdge[1]*coeff[2]-0.49755260400283263*fSkin[0]*coeff[2]+0.1913663861549355*fEdge[0]*coeff[2]; 
  edgeSurf_incr[5] = -(0.3061862178478971*coeff[4]*fSkin[7])-0.5082329989778307*coeff[2]*fSkin[7]+0.3061862178478971*coeff[4]*fEdge[7]-0.15467960838455713*coeff[2]*fEdge[7]-0.3061862178478971*fSkin[4]*coeff[7]+0.3061862178478971*fEdge[4]*coeff[7]-0.26516504294495524*fSkin[2]*coeff[7]-0.26516504294495524*fEdge[2]*coeff[7]-0.26516504294495524*coeff[4]*fSkin[6]-0.49755260400283263*coeff[2]*fSkin[6]-0.26516504294495524*coeff[4]*fEdge[6]+0.1913663861549355*coeff[2]*fEdge[6]-0.5082329989778307*fSkin[4]*coeff[6]-0.15467960838455713*fEdge[4]*coeff[6]-0.49755260400283263*fSkin[2]*coeff[6]+0.1913663861549355*fEdge[2]*coeff[6]-0.3061862178478971*coeff[1]*fSkin[5]-0.5082329989778307*coeff[0]*fSkin[5]+0.3061862178478971*coeff[1]*fEdge[5]-0.15467960838455713*coeff[0]*fEdge[5]-0.3061862178478971*fSkin[1]*coeff[5]+0.3061862178478971*fEdge[1]*coeff[5]-0.26516504294495524*fSkin[0]*coeff[5]-0.26516504294495524*fEdge[0]*coeff[5]-0.26516504294495524*coeff[1]*fSkin[3]-0.49755260400283263*coeff[0]*fSkin[3]-0.26516504294495524*coeff[1]*fEdge[3]+0.1913663861549355*coeff[0]*fEdge[3]-0.5082329989778307*fSkin[1]*coeff[3]-0.15467960838455713*fEdge[1]*coeff[3]-0.49755260400283263*fSkin[0]*coeff[3]+0.1913663861549355*fEdge[0]*coeff[3]; 
  edgeSurf_incr[6] = -(0.19136638615493565*coeff[0]*fSkin[7])-0.19136638615493565*coeff[0]*fEdge[7]-0.19887378220871635*coeff[0]*fSkin[6]+0.19887378220871635*coeff[0]*fEdge[6]-0.19136638615493565*fSkin[1]*coeff[6]-0.19136638615493565*fEdge[1]*coeff[6]-0.19887378220871635*fSkin[0]*coeff[6]+0.19887378220871635*fEdge[0]*coeff[6]-0.19136638615493565*coeff[2]*fSkin[5]-0.19136638615493565*coeff[2]*fEdge[5]-0.19136638615493565*coeff[3]*fSkin[4]-0.19136638615493565*coeff[3]*fEdge[4]-0.19887378220871635*coeff[2]*fSkin[3]+0.19887378220871635*coeff[2]*fEdge[3]-0.19887378220871635*fSkin[2]*coeff[3]+0.19887378220871635*fEdge[2]*coeff[3]; 
  edgeSurf_incr[7] = -(0.3061862178478971*coeff[1]*fSkin[7])-0.5082329989778307*coeff[0]*fSkin[7]+0.3061862178478971*coeff[1]*fEdge[7]-0.15467960838455713*coeff[0]*fEdge[7]-0.3061862178478971*fSkin[1]*coeff[7]+0.3061862178478971*fEdge[1]*coeff[7]-0.26516504294495524*fSkin[0]*coeff[7]-0.26516504294495524*fEdge[0]*coeff[7]-0.26516504294495524*coeff[1]*fSkin[6]-0.49755260400283263*coeff[0]*fSkin[6]-0.26516504294495524*coeff[1]*fEdge[6]+0.1913663861549355*coeff[0]*fEdge[6]-0.5082329989778307*fSkin[1]*coeff[6]-0.15467960838455713*fEdge[1]*coeff[6]-0.49755260400283263*fSkin[0]*coeff[6]+0.1913663861549355*fEdge[0]*coeff[6]-0.3061862178478971*coeff[4]*fSkin[5]-0.5082329989778307*coeff[2]*fSkin[5]+0.3061862178478971*coeff[4]*fEdge[5]-0.15467960838455713*coeff[2]*fEdge[5]-0.3061862178478971*fSkin[4]*coeff[5]+0.3061862178478971*fEdge[4]*coeff[5]-0.26516504294495524*fSkin[2]*coeff[5]-0.26516504294495524*fEdge[2]*coeff[5]-0.5082329989778307*coeff[3]*fSkin[4]-0.15467960838455713*coeff[3]*fEdge[4]-0.26516504294495524*fSkin[3]*coeff[4]-0.26516504294495524*fEdge[3]*coeff[4]-0.49755260400283263*coeff[2]*fSkin[3]+0.1913663861549355*coeff[2]*fEdge[3]-0.49755260400283263*fSkin[2]*coeff[3]+0.1913663861549355*fEdge[2]*coeff[3]; 

  boundSurf_incr[1] = 0.6123724356957944*coeff[7]*fSkin[7]-0.3535533905932737*coeff[6]*fSkin[7]-0.5303300858899105*fSkin[6]*coeff[7]+0.3061862178478971*coeff[6]*fSkin[6]+0.6123724356957944*coeff[5]*fSkin[5]-0.3535533905932737*coeff[3]*fSkin[5]-0.5303300858899105*fSkin[3]*coeff[5]+0.6123724356957944*coeff[4]*fSkin[4]-0.3535533905932737*coeff[2]*fSkin[4]-0.5303300858899105*fSkin[2]*coeff[4]+0.3061862178478971*coeff[3]*fSkin[3]+0.3061862178478971*coeff[2]*fSkin[2]+0.6123724356957944*coeff[1]*fSkin[1]-0.3535533905932737*coeff[0]*fSkin[1]-0.5303300858899105*fSkin[0]*coeff[1]+0.3061862178478971*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 0.6123724356957944*coeff[5]*fSkin[7]-0.3535533905932737*coeff[3]*fSkin[7]+0.6123724356957944*fSkin[5]*coeff[7]-0.5303300858899105*fSkin[3]*coeff[7]-0.5303300858899105*coeff[5]*fSkin[6]+0.3061862178478971*coeff[3]*fSkin[6]-0.3535533905932737*fSkin[5]*coeff[6]+0.3061862178478971*fSkin[3]*coeff[6]+0.6123724356957944*coeff[1]*fSkin[4]-0.3535533905932737*coeff[0]*fSkin[4]+0.6123724356957944*fSkin[1]*coeff[4]-0.5303300858899105*fSkin[0]*coeff[4]-0.5303300858899105*coeff[1]*fSkin[2]+0.3061862178478971*coeff[0]*fSkin[2]-0.3535533905932737*fSkin[1]*coeff[2]+0.3061862178478971*fSkin[0]*coeff[2]; 
  boundSurf_incr[5] = 0.6123724356957944*coeff[4]*fSkin[7]-0.3535533905932737*coeff[2]*fSkin[7]+0.6123724356957944*fSkin[4]*coeff[7]-0.5303300858899105*fSkin[2]*coeff[7]-0.5303300858899105*coeff[4]*fSkin[6]+0.3061862178478971*coeff[2]*fSkin[6]-0.3535533905932737*fSkin[4]*coeff[6]+0.3061862178478971*fSkin[2]*coeff[6]+0.6123724356957944*coeff[1]*fSkin[5]-0.3535533905932737*coeff[0]*fSkin[5]+0.6123724356957944*fSkin[1]*coeff[5]-0.5303300858899105*fSkin[0]*coeff[5]-0.5303300858899105*coeff[1]*fSkin[3]+0.3061862178478971*coeff[0]*fSkin[3]-0.3535533905932737*fSkin[1]*coeff[3]+0.3061862178478971*fSkin[0]*coeff[3]; 
  boundSurf_incr[7] = 0.6123724356957944*coeff[1]*fSkin[7]-0.3535533905932737*coeff[0]*fSkin[7]+0.6123724356957944*fSkin[1]*coeff[7]-0.5303300858899105*fSkin[0]*coeff[7]-0.5303300858899105*coeff[1]*fSkin[6]+0.3061862178478971*coeff[0]*fSkin[6]-0.3535533905932737*fSkin[1]*coeff[6]+0.3061862178478971*fSkin[0]*coeff[6]+0.6123724356957944*coeff[4]*fSkin[5]-0.3535533905932737*coeff[2]*fSkin[5]+0.6123724356957944*fSkin[4]*coeff[5]-0.5303300858899105*fSkin[2]*coeff[5]-0.3535533905932737*coeff[3]*fSkin[4]-0.5303300858899105*fSkin[3]*coeff[4]+0.3061862178478971*coeff[2]*fSkin[3]+0.3061862178478971*fSkin[2]*coeff[3]; 

  } else { 

  edgeSurf_incr[0] = 0.19136638615493565*coeff[6]*fSkin[7]+0.19136638615493565*coeff[6]*fEdge[7]-0.19887378220871635*coeff[6]*fSkin[6]+0.19887378220871635*coeff[6]*fEdge[6]+0.19136638615493565*coeff[3]*fSkin[5]+0.19136638615493565*coeff[3]*fEdge[5]+0.19136638615493565*coeff[2]*fSkin[4]+0.19136638615493565*coeff[2]*fEdge[4]-0.19887378220871635*coeff[3]*fSkin[3]+0.19887378220871635*coeff[3]*fEdge[3]-0.19887378220871635*coeff[2]*fSkin[2]+0.19887378220871635*coeff[2]*fEdge[2]+0.19136638615493565*coeff[0]*fSkin[1]+0.19136638615493565*coeff[0]*fEdge[1]-0.19887378220871635*coeff[0]*fSkin[0]+0.19887378220871635*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.3061862178478971*coeff[7]*fSkin[7]-0.5082329989778307*coeff[6]*fSkin[7]-0.3061862178478971*coeff[7]*fEdge[7]-0.15467960838455713*coeff[6]*fEdge[7]-0.26516504294495524*fSkin[6]*coeff[7]-0.26516504294495524*fEdge[6]*coeff[7]+0.49755260400283263*coeff[6]*fSkin[6]-0.1913663861549355*coeff[6]*fEdge[6]+0.3061862178478971*coeff[5]*fSkin[5]-0.5082329989778307*coeff[3]*fSkin[5]-0.3061862178478971*coeff[5]*fEdge[5]-0.15467960838455713*coeff[3]*fEdge[5]-0.26516504294495524*fSkin[3]*coeff[5]-0.26516504294495524*fEdge[3]*coeff[5]+0.3061862178478971*coeff[4]*fSkin[4]-0.5082329989778307*coeff[2]*fSkin[4]-0.3061862178478971*coeff[4]*fEdge[4]-0.15467960838455713*coeff[2]*fEdge[4]-0.26516504294495524*fSkin[2]*coeff[4]-0.26516504294495524*fEdge[2]*coeff[4]+0.49755260400283263*coeff[3]*fSkin[3]-0.1913663861549355*coeff[3]*fEdge[3]+0.49755260400283263*coeff[2]*fSkin[2]-0.1913663861549355*coeff[2]*fEdge[2]+0.3061862178478971*coeff[1]*fSkin[1]-0.5082329989778307*coeff[0]*fSkin[1]-0.3061862178478971*coeff[1]*fEdge[1]-0.15467960838455713*coeff[0]*fEdge[1]-0.26516504294495524*fSkin[0]*coeff[1]-0.26516504294495524*fEdge[0]*coeff[1]+0.49755260400283263*coeff[0]*fSkin[0]-0.1913663861549355*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.19136638615493565*coeff[3]*fSkin[7]+0.19136638615493565*coeff[3]*fEdge[7]-0.19887378220871635*coeff[3]*fSkin[6]+0.19887378220871635*coeff[3]*fEdge[6]+0.19136638615493565*fSkin[5]*coeff[6]+0.19136638615493565*fEdge[5]*coeff[6]-0.19887378220871635*fSkin[3]*coeff[6]+0.19887378220871635*fEdge[3]*coeff[6]+0.19136638615493565*coeff[0]*fSkin[4]+0.19136638615493565*coeff[0]*fEdge[4]-0.19887378220871635*coeff[0]*fSkin[2]+0.19887378220871635*coeff[0]*fEdge[2]+0.19136638615493565*fSkin[1]*coeff[2]+0.19136638615493565*fEdge[1]*coeff[2]-0.19887378220871635*fSkin[0]*coeff[2]+0.19887378220871635*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = 0.19136638615493565*coeff[2]*fSkin[7]+0.19136638615493565*coeff[2]*fEdge[7]-0.19887378220871635*coeff[2]*fSkin[6]+0.19887378220871635*coeff[2]*fEdge[6]+0.19136638615493565*fSkin[4]*coeff[6]+0.19136638615493565*fEdge[4]*coeff[6]-0.19887378220871635*fSkin[2]*coeff[6]+0.19887378220871635*fEdge[2]*coeff[6]+0.19136638615493565*coeff[0]*fSkin[5]+0.19136638615493565*coeff[0]*fEdge[5]-0.19887378220871635*coeff[0]*fSkin[3]+0.19887378220871635*coeff[0]*fEdge[3]+0.19136638615493565*fSkin[1]*coeff[3]+0.19136638615493565*fEdge[1]*coeff[3]-0.19887378220871635*fSkin[0]*coeff[3]+0.19887378220871635*fEdge[0]*coeff[3]; 
  edgeSurf_incr[4] = 0.3061862178478971*coeff[5]*fSkin[7]-0.5082329989778307*coeff[3]*fSkin[7]-0.3061862178478971*coeff[5]*fEdge[7]-0.15467960838455713*coeff[3]*fEdge[7]+0.3061862178478971*fSkin[5]*coeff[7]-0.3061862178478971*fEdge[5]*coeff[7]-0.26516504294495524*fSkin[3]*coeff[7]-0.26516504294495524*fEdge[3]*coeff[7]-0.26516504294495524*coeff[5]*fSkin[6]+0.49755260400283263*coeff[3]*fSkin[6]-0.26516504294495524*coeff[5]*fEdge[6]-0.1913663861549355*coeff[3]*fEdge[6]-0.5082329989778307*fSkin[5]*coeff[6]-0.15467960838455713*fEdge[5]*coeff[6]+0.49755260400283263*fSkin[3]*coeff[6]-0.1913663861549355*fEdge[3]*coeff[6]+0.3061862178478971*coeff[1]*fSkin[4]-0.5082329989778307*coeff[0]*fSkin[4]-0.3061862178478971*coeff[1]*fEdge[4]-0.15467960838455713*coeff[0]*fEdge[4]+0.3061862178478971*fSkin[1]*coeff[4]-0.3061862178478971*fEdge[1]*coeff[4]-0.26516504294495524*fSkin[0]*coeff[4]-0.26516504294495524*fEdge[0]*coeff[4]-0.26516504294495524*coeff[1]*fSkin[2]+0.49755260400283263*coeff[0]*fSkin[2]-0.26516504294495524*coeff[1]*fEdge[2]-0.1913663861549355*coeff[0]*fEdge[2]-0.5082329989778307*fSkin[1]*coeff[2]-0.15467960838455713*fEdge[1]*coeff[2]+0.49755260400283263*fSkin[0]*coeff[2]-0.1913663861549355*fEdge[0]*coeff[2]; 
  edgeSurf_incr[5] = 0.3061862178478971*coeff[4]*fSkin[7]-0.5082329989778307*coeff[2]*fSkin[7]-0.3061862178478971*coeff[4]*fEdge[7]-0.15467960838455713*coeff[2]*fEdge[7]+0.3061862178478971*fSkin[4]*coeff[7]-0.3061862178478971*fEdge[4]*coeff[7]-0.26516504294495524*fSkin[2]*coeff[7]-0.26516504294495524*fEdge[2]*coeff[7]-0.26516504294495524*coeff[4]*fSkin[6]+0.49755260400283263*coeff[2]*fSkin[6]-0.26516504294495524*coeff[4]*fEdge[6]-0.1913663861549355*coeff[2]*fEdge[6]-0.5082329989778307*fSkin[4]*coeff[6]-0.15467960838455713*fEdge[4]*coeff[6]+0.49755260400283263*fSkin[2]*coeff[6]-0.1913663861549355*fEdge[2]*coeff[6]+0.3061862178478971*coeff[1]*fSkin[5]-0.5082329989778307*coeff[0]*fSkin[5]-0.3061862178478971*coeff[1]*fEdge[5]-0.15467960838455713*coeff[0]*fEdge[5]+0.3061862178478971*fSkin[1]*coeff[5]-0.3061862178478971*fEdge[1]*coeff[5]-0.26516504294495524*fSkin[0]*coeff[5]-0.26516504294495524*fEdge[0]*coeff[5]-0.26516504294495524*coeff[1]*fSkin[3]+0.49755260400283263*coeff[0]*fSkin[3]-0.26516504294495524*coeff[1]*fEdge[3]-0.1913663861549355*coeff[0]*fEdge[3]-0.5082329989778307*fSkin[1]*coeff[3]-0.15467960838455713*fEdge[1]*coeff[3]+0.49755260400283263*fSkin[0]*coeff[3]-0.1913663861549355*fEdge[0]*coeff[3]; 
  edgeSurf_incr[6] = 0.19136638615493565*coeff[0]*fSkin[7]+0.19136638615493565*coeff[0]*fEdge[7]-0.19887378220871635*coeff[0]*fSkin[6]+0.19887378220871635*coeff[0]*fEdge[6]+0.19136638615493565*fSkin[1]*coeff[6]+0.19136638615493565*fEdge[1]*coeff[6]-0.19887378220871635*fSkin[0]*coeff[6]+0.19887378220871635*fEdge[0]*coeff[6]+0.19136638615493565*coeff[2]*fSkin[5]+0.19136638615493565*coeff[2]*fEdge[5]+0.19136638615493565*coeff[3]*fSkin[4]+0.19136638615493565*coeff[3]*fEdge[4]-0.19887378220871635*coeff[2]*fSkin[3]+0.19887378220871635*coeff[2]*fEdge[3]-0.19887378220871635*fSkin[2]*coeff[3]+0.19887378220871635*fEdge[2]*coeff[3]; 
  edgeSurf_incr[7] = 0.3061862178478971*coeff[1]*fSkin[7]-0.5082329989778307*coeff[0]*fSkin[7]-0.3061862178478971*coeff[1]*fEdge[7]-0.15467960838455713*coeff[0]*fEdge[7]+0.3061862178478971*fSkin[1]*coeff[7]-0.3061862178478971*fEdge[1]*coeff[7]-0.26516504294495524*fSkin[0]*coeff[7]-0.26516504294495524*fEdge[0]*coeff[7]-0.26516504294495524*coeff[1]*fSkin[6]+0.49755260400283263*coeff[0]*fSkin[6]-0.26516504294495524*coeff[1]*fEdge[6]-0.1913663861549355*coeff[0]*fEdge[6]-0.5082329989778307*fSkin[1]*coeff[6]-0.15467960838455713*fEdge[1]*coeff[6]+0.49755260400283263*fSkin[0]*coeff[6]-0.1913663861549355*fEdge[0]*coeff[6]+0.3061862178478971*coeff[4]*fSkin[5]-0.5082329989778307*coeff[2]*fSkin[5]-0.3061862178478971*coeff[4]*fEdge[5]-0.15467960838455713*coeff[2]*fEdge[5]+0.3061862178478971*fSkin[4]*coeff[5]-0.3061862178478971*fEdge[4]*coeff[5]-0.26516504294495524*fSkin[2]*coeff[5]-0.26516504294495524*fEdge[2]*coeff[5]-0.5082329989778307*coeff[3]*fSkin[4]-0.15467960838455713*coeff[3]*fEdge[4]-0.26516504294495524*fSkin[3]*coeff[4]-0.26516504294495524*fEdge[3]*coeff[4]+0.49755260400283263*coeff[2]*fSkin[3]-0.1913663861549355*coeff[2]*fEdge[3]+0.49755260400283263*fSkin[2]*coeff[3]-0.1913663861549355*fEdge[2]*coeff[3]; 

  boundSurf_incr[1] = -(0.6123724356957944*coeff[7]*fSkin[7])-0.3535533905932737*coeff[6]*fSkin[7]-0.5303300858899105*fSkin[6]*coeff[7]-0.3061862178478971*coeff[6]*fSkin[6]-0.6123724356957944*coeff[5]*fSkin[5]-0.3535533905932737*coeff[3]*fSkin[5]-0.5303300858899105*fSkin[3]*coeff[5]-0.6123724356957944*coeff[4]*fSkin[4]-0.3535533905932737*coeff[2]*fSkin[4]-0.5303300858899105*fSkin[2]*coeff[4]-0.3061862178478971*coeff[3]*fSkin[3]-0.3061862178478971*coeff[2]*fSkin[2]-0.6123724356957944*coeff[1]*fSkin[1]-0.3535533905932737*coeff[0]*fSkin[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.3061862178478971*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = -(0.6123724356957944*coeff[5]*fSkin[7])-0.3535533905932737*coeff[3]*fSkin[7]-0.6123724356957944*fSkin[5]*coeff[7]-0.5303300858899105*fSkin[3]*coeff[7]-0.5303300858899105*coeff[5]*fSkin[6]-0.3061862178478971*coeff[3]*fSkin[6]-0.3535533905932737*fSkin[5]*coeff[6]-0.3061862178478971*fSkin[3]*coeff[6]-0.6123724356957944*coeff[1]*fSkin[4]-0.3535533905932737*coeff[0]*fSkin[4]-0.6123724356957944*fSkin[1]*coeff[4]-0.5303300858899105*fSkin[0]*coeff[4]-0.5303300858899105*coeff[1]*fSkin[2]-0.3061862178478971*coeff[0]*fSkin[2]-0.3535533905932737*fSkin[1]*coeff[2]-0.3061862178478971*fSkin[0]*coeff[2]; 
  boundSurf_incr[5] = -(0.6123724356957944*coeff[4]*fSkin[7])-0.3535533905932737*coeff[2]*fSkin[7]-0.6123724356957944*fSkin[4]*coeff[7]-0.5303300858899105*fSkin[2]*coeff[7]-0.5303300858899105*coeff[4]*fSkin[6]-0.3061862178478971*coeff[2]*fSkin[6]-0.3535533905932737*fSkin[4]*coeff[6]-0.3061862178478971*fSkin[2]*coeff[6]-0.6123724356957944*coeff[1]*fSkin[5]-0.3535533905932737*coeff[0]*fSkin[5]-0.6123724356957944*fSkin[1]*coeff[5]-0.5303300858899105*fSkin[0]*coeff[5]-0.5303300858899105*coeff[1]*fSkin[3]-0.3061862178478971*coeff[0]*fSkin[3]-0.3535533905932737*fSkin[1]*coeff[3]-0.3061862178478971*fSkin[0]*coeff[3]; 
  boundSurf_incr[7] = -(0.6123724356957944*coeff[1]*fSkin[7])-0.3535533905932737*coeff[0]*fSkin[7]-0.6123724356957944*fSkin[1]*coeff[7]-0.5303300858899105*fSkin[0]*coeff[7]-0.5303300858899105*coeff[1]*fSkin[6]-0.3061862178478971*coeff[0]*fSkin[6]-0.3535533905932737*fSkin[1]*coeff[6]-0.3061862178478971*fSkin[0]*coeff[6]-0.6123724356957944*coeff[4]*fSkin[5]-0.3535533905932737*coeff[2]*fSkin[5]-0.6123724356957944*fSkin[4]*coeff[5]-0.5303300858899105*fSkin[2]*coeff[5]-0.3535533905932737*coeff[3]*fSkin[4]-0.5303300858899105*fSkin[3]*coeff[4]-0.3061862178478971*coeff[2]*fSkin[3]-0.3061862178478971*fSkin[2]*coeff[3]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  return 0.;
}

