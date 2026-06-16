#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order4_boundary_diagx_lower_bound_local_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],4);

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double boundSurf_incr[6] = {0.0}; 


  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 

  return 0.;
}

#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order4_boundary_diagx_lower_bound_recovery_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],4);

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double fGhost[6];
  fGhost[0] = 0.7071067811865476*(JfGhost[1]*jacobgeo_invGhost[1]+JfGhost[0]*jacobgeo_invGhost[0]); 
  fGhost[1] = 0.7071067811865476*(JfGhost[0]*jacobgeo_invGhost[1]+jacobgeo_invGhost[0]*JfGhost[1]); 
  fGhost[2] = 0.7071067811865476*(jacobgeo_invGhost[1]*JfGhost[3]+jacobgeo_invGhost[0]*JfGhost[2]); 
  fGhost[3] = 0.7071067811865476*(jacobgeo_invGhost[0]*JfGhost[3]+jacobgeo_invGhost[1]*JfGhost[2]); 
  fGhost[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invGhost[1]*JfGhost[5]+21.213203435596427*jacobgeo_invGhost[0]*JfGhost[4]); 
  fGhost[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invGhost[0]*JfGhost[5]+21.21320343559643*jacobgeo_invGhost[1]*JfGhost[4]); 

  double boundSurf_incr[6] = {0.0}; 

  boundSurf_incr[0] = 0.9943689110435817*fSkin[1]*nuSkin[1]+0.9943689110435817*fGhost[1]*nuSkin[1]-0.5740991584648069*fSkin[0]*nuSkin[1]+0.5740991584648069*fGhost[0]*nuSkin[1]-0.9943689110435817*fSkin[1]*nuGhost[1]-0.9943689110435817*fGhost[1]*nuGhost[1]+0.5740991584648069*fSkin[0]*nuGhost[1]-0.5740991584648069*fGhost[0]*nuGhost[1]-0.5740991584648069*nuSkin[0]*fSkin[1]-0.5740991584648069*nuGhost[0]*fSkin[1]-0.5740991584648069*nuSkin[0]*fGhost[1]-0.5740991584648069*nuGhost[0]*fGhost[1]+0.33145630368119394*fSkin[0]*nuSkin[0]-0.33145630368119394*fGhost[0]*nuSkin[0]+0.33145630368119394*fSkin[0]*nuGhost[0]-0.33145630368119394*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[1] = -(1.7222974753944202*fSkin[1]*nuSkin[1])-1.7222974753944202*fGhost[1]*nuSkin[1]+0.9943689110435817*fSkin[0]*nuSkin[1]-0.9943689110435817*fGhost[0]*nuSkin[1]+1.7222974753944202*fSkin[1]*nuGhost[1]+1.7222974753944202*fGhost[1]*nuGhost[1]-0.9943689110435817*fSkin[0]*nuGhost[1]+0.9943689110435817*fGhost[0]*nuGhost[1]+0.9943689110435817*nuSkin[0]*fSkin[1]+0.9943689110435817*nuGhost[0]*fSkin[1]+0.9943689110435817*nuSkin[0]*fGhost[1]+0.9943689110435817*nuGhost[0]*fGhost[1]-0.5740991584648069*fSkin[0]*nuSkin[0]+0.5740991584648069*fGhost[0]*nuSkin[0]-0.5740991584648069*fSkin[0]*nuGhost[0]+0.5740991584648069*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[2] = 0.9943689110435817*nuSkin[1]*fSkin[3]-0.9943689110435817*nuGhost[1]*fSkin[3]-0.5740991584648069*nuSkin[0]*fSkin[3]-0.5740991584648069*nuGhost[0]*fSkin[3]+0.9943689110435817*nuSkin[1]*fGhost[3]-0.9943689110435817*nuGhost[1]*fGhost[3]-0.5740991584648069*nuSkin[0]*fGhost[3]-0.5740991584648069*nuGhost[0]*fGhost[3]-0.5740991584648069*nuSkin[1]*fSkin[2]+0.5740991584648069*nuGhost[1]*fSkin[2]+0.33145630368119394*nuSkin[0]*fSkin[2]+0.33145630368119394*nuGhost[0]*fSkin[2]+0.5740991584648069*nuSkin[1]*fGhost[2]-0.5740991584648069*nuGhost[1]*fGhost[2]-0.33145630368119394*nuSkin[0]*fGhost[2]-0.33145630368119394*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[3] = -(1.7222974753944202*nuSkin[1]*fSkin[3])+1.7222974753944202*nuGhost[1]*fSkin[3]+0.9943689110435817*nuSkin[0]*fSkin[3]+0.9943689110435817*nuGhost[0]*fSkin[3]-1.7222974753944202*nuSkin[1]*fGhost[3]+1.7222974753944202*nuGhost[1]*fGhost[3]+0.9943689110435817*nuSkin[0]*fGhost[3]+0.9943689110435817*nuGhost[0]*fGhost[3]+0.9943689110435817*nuSkin[1]*fSkin[2]-0.9943689110435817*nuGhost[1]*fSkin[2]-0.5740991584648069*nuSkin[0]*fSkin[2]-0.5740991584648069*nuGhost[0]*fSkin[2]-0.9943689110435817*nuSkin[1]*fGhost[2]+0.9943689110435817*nuGhost[1]*fGhost[2]+0.5740991584648069*nuSkin[0]*fGhost[2]+0.5740991584648069*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[4] = 0.9943689110435818*nuSkin[1]*fSkin[5]-0.9943689110435818*nuGhost[1]*fSkin[5]-0.574099158464807*nuSkin[0]*fSkin[5]-0.574099158464807*nuGhost[0]*fSkin[5]+0.9943689110435818*nuSkin[1]*fGhost[5]-0.9943689110435818*nuGhost[1]*fGhost[5]-0.574099158464807*nuSkin[0]*fGhost[5]-0.574099158464807*nuGhost[0]*fGhost[5]-0.5740991584648069*nuSkin[1]*fSkin[4]+0.5740991584648069*nuGhost[1]*fSkin[4]+0.33145630368119394*nuSkin[0]*fSkin[4]+0.33145630368119394*nuGhost[0]*fSkin[4]+0.5740991584648069*nuSkin[1]*fGhost[4]-0.5740991584648069*nuGhost[1]*fGhost[4]-0.33145630368119394*nuSkin[0]*fGhost[4]-0.33145630368119394*nuGhost[0]*fGhost[4]; 
  boundSurf_incr[5] = -(1.7222974753944202*nuSkin[1]*fSkin[5])+1.7222974753944202*nuGhost[1]*fSkin[5]+0.9943689110435817*nuSkin[0]*fSkin[5]+0.9943689110435817*nuGhost[0]*fSkin[5]-1.7222974753944202*nuSkin[1]*fGhost[5]+1.7222974753944202*nuGhost[1]*fGhost[5]+0.9943689110435817*nuSkin[0]*fGhost[5]+0.9943689110435817*nuGhost[0]*fGhost[5]+0.9943689110435818*nuSkin[1]*fSkin[4]-0.9943689110435818*nuGhost[1]*fSkin[4]-0.574099158464807*nuSkin[0]*fSkin[4]-0.574099158464807*nuGhost[0]*fSkin[4]-0.9943689110435818*nuSkin[1]*fGhost[4]+0.9943689110435818*nuGhost[1]*fGhost[4]+0.574099158464807*nuSkin[0]*fGhost[4]+0.574099158464807*nuGhost[0]*fGhost[4]; 

  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 

  return 0.;
}

#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order4_boundary_diagx_upper_bound_local_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],4);

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double boundSurf_incr[6] = {0.0}; 


  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 

  return 0.;
}

#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order4_boundary_diagx_upper_bound_recovery_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],4);

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double fGhost[6];
  fGhost[0] = 0.7071067811865476*(JfGhost[1]*jacobgeo_invGhost[1]+JfGhost[0]*jacobgeo_invGhost[0]); 
  fGhost[1] = 0.7071067811865476*(JfGhost[0]*jacobgeo_invGhost[1]+jacobgeo_invGhost[0]*JfGhost[1]); 
  fGhost[2] = 0.7071067811865476*(jacobgeo_invGhost[1]*JfGhost[3]+jacobgeo_invGhost[0]*JfGhost[2]); 
  fGhost[3] = 0.7071067811865476*(jacobgeo_invGhost[0]*JfGhost[3]+jacobgeo_invGhost[1]*JfGhost[2]); 
  fGhost[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invGhost[1]*JfGhost[5]+21.213203435596427*jacobgeo_invGhost[0]*JfGhost[4]); 
  fGhost[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invGhost[0]*JfGhost[5]+21.21320343559643*jacobgeo_invGhost[1]*JfGhost[4]); 

  double boundSurf_incr[6] = {0.0}; 

  boundSurf_incr[0] = 0.9943689110435817*fSkin[1]*nuSkin[1]+0.9943689110435817*fGhost[1]*nuSkin[1]+0.5740991584648069*fSkin[0]*nuSkin[1]-0.5740991584648069*fGhost[0]*nuSkin[1]-0.9943689110435817*fSkin[1]*nuGhost[1]-0.9943689110435817*fGhost[1]*nuGhost[1]-0.5740991584648069*fSkin[0]*nuGhost[1]+0.5740991584648069*fGhost[0]*nuGhost[1]+0.5740991584648069*nuSkin[0]*fSkin[1]+0.5740991584648069*nuGhost[0]*fSkin[1]+0.5740991584648069*nuSkin[0]*fGhost[1]+0.5740991584648069*nuGhost[0]*fGhost[1]+0.33145630368119394*fSkin[0]*nuSkin[0]-0.33145630368119394*fGhost[0]*nuSkin[0]+0.33145630368119394*fSkin[0]*nuGhost[0]-0.33145630368119394*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[1] = 1.7222974753944202*fSkin[1]*nuSkin[1]+1.7222974753944202*fGhost[1]*nuSkin[1]+0.9943689110435817*fSkin[0]*nuSkin[1]-0.9943689110435817*fGhost[0]*nuSkin[1]-1.7222974753944202*fSkin[1]*nuGhost[1]-1.7222974753944202*fGhost[1]*nuGhost[1]-0.9943689110435817*fSkin[0]*nuGhost[1]+0.9943689110435817*fGhost[0]*nuGhost[1]+0.9943689110435817*nuSkin[0]*fSkin[1]+0.9943689110435817*nuGhost[0]*fSkin[1]+0.9943689110435817*nuSkin[0]*fGhost[1]+0.9943689110435817*nuGhost[0]*fGhost[1]+0.5740991584648069*fSkin[0]*nuSkin[0]-0.5740991584648069*fGhost[0]*nuSkin[0]+0.5740991584648069*fSkin[0]*nuGhost[0]-0.5740991584648069*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[2] = 0.9943689110435817*nuSkin[1]*fSkin[3]-0.9943689110435817*nuGhost[1]*fSkin[3]+0.5740991584648069*nuSkin[0]*fSkin[3]+0.5740991584648069*nuGhost[0]*fSkin[3]+0.9943689110435817*nuSkin[1]*fGhost[3]-0.9943689110435817*nuGhost[1]*fGhost[3]+0.5740991584648069*nuSkin[0]*fGhost[3]+0.5740991584648069*nuGhost[0]*fGhost[3]+0.5740991584648069*nuSkin[1]*fSkin[2]-0.5740991584648069*nuGhost[1]*fSkin[2]+0.33145630368119394*nuSkin[0]*fSkin[2]+0.33145630368119394*nuGhost[0]*fSkin[2]-0.5740991584648069*nuSkin[1]*fGhost[2]+0.5740991584648069*nuGhost[1]*fGhost[2]-0.33145630368119394*nuSkin[0]*fGhost[2]-0.33145630368119394*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[3] = 1.7222974753944202*nuSkin[1]*fSkin[3]-1.7222974753944202*nuGhost[1]*fSkin[3]+0.9943689110435817*nuSkin[0]*fSkin[3]+0.9943689110435817*nuGhost[0]*fSkin[3]+1.7222974753944202*nuSkin[1]*fGhost[3]-1.7222974753944202*nuGhost[1]*fGhost[3]+0.9943689110435817*nuSkin[0]*fGhost[3]+0.9943689110435817*nuGhost[0]*fGhost[3]+0.9943689110435817*nuSkin[1]*fSkin[2]-0.9943689110435817*nuGhost[1]*fSkin[2]+0.5740991584648069*nuSkin[0]*fSkin[2]+0.5740991584648069*nuGhost[0]*fSkin[2]-0.9943689110435817*nuSkin[1]*fGhost[2]+0.9943689110435817*nuGhost[1]*fGhost[2]-0.5740991584648069*nuSkin[0]*fGhost[2]-0.5740991584648069*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[4] = 0.9943689110435818*nuSkin[1]*fSkin[5]-0.9943689110435818*nuGhost[1]*fSkin[5]+0.574099158464807*nuSkin[0]*fSkin[5]+0.574099158464807*nuGhost[0]*fSkin[5]+0.9943689110435818*nuSkin[1]*fGhost[5]-0.9943689110435818*nuGhost[1]*fGhost[5]+0.574099158464807*nuSkin[0]*fGhost[5]+0.574099158464807*nuGhost[0]*fGhost[5]+0.5740991584648069*nuSkin[1]*fSkin[4]-0.5740991584648069*nuGhost[1]*fSkin[4]+0.33145630368119394*nuSkin[0]*fSkin[4]+0.33145630368119394*nuGhost[0]*fSkin[4]-0.5740991584648069*nuSkin[1]*fGhost[4]+0.5740991584648069*nuGhost[1]*fGhost[4]-0.33145630368119394*nuSkin[0]*fGhost[4]-0.33145630368119394*nuGhost[0]*fGhost[4]; 
  boundSurf_incr[5] = 1.7222974753944202*nuSkin[1]*fSkin[5]-1.7222974753944202*nuGhost[1]*fSkin[5]+0.9943689110435817*nuSkin[0]*fSkin[5]+0.9943689110435817*nuGhost[0]*fSkin[5]+1.7222974753944202*nuSkin[1]*fGhost[5]-1.7222974753944202*nuGhost[1]*fGhost[5]+0.9943689110435817*nuSkin[0]*fGhost[5]+0.9943689110435817*nuGhost[0]*fGhost[5]+0.9943689110435818*nuSkin[1]*fSkin[4]-0.9943689110435818*nuGhost[1]*fSkin[4]+0.574099158464807*nuSkin[0]*fSkin[4]+0.574099158464807*nuGhost[0]*fSkin[4]-0.9943689110435818*nuSkin[1]*fGhost[4]+0.9943689110435818*nuGhost[1]*fGhost[4]-0.574099158464807*nuSkin[0]*fGhost[4]-0.574099158464807*nuGhost[0]*fGhost[4]; 

  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 

  return 0.;
}

