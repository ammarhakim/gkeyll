#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_diagx_lower_bound_local_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fSkin[12];
  fSkin[0] = 0.7071067811865475*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865475*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[4]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[5]+jacobgeo_invSkin[0]*JfSkin[3]); 
  fSkin[4] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[4]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[5] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[5]+jacobgeo_invSkin[1]*JfSkin[3]); 
  fSkin[6] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[7]+jacobgeo_invSkin[0]*JfSkin[6]); 
  fSkin[7] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[7]+jacobgeo_invSkin[1]*JfSkin[6]); 
  fSkin[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[9]+15.0*jacobgeo_invSkin[0]*JfSkin[8]); 
  fSkin[9] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[9]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[8]); 
  fSkin[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[11]+15.0*jacobgeo_invSkin[0]*JfSkin[10]); 
  fSkin[11] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[11]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[10]); 

  double boundSurf_incr[12] = {0.0}; 

  boundSurf_incr[0] = 0.6123724356957944*nuSkin[0]*fSkin[1]-1.060660171779821*fSkin[1]*nuSkin[1]; 
  boundSurf_incr[1] = 1.837117307087383*fSkin[1]*nuSkin[1]-1.060660171779821*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[2] = 0.6123724356957944*nuSkin[0]*fSkin[4]-1.060660171779821*nuSkin[1]*fSkin[4]; 
  boundSurf_incr[3] = 0.6123724356957944*nuSkin[0]*fSkin[5]-1.060660171779821*nuSkin[1]*fSkin[5]; 
  boundSurf_incr[4] = 1.837117307087383*nuSkin[1]*fSkin[4]-1.060660171779821*nuSkin[0]*fSkin[4]; 
  boundSurf_incr[5] = 1.837117307087383*nuSkin[1]*fSkin[5]-1.060660171779821*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[6] = 0.6123724356957944*nuSkin[0]*fSkin[7]-1.060660171779821*nuSkin[1]*fSkin[7]; 
  boundSurf_incr[7] = 1.837117307087383*nuSkin[1]*fSkin[7]-1.060660171779821*nuSkin[0]*fSkin[7]; 
  boundSurf_incr[8] = 0.6123724356957944*nuSkin[0]*fSkin[9]-1.060660171779821*nuSkin[1]*fSkin[9]; 
  boundSurf_incr[9] = 1.837117307087383*nuSkin[1]*fSkin[9]-1.060660171779821*nuSkin[0]*fSkin[9]; 
  boundSurf_incr[10] = 0.6123724356957944*nuSkin[0]*fSkin[11]-1.060660171779821*nuSkin[1]*fSkin[11]; 
  boundSurf_incr[11] = 1.837117307087383*nuSkin[1]*fSkin[11]-1.060660171779821*nuSkin[0]*fSkin[11]; 

  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 
  out[6] += boundSurf_incr[6]*rdx2fac; 
  out[7] += boundSurf_incr[7]*rdx2fac; 
  out[8] += boundSurf_incr[8]*rdx2fac; 
  out[9] += boundSurf_incr[9]*rdx2fac; 
  out[10] += boundSurf_incr[10]*rdx2fac; 
  out[11] += boundSurf_incr[11]*rdx2fac; 

  return 0.;
}

#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_diagx_lower_bound_recovery_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fSkin[12];
  fSkin[0] = 0.7071067811865475*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865475*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[4]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[5]+jacobgeo_invSkin[0]*JfSkin[3]); 
  fSkin[4] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[4]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[5] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[5]+jacobgeo_invSkin[1]*JfSkin[3]); 
  fSkin[6] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[7]+jacobgeo_invSkin[0]*JfSkin[6]); 
  fSkin[7] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[7]+jacobgeo_invSkin[1]*JfSkin[6]); 
  fSkin[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[9]+15.0*jacobgeo_invSkin[0]*JfSkin[8]); 
  fSkin[9] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[9]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[8]); 
  fSkin[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[11]+15.0*jacobgeo_invSkin[0]*JfSkin[10]); 
  fSkin[11] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[11]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[10]); 

  double fGhost[12];
  fGhost[0] = 0.7071067811865475*(JfGhost[1]*jacobgeo_invGhost[1]+JfGhost[0]*jacobgeo_invGhost[0]); 
  fGhost[1] = 0.7071067811865475*(JfGhost[0]*jacobgeo_invGhost[1]+jacobgeo_invGhost[0]*JfGhost[1]); 
  fGhost[2] = 0.7071067811865475*(jacobgeo_invGhost[1]*JfGhost[4]+jacobgeo_invGhost[0]*JfGhost[2]); 
  fGhost[3] = 0.7071067811865475*(jacobgeo_invGhost[1]*JfGhost[5]+jacobgeo_invGhost[0]*JfGhost[3]); 
  fGhost[4] = 0.7071067811865475*(jacobgeo_invGhost[0]*JfGhost[4]+jacobgeo_invGhost[1]*JfGhost[2]); 
  fGhost[5] = 0.7071067811865475*(jacobgeo_invGhost[0]*JfGhost[5]+jacobgeo_invGhost[1]*JfGhost[3]); 
  fGhost[6] = 0.7071067811865475*(jacobgeo_invGhost[1]*JfGhost[7]+jacobgeo_invGhost[0]*JfGhost[6]); 
  fGhost[7] = 0.7071067811865475*(jacobgeo_invGhost[0]*JfGhost[7]+jacobgeo_invGhost[1]*JfGhost[6]); 
  fGhost[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invGhost[1]*JfGhost[9]+15.0*jacobgeo_invGhost[0]*JfGhost[8]); 
  fGhost[9] = 0.04714045207910316*(15.0*jacobgeo_invGhost[0]*JfGhost[9]+15.000000000000002*jacobgeo_invGhost[1]*JfGhost[8]); 
  fGhost[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invGhost[1]*JfGhost[11]+15.0*jacobgeo_invGhost[0]*JfGhost[10]); 
  fGhost[11] = 0.04714045207910316*(15.0*jacobgeo_invGhost[0]*JfGhost[11]+15.000000000000002*jacobgeo_invGhost[1]*JfGhost[10]); 

  double boundSurf_incr[12] = {0.0}; 

  boundSurf_incr[0] = 0.33145630368119394*fSkin[1]*nuSkin[1]+0.33145630368119394*fGhost[1]*nuSkin[1]-0.34445949507888407*fSkin[0]*nuSkin[1]+0.34445949507888407*fGhost[0]*nuSkin[1]-0.33145630368119394*fSkin[1]*nuGhost[1]-0.33145630368119394*fGhost[1]*nuGhost[1]+0.34445949507888407*fSkin[0]*nuGhost[1]-0.34445949507888407*fGhost[0]*nuGhost[1]-0.19136638615493565*nuSkin[0]*fSkin[1]-0.19136638615493565*nuGhost[0]*fSkin[1]-0.19136638615493565*nuSkin[0]*fGhost[1]-0.19136638615493565*nuGhost[0]*fGhost[1]+0.19887378220871635*fSkin[0]*nuSkin[0]-0.19887378220871635*fGhost[0]*nuSkin[0]+0.19887378220871635*fSkin[0]*nuGhost[0]-0.19887378220871635*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[1] = -(0.5740991584648069*fSkin[1]*nuSkin[1])-0.5740991584648069*fGhost[1]*nuSkin[1]+0.5966213466261491*fSkin[0]*nuSkin[1]-0.5966213466261491*fGhost[0]*nuSkin[1]+0.5740991584648069*fSkin[1]*nuGhost[1]+0.5740991584648069*fGhost[1]*nuGhost[1]-0.5966213466261491*fSkin[0]*nuGhost[1]+0.5966213466261491*fGhost[0]*nuGhost[1]+0.33145630368119394*nuSkin[0]*fSkin[1]+0.33145630368119394*nuGhost[0]*fSkin[1]+0.33145630368119394*nuSkin[0]*fGhost[1]+0.33145630368119394*nuGhost[0]*fGhost[1]-0.34445949507888407*fSkin[0]*nuSkin[0]+0.34445949507888407*fGhost[0]*nuSkin[0]-0.34445949507888407*fSkin[0]*nuGhost[0]+0.34445949507888407*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[2] = 0.33145630368119394*nuSkin[1]*fSkin[4]-0.33145630368119394*nuGhost[1]*fSkin[4]-0.19136638615493565*nuSkin[0]*fSkin[4]-0.19136638615493565*nuGhost[0]*fSkin[4]+0.33145630368119394*nuSkin[1]*fGhost[4]-0.33145630368119394*nuGhost[1]*fGhost[4]-0.19136638615493565*nuSkin[0]*fGhost[4]-0.19136638615493565*nuGhost[0]*fGhost[4]-0.34445949507888407*nuSkin[1]*fSkin[2]+0.34445949507888407*nuGhost[1]*fSkin[2]+0.19887378220871635*nuSkin[0]*fSkin[2]+0.19887378220871635*nuGhost[0]*fSkin[2]+0.34445949507888407*nuSkin[1]*fGhost[2]-0.34445949507888407*nuGhost[1]*fGhost[2]-0.19887378220871635*nuSkin[0]*fGhost[2]-0.19887378220871635*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[3] = 0.33145630368119394*nuSkin[1]*fSkin[5]-0.33145630368119394*nuGhost[1]*fSkin[5]-0.19136638615493565*nuSkin[0]*fSkin[5]-0.19136638615493565*nuGhost[0]*fSkin[5]+0.33145630368119394*nuSkin[1]*fGhost[5]-0.33145630368119394*nuGhost[1]*fGhost[5]-0.19136638615493565*nuSkin[0]*fGhost[5]-0.19136638615493565*nuGhost[0]*fGhost[5]-0.34445949507888407*nuSkin[1]*fSkin[3]+0.34445949507888407*nuGhost[1]*fSkin[3]+0.19887378220871635*nuSkin[0]*fSkin[3]+0.19887378220871635*nuGhost[0]*fSkin[3]+0.34445949507888407*nuSkin[1]*fGhost[3]-0.34445949507888407*nuGhost[1]*fGhost[3]-0.19887378220871635*nuSkin[0]*fGhost[3]-0.19887378220871635*nuGhost[0]*fGhost[3]; 
  boundSurf_incr[4] = -(0.5740991584648069*nuSkin[1]*fSkin[4])+0.5740991584648069*nuGhost[1]*fSkin[4]+0.33145630368119394*nuSkin[0]*fSkin[4]+0.33145630368119394*nuGhost[0]*fSkin[4]-0.5740991584648069*nuSkin[1]*fGhost[4]+0.5740991584648069*nuGhost[1]*fGhost[4]+0.33145630368119394*nuSkin[0]*fGhost[4]+0.33145630368119394*nuGhost[0]*fGhost[4]+0.5966213466261491*nuSkin[1]*fSkin[2]-0.5966213466261491*nuGhost[1]*fSkin[2]-0.34445949507888407*nuSkin[0]*fSkin[2]-0.34445949507888407*nuGhost[0]*fSkin[2]-0.5966213466261491*nuSkin[1]*fGhost[2]+0.5966213466261491*nuGhost[1]*fGhost[2]+0.34445949507888407*nuSkin[0]*fGhost[2]+0.34445949507888407*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[5] = -(0.5740991584648069*nuSkin[1]*fSkin[5])+0.5740991584648069*nuGhost[1]*fSkin[5]+0.33145630368119394*nuSkin[0]*fSkin[5]+0.33145630368119394*nuGhost[0]*fSkin[5]-0.5740991584648069*nuSkin[1]*fGhost[5]+0.5740991584648069*nuGhost[1]*fGhost[5]+0.33145630368119394*nuSkin[0]*fGhost[5]+0.33145630368119394*nuGhost[0]*fGhost[5]+0.5966213466261491*nuSkin[1]*fSkin[3]-0.5966213466261491*nuGhost[1]*fSkin[3]-0.34445949507888407*nuSkin[0]*fSkin[3]-0.34445949507888407*nuGhost[0]*fSkin[3]-0.5966213466261491*nuSkin[1]*fGhost[3]+0.5966213466261491*nuGhost[1]*fGhost[3]+0.34445949507888407*nuSkin[0]*fGhost[3]+0.34445949507888407*nuGhost[0]*fGhost[3]; 
  boundSurf_incr[6] = 0.33145630368119394*nuSkin[1]*fSkin[7]-0.33145630368119394*nuGhost[1]*fSkin[7]-0.19136638615493565*nuSkin[0]*fSkin[7]-0.19136638615493565*nuGhost[0]*fSkin[7]+0.33145630368119394*nuSkin[1]*fGhost[7]-0.33145630368119394*nuGhost[1]*fGhost[7]-0.19136638615493565*nuSkin[0]*fGhost[7]-0.19136638615493565*nuGhost[0]*fGhost[7]-0.34445949507888407*nuSkin[1]*fSkin[6]+0.34445949507888407*nuGhost[1]*fSkin[6]+0.19887378220871635*nuSkin[0]*fSkin[6]+0.19887378220871635*nuGhost[0]*fSkin[6]+0.34445949507888407*nuSkin[1]*fGhost[6]-0.34445949507888407*nuGhost[1]*fGhost[6]-0.19887378220871635*nuSkin[0]*fGhost[6]-0.19887378220871635*nuGhost[0]*fGhost[6]; 
  boundSurf_incr[7] = -(0.5740991584648069*nuSkin[1]*fSkin[7])+0.5740991584648069*nuGhost[1]*fSkin[7]+0.33145630368119394*nuSkin[0]*fSkin[7]+0.33145630368119394*nuGhost[0]*fSkin[7]-0.5740991584648069*nuSkin[1]*fGhost[7]+0.5740991584648069*nuGhost[1]*fGhost[7]+0.33145630368119394*nuSkin[0]*fGhost[7]+0.33145630368119394*nuGhost[0]*fGhost[7]+0.5966213466261491*nuSkin[1]*fSkin[6]-0.5966213466261491*nuGhost[1]*fSkin[6]-0.34445949507888407*nuSkin[0]*fSkin[6]-0.34445949507888407*nuGhost[0]*fSkin[6]-0.5966213466261491*nuSkin[1]*fGhost[6]+0.5966213466261491*nuGhost[1]*fGhost[6]+0.34445949507888407*nuSkin[0]*fGhost[6]+0.34445949507888407*nuGhost[0]*fGhost[6]; 
  boundSurf_incr[8] = 0.3314563036811939*nuSkin[1]*fSkin[9]-0.3314563036811939*nuGhost[1]*fSkin[9]-0.19136638615493567*nuSkin[0]*fSkin[9]-0.19136638615493567*nuGhost[0]*fSkin[9]+0.3314563036811939*nuSkin[1]*fGhost[9]-0.3314563036811939*nuGhost[1]*fGhost[9]-0.19136638615493567*nuSkin[0]*fGhost[9]-0.19136638615493567*nuGhost[0]*fGhost[9]-0.34445949507888407*nuSkin[1]*fSkin[8]+0.34445949507888407*nuGhost[1]*fSkin[8]+0.19887378220871635*nuSkin[0]*fSkin[8]+0.19887378220871635*nuGhost[0]*fSkin[8]+0.34445949507888407*nuSkin[1]*fGhost[8]-0.34445949507888407*nuGhost[1]*fGhost[8]-0.19887378220871635*nuSkin[0]*fGhost[8]-0.19887378220871635*nuGhost[0]*fGhost[8]; 
  boundSurf_incr[9] = -(0.5740991584648069*nuSkin[1]*fSkin[9])+0.5740991584648069*nuGhost[1]*fSkin[9]+0.33145630368119394*nuSkin[0]*fSkin[9]+0.33145630368119394*nuGhost[0]*fSkin[9]-0.5740991584648069*nuSkin[1]*fGhost[9]+0.5740991584648069*nuGhost[1]*fGhost[9]+0.33145630368119394*nuSkin[0]*fGhost[9]+0.33145630368119394*nuGhost[0]*fGhost[9]+0.5966213466261489*nuSkin[1]*fSkin[8]-0.5966213466261489*nuGhost[1]*fSkin[8]-0.3444594950788842*nuSkin[0]*fSkin[8]-0.3444594950788842*nuGhost[0]*fSkin[8]-0.5966213466261489*nuSkin[1]*fGhost[8]+0.5966213466261489*nuGhost[1]*fGhost[8]+0.3444594950788842*nuSkin[0]*fGhost[8]+0.3444594950788842*nuGhost[0]*fGhost[8]; 
  boundSurf_incr[10] = 0.3314563036811939*nuSkin[1]*fSkin[11]-0.3314563036811939*nuGhost[1]*fSkin[11]-0.19136638615493567*nuSkin[0]*fSkin[11]-0.19136638615493567*nuGhost[0]*fSkin[11]+0.3314563036811939*nuSkin[1]*fGhost[11]-0.3314563036811939*nuGhost[1]*fGhost[11]-0.19136638615493567*nuSkin[0]*fGhost[11]-0.19136638615493567*nuGhost[0]*fGhost[11]-0.34445949507888407*nuSkin[1]*fSkin[10]+0.34445949507888407*nuGhost[1]*fSkin[10]+0.19887378220871635*nuSkin[0]*fSkin[10]+0.19887378220871635*nuGhost[0]*fSkin[10]+0.34445949507888407*nuSkin[1]*fGhost[10]-0.34445949507888407*nuGhost[1]*fGhost[10]-0.19887378220871635*nuSkin[0]*fGhost[10]-0.19887378220871635*nuGhost[0]*fGhost[10]; 
  boundSurf_incr[11] = -(0.5740991584648069*nuSkin[1]*fSkin[11])+0.5740991584648069*nuGhost[1]*fSkin[11]+0.33145630368119394*nuSkin[0]*fSkin[11]+0.33145630368119394*nuGhost[0]*fSkin[11]-0.5740991584648069*nuSkin[1]*fGhost[11]+0.5740991584648069*nuGhost[1]*fGhost[11]+0.33145630368119394*nuSkin[0]*fGhost[11]+0.33145630368119394*nuGhost[0]*fGhost[11]+0.5966213466261489*nuSkin[1]*fSkin[10]-0.5966213466261489*nuGhost[1]*fSkin[10]-0.3444594950788842*nuSkin[0]*fSkin[10]-0.3444594950788842*nuGhost[0]*fSkin[10]-0.5966213466261489*nuSkin[1]*fGhost[10]+0.5966213466261489*nuGhost[1]*fGhost[10]+0.3444594950788842*nuSkin[0]*fGhost[10]+0.3444594950788842*nuGhost[0]*fGhost[10]; 

  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 
  out[6] += boundSurf_incr[6]*rdx2fac; 
  out[7] += boundSurf_incr[7]*rdx2fac; 
  out[8] += boundSurf_incr[8]*rdx2fac; 
  out[9] += boundSurf_incr[9]*rdx2fac; 
  out[10] += boundSurf_incr[10]*rdx2fac; 
  out[11] += boundSurf_incr[11]*rdx2fac; 

  return 0.;
}

#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_diagx_upper_bound_local_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fSkin[12];
  fSkin[0] = 0.7071067811865475*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865475*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[4]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[5]+jacobgeo_invSkin[0]*JfSkin[3]); 
  fSkin[4] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[4]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[5] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[5]+jacobgeo_invSkin[1]*JfSkin[3]); 
  fSkin[6] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[7]+jacobgeo_invSkin[0]*JfSkin[6]); 
  fSkin[7] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[7]+jacobgeo_invSkin[1]*JfSkin[6]); 
  fSkin[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[9]+15.0*jacobgeo_invSkin[0]*JfSkin[8]); 
  fSkin[9] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[9]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[8]); 
  fSkin[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[11]+15.0*jacobgeo_invSkin[0]*JfSkin[10]); 
  fSkin[11] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[11]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[10]); 

  double boundSurf_incr[12] = {0.0}; 

  boundSurf_incr[0] = -(1.060660171779821*fSkin[1]*nuSkin[1])-0.6123724356957944*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[1] = -(1.837117307087383*fSkin[1]*nuSkin[1])-1.060660171779821*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[2] = -(1.060660171779821*nuSkin[1]*fSkin[4])-0.6123724356957944*nuSkin[0]*fSkin[4]; 
  boundSurf_incr[3] = -(1.060660171779821*nuSkin[1]*fSkin[5])-0.6123724356957944*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[4] = -(1.837117307087383*nuSkin[1]*fSkin[4])-1.060660171779821*nuSkin[0]*fSkin[4]; 
  boundSurf_incr[5] = -(1.837117307087383*nuSkin[1]*fSkin[5])-1.060660171779821*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[6] = -(1.060660171779821*nuSkin[1]*fSkin[7])-0.6123724356957944*nuSkin[0]*fSkin[7]; 
  boundSurf_incr[7] = -(1.837117307087383*nuSkin[1]*fSkin[7])-1.060660171779821*nuSkin[0]*fSkin[7]; 
  boundSurf_incr[8] = -(1.060660171779821*nuSkin[1]*fSkin[9])-0.6123724356957944*nuSkin[0]*fSkin[9]; 
  boundSurf_incr[9] = -(1.837117307087383*nuSkin[1]*fSkin[9])-1.060660171779821*nuSkin[0]*fSkin[9]; 
  boundSurf_incr[10] = -(1.060660171779821*nuSkin[1]*fSkin[11])-0.6123724356957944*nuSkin[0]*fSkin[11]; 
  boundSurf_incr[11] = -(1.837117307087383*nuSkin[1]*fSkin[11])-1.060660171779821*nuSkin[0]*fSkin[11]; 

  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 
  out[6] += boundSurf_incr[6]*rdx2fac; 
  out[7] += boundSurf_incr[7]*rdx2fac; 
  out[8] += boundSurf_incr[8]*rdx2fac; 
  out[9] += boundSurf_incr[9]*rdx2fac; 
  out[10] += boundSurf_incr[10]*rdx2fac; 
  out[11] += boundSurf_incr[11]*rdx2fac; 

  return 0.;
}

#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_diagx_upper_bound_recovery_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_invSkin/jacobgeo_invGhost: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fSkin[12];
  fSkin[0] = 0.7071067811865475*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865475*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[4]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[5]+jacobgeo_invSkin[0]*JfSkin[3]); 
  fSkin[4] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[4]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[5] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[5]+jacobgeo_invSkin[1]*JfSkin[3]); 
  fSkin[6] = 0.7071067811865475*(jacobgeo_invSkin[1]*JfSkin[7]+jacobgeo_invSkin[0]*JfSkin[6]); 
  fSkin[7] = 0.7071067811865475*(jacobgeo_invSkin[0]*JfSkin[7]+jacobgeo_invSkin[1]*JfSkin[6]); 
  fSkin[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[9]+15.0*jacobgeo_invSkin[0]*JfSkin[8]); 
  fSkin[9] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[9]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[8]); 
  fSkin[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invSkin[1]*JfSkin[11]+15.0*jacobgeo_invSkin[0]*JfSkin[10]); 
  fSkin[11] = 0.04714045207910316*(15.0*jacobgeo_invSkin[0]*JfSkin[11]+15.000000000000002*jacobgeo_invSkin[1]*JfSkin[10]); 

  double fGhost[12];
  fGhost[0] = 0.7071067811865475*(JfGhost[1]*jacobgeo_invGhost[1]+JfGhost[0]*jacobgeo_invGhost[0]); 
  fGhost[1] = 0.7071067811865475*(JfGhost[0]*jacobgeo_invGhost[1]+jacobgeo_invGhost[0]*JfGhost[1]); 
  fGhost[2] = 0.7071067811865475*(jacobgeo_invGhost[1]*JfGhost[4]+jacobgeo_invGhost[0]*JfGhost[2]); 
  fGhost[3] = 0.7071067811865475*(jacobgeo_invGhost[1]*JfGhost[5]+jacobgeo_invGhost[0]*JfGhost[3]); 
  fGhost[4] = 0.7071067811865475*(jacobgeo_invGhost[0]*JfGhost[4]+jacobgeo_invGhost[1]*JfGhost[2]); 
  fGhost[5] = 0.7071067811865475*(jacobgeo_invGhost[0]*JfGhost[5]+jacobgeo_invGhost[1]*JfGhost[3]); 
  fGhost[6] = 0.7071067811865475*(jacobgeo_invGhost[1]*JfGhost[7]+jacobgeo_invGhost[0]*JfGhost[6]); 
  fGhost[7] = 0.7071067811865475*(jacobgeo_invGhost[0]*JfGhost[7]+jacobgeo_invGhost[1]*JfGhost[6]); 
  fGhost[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invGhost[1]*JfGhost[9]+15.0*jacobgeo_invGhost[0]*JfGhost[8]); 
  fGhost[9] = 0.04714045207910316*(15.0*jacobgeo_invGhost[0]*JfGhost[9]+15.000000000000002*jacobgeo_invGhost[1]*JfGhost[8]); 
  fGhost[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invGhost[1]*JfGhost[11]+15.0*jacobgeo_invGhost[0]*JfGhost[10]); 
  fGhost[11] = 0.04714045207910316*(15.0*jacobgeo_invGhost[0]*JfGhost[11]+15.000000000000002*jacobgeo_invGhost[1]*JfGhost[10]); 

  double boundSurf_incr[12] = {0.0}; 

  boundSurf_incr[0] = 0.33145630368119394*fSkin[1]*nuSkin[1]+0.33145630368119394*fGhost[1]*nuSkin[1]+0.34445949507888407*fSkin[0]*nuSkin[1]-0.34445949507888407*fGhost[0]*nuSkin[1]-0.33145630368119394*fSkin[1]*nuGhost[1]-0.33145630368119394*fGhost[1]*nuGhost[1]-0.34445949507888407*fSkin[0]*nuGhost[1]+0.34445949507888407*fGhost[0]*nuGhost[1]+0.19136638615493565*nuSkin[0]*fSkin[1]+0.19136638615493565*nuGhost[0]*fSkin[1]+0.19136638615493565*nuSkin[0]*fGhost[1]+0.19136638615493565*nuGhost[0]*fGhost[1]+0.19887378220871635*fSkin[0]*nuSkin[0]-0.19887378220871635*fGhost[0]*nuSkin[0]+0.19887378220871635*fSkin[0]*nuGhost[0]-0.19887378220871635*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[1] = 0.5740991584648069*fSkin[1]*nuSkin[1]+0.5740991584648069*fGhost[1]*nuSkin[1]+0.5966213466261491*fSkin[0]*nuSkin[1]-0.5966213466261491*fGhost[0]*nuSkin[1]-0.5740991584648069*fSkin[1]*nuGhost[1]-0.5740991584648069*fGhost[1]*nuGhost[1]-0.5966213466261491*fSkin[0]*nuGhost[1]+0.5966213466261491*fGhost[0]*nuGhost[1]+0.33145630368119394*nuSkin[0]*fSkin[1]+0.33145630368119394*nuGhost[0]*fSkin[1]+0.33145630368119394*nuSkin[0]*fGhost[1]+0.33145630368119394*nuGhost[0]*fGhost[1]+0.34445949507888407*fSkin[0]*nuSkin[0]-0.34445949507888407*fGhost[0]*nuSkin[0]+0.34445949507888407*fSkin[0]*nuGhost[0]-0.34445949507888407*fGhost[0]*nuGhost[0]; 
  boundSurf_incr[2] = 0.33145630368119394*nuSkin[1]*fSkin[4]-0.33145630368119394*nuGhost[1]*fSkin[4]+0.19136638615493565*nuSkin[0]*fSkin[4]+0.19136638615493565*nuGhost[0]*fSkin[4]+0.33145630368119394*nuSkin[1]*fGhost[4]-0.33145630368119394*nuGhost[1]*fGhost[4]+0.19136638615493565*nuSkin[0]*fGhost[4]+0.19136638615493565*nuGhost[0]*fGhost[4]+0.34445949507888407*nuSkin[1]*fSkin[2]-0.34445949507888407*nuGhost[1]*fSkin[2]+0.19887378220871635*nuSkin[0]*fSkin[2]+0.19887378220871635*nuGhost[0]*fSkin[2]-0.34445949507888407*nuSkin[1]*fGhost[2]+0.34445949507888407*nuGhost[1]*fGhost[2]-0.19887378220871635*nuSkin[0]*fGhost[2]-0.19887378220871635*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[3] = 0.33145630368119394*nuSkin[1]*fSkin[5]-0.33145630368119394*nuGhost[1]*fSkin[5]+0.19136638615493565*nuSkin[0]*fSkin[5]+0.19136638615493565*nuGhost[0]*fSkin[5]+0.33145630368119394*nuSkin[1]*fGhost[5]-0.33145630368119394*nuGhost[1]*fGhost[5]+0.19136638615493565*nuSkin[0]*fGhost[5]+0.19136638615493565*nuGhost[0]*fGhost[5]+0.34445949507888407*nuSkin[1]*fSkin[3]-0.34445949507888407*nuGhost[1]*fSkin[3]+0.19887378220871635*nuSkin[0]*fSkin[3]+0.19887378220871635*nuGhost[0]*fSkin[3]-0.34445949507888407*nuSkin[1]*fGhost[3]+0.34445949507888407*nuGhost[1]*fGhost[3]-0.19887378220871635*nuSkin[0]*fGhost[3]-0.19887378220871635*nuGhost[0]*fGhost[3]; 
  boundSurf_incr[4] = 0.5740991584648069*nuSkin[1]*fSkin[4]-0.5740991584648069*nuGhost[1]*fSkin[4]+0.33145630368119394*nuSkin[0]*fSkin[4]+0.33145630368119394*nuGhost[0]*fSkin[4]+0.5740991584648069*nuSkin[1]*fGhost[4]-0.5740991584648069*nuGhost[1]*fGhost[4]+0.33145630368119394*nuSkin[0]*fGhost[4]+0.33145630368119394*nuGhost[0]*fGhost[4]+0.5966213466261491*nuSkin[1]*fSkin[2]-0.5966213466261491*nuGhost[1]*fSkin[2]+0.34445949507888407*nuSkin[0]*fSkin[2]+0.34445949507888407*nuGhost[0]*fSkin[2]-0.5966213466261491*nuSkin[1]*fGhost[2]+0.5966213466261491*nuGhost[1]*fGhost[2]-0.34445949507888407*nuSkin[0]*fGhost[2]-0.34445949507888407*nuGhost[0]*fGhost[2]; 
  boundSurf_incr[5] = 0.5740991584648069*nuSkin[1]*fSkin[5]-0.5740991584648069*nuGhost[1]*fSkin[5]+0.33145630368119394*nuSkin[0]*fSkin[5]+0.33145630368119394*nuGhost[0]*fSkin[5]+0.5740991584648069*nuSkin[1]*fGhost[5]-0.5740991584648069*nuGhost[1]*fGhost[5]+0.33145630368119394*nuSkin[0]*fGhost[5]+0.33145630368119394*nuGhost[0]*fGhost[5]+0.5966213466261491*nuSkin[1]*fSkin[3]-0.5966213466261491*nuGhost[1]*fSkin[3]+0.34445949507888407*nuSkin[0]*fSkin[3]+0.34445949507888407*nuGhost[0]*fSkin[3]-0.5966213466261491*nuSkin[1]*fGhost[3]+0.5966213466261491*nuGhost[1]*fGhost[3]-0.34445949507888407*nuSkin[0]*fGhost[3]-0.34445949507888407*nuGhost[0]*fGhost[3]; 
  boundSurf_incr[6] = 0.33145630368119394*nuSkin[1]*fSkin[7]-0.33145630368119394*nuGhost[1]*fSkin[7]+0.19136638615493565*nuSkin[0]*fSkin[7]+0.19136638615493565*nuGhost[0]*fSkin[7]+0.33145630368119394*nuSkin[1]*fGhost[7]-0.33145630368119394*nuGhost[1]*fGhost[7]+0.19136638615493565*nuSkin[0]*fGhost[7]+0.19136638615493565*nuGhost[0]*fGhost[7]+0.34445949507888407*nuSkin[1]*fSkin[6]-0.34445949507888407*nuGhost[1]*fSkin[6]+0.19887378220871635*nuSkin[0]*fSkin[6]+0.19887378220871635*nuGhost[0]*fSkin[6]-0.34445949507888407*nuSkin[1]*fGhost[6]+0.34445949507888407*nuGhost[1]*fGhost[6]-0.19887378220871635*nuSkin[0]*fGhost[6]-0.19887378220871635*nuGhost[0]*fGhost[6]; 
  boundSurf_incr[7] = 0.5740991584648069*nuSkin[1]*fSkin[7]-0.5740991584648069*nuGhost[1]*fSkin[7]+0.33145630368119394*nuSkin[0]*fSkin[7]+0.33145630368119394*nuGhost[0]*fSkin[7]+0.5740991584648069*nuSkin[1]*fGhost[7]-0.5740991584648069*nuGhost[1]*fGhost[7]+0.33145630368119394*nuSkin[0]*fGhost[7]+0.33145630368119394*nuGhost[0]*fGhost[7]+0.5966213466261491*nuSkin[1]*fSkin[6]-0.5966213466261491*nuGhost[1]*fSkin[6]+0.34445949507888407*nuSkin[0]*fSkin[6]+0.34445949507888407*nuGhost[0]*fSkin[6]-0.5966213466261491*nuSkin[1]*fGhost[6]+0.5966213466261491*nuGhost[1]*fGhost[6]-0.34445949507888407*nuSkin[0]*fGhost[6]-0.34445949507888407*nuGhost[0]*fGhost[6]; 
  boundSurf_incr[8] = 0.3314563036811939*nuSkin[1]*fSkin[9]-0.3314563036811939*nuGhost[1]*fSkin[9]+0.19136638615493567*nuSkin[0]*fSkin[9]+0.19136638615493567*nuGhost[0]*fSkin[9]+0.3314563036811939*nuSkin[1]*fGhost[9]-0.3314563036811939*nuGhost[1]*fGhost[9]+0.19136638615493567*nuSkin[0]*fGhost[9]+0.19136638615493567*nuGhost[0]*fGhost[9]+0.34445949507888407*nuSkin[1]*fSkin[8]-0.34445949507888407*nuGhost[1]*fSkin[8]+0.19887378220871635*nuSkin[0]*fSkin[8]+0.19887378220871635*nuGhost[0]*fSkin[8]-0.34445949507888407*nuSkin[1]*fGhost[8]+0.34445949507888407*nuGhost[1]*fGhost[8]-0.19887378220871635*nuSkin[0]*fGhost[8]-0.19887378220871635*nuGhost[0]*fGhost[8]; 
  boundSurf_incr[9] = 0.5740991584648069*nuSkin[1]*fSkin[9]-0.5740991584648069*nuGhost[1]*fSkin[9]+0.33145630368119394*nuSkin[0]*fSkin[9]+0.33145630368119394*nuGhost[0]*fSkin[9]+0.5740991584648069*nuSkin[1]*fGhost[9]-0.5740991584648069*nuGhost[1]*fGhost[9]+0.33145630368119394*nuSkin[0]*fGhost[9]+0.33145630368119394*nuGhost[0]*fGhost[9]+0.5966213466261489*nuSkin[1]*fSkin[8]-0.5966213466261489*nuGhost[1]*fSkin[8]+0.3444594950788842*nuSkin[0]*fSkin[8]+0.3444594950788842*nuGhost[0]*fSkin[8]-0.5966213466261489*nuSkin[1]*fGhost[8]+0.5966213466261489*nuGhost[1]*fGhost[8]-0.3444594950788842*nuSkin[0]*fGhost[8]-0.3444594950788842*nuGhost[0]*fGhost[8]; 
  boundSurf_incr[10] = 0.3314563036811939*nuSkin[1]*fSkin[11]-0.3314563036811939*nuGhost[1]*fSkin[11]+0.19136638615493567*nuSkin[0]*fSkin[11]+0.19136638615493567*nuGhost[0]*fSkin[11]+0.3314563036811939*nuSkin[1]*fGhost[11]-0.3314563036811939*nuGhost[1]*fGhost[11]+0.19136638615493567*nuSkin[0]*fGhost[11]+0.19136638615493567*nuGhost[0]*fGhost[11]+0.34445949507888407*nuSkin[1]*fSkin[10]-0.34445949507888407*nuGhost[1]*fSkin[10]+0.19887378220871635*nuSkin[0]*fSkin[10]+0.19887378220871635*nuGhost[0]*fSkin[10]-0.34445949507888407*nuSkin[1]*fGhost[10]+0.34445949507888407*nuGhost[1]*fGhost[10]-0.19887378220871635*nuSkin[0]*fGhost[10]-0.19887378220871635*nuGhost[0]*fGhost[10]; 
  boundSurf_incr[11] = 0.5740991584648069*nuSkin[1]*fSkin[11]-0.5740991584648069*nuGhost[1]*fSkin[11]+0.33145630368119394*nuSkin[0]*fSkin[11]+0.33145630368119394*nuGhost[0]*fSkin[11]+0.5740991584648069*nuSkin[1]*fGhost[11]-0.5740991584648069*nuGhost[1]*fGhost[11]+0.33145630368119394*nuSkin[0]*fGhost[11]+0.33145630368119394*nuGhost[0]*fGhost[11]+0.5966213466261489*nuSkin[1]*fSkin[10]-0.5966213466261489*nuGhost[1]*fSkin[10]+0.3444594950788842*nuSkin[0]*fSkin[10]+0.3444594950788842*nuGhost[0]*fSkin[10]-0.5966213466261489*nuSkin[1]*fGhost[10]+0.5966213466261489*nuGhost[1]*fGhost[10]-0.3444594950788842*nuSkin[0]*fGhost[10]-0.3444594950788842*nuGhost[0]*fGhost[10]; 

  out[0] += boundSurf_incr[0]*rdx2fac; 
  out[1] += boundSurf_incr[1]*rdx2fac; 
  out[2] += boundSurf_incr[2]*rdx2fac; 
  out[3] += boundSurf_incr[3]*rdx2fac; 
  out[4] += boundSurf_incr[4]*rdx2fac; 
  out[5] += boundSurf_incr[5]*rdx2fac; 
  out[6] += boundSurf_incr[6]*rdx2fac; 
  out[7] += boundSurf_incr[7]*rdx2fac; 
  out[8] += boundSurf_incr[8]*rdx2fac; 
  out[9] += boundSurf_incr[9]*rdx2fac; 
  out[10] += boundSurf_incr[10]*rdx2fac; 
  out[11] += boundSurf_incr[11]*rdx2fac; 

  return 0.;
}

