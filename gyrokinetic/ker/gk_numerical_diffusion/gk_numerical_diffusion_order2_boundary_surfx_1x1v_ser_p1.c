#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_lower_zero_flux_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[6];
  fEdge[0] = 0.7071067811865476*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865476*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865476*(jacobgeo_invEdge[1]*JfEdge[3]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865476*(jacobgeo_invEdge[0]*JfEdge[3]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invEdge[1]*JfEdge[5]+21.213203435596427*jacobgeo_invEdge[0]*JfEdge[4]); 
  fEdge[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invEdge[0]*JfEdge[5]+21.21320343559643*jacobgeo_invEdge[1]*JfEdge[4]); 

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double edgeSurf_incr[6] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]-0.34445949507888407*fSkin[0]*nuSkin[1]+0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]+0.34445949507888407*fSkin[0]*nuEdge[1]-0.34445949507888407*fEdge[0]*nuEdge[1]-0.19136638615493565*nuSkin[0]*fSkin[1]-0.19136638615493565*nuEdge[0]*fSkin[1]-0.19136638615493565*nuSkin[0]*fEdge[1]-0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = -(0.5740991584648069*fSkin[1]*nuSkin[1])-0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]+0.5740991584648069*fSkin[1]*nuEdge[1]+0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]-0.34445949507888407*fSkin[0]*nuSkin[0]+0.34445949507888407*fEdge[0]*nuSkin[0]-0.34445949507888407*fSkin[0]*nuEdge[0]+0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[3])+0.33145630368119394*nuEdge[1]*fSkin[3]-0.19136638615493565*nuSkin[0]*fSkin[3]-0.19136638615493565*nuEdge[0]*fSkin[3]-0.33145630368119394*nuSkin[1]*fEdge[3]+0.33145630368119394*nuEdge[1]*fEdge[3]-0.19136638615493565*nuSkin[0]*fEdge[3]-0.19136638615493565*nuEdge[0]*fEdge[3]-0.34445949507888407*nuSkin[1]*fSkin[2]+0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]+0.34445949507888407*nuSkin[1]*fEdge[2]-0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.5740991584648069*nuSkin[1]*fSkin[3])+0.5740991584648069*nuEdge[1]*fSkin[3]-0.33145630368119394*nuSkin[0]*fSkin[3]-0.33145630368119394*nuEdge[0]*fSkin[3]-0.5740991584648069*nuSkin[1]*fEdge[3]+0.5740991584648069*nuEdge[1]*fEdge[3]-0.33145630368119394*nuSkin[0]*fEdge[3]-0.33145630368119394*nuEdge[0]*fEdge[3]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]-0.34445949507888407*nuSkin[0]*fSkin[2]-0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]+0.34445949507888407*nuSkin[0]*fEdge[2]+0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(0.3314563036811939*nuSkin[1]*fSkin[5])+0.3314563036811939*nuEdge[1]*fSkin[5]-0.19136638615493567*nuSkin[0]*fSkin[5]-0.19136638615493567*nuEdge[0]*fSkin[5]-0.3314563036811939*nuSkin[1]*fEdge[5]+0.3314563036811939*nuEdge[1]*fEdge[5]-0.19136638615493567*nuSkin[0]*fEdge[5]-0.19136638615493567*nuEdge[0]*fEdge[5]-0.34445949507888407*nuSkin[1]*fSkin[4]+0.34445949507888407*nuEdge[1]*fSkin[4]-0.19887378220871635*nuSkin[0]*fSkin[4]-0.19887378220871635*nuEdge[0]*fSkin[4]+0.34445949507888407*nuSkin[1]*fEdge[4]-0.34445949507888407*nuEdge[1]*fEdge[4]+0.19887378220871635*nuSkin[0]*fEdge[4]+0.19887378220871635*nuEdge[0]*fEdge[4]; 
  edgeSurf_incr[5] = -(0.5740991584648069*nuSkin[1]*fSkin[5])+0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]-0.5740991584648069*nuSkin[1]*fEdge[5]+0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261489*nuSkin[1]*fSkin[4]+0.5966213466261489*nuEdge[1]*fSkin[4]-0.3444594950788842*nuSkin[0]*fSkin[4]-0.3444594950788842*nuEdge[0]*fSkin[4]+0.5966213466261489*nuSkin[1]*fEdge[4]-0.5966213466261489*nuEdge[1]*fEdge[4]+0.3444594950788842*nuSkin[0]*fEdge[4]+0.3444594950788842*nuEdge[0]*fEdge[4]; 


  out[0] += edgeSurf_incr[0]*rdx2fac; 
  out[1] += edgeSurf_incr[1]*rdx2fac; 
  out[2] += edgeSurf_incr[2]*rdx2fac; 
  out[3] += edgeSurf_incr[3]*rdx2fac; 
  out[4] += edgeSurf_incr[4]*rdx2fac; 
  out[5] += edgeSurf_incr[5]*rdx2fac; 

  return 0.;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_lower_bound_local_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[6];
  fEdge[0] = 0.7071067811865476*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865476*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865476*(jacobgeo_invEdge[1]*JfEdge[3]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865476*(jacobgeo_invEdge[0]*JfEdge[3]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invEdge[1]*JfEdge[5]+21.213203435596427*jacobgeo_invEdge[0]*JfEdge[4]); 
  fEdge[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invEdge[0]*JfEdge[5]+21.21320343559643*jacobgeo_invEdge[1]*JfEdge[4]); 

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double edgeSurf_incr[6] = {0.0}; 
  double boundSurf_incr[6] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]-0.34445949507888407*fSkin[0]*nuSkin[1]+0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]+0.34445949507888407*fSkin[0]*nuEdge[1]-0.34445949507888407*fEdge[0]*nuEdge[1]-0.19136638615493565*nuSkin[0]*fSkin[1]-0.19136638615493565*nuEdge[0]*fSkin[1]-0.19136638615493565*nuSkin[0]*fEdge[1]-0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = -(0.5740991584648069*fSkin[1]*nuSkin[1])-0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]+0.5740991584648069*fSkin[1]*nuEdge[1]+0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]-0.34445949507888407*fSkin[0]*nuSkin[0]+0.34445949507888407*fEdge[0]*nuSkin[0]-0.34445949507888407*fSkin[0]*nuEdge[0]+0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[3])+0.33145630368119394*nuEdge[1]*fSkin[3]-0.19136638615493565*nuSkin[0]*fSkin[3]-0.19136638615493565*nuEdge[0]*fSkin[3]-0.33145630368119394*nuSkin[1]*fEdge[3]+0.33145630368119394*nuEdge[1]*fEdge[3]-0.19136638615493565*nuSkin[0]*fEdge[3]-0.19136638615493565*nuEdge[0]*fEdge[3]-0.34445949507888407*nuSkin[1]*fSkin[2]+0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]+0.34445949507888407*nuSkin[1]*fEdge[2]-0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.5740991584648069*nuSkin[1]*fSkin[3])+0.5740991584648069*nuEdge[1]*fSkin[3]-0.33145630368119394*nuSkin[0]*fSkin[3]-0.33145630368119394*nuEdge[0]*fSkin[3]-0.5740991584648069*nuSkin[1]*fEdge[3]+0.5740991584648069*nuEdge[1]*fEdge[3]-0.33145630368119394*nuSkin[0]*fEdge[3]-0.33145630368119394*nuEdge[0]*fEdge[3]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]-0.34445949507888407*nuSkin[0]*fSkin[2]-0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]+0.34445949507888407*nuSkin[0]*fEdge[2]+0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(0.3314563036811939*nuSkin[1]*fSkin[5])+0.3314563036811939*nuEdge[1]*fSkin[5]-0.19136638615493567*nuSkin[0]*fSkin[5]-0.19136638615493567*nuEdge[0]*fSkin[5]-0.3314563036811939*nuSkin[1]*fEdge[5]+0.3314563036811939*nuEdge[1]*fEdge[5]-0.19136638615493567*nuSkin[0]*fEdge[5]-0.19136638615493567*nuEdge[0]*fEdge[5]-0.34445949507888407*nuSkin[1]*fSkin[4]+0.34445949507888407*nuEdge[1]*fSkin[4]-0.19887378220871635*nuSkin[0]*fSkin[4]-0.19887378220871635*nuEdge[0]*fSkin[4]+0.34445949507888407*nuSkin[1]*fEdge[4]-0.34445949507888407*nuEdge[1]*fEdge[4]+0.19887378220871635*nuSkin[0]*fEdge[4]+0.19887378220871635*nuEdge[0]*fEdge[4]; 
  edgeSurf_incr[5] = -(0.5740991584648069*nuSkin[1]*fSkin[5])+0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]-0.5740991584648069*nuSkin[1]*fEdge[5]+0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261489*nuSkin[1]*fSkin[4]+0.5966213466261489*nuEdge[1]*fSkin[4]-0.3444594950788842*nuSkin[0]*fSkin[4]-0.3444594950788842*nuEdge[0]*fSkin[4]+0.5966213466261489*nuSkin[1]*fEdge[4]-0.5966213466261489*nuEdge[1]*fEdge[4]+0.3444594950788842*nuSkin[0]*fEdge[4]+0.3444594950788842*nuEdge[0]*fEdge[4]; 

  boundSurf_incr[0] = 1.060660171779821*fSkin[1]*nuSkin[1]-0.6123724356957944*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[1] = 1.060660171779821*nuSkin[0]*fSkin[1]-1.837117307087383*fSkin[1]*nuSkin[1]; 
  boundSurf_incr[2] = 1.060660171779821*nuSkin[1]*fSkin[3]-0.6123724356957944*nuSkin[0]*fSkin[3]; 
  boundSurf_incr[3] = 1.060660171779821*nuSkin[0]*fSkin[3]-1.837117307087383*nuSkin[1]*fSkin[3]; 
  boundSurf_incr[4] = 1.060660171779821*nuSkin[1]*fSkin[5]-0.6123724356957944*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[5] = 1.060660171779821*nuSkin[0]*fSkin[5]-1.837117307087383*nuSkin[1]*fSkin[5]; 

  out[0] += (edgeSurf_incr[0]+boundSurf_incr[0])*rdx2fac; 
  out[1] += (edgeSurf_incr[1]+boundSurf_incr[1])*rdx2fac; 
  out[2] += (edgeSurf_incr[2]+boundSurf_incr[2])*rdx2fac; 
  out[3] += (edgeSurf_incr[3]+boundSurf_incr[3])*rdx2fac; 
  out[4] += (edgeSurf_incr[4]+boundSurf_incr[4])*rdx2fac; 
  out[5] += (edgeSurf_incr[5]+boundSurf_incr[5])*rdx2fac; 

  return 0.;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_upper_zero_flux_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[6];
  fEdge[0] = 0.7071067811865476*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865476*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865476*(jacobgeo_invEdge[1]*JfEdge[3]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865476*(jacobgeo_invEdge[0]*JfEdge[3]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invEdge[1]*JfEdge[5]+21.213203435596427*jacobgeo_invEdge[0]*JfEdge[4]); 
  fEdge[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invEdge[0]*JfEdge[5]+21.21320343559643*jacobgeo_invEdge[1]*JfEdge[4]); 

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double edgeSurf_incr[6] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]+0.34445949507888407*fSkin[0]*nuSkin[1]-0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]-0.34445949507888407*fSkin[0]*nuEdge[1]+0.34445949507888407*fEdge[0]*nuEdge[1]+0.19136638615493565*nuSkin[0]*fSkin[1]+0.19136638615493565*nuEdge[0]*fSkin[1]+0.19136638615493565*nuSkin[0]*fEdge[1]+0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = 0.5740991584648069*fSkin[1]*nuSkin[1]+0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]-0.5740991584648069*fSkin[1]*nuEdge[1]-0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]+0.34445949507888407*fSkin[0]*nuSkin[0]-0.34445949507888407*fEdge[0]*nuSkin[0]+0.34445949507888407*fSkin[0]*nuEdge[0]-0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[3])+0.33145630368119394*nuEdge[1]*fSkin[3]+0.19136638615493565*nuSkin[0]*fSkin[3]+0.19136638615493565*nuEdge[0]*fSkin[3]-0.33145630368119394*nuSkin[1]*fEdge[3]+0.33145630368119394*nuEdge[1]*fEdge[3]+0.19136638615493565*nuSkin[0]*fEdge[3]+0.19136638615493565*nuEdge[0]*fEdge[3]+0.34445949507888407*nuSkin[1]*fSkin[2]-0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]-0.34445949507888407*nuSkin[1]*fEdge[2]+0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5740991584648069*nuSkin[1]*fSkin[3]-0.5740991584648069*nuEdge[1]*fSkin[3]-0.33145630368119394*nuSkin[0]*fSkin[3]-0.33145630368119394*nuEdge[0]*fSkin[3]+0.5740991584648069*nuSkin[1]*fEdge[3]-0.5740991584648069*nuEdge[1]*fEdge[3]-0.33145630368119394*nuSkin[0]*fEdge[3]-0.33145630368119394*nuEdge[0]*fEdge[3]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]+0.34445949507888407*nuSkin[0]*fSkin[2]+0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]-0.34445949507888407*nuSkin[0]*fEdge[2]-0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(0.3314563036811939*nuSkin[1]*fSkin[5])+0.3314563036811939*nuEdge[1]*fSkin[5]+0.19136638615493567*nuSkin[0]*fSkin[5]+0.19136638615493567*nuEdge[0]*fSkin[5]-0.3314563036811939*nuSkin[1]*fEdge[5]+0.3314563036811939*nuEdge[1]*fEdge[5]+0.19136638615493567*nuSkin[0]*fEdge[5]+0.19136638615493567*nuEdge[0]*fEdge[5]+0.34445949507888407*nuSkin[1]*fSkin[4]-0.34445949507888407*nuEdge[1]*fSkin[4]-0.19887378220871635*nuSkin[0]*fSkin[4]-0.19887378220871635*nuEdge[0]*fSkin[4]-0.34445949507888407*nuSkin[1]*fEdge[4]+0.34445949507888407*nuEdge[1]*fEdge[4]+0.19887378220871635*nuSkin[0]*fEdge[4]+0.19887378220871635*nuEdge[0]*fEdge[4]; 
  edgeSurf_incr[5] = 0.5740991584648069*nuSkin[1]*fSkin[5]-0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]+0.5740991584648069*nuSkin[1]*fEdge[5]-0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261489*nuSkin[1]*fSkin[4]+0.5966213466261489*nuEdge[1]*fSkin[4]+0.3444594950788842*nuSkin[0]*fSkin[4]+0.3444594950788842*nuEdge[0]*fSkin[4]+0.5966213466261489*nuSkin[1]*fEdge[4]-0.5966213466261489*nuEdge[1]*fEdge[4]-0.3444594950788842*nuSkin[0]*fEdge[4]-0.3444594950788842*nuEdge[0]*fEdge[4]; 


  out[0] += edgeSurf_incr[0]*rdx2fac; 
  out[1] += edgeSurf_incr[1]*rdx2fac; 
  out[2] += edgeSurf_incr[2]*rdx2fac; 
  out[3] += edgeSurf_incr[3]*rdx2fac; 
  out[4] += edgeSurf_incr[4]*rdx2fac; 
  out[5] += edgeSurf_incr[5]*rdx2fac; 

  return 0.;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_upper_bound_local_1x1v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[6];
  fEdge[0] = 0.7071067811865476*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865476*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865476*(jacobgeo_invEdge[1]*JfEdge[3]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865476*(jacobgeo_invEdge[0]*JfEdge[3]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invEdge[1]*JfEdge[5]+21.213203435596427*jacobgeo_invEdge[0]*JfEdge[4]); 
  fEdge[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invEdge[0]*JfEdge[5]+21.21320343559643*jacobgeo_invEdge[1]*JfEdge[4]); 

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(JfSkin[1]*jacobgeo_invSkin[1]+JfSkin[0]*jacobgeo_invSkin[0]); 
  fSkin[1] = 0.7071067811865476*(JfSkin[0]*jacobgeo_invSkin[1]+jacobgeo_invSkin[0]*JfSkin[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_invSkin[1]*JfSkin[3]+jacobgeo_invSkin[0]*JfSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_invSkin[0]*JfSkin[3]+jacobgeo_invSkin[1]*JfSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_invSkin[1]*JfSkin[5]+21.213203435596427*jacobgeo_invSkin[0]*JfSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_invSkin[0]*JfSkin[5]+21.21320343559643*jacobgeo_invSkin[1]*JfSkin[4]); 

  double edgeSurf_incr[6] = {0.0}; 
  double boundSurf_incr[6] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]+0.34445949507888407*fSkin[0]*nuSkin[1]-0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]-0.34445949507888407*fSkin[0]*nuEdge[1]+0.34445949507888407*fEdge[0]*nuEdge[1]+0.19136638615493565*nuSkin[0]*fSkin[1]+0.19136638615493565*nuEdge[0]*fSkin[1]+0.19136638615493565*nuSkin[0]*fEdge[1]+0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = 0.5740991584648069*fSkin[1]*nuSkin[1]+0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]-0.5740991584648069*fSkin[1]*nuEdge[1]-0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]+0.34445949507888407*fSkin[0]*nuSkin[0]-0.34445949507888407*fEdge[0]*nuSkin[0]+0.34445949507888407*fSkin[0]*nuEdge[0]-0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[3])+0.33145630368119394*nuEdge[1]*fSkin[3]+0.19136638615493565*nuSkin[0]*fSkin[3]+0.19136638615493565*nuEdge[0]*fSkin[3]-0.33145630368119394*nuSkin[1]*fEdge[3]+0.33145630368119394*nuEdge[1]*fEdge[3]+0.19136638615493565*nuSkin[0]*fEdge[3]+0.19136638615493565*nuEdge[0]*fEdge[3]+0.34445949507888407*nuSkin[1]*fSkin[2]-0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]-0.34445949507888407*nuSkin[1]*fEdge[2]+0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5740991584648069*nuSkin[1]*fSkin[3]-0.5740991584648069*nuEdge[1]*fSkin[3]-0.33145630368119394*nuSkin[0]*fSkin[3]-0.33145630368119394*nuEdge[0]*fSkin[3]+0.5740991584648069*nuSkin[1]*fEdge[3]-0.5740991584648069*nuEdge[1]*fEdge[3]-0.33145630368119394*nuSkin[0]*fEdge[3]-0.33145630368119394*nuEdge[0]*fEdge[3]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]+0.34445949507888407*nuSkin[0]*fSkin[2]+0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]-0.34445949507888407*nuSkin[0]*fEdge[2]-0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(0.3314563036811939*nuSkin[1]*fSkin[5])+0.3314563036811939*nuEdge[1]*fSkin[5]+0.19136638615493567*nuSkin[0]*fSkin[5]+0.19136638615493567*nuEdge[0]*fSkin[5]-0.3314563036811939*nuSkin[1]*fEdge[5]+0.3314563036811939*nuEdge[1]*fEdge[5]+0.19136638615493567*nuSkin[0]*fEdge[5]+0.19136638615493567*nuEdge[0]*fEdge[5]+0.34445949507888407*nuSkin[1]*fSkin[4]-0.34445949507888407*nuEdge[1]*fSkin[4]-0.19887378220871635*nuSkin[0]*fSkin[4]-0.19887378220871635*nuEdge[0]*fSkin[4]-0.34445949507888407*nuSkin[1]*fEdge[4]+0.34445949507888407*nuEdge[1]*fEdge[4]+0.19887378220871635*nuSkin[0]*fEdge[4]+0.19887378220871635*nuEdge[0]*fEdge[4]; 
  edgeSurf_incr[5] = 0.5740991584648069*nuSkin[1]*fSkin[5]-0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]+0.5740991584648069*nuSkin[1]*fEdge[5]-0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261489*nuSkin[1]*fSkin[4]+0.5966213466261489*nuEdge[1]*fSkin[4]+0.3444594950788842*nuSkin[0]*fSkin[4]+0.3444594950788842*nuEdge[0]*fSkin[4]+0.5966213466261489*nuSkin[1]*fEdge[4]-0.5966213466261489*nuEdge[1]*fEdge[4]-0.3444594950788842*nuSkin[0]*fEdge[4]-0.3444594950788842*nuEdge[0]*fEdge[4]; 

  boundSurf_incr[0] = 1.060660171779821*fSkin[1]*nuSkin[1]+0.6123724356957944*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[1] = 1.837117307087383*fSkin[1]*nuSkin[1]+1.060660171779821*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[2] = 1.060660171779821*nuSkin[1]*fSkin[3]+0.6123724356957944*nuSkin[0]*fSkin[3]; 
  boundSurf_incr[3] = 1.837117307087383*nuSkin[1]*fSkin[3]+1.060660171779821*nuSkin[0]*fSkin[3]; 
  boundSurf_incr[4] = 1.060660171779821*nuSkin[1]*fSkin[5]+0.6123724356957944*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[5] = 1.837117307087383*nuSkin[1]*fSkin[5]+1.060660171779821*nuSkin[0]*fSkin[5]; 

  out[0] += (edgeSurf_incr[0]+boundSurf_incr[0])*rdx2fac; 
  out[1] += (edgeSurf_incr[1]+boundSurf_incr[1])*rdx2fac; 
  out[2] += (edgeSurf_incr[2]+boundSurf_incr[2])*rdx2fac; 
  out[3] += (edgeSurf_incr[3]+boundSurf_incr[3])*rdx2fac; 
  out[4] += (edgeSurf_incr[4]+boundSurf_incr[4])*rdx2fac; 
  out[5] += (edgeSurf_incr[5]+boundSurf_incr[5])*rdx2fac; 

  return 0.;
}

