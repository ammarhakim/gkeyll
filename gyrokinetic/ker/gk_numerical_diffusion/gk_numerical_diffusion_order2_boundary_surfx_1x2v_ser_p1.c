#include <gkyl_gk_numerical_diffusion_kernels.h>

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_lower_zero_flux_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[12];
  fEdge[0] = 0.7071067811865475*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865475*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[4]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[5]+jacobgeo_invEdge[0]*JfEdge[3]); 
  fEdge[4] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[4]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[5] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[5]+jacobgeo_invEdge[1]*JfEdge[3]); 
  fEdge[6] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[7]+jacobgeo_invEdge[0]*JfEdge[6]); 
  fEdge[7] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[7]+jacobgeo_invEdge[1]*JfEdge[6]); 
  fEdge[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[9]+15.0*jacobgeo_invEdge[0]*JfEdge[8]); 
  fEdge[9] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[9]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[8]); 
  fEdge[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[11]+15.0*jacobgeo_invEdge[0]*JfEdge[10]); 
  fEdge[11] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[11]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[10]); 

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

  double edgeSurf_incr[12] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]-0.34445949507888407*fSkin[0]*nuSkin[1]+0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]+0.34445949507888407*fSkin[0]*nuEdge[1]-0.34445949507888407*fEdge[0]*nuEdge[1]-0.19136638615493565*nuSkin[0]*fSkin[1]-0.19136638615493565*nuEdge[0]*fSkin[1]-0.19136638615493565*nuSkin[0]*fEdge[1]-0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = -(0.5740991584648069*fSkin[1]*nuSkin[1])-0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]+0.5740991584648069*fSkin[1]*nuEdge[1]+0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]-0.34445949507888407*fSkin[0]*nuSkin[0]+0.34445949507888407*fEdge[0]*nuSkin[0]-0.34445949507888407*fSkin[0]*nuEdge[0]+0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[4])+0.33145630368119394*nuEdge[1]*fSkin[4]-0.19136638615493565*nuSkin[0]*fSkin[4]-0.19136638615493565*nuEdge[0]*fSkin[4]-0.33145630368119394*nuSkin[1]*fEdge[4]+0.33145630368119394*nuEdge[1]*fEdge[4]-0.19136638615493565*nuSkin[0]*fEdge[4]-0.19136638615493565*nuEdge[0]*fEdge[4]-0.34445949507888407*nuSkin[1]*fSkin[2]+0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]+0.34445949507888407*nuSkin[1]*fEdge[2]-0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.33145630368119394*nuSkin[1]*fSkin[5])+0.33145630368119394*nuEdge[1]*fSkin[5]-0.19136638615493565*nuSkin[0]*fSkin[5]-0.19136638615493565*nuEdge[0]*fSkin[5]-0.33145630368119394*nuSkin[1]*fEdge[5]+0.33145630368119394*nuEdge[1]*fEdge[5]-0.19136638615493565*nuSkin[0]*fEdge[5]-0.19136638615493565*nuEdge[0]*fEdge[5]-0.34445949507888407*nuSkin[1]*fSkin[3]+0.34445949507888407*nuEdge[1]*fSkin[3]-0.19887378220871635*nuSkin[0]*fSkin[3]-0.19887378220871635*nuEdge[0]*fSkin[3]+0.34445949507888407*nuSkin[1]*fEdge[3]-0.34445949507888407*nuEdge[1]*fEdge[3]+0.19887378220871635*nuSkin[0]*fEdge[3]+0.19887378220871635*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(0.5740991584648069*nuSkin[1]*fSkin[4])+0.5740991584648069*nuEdge[1]*fSkin[4]-0.33145630368119394*nuSkin[0]*fSkin[4]-0.33145630368119394*nuEdge[0]*fSkin[4]-0.5740991584648069*nuSkin[1]*fEdge[4]+0.5740991584648069*nuEdge[1]*fEdge[4]-0.33145630368119394*nuSkin[0]*fEdge[4]-0.33145630368119394*nuEdge[0]*fEdge[4]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]-0.34445949507888407*nuSkin[0]*fSkin[2]-0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]+0.34445949507888407*nuSkin[0]*fEdge[2]+0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[5] = -(0.5740991584648069*nuSkin[1]*fSkin[5])+0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]-0.5740991584648069*nuSkin[1]*fEdge[5]+0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261491*nuSkin[1]*fSkin[3]+0.5966213466261491*nuEdge[1]*fSkin[3]-0.34445949507888407*nuSkin[0]*fSkin[3]-0.34445949507888407*nuEdge[0]*fSkin[3]+0.5966213466261491*nuSkin[1]*fEdge[3]-0.5966213466261491*nuEdge[1]*fEdge[3]+0.34445949507888407*nuSkin[0]*fEdge[3]+0.34445949507888407*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.33145630368119394*nuSkin[1]*fSkin[7])+0.33145630368119394*nuEdge[1]*fSkin[7]-0.19136638615493565*nuSkin[0]*fSkin[7]-0.19136638615493565*nuEdge[0]*fSkin[7]-0.33145630368119394*nuSkin[1]*fEdge[7]+0.33145630368119394*nuEdge[1]*fEdge[7]-0.19136638615493565*nuSkin[0]*fEdge[7]-0.19136638615493565*nuEdge[0]*fEdge[7]-0.34445949507888407*nuSkin[1]*fSkin[6]+0.34445949507888407*nuEdge[1]*fSkin[6]-0.19887378220871635*nuSkin[0]*fSkin[6]-0.19887378220871635*nuEdge[0]*fSkin[6]+0.34445949507888407*nuSkin[1]*fEdge[6]-0.34445949507888407*nuEdge[1]*fEdge[6]+0.19887378220871635*nuSkin[0]*fEdge[6]+0.19887378220871635*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[7] = -(0.5740991584648069*nuSkin[1]*fSkin[7])+0.5740991584648069*nuEdge[1]*fSkin[7]-0.33145630368119394*nuSkin[0]*fSkin[7]-0.33145630368119394*nuEdge[0]*fSkin[7]-0.5740991584648069*nuSkin[1]*fEdge[7]+0.5740991584648069*nuEdge[1]*fEdge[7]-0.33145630368119394*nuSkin[0]*fEdge[7]-0.33145630368119394*nuEdge[0]*fEdge[7]-0.5966213466261491*nuSkin[1]*fSkin[6]+0.5966213466261491*nuEdge[1]*fSkin[6]-0.34445949507888407*nuSkin[0]*fSkin[6]-0.34445949507888407*nuEdge[0]*fSkin[6]+0.5966213466261491*nuSkin[1]*fEdge[6]-0.5966213466261491*nuEdge[1]*fEdge[6]+0.34445949507888407*nuSkin[0]*fEdge[6]+0.34445949507888407*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[8] = -(0.3314563036811939*nuSkin[1]*fSkin[9])+0.3314563036811939*nuEdge[1]*fSkin[9]-0.19136638615493567*nuSkin[0]*fSkin[9]-0.19136638615493567*nuEdge[0]*fSkin[9]-0.3314563036811939*nuSkin[1]*fEdge[9]+0.3314563036811939*nuEdge[1]*fEdge[9]-0.19136638615493567*nuSkin[0]*fEdge[9]-0.19136638615493567*nuEdge[0]*fEdge[9]-0.34445949507888407*nuSkin[1]*fSkin[8]+0.34445949507888407*nuEdge[1]*fSkin[8]-0.19887378220871635*nuSkin[0]*fSkin[8]-0.19887378220871635*nuEdge[0]*fSkin[8]+0.34445949507888407*nuSkin[1]*fEdge[8]-0.34445949507888407*nuEdge[1]*fEdge[8]+0.19887378220871635*nuSkin[0]*fEdge[8]+0.19887378220871635*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[9] = -(0.5740991584648069*nuSkin[1]*fSkin[9])+0.5740991584648069*nuEdge[1]*fSkin[9]-0.33145630368119394*nuSkin[0]*fSkin[9]-0.33145630368119394*nuEdge[0]*fSkin[9]-0.5740991584648069*nuSkin[1]*fEdge[9]+0.5740991584648069*nuEdge[1]*fEdge[9]-0.33145630368119394*nuSkin[0]*fEdge[9]-0.33145630368119394*nuEdge[0]*fEdge[9]-0.5966213466261489*nuSkin[1]*fSkin[8]+0.5966213466261489*nuEdge[1]*fSkin[8]-0.3444594950788842*nuSkin[0]*fSkin[8]-0.3444594950788842*nuEdge[0]*fSkin[8]+0.5966213466261489*nuSkin[1]*fEdge[8]-0.5966213466261489*nuEdge[1]*fEdge[8]+0.3444594950788842*nuSkin[0]*fEdge[8]+0.3444594950788842*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[10] = -(0.3314563036811939*nuSkin[1]*fSkin[11])+0.3314563036811939*nuEdge[1]*fSkin[11]-0.19136638615493567*nuSkin[0]*fSkin[11]-0.19136638615493567*nuEdge[0]*fSkin[11]-0.3314563036811939*nuSkin[1]*fEdge[11]+0.3314563036811939*nuEdge[1]*fEdge[11]-0.19136638615493567*nuSkin[0]*fEdge[11]-0.19136638615493567*nuEdge[0]*fEdge[11]-0.34445949507888407*nuSkin[1]*fSkin[10]+0.34445949507888407*nuEdge[1]*fSkin[10]-0.19887378220871635*nuSkin[0]*fSkin[10]-0.19887378220871635*nuEdge[0]*fSkin[10]+0.34445949507888407*nuSkin[1]*fEdge[10]-0.34445949507888407*nuEdge[1]*fEdge[10]+0.19887378220871635*nuSkin[0]*fEdge[10]+0.19887378220871635*nuEdge[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.5740991584648069*nuSkin[1]*fSkin[11])+0.5740991584648069*nuEdge[1]*fSkin[11]-0.33145630368119394*nuSkin[0]*fSkin[11]-0.33145630368119394*nuEdge[0]*fSkin[11]-0.5740991584648069*nuSkin[1]*fEdge[11]+0.5740991584648069*nuEdge[1]*fEdge[11]-0.33145630368119394*nuSkin[0]*fEdge[11]-0.33145630368119394*nuEdge[0]*fEdge[11]-0.5966213466261489*nuSkin[1]*fSkin[10]+0.5966213466261489*nuEdge[1]*fSkin[10]-0.3444594950788842*nuSkin[0]*fSkin[10]-0.3444594950788842*nuEdge[0]*fSkin[10]+0.5966213466261489*nuSkin[1]*fEdge[10]-0.5966213466261489*nuEdge[1]*fEdge[10]+0.3444594950788842*nuSkin[0]*fEdge[10]+0.3444594950788842*nuEdge[0]*fEdge[10]; 


  out[0] += edgeSurf_incr[0]*rdx2fac; 
  out[1] += edgeSurf_incr[1]*rdx2fac; 
  out[2] += edgeSurf_incr[2]*rdx2fac; 
  out[3] += edgeSurf_incr[3]*rdx2fac; 
  out[4] += edgeSurf_incr[4]*rdx2fac; 
  out[5] += edgeSurf_incr[5]*rdx2fac; 
  out[6] += edgeSurf_incr[6]*rdx2fac; 
  out[7] += edgeSurf_incr[7]*rdx2fac; 
  out[8] += edgeSurf_incr[8]*rdx2fac; 
  out[9] += edgeSurf_incr[9]*rdx2fac; 
  out[10] += edgeSurf_incr[10]*rdx2fac; 
  out[11] += edgeSurf_incr[11]*rdx2fac; 

  return 0.;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_lower_bound_local_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[12];
  fEdge[0] = 0.7071067811865475*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865475*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[4]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[5]+jacobgeo_invEdge[0]*JfEdge[3]); 
  fEdge[4] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[4]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[5] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[5]+jacobgeo_invEdge[1]*JfEdge[3]); 
  fEdge[6] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[7]+jacobgeo_invEdge[0]*JfEdge[6]); 
  fEdge[7] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[7]+jacobgeo_invEdge[1]*JfEdge[6]); 
  fEdge[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[9]+15.0*jacobgeo_invEdge[0]*JfEdge[8]); 
  fEdge[9] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[9]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[8]); 
  fEdge[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[11]+15.0*jacobgeo_invEdge[0]*JfEdge[10]); 
  fEdge[11] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[11]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[10]); 

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

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]-0.34445949507888407*fSkin[0]*nuSkin[1]+0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]+0.34445949507888407*fSkin[0]*nuEdge[1]-0.34445949507888407*fEdge[0]*nuEdge[1]-0.19136638615493565*nuSkin[0]*fSkin[1]-0.19136638615493565*nuEdge[0]*fSkin[1]-0.19136638615493565*nuSkin[0]*fEdge[1]-0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = -(0.5740991584648069*fSkin[1]*nuSkin[1])-0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]+0.5740991584648069*fSkin[1]*nuEdge[1]+0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]-0.34445949507888407*fSkin[0]*nuSkin[0]+0.34445949507888407*fEdge[0]*nuSkin[0]-0.34445949507888407*fSkin[0]*nuEdge[0]+0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[4])+0.33145630368119394*nuEdge[1]*fSkin[4]-0.19136638615493565*nuSkin[0]*fSkin[4]-0.19136638615493565*nuEdge[0]*fSkin[4]-0.33145630368119394*nuSkin[1]*fEdge[4]+0.33145630368119394*nuEdge[1]*fEdge[4]-0.19136638615493565*nuSkin[0]*fEdge[4]-0.19136638615493565*nuEdge[0]*fEdge[4]-0.34445949507888407*nuSkin[1]*fSkin[2]+0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]+0.34445949507888407*nuSkin[1]*fEdge[2]-0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.33145630368119394*nuSkin[1]*fSkin[5])+0.33145630368119394*nuEdge[1]*fSkin[5]-0.19136638615493565*nuSkin[0]*fSkin[5]-0.19136638615493565*nuEdge[0]*fSkin[5]-0.33145630368119394*nuSkin[1]*fEdge[5]+0.33145630368119394*nuEdge[1]*fEdge[5]-0.19136638615493565*nuSkin[0]*fEdge[5]-0.19136638615493565*nuEdge[0]*fEdge[5]-0.34445949507888407*nuSkin[1]*fSkin[3]+0.34445949507888407*nuEdge[1]*fSkin[3]-0.19887378220871635*nuSkin[0]*fSkin[3]-0.19887378220871635*nuEdge[0]*fSkin[3]+0.34445949507888407*nuSkin[1]*fEdge[3]-0.34445949507888407*nuEdge[1]*fEdge[3]+0.19887378220871635*nuSkin[0]*fEdge[3]+0.19887378220871635*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(0.5740991584648069*nuSkin[1]*fSkin[4])+0.5740991584648069*nuEdge[1]*fSkin[4]-0.33145630368119394*nuSkin[0]*fSkin[4]-0.33145630368119394*nuEdge[0]*fSkin[4]-0.5740991584648069*nuSkin[1]*fEdge[4]+0.5740991584648069*nuEdge[1]*fEdge[4]-0.33145630368119394*nuSkin[0]*fEdge[4]-0.33145630368119394*nuEdge[0]*fEdge[4]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]-0.34445949507888407*nuSkin[0]*fSkin[2]-0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]+0.34445949507888407*nuSkin[0]*fEdge[2]+0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[5] = -(0.5740991584648069*nuSkin[1]*fSkin[5])+0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]-0.5740991584648069*nuSkin[1]*fEdge[5]+0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261491*nuSkin[1]*fSkin[3]+0.5966213466261491*nuEdge[1]*fSkin[3]-0.34445949507888407*nuSkin[0]*fSkin[3]-0.34445949507888407*nuEdge[0]*fSkin[3]+0.5966213466261491*nuSkin[1]*fEdge[3]-0.5966213466261491*nuEdge[1]*fEdge[3]+0.34445949507888407*nuSkin[0]*fEdge[3]+0.34445949507888407*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.33145630368119394*nuSkin[1]*fSkin[7])+0.33145630368119394*nuEdge[1]*fSkin[7]-0.19136638615493565*nuSkin[0]*fSkin[7]-0.19136638615493565*nuEdge[0]*fSkin[7]-0.33145630368119394*nuSkin[1]*fEdge[7]+0.33145630368119394*nuEdge[1]*fEdge[7]-0.19136638615493565*nuSkin[0]*fEdge[7]-0.19136638615493565*nuEdge[0]*fEdge[7]-0.34445949507888407*nuSkin[1]*fSkin[6]+0.34445949507888407*nuEdge[1]*fSkin[6]-0.19887378220871635*nuSkin[0]*fSkin[6]-0.19887378220871635*nuEdge[0]*fSkin[6]+0.34445949507888407*nuSkin[1]*fEdge[6]-0.34445949507888407*nuEdge[1]*fEdge[6]+0.19887378220871635*nuSkin[0]*fEdge[6]+0.19887378220871635*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[7] = -(0.5740991584648069*nuSkin[1]*fSkin[7])+0.5740991584648069*nuEdge[1]*fSkin[7]-0.33145630368119394*nuSkin[0]*fSkin[7]-0.33145630368119394*nuEdge[0]*fSkin[7]-0.5740991584648069*nuSkin[1]*fEdge[7]+0.5740991584648069*nuEdge[1]*fEdge[7]-0.33145630368119394*nuSkin[0]*fEdge[7]-0.33145630368119394*nuEdge[0]*fEdge[7]-0.5966213466261491*nuSkin[1]*fSkin[6]+0.5966213466261491*nuEdge[1]*fSkin[6]-0.34445949507888407*nuSkin[0]*fSkin[6]-0.34445949507888407*nuEdge[0]*fSkin[6]+0.5966213466261491*nuSkin[1]*fEdge[6]-0.5966213466261491*nuEdge[1]*fEdge[6]+0.34445949507888407*nuSkin[0]*fEdge[6]+0.34445949507888407*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[8] = -(0.3314563036811939*nuSkin[1]*fSkin[9])+0.3314563036811939*nuEdge[1]*fSkin[9]-0.19136638615493567*nuSkin[0]*fSkin[9]-0.19136638615493567*nuEdge[0]*fSkin[9]-0.3314563036811939*nuSkin[1]*fEdge[9]+0.3314563036811939*nuEdge[1]*fEdge[9]-0.19136638615493567*nuSkin[0]*fEdge[9]-0.19136638615493567*nuEdge[0]*fEdge[9]-0.34445949507888407*nuSkin[1]*fSkin[8]+0.34445949507888407*nuEdge[1]*fSkin[8]-0.19887378220871635*nuSkin[0]*fSkin[8]-0.19887378220871635*nuEdge[0]*fSkin[8]+0.34445949507888407*nuSkin[1]*fEdge[8]-0.34445949507888407*nuEdge[1]*fEdge[8]+0.19887378220871635*nuSkin[0]*fEdge[8]+0.19887378220871635*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[9] = -(0.5740991584648069*nuSkin[1]*fSkin[9])+0.5740991584648069*nuEdge[1]*fSkin[9]-0.33145630368119394*nuSkin[0]*fSkin[9]-0.33145630368119394*nuEdge[0]*fSkin[9]-0.5740991584648069*nuSkin[1]*fEdge[9]+0.5740991584648069*nuEdge[1]*fEdge[9]-0.33145630368119394*nuSkin[0]*fEdge[9]-0.33145630368119394*nuEdge[0]*fEdge[9]-0.5966213466261489*nuSkin[1]*fSkin[8]+0.5966213466261489*nuEdge[1]*fSkin[8]-0.3444594950788842*nuSkin[0]*fSkin[8]-0.3444594950788842*nuEdge[0]*fSkin[8]+0.5966213466261489*nuSkin[1]*fEdge[8]-0.5966213466261489*nuEdge[1]*fEdge[8]+0.3444594950788842*nuSkin[0]*fEdge[8]+0.3444594950788842*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[10] = -(0.3314563036811939*nuSkin[1]*fSkin[11])+0.3314563036811939*nuEdge[1]*fSkin[11]-0.19136638615493567*nuSkin[0]*fSkin[11]-0.19136638615493567*nuEdge[0]*fSkin[11]-0.3314563036811939*nuSkin[1]*fEdge[11]+0.3314563036811939*nuEdge[1]*fEdge[11]-0.19136638615493567*nuSkin[0]*fEdge[11]-0.19136638615493567*nuEdge[0]*fEdge[11]-0.34445949507888407*nuSkin[1]*fSkin[10]+0.34445949507888407*nuEdge[1]*fSkin[10]-0.19887378220871635*nuSkin[0]*fSkin[10]-0.19887378220871635*nuEdge[0]*fSkin[10]+0.34445949507888407*nuSkin[1]*fEdge[10]-0.34445949507888407*nuEdge[1]*fEdge[10]+0.19887378220871635*nuSkin[0]*fEdge[10]+0.19887378220871635*nuEdge[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.5740991584648069*nuSkin[1]*fSkin[11])+0.5740991584648069*nuEdge[1]*fSkin[11]-0.33145630368119394*nuSkin[0]*fSkin[11]-0.33145630368119394*nuEdge[0]*fSkin[11]-0.5740991584648069*nuSkin[1]*fEdge[11]+0.5740991584648069*nuEdge[1]*fEdge[11]-0.33145630368119394*nuSkin[0]*fEdge[11]-0.33145630368119394*nuEdge[0]*fEdge[11]-0.5966213466261489*nuSkin[1]*fSkin[10]+0.5966213466261489*nuEdge[1]*fSkin[10]-0.3444594950788842*nuSkin[0]*fSkin[10]-0.3444594950788842*nuEdge[0]*fSkin[10]+0.5966213466261489*nuSkin[1]*fEdge[10]-0.5966213466261489*nuEdge[1]*fEdge[10]+0.3444594950788842*nuSkin[0]*fEdge[10]+0.3444594950788842*nuEdge[0]*fEdge[10]; 

  boundSurf_incr[0] = 1.060660171779821*fSkin[1]*nuSkin[1]-0.6123724356957944*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[1] = 1.060660171779821*nuSkin[0]*fSkin[1]-1.837117307087383*fSkin[1]*nuSkin[1]; 
  boundSurf_incr[2] = 1.060660171779821*nuSkin[1]*fSkin[4]-0.6123724356957944*nuSkin[0]*fSkin[4]; 
  boundSurf_incr[3] = 1.060660171779821*nuSkin[1]*fSkin[5]-0.6123724356957944*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[4] = 1.060660171779821*nuSkin[0]*fSkin[4]-1.837117307087383*nuSkin[1]*fSkin[4]; 
  boundSurf_incr[5] = 1.060660171779821*nuSkin[0]*fSkin[5]-1.837117307087383*nuSkin[1]*fSkin[5]; 
  boundSurf_incr[6] = 1.060660171779821*nuSkin[1]*fSkin[7]-0.6123724356957944*nuSkin[0]*fSkin[7]; 
  boundSurf_incr[7] = 1.060660171779821*nuSkin[0]*fSkin[7]-1.837117307087383*nuSkin[1]*fSkin[7]; 
  boundSurf_incr[8] = 1.060660171779821*nuSkin[1]*fSkin[9]-0.6123724356957944*nuSkin[0]*fSkin[9]; 
  boundSurf_incr[9] = 1.060660171779821*nuSkin[0]*fSkin[9]-1.837117307087383*nuSkin[1]*fSkin[9]; 
  boundSurf_incr[10] = 1.060660171779821*nuSkin[1]*fSkin[11]-0.6123724356957944*nuSkin[0]*fSkin[11]; 
  boundSurf_incr[11] = 1.060660171779821*nuSkin[0]*fSkin[11]-1.837117307087383*nuSkin[1]*fSkin[11]; 

  out[0] += (edgeSurf_incr[0]+boundSurf_incr[0])*rdx2fac; 
  out[1] += (edgeSurf_incr[1]+boundSurf_incr[1])*rdx2fac; 
  out[2] += (edgeSurf_incr[2]+boundSurf_incr[2])*rdx2fac; 
  out[3] += (edgeSurf_incr[3]+boundSurf_incr[3])*rdx2fac; 
  out[4] += (edgeSurf_incr[4]+boundSurf_incr[4])*rdx2fac; 
  out[5] += (edgeSurf_incr[5]+boundSurf_incr[5])*rdx2fac; 
  out[6] += (edgeSurf_incr[6]+boundSurf_incr[6])*rdx2fac; 
  out[7] += (edgeSurf_incr[7]+boundSurf_incr[7])*rdx2fac; 
  out[8] += (edgeSurf_incr[8]+boundSurf_incr[8])*rdx2fac; 
  out[9] += (edgeSurf_incr[9]+boundSurf_incr[9])*rdx2fac; 
  out[10] += (edgeSurf_incr[10]+boundSurf_incr[10])*rdx2fac; 
  out[11] += (edgeSurf_incr[11]+boundSurf_incr[11])*rdx2fac; 

  return 0.;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_upper_zero_flux_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[12];
  fEdge[0] = 0.7071067811865475*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865475*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[4]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[5]+jacobgeo_invEdge[0]*JfEdge[3]); 
  fEdge[4] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[4]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[5] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[5]+jacobgeo_invEdge[1]*JfEdge[3]); 
  fEdge[6] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[7]+jacobgeo_invEdge[0]*JfEdge[6]); 
  fEdge[7] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[7]+jacobgeo_invEdge[1]*JfEdge[6]); 
  fEdge[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[9]+15.0*jacobgeo_invEdge[0]*JfEdge[8]); 
  fEdge[9] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[9]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[8]); 
  fEdge[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[11]+15.0*jacobgeo_invEdge[0]*JfEdge[10]); 
  fEdge[11] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[11]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[10]); 

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

  double edgeSurf_incr[12] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]+0.34445949507888407*fSkin[0]*nuSkin[1]-0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]-0.34445949507888407*fSkin[0]*nuEdge[1]+0.34445949507888407*fEdge[0]*nuEdge[1]+0.19136638615493565*nuSkin[0]*fSkin[1]+0.19136638615493565*nuEdge[0]*fSkin[1]+0.19136638615493565*nuSkin[0]*fEdge[1]+0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = 0.5740991584648069*fSkin[1]*nuSkin[1]+0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]-0.5740991584648069*fSkin[1]*nuEdge[1]-0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]+0.34445949507888407*fSkin[0]*nuSkin[0]-0.34445949507888407*fEdge[0]*nuSkin[0]+0.34445949507888407*fSkin[0]*nuEdge[0]-0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[4])+0.33145630368119394*nuEdge[1]*fSkin[4]+0.19136638615493565*nuSkin[0]*fSkin[4]+0.19136638615493565*nuEdge[0]*fSkin[4]-0.33145630368119394*nuSkin[1]*fEdge[4]+0.33145630368119394*nuEdge[1]*fEdge[4]+0.19136638615493565*nuSkin[0]*fEdge[4]+0.19136638615493565*nuEdge[0]*fEdge[4]+0.34445949507888407*nuSkin[1]*fSkin[2]-0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]-0.34445949507888407*nuSkin[1]*fEdge[2]+0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.33145630368119394*nuSkin[1]*fSkin[5])+0.33145630368119394*nuEdge[1]*fSkin[5]+0.19136638615493565*nuSkin[0]*fSkin[5]+0.19136638615493565*nuEdge[0]*fSkin[5]-0.33145630368119394*nuSkin[1]*fEdge[5]+0.33145630368119394*nuEdge[1]*fEdge[5]+0.19136638615493565*nuSkin[0]*fEdge[5]+0.19136638615493565*nuEdge[0]*fEdge[5]+0.34445949507888407*nuSkin[1]*fSkin[3]-0.34445949507888407*nuEdge[1]*fSkin[3]-0.19887378220871635*nuSkin[0]*fSkin[3]-0.19887378220871635*nuEdge[0]*fSkin[3]-0.34445949507888407*nuSkin[1]*fEdge[3]+0.34445949507888407*nuEdge[1]*fEdge[3]+0.19887378220871635*nuSkin[0]*fEdge[3]+0.19887378220871635*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.5740991584648069*nuSkin[1]*fSkin[4]-0.5740991584648069*nuEdge[1]*fSkin[4]-0.33145630368119394*nuSkin[0]*fSkin[4]-0.33145630368119394*nuEdge[0]*fSkin[4]+0.5740991584648069*nuSkin[1]*fEdge[4]-0.5740991584648069*nuEdge[1]*fEdge[4]-0.33145630368119394*nuSkin[0]*fEdge[4]-0.33145630368119394*nuEdge[0]*fEdge[4]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]+0.34445949507888407*nuSkin[0]*fSkin[2]+0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]-0.34445949507888407*nuSkin[0]*fEdge[2]-0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[5] = 0.5740991584648069*nuSkin[1]*fSkin[5]-0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]+0.5740991584648069*nuSkin[1]*fEdge[5]-0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261491*nuSkin[1]*fSkin[3]+0.5966213466261491*nuEdge[1]*fSkin[3]+0.34445949507888407*nuSkin[0]*fSkin[3]+0.34445949507888407*nuEdge[0]*fSkin[3]+0.5966213466261491*nuSkin[1]*fEdge[3]-0.5966213466261491*nuEdge[1]*fEdge[3]-0.34445949507888407*nuSkin[0]*fEdge[3]-0.34445949507888407*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.33145630368119394*nuSkin[1]*fSkin[7])+0.33145630368119394*nuEdge[1]*fSkin[7]+0.19136638615493565*nuSkin[0]*fSkin[7]+0.19136638615493565*nuEdge[0]*fSkin[7]-0.33145630368119394*nuSkin[1]*fEdge[7]+0.33145630368119394*nuEdge[1]*fEdge[7]+0.19136638615493565*nuSkin[0]*fEdge[7]+0.19136638615493565*nuEdge[0]*fEdge[7]+0.34445949507888407*nuSkin[1]*fSkin[6]-0.34445949507888407*nuEdge[1]*fSkin[6]-0.19887378220871635*nuSkin[0]*fSkin[6]-0.19887378220871635*nuEdge[0]*fSkin[6]-0.34445949507888407*nuSkin[1]*fEdge[6]+0.34445949507888407*nuEdge[1]*fEdge[6]+0.19887378220871635*nuSkin[0]*fEdge[6]+0.19887378220871635*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[7] = 0.5740991584648069*nuSkin[1]*fSkin[7]-0.5740991584648069*nuEdge[1]*fSkin[7]-0.33145630368119394*nuSkin[0]*fSkin[7]-0.33145630368119394*nuEdge[0]*fSkin[7]+0.5740991584648069*nuSkin[1]*fEdge[7]-0.5740991584648069*nuEdge[1]*fEdge[7]-0.33145630368119394*nuSkin[0]*fEdge[7]-0.33145630368119394*nuEdge[0]*fEdge[7]-0.5966213466261491*nuSkin[1]*fSkin[6]+0.5966213466261491*nuEdge[1]*fSkin[6]+0.34445949507888407*nuSkin[0]*fSkin[6]+0.34445949507888407*nuEdge[0]*fSkin[6]+0.5966213466261491*nuSkin[1]*fEdge[6]-0.5966213466261491*nuEdge[1]*fEdge[6]-0.34445949507888407*nuSkin[0]*fEdge[6]-0.34445949507888407*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[8] = -(0.3314563036811939*nuSkin[1]*fSkin[9])+0.3314563036811939*nuEdge[1]*fSkin[9]+0.19136638615493567*nuSkin[0]*fSkin[9]+0.19136638615493567*nuEdge[0]*fSkin[9]-0.3314563036811939*nuSkin[1]*fEdge[9]+0.3314563036811939*nuEdge[1]*fEdge[9]+0.19136638615493567*nuSkin[0]*fEdge[9]+0.19136638615493567*nuEdge[0]*fEdge[9]+0.34445949507888407*nuSkin[1]*fSkin[8]-0.34445949507888407*nuEdge[1]*fSkin[8]-0.19887378220871635*nuSkin[0]*fSkin[8]-0.19887378220871635*nuEdge[0]*fSkin[8]-0.34445949507888407*nuSkin[1]*fEdge[8]+0.34445949507888407*nuEdge[1]*fEdge[8]+0.19887378220871635*nuSkin[0]*fEdge[8]+0.19887378220871635*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[9] = 0.5740991584648069*nuSkin[1]*fSkin[9]-0.5740991584648069*nuEdge[1]*fSkin[9]-0.33145630368119394*nuSkin[0]*fSkin[9]-0.33145630368119394*nuEdge[0]*fSkin[9]+0.5740991584648069*nuSkin[1]*fEdge[9]-0.5740991584648069*nuEdge[1]*fEdge[9]-0.33145630368119394*nuSkin[0]*fEdge[9]-0.33145630368119394*nuEdge[0]*fEdge[9]-0.5966213466261489*nuSkin[1]*fSkin[8]+0.5966213466261489*nuEdge[1]*fSkin[8]+0.3444594950788842*nuSkin[0]*fSkin[8]+0.3444594950788842*nuEdge[0]*fSkin[8]+0.5966213466261489*nuSkin[1]*fEdge[8]-0.5966213466261489*nuEdge[1]*fEdge[8]-0.3444594950788842*nuSkin[0]*fEdge[8]-0.3444594950788842*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[10] = -(0.3314563036811939*nuSkin[1]*fSkin[11])+0.3314563036811939*nuEdge[1]*fSkin[11]+0.19136638615493567*nuSkin[0]*fSkin[11]+0.19136638615493567*nuEdge[0]*fSkin[11]-0.3314563036811939*nuSkin[1]*fEdge[11]+0.3314563036811939*nuEdge[1]*fEdge[11]+0.19136638615493567*nuSkin[0]*fEdge[11]+0.19136638615493567*nuEdge[0]*fEdge[11]+0.34445949507888407*nuSkin[1]*fSkin[10]-0.34445949507888407*nuEdge[1]*fSkin[10]-0.19887378220871635*nuSkin[0]*fSkin[10]-0.19887378220871635*nuEdge[0]*fSkin[10]-0.34445949507888407*nuSkin[1]*fEdge[10]+0.34445949507888407*nuEdge[1]*fEdge[10]+0.19887378220871635*nuSkin[0]*fEdge[10]+0.19887378220871635*nuEdge[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.5740991584648069*nuSkin[1]*fSkin[11]-0.5740991584648069*nuEdge[1]*fSkin[11]-0.33145630368119394*nuSkin[0]*fSkin[11]-0.33145630368119394*nuEdge[0]*fSkin[11]+0.5740991584648069*nuSkin[1]*fEdge[11]-0.5740991584648069*nuEdge[1]*fEdge[11]-0.33145630368119394*nuSkin[0]*fEdge[11]-0.33145630368119394*nuEdge[0]*fEdge[11]-0.5966213466261489*nuSkin[1]*fSkin[10]+0.5966213466261489*nuEdge[1]*fSkin[10]+0.3444594950788842*nuSkin[0]*fSkin[10]+0.3444594950788842*nuEdge[0]*fSkin[10]+0.5966213466261489*nuSkin[1]*fEdge[10]-0.5966213466261489*nuEdge[1]*fEdge[10]-0.3444594950788842*nuSkin[0]*fEdge[10]-0.3444594950788842*nuEdge[0]*fEdge[10]; 


  out[0] += edgeSurf_incr[0]*rdx2fac; 
  out[1] += edgeSurf_incr[1]*rdx2fac; 
  out[2] += edgeSurf_incr[2]*rdx2fac; 
  out[3] += edgeSurf_incr[3]*rdx2fac; 
  out[4] += edgeSurf_incr[4]*rdx2fac; 
  out[5] += edgeSurf_incr[5]*rdx2fac; 
  out[6] += edgeSurf_incr[6]*rdx2fac; 
  out[7] += edgeSurf_incr[7]*rdx2fac; 
  out[8] += edgeSurf_incr[8]*rdx2fac; 
  out[9] += edgeSurf_incr[9]*rdx2fac; 
  out[10] += edgeSurf_incr[10]*rdx2fac; 
  out[11] += edgeSurf_incr[11]*rdx2fac; 

  return 0.;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_boundary_surfx_upper_bound_local_1x2v_ser_p1_varnu(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuEdge/nuSkin: Diffusivity in edge and skin cells.
  // jacobgeo_invEdge/jacobgeo_invSkin: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfEdge/JfSkin: distribution times conf-space Jacobian in egde and skin cells.
  // out: Incremented output.

  const double rdx2fac = pow(2./dxSkin[0],2);

  double fEdge[12];
  fEdge[0] = 0.7071067811865475*(JfEdge[1]*jacobgeo_invEdge[1]+JfEdge[0]*jacobgeo_invEdge[0]); 
  fEdge[1] = 0.7071067811865475*(JfEdge[0]*jacobgeo_invEdge[1]+jacobgeo_invEdge[0]*JfEdge[1]); 
  fEdge[2] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[4]+jacobgeo_invEdge[0]*JfEdge[2]); 
  fEdge[3] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[5]+jacobgeo_invEdge[0]*JfEdge[3]); 
  fEdge[4] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[4]+jacobgeo_invEdge[1]*JfEdge[2]); 
  fEdge[5] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[5]+jacobgeo_invEdge[1]*JfEdge[3]); 
  fEdge[6] = 0.7071067811865475*(jacobgeo_invEdge[1]*JfEdge[7]+jacobgeo_invEdge[0]*JfEdge[6]); 
  fEdge[7] = 0.7071067811865475*(jacobgeo_invEdge[0]*JfEdge[7]+jacobgeo_invEdge[1]*JfEdge[6]); 
  fEdge[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[9]+15.0*jacobgeo_invEdge[0]*JfEdge[8]); 
  fEdge[9] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[9]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[8]); 
  fEdge[10] = 0.04714045207910316*(15.000000000000002*jacobgeo_invEdge[1]*JfEdge[11]+15.0*jacobgeo_invEdge[0]*JfEdge[10]); 
  fEdge[11] = 0.04714045207910316*(15.0*jacobgeo_invEdge[0]*JfEdge[11]+15.000000000000002*jacobgeo_invEdge[1]*JfEdge[10]); 

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

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  edgeSurf_incr[0] = -(0.33145630368119394*fSkin[1]*nuSkin[1])-0.33145630368119394*fEdge[1]*nuSkin[1]+0.34445949507888407*fSkin[0]*nuSkin[1]-0.34445949507888407*fEdge[0]*nuSkin[1]+0.33145630368119394*fSkin[1]*nuEdge[1]+0.33145630368119394*fEdge[1]*nuEdge[1]-0.34445949507888407*fSkin[0]*nuEdge[1]+0.34445949507888407*fEdge[0]*nuEdge[1]+0.19136638615493565*nuSkin[0]*fSkin[1]+0.19136638615493565*nuEdge[0]*fSkin[1]+0.19136638615493565*nuSkin[0]*fEdge[1]+0.19136638615493565*nuEdge[0]*fEdge[1]-0.19887378220871635*fSkin[0]*nuSkin[0]+0.19887378220871635*fEdge[0]*nuSkin[0]-0.19887378220871635*fSkin[0]*nuEdge[0]+0.19887378220871635*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[1] = 0.5740991584648069*fSkin[1]*nuSkin[1]+0.5740991584648069*fEdge[1]*nuSkin[1]-0.5966213466261491*fSkin[0]*nuSkin[1]+0.5966213466261491*fEdge[0]*nuSkin[1]-0.5740991584648069*fSkin[1]*nuEdge[1]-0.5740991584648069*fEdge[1]*nuEdge[1]+0.5966213466261491*fSkin[0]*nuEdge[1]-0.5966213466261491*fEdge[0]*nuEdge[1]-0.33145630368119394*nuSkin[0]*fSkin[1]-0.33145630368119394*nuEdge[0]*fSkin[1]-0.33145630368119394*nuSkin[0]*fEdge[1]-0.33145630368119394*nuEdge[0]*fEdge[1]+0.34445949507888407*fSkin[0]*nuSkin[0]-0.34445949507888407*fEdge[0]*nuSkin[0]+0.34445949507888407*fSkin[0]*nuEdge[0]-0.34445949507888407*fEdge[0]*nuEdge[0]; 
  edgeSurf_incr[2] = -(0.33145630368119394*nuSkin[1]*fSkin[4])+0.33145630368119394*nuEdge[1]*fSkin[4]+0.19136638615493565*nuSkin[0]*fSkin[4]+0.19136638615493565*nuEdge[0]*fSkin[4]-0.33145630368119394*nuSkin[1]*fEdge[4]+0.33145630368119394*nuEdge[1]*fEdge[4]+0.19136638615493565*nuSkin[0]*fEdge[4]+0.19136638615493565*nuEdge[0]*fEdge[4]+0.34445949507888407*nuSkin[1]*fSkin[2]-0.34445949507888407*nuEdge[1]*fSkin[2]-0.19887378220871635*nuSkin[0]*fSkin[2]-0.19887378220871635*nuEdge[0]*fSkin[2]-0.34445949507888407*nuSkin[1]*fEdge[2]+0.34445949507888407*nuEdge[1]*fEdge[2]+0.19887378220871635*nuSkin[0]*fEdge[2]+0.19887378220871635*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.33145630368119394*nuSkin[1]*fSkin[5])+0.33145630368119394*nuEdge[1]*fSkin[5]+0.19136638615493565*nuSkin[0]*fSkin[5]+0.19136638615493565*nuEdge[0]*fSkin[5]-0.33145630368119394*nuSkin[1]*fEdge[5]+0.33145630368119394*nuEdge[1]*fEdge[5]+0.19136638615493565*nuSkin[0]*fEdge[5]+0.19136638615493565*nuEdge[0]*fEdge[5]+0.34445949507888407*nuSkin[1]*fSkin[3]-0.34445949507888407*nuEdge[1]*fSkin[3]-0.19887378220871635*nuSkin[0]*fSkin[3]-0.19887378220871635*nuEdge[0]*fSkin[3]-0.34445949507888407*nuSkin[1]*fEdge[3]+0.34445949507888407*nuEdge[1]*fEdge[3]+0.19887378220871635*nuSkin[0]*fEdge[3]+0.19887378220871635*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.5740991584648069*nuSkin[1]*fSkin[4]-0.5740991584648069*nuEdge[1]*fSkin[4]-0.33145630368119394*nuSkin[0]*fSkin[4]-0.33145630368119394*nuEdge[0]*fSkin[4]+0.5740991584648069*nuSkin[1]*fEdge[4]-0.5740991584648069*nuEdge[1]*fEdge[4]-0.33145630368119394*nuSkin[0]*fEdge[4]-0.33145630368119394*nuEdge[0]*fEdge[4]-0.5966213466261491*nuSkin[1]*fSkin[2]+0.5966213466261491*nuEdge[1]*fSkin[2]+0.34445949507888407*nuSkin[0]*fSkin[2]+0.34445949507888407*nuEdge[0]*fSkin[2]+0.5966213466261491*nuSkin[1]*fEdge[2]-0.5966213466261491*nuEdge[1]*fEdge[2]-0.34445949507888407*nuSkin[0]*fEdge[2]-0.34445949507888407*nuEdge[0]*fEdge[2]; 
  edgeSurf_incr[5] = 0.5740991584648069*nuSkin[1]*fSkin[5]-0.5740991584648069*nuEdge[1]*fSkin[5]-0.33145630368119394*nuSkin[0]*fSkin[5]-0.33145630368119394*nuEdge[0]*fSkin[5]+0.5740991584648069*nuSkin[1]*fEdge[5]-0.5740991584648069*nuEdge[1]*fEdge[5]-0.33145630368119394*nuSkin[0]*fEdge[5]-0.33145630368119394*nuEdge[0]*fEdge[5]-0.5966213466261491*nuSkin[1]*fSkin[3]+0.5966213466261491*nuEdge[1]*fSkin[3]+0.34445949507888407*nuSkin[0]*fSkin[3]+0.34445949507888407*nuEdge[0]*fSkin[3]+0.5966213466261491*nuSkin[1]*fEdge[3]-0.5966213466261491*nuEdge[1]*fEdge[3]-0.34445949507888407*nuSkin[0]*fEdge[3]-0.34445949507888407*nuEdge[0]*fEdge[3]; 
  edgeSurf_incr[6] = -(0.33145630368119394*nuSkin[1]*fSkin[7])+0.33145630368119394*nuEdge[1]*fSkin[7]+0.19136638615493565*nuSkin[0]*fSkin[7]+0.19136638615493565*nuEdge[0]*fSkin[7]-0.33145630368119394*nuSkin[1]*fEdge[7]+0.33145630368119394*nuEdge[1]*fEdge[7]+0.19136638615493565*nuSkin[0]*fEdge[7]+0.19136638615493565*nuEdge[0]*fEdge[7]+0.34445949507888407*nuSkin[1]*fSkin[6]-0.34445949507888407*nuEdge[1]*fSkin[6]-0.19887378220871635*nuSkin[0]*fSkin[6]-0.19887378220871635*nuEdge[0]*fSkin[6]-0.34445949507888407*nuSkin[1]*fEdge[6]+0.34445949507888407*nuEdge[1]*fEdge[6]+0.19887378220871635*nuSkin[0]*fEdge[6]+0.19887378220871635*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[7] = 0.5740991584648069*nuSkin[1]*fSkin[7]-0.5740991584648069*nuEdge[1]*fSkin[7]-0.33145630368119394*nuSkin[0]*fSkin[7]-0.33145630368119394*nuEdge[0]*fSkin[7]+0.5740991584648069*nuSkin[1]*fEdge[7]-0.5740991584648069*nuEdge[1]*fEdge[7]-0.33145630368119394*nuSkin[0]*fEdge[7]-0.33145630368119394*nuEdge[0]*fEdge[7]-0.5966213466261491*nuSkin[1]*fSkin[6]+0.5966213466261491*nuEdge[1]*fSkin[6]+0.34445949507888407*nuSkin[0]*fSkin[6]+0.34445949507888407*nuEdge[0]*fSkin[6]+0.5966213466261491*nuSkin[1]*fEdge[6]-0.5966213466261491*nuEdge[1]*fEdge[6]-0.34445949507888407*nuSkin[0]*fEdge[6]-0.34445949507888407*nuEdge[0]*fEdge[6]; 
  edgeSurf_incr[8] = -(0.3314563036811939*nuSkin[1]*fSkin[9])+0.3314563036811939*nuEdge[1]*fSkin[9]+0.19136638615493567*nuSkin[0]*fSkin[9]+0.19136638615493567*nuEdge[0]*fSkin[9]-0.3314563036811939*nuSkin[1]*fEdge[9]+0.3314563036811939*nuEdge[1]*fEdge[9]+0.19136638615493567*nuSkin[0]*fEdge[9]+0.19136638615493567*nuEdge[0]*fEdge[9]+0.34445949507888407*nuSkin[1]*fSkin[8]-0.34445949507888407*nuEdge[1]*fSkin[8]-0.19887378220871635*nuSkin[0]*fSkin[8]-0.19887378220871635*nuEdge[0]*fSkin[8]-0.34445949507888407*nuSkin[1]*fEdge[8]+0.34445949507888407*nuEdge[1]*fEdge[8]+0.19887378220871635*nuSkin[0]*fEdge[8]+0.19887378220871635*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[9] = 0.5740991584648069*nuSkin[1]*fSkin[9]-0.5740991584648069*nuEdge[1]*fSkin[9]-0.33145630368119394*nuSkin[0]*fSkin[9]-0.33145630368119394*nuEdge[0]*fSkin[9]+0.5740991584648069*nuSkin[1]*fEdge[9]-0.5740991584648069*nuEdge[1]*fEdge[9]-0.33145630368119394*nuSkin[0]*fEdge[9]-0.33145630368119394*nuEdge[0]*fEdge[9]-0.5966213466261489*nuSkin[1]*fSkin[8]+0.5966213466261489*nuEdge[1]*fSkin[8]+0.3444594950788842*nuSkin[0]*fSkin[8]+0.3444594950788842*nuEdge[0]*fSkin[8]+0.5966213466261489*nuSkin[1]*fEdge[8]-0.5966213466261489*nuEdge[1]*fEdge[8]-0.3444594950788842*nuSkin[0]*fEdge[8]-0.3444594950788842*nuEdge[0]*fEdge[8]; 
  edgeSurf_incr[10] = -(0.3314563036811939*nuSkin[1]*fSkin[11])+0.3314563036811939*nuEdge[1]*fSkin[11]+0.19136638615493567*nuSkin[0]*fSkin[11]+0.19136638615493567*nuEdge[0]*fSkin[11]-0.3314563036811939*nuSkin[1]*fEdge[11]+0.3314563036811939*nuEdge[1]*fEdge[11]+0.19136638615493567*nuSkin[0]*fEdge[11]+0.19136638615493567*nuEdge[0]*fEdge[11]+0.34445949507888407*nuSkin[1]*fSkin[10]-0.34445949507888407*nuEdge[1]*fSkin[10]-0.19887378220871635*nuSkin[0]*fSkin[10]-0.19887378220871635*nuEdge[0]*fSkin[10]-0.34445949507888407*nuSkin[1]*fEdge[10]+0.34445949507888407*nuEdge[1]*fEdge[10]+0.19887378220871635*nuSkin[0]*fEdge[10]+0.19887378220871635*nuEdge[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.5740991584648069*nuSkin[1]*fSkin[11]-0.5740991584648069*nuEdge[1]*fSkin[11]-0.33145630368119394*nuSkin[0]*fSkin[11]-0.33145630368119394*nuEdge[0]*fSkin[11]+0.5740991584648069*nuSkin[1]*fEdge[11]-0.5740991584648069*nuEdge[1]*fEdge[11]-0.33145630368119394*nuSkin[0]*fEdge[11]-0.33145630368119394*nuEdge[0]*fEdge[11]-0.5966213466261489*nuSkin[1]*fSkin[10]+0.5966213466261489*nuEdge[1]*fSkin[10]+0.3444594950788842*nuSkin[0]*fSkin[10]+0.3444594950788842*nuEdge[0]*fSkin[10]+0.5966213466261489*nuSkin[1]*fEdge[10]-0.5966213466261489*nuEdge[1]*fEdge[10]-0.3444594950788842*nuSkin[0]*fEdge[10]-0.3444594950788842*nuEdge[0]*fEdge[10]; 

  boundSurf_incr[0] = 1.060660171779821*fSkin[1]*nuSkin[1]+0.6123724356957944*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[1] = 1.837117307087383*fSkin[1]*nuSkin[1]+1.060660171779821*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[2] = 1.060660171779821*nuSkin[1]*fSkin[4]+0.6123724356957944*nuSkin[0]*fSkin[4]; 
  boundSurf_incr[3] = 1.060660171779821*nuSkin[1]*fSkin[5]+0.6123724356957944*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[4] = 1.837117307087383*nuSkin[1]*fSkin[4]+1.060660171779821*nuSkin[0]*fSkin[4]; 
  boundSurf_incr[5] = 1.837117307087383*nuSkin[1]*fSkin[5]+1.060660171779821*nuSkin[0]*fSkin[5]; 
  boundSurf_incr[6] = 1.060660171779821*nuSkin[1]*fSkin[7]+0.6123724356957944*nuSkin[0]*fSkin[7]; 
  boundSurf_incr[7] = 1.837117307087383*nuSkin[1]*fSkin[7]+1.060660171779821*nuSkin[0]*fSkin[7]; 
  boundSurf_incr[8] = 1.060660171779821*nuSkin[1]*fSkin[9]+0.6123724356957944*nuSkin[0]*fSkin[9]; 
  boundSurf_incr[9] = 1.837117307087383*nuSkin[1]*fSkin[9]+1.060660171779821*nuSkin[0]*fSkin[9]; 
  boundSurf_incr[10] = 1.060660171779821*nuSkin[1]*fSkin[11]+0.6123724356957944*nuSkin[0]*fSkin[11]; 
  boundSurf_incr[11] = 1.837117307087383*nuSkin[1]*fSkin[11]+1.060660171779821*nuSkin[0]*fSkin[11]; 

  out[0] += (edgeSurf_incr[0]+boundSurf_incr[0])*rdx2fac; 
  out[1] += (edgeSurf_incr[1]+boundSurf_incr[1])*rdx2fac; 
  out[2] += (edgeSurf_incr[2]+boundSurf_incr[2])*rdx2fac; 
  out[3] += (edgeSurf_incr[3]+boundSurf_incr[3])*rdx2fac; 
  out[4] += (edgeSurf_incr[4]+boundSurf_incr[4])*rdx2fac; 
  out[5] += (edgeSurf_incr[5]+boundSurf_incr[5])*rdx2fac; 
  out[6] += (edgeSurf_incr[6]+boundSurf_incr[6])*rdx2fac; 
  out[7] += (edgeSurf_incr[7]+boundSurf_incr[7])*rdx2fac; 
  out[8] += (edgeSurf_incr[8]+boundSurf_incr[8])*rdx2fac; 
  out[9] += (edgeSurf_incr[9]+boundSurf_incr[9])*rdx2fac; 
  out[10] += (edgeSurf_incr[10]+boundSurf_incr[10])*rdx2fac; 
  out[11] += (edgeSurf_incr[11]+boundSurf_incr[11])*rdx2fac; 

  return 0.;
}

