#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2.0/dxv[1]; 
  const double *E0 = &qmem[0]; 
  double alpha_l[2] = {0.0}; 
  double alpha_r[2] = {0.0}; 

  alpha_l[0] = E0[0]; 
  alpha_l[1] = E0[1]; 

  alpha_r[0] = E0[0]; 
  alpha_r[1] = E0[1]; 

  double fUpwindQuad_l[2] = {0.0};
  double fUpwindQuad_r[2] = {0.0};
  double fUpwind_l[2] = {0.0};;
  double fUpwind_r[2] = {0.0};
  double Ghat_l[2] = {0.0}; 
  double Ghat_r[2] = {0.0}; 

  if (0.7071067811865475*alpha_l[0]-0.7071067811865475*alpha_l[1] > 0) { 
    fUpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fc); 
  } 
  if (0.7071067811865475*alpha_r[0]-0.7071067811865475*alpha_r[1] > 0) { 
    fUpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.7071067811865475*(alpha_l[1]+alpha_l[0]) > 0) { 
    fUpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fc); 
  } 
  if (0.7071067811865475*(alpha_r[1]+alpha_r[0]) > 0) { 
    fUpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 
  Ghat_l[0] = 0.7071067811865475*alpha_l[1]*fUpwind_l[1]+0.7071067811865475*alpha_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.7071067811865475*alpha_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alpha_l[1]; 

  Ghat_r[0] = 0.7071067811865475*alpha_r[1]*fUpwind_r[1]+0.7071067811865475*alpha_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.7071067811865475*alpha_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alpha_r[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[5] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 

  return 0.;

} 
