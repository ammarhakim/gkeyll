#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_1x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 

  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dx10 = 2/dxv[0]; 

  double alpha[6] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ax = &pots[2]; 
  const double *Ay = &pots[4]; 

  alpha[0] = -(2.4494897427831783*Ay[1]*dx10*wv1); 
  alpha[2] = -(0.7071067811865475*Ay[1]*dv1*dx10); 

  double fUpwindQuad_l[6] = {0.0};
  double fUpwindQuad_r[6] = {0.0};
  double fUpwind_l[6] = {0.0};;
  double fUpwind_r[6] = {0.0};
  double Ghat_l[6] = {0.0}; 
  double Ghat_r[6] = {0.0}; 

  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 
    fUpwindQuad_l[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (0.5*alpha[0] > 0) { 
    fUpwindQuad_l[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fr); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 
    fUpwindQuad_l[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.5*alpha[0] > 0) { 
    fUpwindQuad_l[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alpha[2]*fUpwind_l[2]+0.5*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.5*alpha[2]*fUpwind_l[3]+0.5*alpha[0]*fUpwind_l[1]; 
  Ghat_l[2] = 0.4472135954999579*alpha[2]*fUpwind_l[4]+0.5*alpha[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.44721359549995804*alpha[2]*fUpwind_l[5]+0.5*alpha[0]*fUpwind_l[3]+0.5*fUpwind_l[1]*alpha[2]; 
  Ghat_l[4] = 0.5*alpha[0]*fUpwind_l[4]+0.4472135954999579*alpha[2]*fUpwind_l[2]; 
  Ghat_l[5] = 0.5*alpha[0]*fUpwind_l[5]+0.44721359549995804*alpha[2]*fUpwind_l[3]; 

  Ghat_r[0] = 0.5*alpha[2]*fUpwind_r[2]+0.5*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.5*alpha[2]*fUpwind_r[3]+0.5*alpha[0]*fUpwind_r[1]; 
  Ghat_r[2] = 0.4472135954999579*alpha[2]*fUpwind_r[4]+0.5*alpha[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.44721359549995804*alpha[2]*fUpwind_r[5]+0.5*alpha[0]*fUpwind_r[3]+0.5*fUpwind_r[1]*alpha[2]; 
  Ghat_r[4] = 0.5*alpha[0]*fUpwind_r[4]+0.4472135954999579*alpha[2]*fUpwind_r[2]; 
  Ghat_r[5] = 0.5*alpha[0]*fUpwind_r[5]+0.44721359549995804*alpha[2]*fUpwind_r[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11); 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11); 
  out[6] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11); 
  out[7] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11); 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[10] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11); 
  out[11] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11); 
  out[12] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv11; 
  out[13] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv11; 
  out[14] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv11; 
  out[15] += (1.5811388300841895*Ghat_l[3]-1.5811388300841895*Ghat_r[3])*dv11; 

  return 0.;

} 
