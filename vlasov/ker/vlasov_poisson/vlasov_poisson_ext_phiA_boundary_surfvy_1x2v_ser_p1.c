#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_1x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // pots:        potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:       external E and B fields (scaled by q/m).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 

  const double dv11 = 2/dxv[2]; 
  const double dx10 = 2/dxv[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  double alpha[6] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ax = &pots[2]; 
  const double *Ay = &pots[4]; 

  alpha[0] = -(2.4494897427831783*Ay[1]*dx10*wv1); 
  alpha[2] = -(0.7071067811865475*Ay[1]*dv1*dx10); 

  double fUpwindQuad[6] = {0.0};
  double fUpwind[6] = {0.0};
  double Ghat[6] = {0.0}; 

  if (edge == -1) { 

  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if (0.5*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[2]*fUpwind[2]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[2]*fUpwind[3]+0.5*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.4472135954999579*alpha[2]*fUpwind[4]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.44721359549995804*alpha[2]*fUpwind[5]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.5*alpha[0]*fUpwind[4]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[5] = 0.5*alpha[0]*fUpwind[5]+0.44721359549995804*alpha[2]*fUpwind[3]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv11); 
  out[3] += -(1.224744871391589*Ghat[0]*dv11); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv11); 
  out[5] += -(1.224744871391589*Ghat[1]*dv11); 
  out[6] += -(1.224744871391589*Ghat[2]*dv11); 
  out[7] += -(1.224744871391589*Ghat[3]*dv11); 
  out[8] += -(0.7071067811865475*Ghat[4]*dv11); 
  out[9] += -(0.7071067811865475*Ghat[5]*dv11); 
  out[10] += -(1.224744871391589*Ghat[4]*dv11); 
  out[11] += -(1.224744871391589*Ghat[5]*dv11); 
  out[12] += -(1.5811388300841895*Ghat[0]*dv11); 
  out[13] += -(1.5811388300841898*Ghat[1]*dv11); 
  out[14] += -(1.5811388300841898*Ghat[2]*dv11); 
  out[15] += -(1.5811388300841895*Ghat[3]*dv11); 

  } else { 

  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if (0.5*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[2]*fUpwind[2]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[2]*fUpwind[3]+0.5*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.4472135954999579*alpha[2]*fUpwind[4]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.44721359549995804*alpha[2]*fUpwind[5]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.5*alpha[0]*fUpwind[4]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[5] = 0.5*alpha[0]*fUpwind[5]+0.44721359549995804*alpha[2]*fUpwind[3]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -(1.224744871391589*Ghat[0]*dv11); 
  out[4] += 0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -(1.224744871391589*Ghat[1]*dv11); 
  out[6] += -(1.224744871391589*Ghat[2]*dv11); 
  out[7] += -(1.224744871391589*Ghat[3]*dv11); 
  out[8] += 0.7071067811865475*Ghat[4]*dv11; 
  out[9] += 0.7071067811865475*Ghat[5]*dv11; 
  out[10] += -(1.224744871391589*Ghat[4]*dv11); 
  out[11] += -(1.224744871391589*Ghat[5]*dv11); 
  out[12] += 1.5811388300841895*Ghat[0]*dv11; 
  out[13] += 1.5811388300841898*Ghat[1]*dv11; 
  out[14] += 1.5811388300841898*Ghat[2]*dv11; 
  out[15] += 1.5811388300841895*Ghat[3]*dv11; 

  } 
  return 0.;

} 
