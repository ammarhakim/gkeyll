#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x2v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // field:       q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E0 = &field[0]; 
  const double *B2 = &field[10]; 

  double alpha[6] = {0.0}; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 

  double fUpwindQuad[6] = {0.0};
  double fUpwind[6] = {0.0};
  double Ghat[6] = {0.0}; 

  if (edge == -1) { 

  if (0.6708203932499369*alpha[3]-0.6708203932499369*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.5*alpha[0]-0.5*alpha[1] > 0) { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if ((-0.6708203932499369*alpha[3])+0.6708203932499369*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_l(fEdge); 
  } 
  if (0.5*(alpha[1]+alpha[0])-0.6708203932499369*(alpha[3]+alpha[2]) > 0) { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_l(fEdge); 
  } 
  if (0.5*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_l(fEdge); 
  } 
  if (0.6708203932499369*(alpha[3]+alpha[2])+0.5*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.447213595499958*alpha[3]*fUpwind[5]+0.4472135954999579*alpha[2]*fUpwind[4]+0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.447213595499958*alpha[2]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[4]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.5000000000000001*alpha[1]*fUpwind[5]+0.5*alpha[0]*fUpwind[4]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[5] = 0.5*alpha[0]*fUpwind[5]+0.5000000000000001*alpha[1]*fUpwind[4]+0.447213595499958*alpha[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alpha[3]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += -0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 
  out[8] += -1.58113883008419*Ghat[0]*dv10; 
  out[9] += -1.58113883008419*Ghat[1]*dv10; 
  out[10] += -1.58113883008419*Ghat[2]*dv10; 
  out[11] += -1.58113883008419*Ghat[3]*dv10; 
  out[12] += -0.7071067811865475*Ghat[4]*dv10; 
  out[13] += -0.7071067811865475*Ghat[5]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += -1.224744871391589*Ghat[5]*dv10; 

  } else { 

  if (0.6708203932499369*alpha[3]-0.6708203932499369*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.5*alpha[0]-0.5*alpha[1] > 0) { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if ((-0.6708203932499369*alpha[3])+0.6708203932499369*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_l(fSkin); 
  } 
  if (0.5*(alpha[1]+alpha[0])-0.6708203932499369*(alpha[3]+alpha[2]) > 0) { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_l(fSkin); 
  } 
  if (0.5*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_l(fSkin); 
  } 
  if (0.6708203932499369*(alpha[3]+alpha[2])+0.5*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.447213595499958*alpha[3]*fUpwind[5]+0.4472135954999579*alpha[2]*fUpwind[4]+0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.447213595499958*alpha[2]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[4]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.5000000000000001*alpha[1]*fUpwind[5]+0.5*alpha[0]*fUpwind[4]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[5] = 0.5*alpha[0]*fUpwind[5]+0.5000000000000001*alpha[1]*fUpwind[4]+0.447213595499958*alpha[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alpha[3]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 
  out[8] += 1.58113883008419*Ghat[0]*dv10; 
  out[9] += 1.58113883008419*Ghat[1]*dv10; 
  out[10] += 1.58113883008419*Ghat[2]*dv10; 
  out[11] += 1.58113883008419*Ghat[3]*dv10; 
  out[12] += 0.7071067811865475*Ghat[4]*dv10; 
  out[13] += 0.7071067811865475*Ghat[5]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += -1.224744871391589*Ghat[5]*dv10; 

  } 
  return 0.;

} 
