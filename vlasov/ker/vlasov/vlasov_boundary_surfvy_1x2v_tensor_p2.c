#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // field:       q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E1 = &field[3]; 
  const double *B2 = &field[15]; 

  double alpha[9] = {0.0}; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 
  alpha[4] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 

  double fUpwindQuad[9] = {0.0};
  double fUpwind[9] = {0.0};
  double Ghat[9] = {0.0}; 

  if (edge == -1) { 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if (0.5*alpha[0]-0.5590169943749475*alpha[4] > 0) { 
    fUpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(fEdge); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[6]*fUpwind[6]+0.5*alpha[4]*fUpwind[4]+0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.447213595499958*alpha[3]*fUpwind[6]+0.447213595499958*fUpwind[3]*alpha[6]+0.4472135954999579*alpha[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alpha[4]+0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.447213595499958*alpha[6]*fUpwind[8]+0.447213595499958*alpha[3]*fUpwind[7]+0.5000000000000001*alpha[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alpha[6]+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.4*alpha[3]*fUpwind[8]+0.4*alpha[6]*fUpwind[7]+0.447213595499958*alpha[2]*fUpwind[7]+0.447213595499958*alpha[1]*fUpwind[6]+0.447213595499958*fUpwind[1]*alpha[6]+0.4472135954999579*alpha[3]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alpha[4]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.31943828249997*alpha[6]*fUpwind[6]+0.5000000000000001*alpha[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alpha[6]+0.31943828249997*alpha[4]*fUpwind[4]+0.5*alpha[0]*fUpwind[4]+0.5*fUpwind[0]*alpha[4]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[1]*fUpwind[1]; 
  Ghat[5] = 0.5*alpha[4]*fUpwind[8]+0.5000000000000001*alpha[1]*fUpwind[7]+0.4472135954999579*alpha[6]*fUpwind[6]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[6] = 0.2857142857142857*alpha[6]*fUpwind[8]+0.447213595499958*alpha[2]*fUpwind[8]+0.4*alpha[3]*fUpwind[7]+0.31943828249997*alpha[4]*fUpwind[6]+0.5*alpha[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alpha[6]+0.31943828249997*fUpwind[4]*alpha[6]+0.5*fUpwind[0]*alpha[6]+0.5000000000000001*alpha[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alpha[4]+0.447213595499958*alpha[1]*fUpwind[3]+0.447213595499958*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.447213595499958*alpha[1]*fUpwind[8]+0.4472135954999579*alpha[4]*fUpwind[7]+0.5*alpha[0]*fUpwind[7]+0.4*alpha[3]*fUpwind[6]+0.4*fUpwind[3]*alpha[6]+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*alpha[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.31943828249997*alpha[4]*fUpwind[8]+0.5*alpha[0]*fUpwind[8]+0.447213595499958*alpha[1]*fUpwind[7]+0.2857142857142857*alpha[6]*fUpwind[6]+0.447213595499958*alpha[2]*fUpwind[6]+0.447213595499958*fUpwind[2]*alpha[6]+0.5*alpha[4]*fUpwind[5]+0.4*alpha[3]*fUpwind[3]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += -0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += -0.7071067811865475*Ghat[4]*dv11; 
  out[8] += -0.7071067811865475*Ghat[5]*dv11; 
  out[9] += -1.58113883008419*Ghat[0]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += -0.7071067811865475*Ghat[6]*dv11; 
  out[12] += -0.7071067811865475*Ghat[7]*dv11; 
  out[13] += -1.224744871391589*Ghat[4]*dv11; 
  out[14] += -1.224744871391589*Ghat[5]*dv11; 
  out[15] += -1.58113883008419*Ghat[1]*dv11; 
  out[16] += -1.58113883008419*Ghat[2]*dv11; 
  out[17] += -1.224744871391589*Ghat[6]*dv11; 
  out[18] += -1.224744871391589*Ghat[7]*dv11; 
  out[19] += -1.58113883008419*Ghat[3]*dv11; 
  out[20] += -0.7071067811865475*Ghat[8]*dv11; 
  out[21] += -1.58113883008419*Ghat[4]*dv11; 
  out[22] += -1.58113883008419*Ghat[5]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += -1.58113883008419*Ghat[6]*dv11; 
  out[25] += -1.58113883008419*Ghat[7]*dv11; 
  out[26] += -1.58113883008419*Ghat[8]*dv11; 

  } else { 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if (0.5*alpha[0]-0.5590169943749475*alpha[4] > 0) { 
    fUpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(fSkin); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[6]*fUpwind[6]+0.5*alpha[4]*fUpwind[4]+0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.447213595499958*alpha[3]*fUpwind[6]+0.447213595499958*fUpwind[3]*alpha[6]+0.4472135954999579*alpha[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alpha[4]+0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.447213595499958*alpha[6]*fUpwind[8]+0.447213595499958*alpha[3]*fUpwind[7]+0.5000000000000001*alpha[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alpha[6]+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.4*alpha[3]*fUpwind[8]+0.4*alpha[6]*fUpwind[7]+0.447213595499958*alpha[2]*fUpwind[7]+0.447213595499958*alpha[1]*fUpwind[6]+0.447213595499958*fUpwind[1]*alpha[6]+0.4472135954999579*alpha[3]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alpha[4]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.31943828249997*alpha[6]*fUpwind[6]+0.5000000000000001*alpha[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alpha[6]+0.31943828249997*alpha[4]*fUpwind[4]+0.5*alpha[0]*fUpwind[4]+0.5*fUpwind[0]*alpha[4]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[1]*fUpwind[1]; 
  Ghat[5] = 0.5*alpha[4]*fUpwind[8]+0.5000000000000001*alpha[1]*fUpwind[7]+0.4472135954999579*alpha[6]*fUpwind[6]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[6] = 0.2857142857142857*alpha[6]*fUpwind[8]+0.447213595499958*alpha[2]*fUpwind[8]+0.4*alpha[3]*fUpwind[7]+0.31943828249997*alpha[4]*fUpwind[6]+0.5*alpha[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alpha[6]+0.31943828249997*fUpwind[4]*alpha[6]+0.5*fUpwind[0]*alpha[6]+0.5000000000000001*alpha[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alpha[4]+0.447213595499958*alpha[1]*fUpwind[3]+0.447213595499958*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.447213595499958*alpha[1]*fUpwind[8]+0.4472135954999579*alpha[4]*fUpwind[7]+0.5*alpha[0]*fUpwind[7]+0.4*alpha[3]*fUpwind[6]+0.4*fUpwind[3]*alpha[6]+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*alpha[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.31943828249997*alpha[4]*fUpwind[8]+0.5*alpha[0]*fUpwind[8]+0.447213595499958*alpha[1]*fUpwind[7]+0.2857142857142857*alpha[6]*fUpwind[6]+0.447213595499958*alpha[2]*fUpwind[6]+0.447213595499958*fUpwind[2]*alpha[6]+0.5*alpha[4]*fUpwind[5]+0.4*alpha[3]*fUpwind[3]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += 0.7071067811865475*Ghat[4]*dv11; 
  out[8] += 0.7071067811865475*Ghat[5]*dv11; 
  out[9] += 1.58113883008419*Ghat[0]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += 0.7071067811865475*Ghat[6]*dv11; 
  out[12] += 0.7071067811865475*Ghat[7]*dv11; 
  out[13] += -1.224744871391589*Ghat[4]*dv11; 
  out[14] += -1.224744871391589*Ghat[5]*dv11; 
  out[15] += 1.58113883008419*Ghat[1]*dv11; 
  out[16] += 1.58113883008419*Ghat[2]*dv11; 
  out[17] += -1.224744871391589*Ghat[6]*dv11; 
  out[18] += -1.224744871391589*Ghat[7]*dv11; 
  out[19] += 1.58113883008419*Ghat[3]*dv11; 
  out[20] += 0.7071067811865475*Ghat[8]*dv11; 
  out[21] += 1.58113883008419*Ghat[4]*dv11; 
  out[22] += 1.58113883008419*Ghat[5]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += 1.58113883008419*Ghat[6]*dv11; 
  out[25] += 1.58113883008419*Ghat[7]*dv11; 
  out[26] += 1.58113883008419*Ghat[8]*dv11; 

  } 
  return 0.;

} 
