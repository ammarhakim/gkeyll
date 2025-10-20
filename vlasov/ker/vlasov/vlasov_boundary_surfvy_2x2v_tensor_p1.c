#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_tensor_4x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_tensor_4x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // field:       q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *E1 = &field[4]; 
  const double *B2 = &field[20]; 

  double alpha[8] = {0.0}; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[3] = -0.408248290463863*B2[0]*dv1; 
  alpha[4] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[5] = -0.408248290463863*B2[1]*dv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[7] = -0.408248290463863*B2[3]*dv1; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if ((-0.3535533905932737*alpha[7])+0.3535533905932737*(alpha[6]+alpha[5]+alpha[4])-0.3535533905932737*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_4x_p1_surfx4_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = tensor_4x_p1_surfx4_eval_quad_node_0_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[7]-0.3535533905932737*(alpha[6]+alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_4x_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = tensor_4x_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[7]-0.3535533905932737*alpha[6]+0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_4x_p1_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = tensor_4x_p1_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if ((-0.3535533905932737*alpha[7])+0.3535533905932737*alpha[6]-0.3535533905932737*(alpha[5]+alpha[4])+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_4x_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = tensor_4x_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (0.3535533905932737*(alpha[7]+alpha[6])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[4] = tensor_4x_p1_surfx4_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = tensor_4x_p1_surfx4_eval_quad_node_4_l(fEdge); 
  } 
  if ((-0.3535533905932737*(alpha[7]+alpha[6]))+0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[5] = tensor_4x_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = tensor_4x_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]))+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[6] = tensor_4x_p1_surfx4_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = tensor_4x_p1_surfx4_eval_quad_node_6_l(fEdge); 
  } 
  if (0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[7] = tensor_4x_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = tensor_4x_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_4x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[7]*fUpwind[7]+0.3535533905932737*alpha[6]*fUpwind[6]+0.3535533905932737*alpha[5]*fUpwind[5]+0.3535533905932737*alpha[4]*fUpwind[4]+0.3535533905932737*alpha[3]*fUpwind[3]+0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.3535533905932737*alpha[5]*fUpwind[7]+0.3535533905932737*fUpwind[5]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3535533905932737*alpha[4]*fUpwind[7]+0.3535533905932737*fUpwind[4]*alpha[7]+0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.3535533905932737*alpha[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3535533905932737*alpha[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alpha[3]; 
  Ghat[6] = 0.3535533905932737*alpha[1]*fUpwind[7]+0.3535533905932737*fUpwind[1]*alpha[7]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alpha[3]; 
  Ghat[7] = 0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alpha[7]+0.3535533905932737*alpha[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alpha[4]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += -0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -0.7071067811865475*Ghat[5]*dv11; 
  out[7] += -0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -1.224744871391589*Ghat[1]*dv11; 
  out[9] += -1.224744871391589*Ghat[2]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += -0.7071067811865475*Ghat[7]*dv11; 
  out[12] += -1.224744871391589*Ghat[4]*dv11; 
  out[13] += -1.224744871391589*Ghat[5]*dv11; 
  out[14] += -1.224744871391589*Ghat[6]*dv11; 
  out[15] += -1.224744871391589*Ghat[7]*dv11; 

  } else { 

  if ((-0.3535533905932737*alpha[7])+0.3535533905932737*(alpha[6]+alpha[5]+alpha[4])-0.3535533905932737*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_4x_p1_surfx4_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = tensor_4x_p1_surfx4_eval_quad_node_0_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[7]-0.3535533905932737*(alpha[6]+alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_4x_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = tensor_4x_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[7]-0.3535533905932737*alpha[6]+0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_4x_p1_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = tensor_4x_p1_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if ((-0.3535533905932737*alpha[7])+0.3535533905932737*alpha[6]-0.3535533905932737*(alpha[5]+alpha[4])+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_4x_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = tensor_4x_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (0.3535533905932737*(alpha[7]+alpha[6])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[4] = tensor_4x_p1_surfx4_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = tensor_4x_p1_surfx4_eval_quad_node_4_l(fSkin); 
  } 
  if ((-0.3535533905932737*(alpha[7]+alpha[6]))+0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[5] = tensor_4x_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = tensor_4x_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]))+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[6] = tensor_4x_p1_surfx4_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = tensor_4x_p1_surfx4_eval_quad_node_6_l(fSkin); 
  } 
  if (0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[7] = tensor_4x_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = tensor_4x_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_4x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[7]*fUpwind[7]+0.3535533905932737*alpha[6]*fUpwind[6]+0.3535533905932737*alpha[5]*fUpwind[5]+0.3535533905932737*alpha[4]*fUpwind[4]+0.3535533905932737*alpha[3]*fUpwind[3]+0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.3535533905932737*alpha[5]*fUpwind[7]+0.3535533905932737*fUpwind[5]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3535533905932737*alpha[4]*fUpwind[7]+0.3535533905932737*fUpwind[4]*alpha[7]+0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.3535533905932737*alpha[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3535533905932737*alpha[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alpha[3]; 
  Ghat[6] = 0.3535533905932737*alpha[1]*fUpwind[7]+0.3535533905932737*fUpwind[1]*alpha[7]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alpha[3]; 
  Ghat[7] = 0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alpha[7]+0.3535533905932737*alpha[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alpha[4]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -1.224744871391589*Ghat[1]*dv11; 
  out[9] += -1.224744871391589*Ghat[2]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += 0.7071067811865475*Ghat[7]*dv11; 
  out[12] += -1.224744871391589*Ghat[4]*dv11; 
  out[13] += -1.224744871391589*Ghat[5]*dv11; 
  out[14] += -1.224744871391589*Ghat[6]*dv11; 
  out[15] += -1.224744871391589*Ghat[7]*dv11; 

  } 
  return 0.;

} 
