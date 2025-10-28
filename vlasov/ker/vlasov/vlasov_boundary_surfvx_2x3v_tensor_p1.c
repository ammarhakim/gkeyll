#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_tensor_5x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_tensor_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // field:       q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E0 = &field[0]; 
  const double *B0 = &field[12]; 
  const double *B1 = &field[16]; 
  const double *B2 = &field[20]; 

  double alpha[16] = {0.0}; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3])-2.0*B1[3]*wv3; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = 0.5773502691896258*B2[3]*dv2; 
  alpha[12] = -0.5773502691896258*B1[3]*dv3; 

  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  if ((-0.25*(alpha[12]+alpha[11]))+0.25*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.25*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_5x_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = tensor_5x_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (0.25*alpha[12]-0.25*(alpha[11]+alpha[9]+alpha[8])+0.25*(alpha[7]+alpha[6]+alpha[5]+alpha[4])-0.25*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_5x_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = tensor_5x_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if ((-0.25*alpha[12])+0.25*(alpha[11]+alpha[9]+alpha[8])-0.25*(alpha[7]+alpha[6])+0.25*alpha[5]-0.25*alpha[4]+0.25*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_5x_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = tensor_5x_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.25*(alpha[12]+alpha[11])-0.25*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_5x_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = tensor_5x_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if (0.25*(alpha[12]+alpha[11])-0.25*alpha[9]+0.25*alpha[8]-0.25*alpha[7]+0.25*alpha[6]-0.25*(alpha[5]+alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = tensor_5x_p1_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = tensor_5x_p1_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if ((-0.25*alpha[12])+0.25*(alpha[11]+alpha[9])-0.25*(alpha[8]+alpha[7])+0.25*alpha[6]-0.25*alpha[5]+0.25*alpha[4]-0.25*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = tensor_5x_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = tensor_5x_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (0.25*alpha[12]-0.25*(alpha[11]+alpha[9])+0.25*(alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5]+alpha[4])+0.25*(alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = tensor_5x_p1_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = tensor_5x_p1_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]))+0.25*alpha[9]-0.25*alpha[8]+0.25*alpha[7]-0.25*(alpha[6]+alpha[5])+0.25*(alpha[4]+alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = tensor_5x_p1_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = tensor_5x_p1_surfx3_eval_quad_node_7_l(fEdge); 
  } 
  if (0.25*(alpha[12]+alpha[11]+alpha[9])-0.25*alpha[8]+0.25*alpha[7]-0.25*(alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[8] = tensor_5x_p1_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = tensor_5x_p1_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if ((-0.25*alpha[12])+0.25*alpha[11]-0.25*alpha[9]+0.25*(alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5])+0.25*alpha[4]-0.25*(alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[9] = tensor_5x_p1_surfx3_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = tensor_5x_p1_surfx3_eval_quad_node_9_l(fEdge); 
  } 
  if (0.25*alpha[12]-0.25*alpha[11]+0.25*alpha[9]-0.25*(alpha[8]+alpha[7])+0.25*alpha[6]-0.25*(alpha[5]+alpha[4])+0.25*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[10] = tensor_5x_p1_surfx3_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = tensor_5x_p1_surfx3_eval_quad_node_10_l(fEdge); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]+alpha[9]))+0.25*alpha[8]-0.25*alpha[7]+0.25*alpha[6]-0.25*alpha[5]+0.25*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[11] = tensor_5x_p1_surfx3_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = tensor_5x_p1_surfx3_eval_quad_node_11_l(fEdge); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.25*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[12] = tensor_5x_p1_surfx3_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = tensor_5x_p1_surfx3_eval_quad_node_12_l(fEdge); 
  } 
  if (0.25*alpha[12]-0.25*alpha[11]+0.25*(alpha[9]+alpha[8])-0.25*(alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4])-0.25*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[13] = tensor_5x_p1_surfx3_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = tensor_5x_p1_surfx3_eval_quad_node_13_l(fEdge); 
  } 
  if ((-0.25*alpha[12])+0.25*alpha[11]-0.25*(alpha[9]+alpha[8])+0.25*(alpha[7]+alpha[6]+alpha[5])-0.25*alpha[4]+0.25*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[14] = tensor_5x_p1_surfx3_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = tensor_5x_p1_surfx3_eval_quad_node_14_l(fEdge); 
  } 
  if (0.25*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[15] = tensor_5x_p1_surfx3_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = tensor_5x_p1_surfx3_eval_quad_node_15_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_5x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[12]*fUpwind[12]+0.25*alpha[11]*fUpwind[11]+0.25*alpha[9]*fUpwind[9]+0.25*alpha[8]*fUpwind[8]+0.25*alpha[7]*fUpwind[7]+0.25*alpha[6]*fUpwind[6]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[4]*fUpwind[4]+0.25*alpha[3]*fUpwind[3]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.25*alpha[9]*fUpwind[12]+0.25*fUpwind[9]*alpha[12]+0.25*alpha[7]*fUpwind[11]+0.25*fUpwind[7]*alpha[11]+0.25*alpha[4]*fUpwind[8]+0.25*fUpwind[4]*alpha[8]+0.25*alpha[3]*fUpwind[6]+0.25*fUpwind[3]*alpha[6]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.25*alpha[8]*fUpwind[12]+0.25*fUpwind[8]*alpha[12]+0.25*alpha[6]*fUpwind[11]+0.25*fUpwind[6]*alpha[11]+0.25*alpha[4]*fUpwind[9]+0.25*fUpwind[4]*alpha[9]+0.25*alpha[3]*fUpwind[7]+0.25*fUpwind[3]*alpha[7]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.25*alpha[12]*fUpwind[15]+0.25*alpha[9]*fUpwind[14]+0.25*alpha[8]*fUpwind[13]+0.25*alpha[5]*fUpwind[11]+0.25*fUpwind[5]*alpha[11]+0.25*alpha[4]*fUpwind[10]+0.25*alpha[2]*fUpwind[7]+0.25*fUpwind[2]*alpha[7]+0.25*alpha[1]*fUpwind[6]+0.25*fUpwind[1]*alpha[6]+0.25*alpha[0]*fUpwind[3]+0.25*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.25*alpha[11]*fUpwind[15]+0.25*alpha[7]*fUpwind[14]+0.25*alpha[6]*fUpwind[13]+0.25*alpha[5]*fUpwind[12]+0.25*fUpwind[5]*alpha[12]+0.25*alpha[3]*fUpwind[10]+0.25*alpha[2]*fUpwind[9]+0.25*fUpwind[2]*alpha[9]+0.25*alpha[1]*fUpwind[8]+0.25*fUpwind[1]*alpha[8]+0.25*alpha[0]*fUpwind[4]+0.25*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.25*alpha[4]*fUpwind[12]+0.25*fUpwind[4]*alpha[12]+0.25*alpha[3]*fUpwind[11]+0.25*fUpwind[3]*alpha[11]+0.25*alpha[8]*fUpwind[9]+0.25*fUpwind[8]*alpha[9]+0.25*alpha[6]*fUpwind[7]+0.25*fUpwind[6]*alpha[7]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.25*alpha[9]*fUpwind[15]+0.25*alpha[12]*fUpwind[14]+0.25*alpha[4]*fUpwind[13]+0.25*alpha[2]*fUpwind[11]+0.25*fUpwind[2]*alpha[11]+0.25*alpha[8]*fUpwind[10]+0.25*alpha[5]*fUpwind[7]+0.25*fUpwind[5]*alpha[7]+0.25*alpha[0]*fUpwind[6]+0.25*fUpwind[0]*alpha[6]+0.25*alpha[1]*fUpwind[3]+0.25*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.25*alpha[8]*fUpwind[15]+0.25*alpha[4]*fUpwind[14]+0.25*alpha[12]*fUpwind[13]+0.25*alpha[1]*fUpwind[11]+0.25*fUpwind[1]*alpha[11]+0.25*alpha[9]*fUpwind[10]+0.25*alpha[0]*fUpwind[7]+0.25*fUpwind[0]*alpha[7]+0.25*alpha[5]*fUpwind[6]+0.25*fUpwind[5]*alpha[6]+0.25*alpha[2]*fUpwind[3]+0.25*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.25*alpha[7]*fUpwind[15]+0.25*alpha[11]*fUpwind[14]+0.25*alpha[3]*fUpwind[13]+0.25*alpha[2]*fUpwind[12]+0.25*fUpwind[2]*alpha[12]+0.25*alpha[6]*fUpwind[10]+0.25*alpha[5]*fUpwind[9]+0.25*fUpwind[5]*alpha[9]+0.25*alpha[0]*fUpwind[8]+0.25*fUpwind[0]*alpha[8]+0.25*alpha[1]*fUpwind[4]+0.25*fUpwind[1]*alpha[4]; 
  Ghat[9] = 0.25*alpha[6]*fUpwind[15]+0.25*alpha[3]*fUpwind[14]+0.25*alpha[11]*fUpwind[13]+0.25*alpha[1]*fUpwind[12]+0.25*fUpwind[1]*alpha[12]+0.25*alpha[7]*fUpwind[10]+0.25*alpha[0]*fUpwind[9]+0.25*fUpwind[0]*alpha[9]+0.25*alpha[5]*fUpwind[8]+0.25*fUpwind[5]*alpha[8]+0.25*alpha[2]*fUpwind[4]+0.25*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.25*alpha[5]*fUpwind[15]+0.25*alpha[2]*fUpwind[14]+0.25*alpha[1]*fUpwind[13]+0.25*alpha[11]*fUpwind[12]+0.25*fUpwind[11]*alpha[12]+0.25*alpha[0]*fUpwind[10]+0.25*alpha[7]*fUpwind[9]+0.25*fUpwind[7]*alpha[9]+0.25*alpha[6]*fUpwind[8]+0.25*fUpwind[6]*alpha[8]+0.25*alpha[3]*fUpwind[4]+0.25*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.25*alpha[4]*fUpwind[15]+0.25*alpha[8]*fUpwind[14]+0.25*alpha[9]*fUpwind[13]+0.25*fUpwind[10]*alpha[12]+0.25*alpha[0]*fUpwind[11]+0.25*fUpwind[0]*alpha[11]+0.25*alpha[1]*fUpwind[7]+0.25*fUpwind[1]*alpha[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[2]*alpha[6]+0.25*alpha[3]*fUpwind[5]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[12] = 0.25*alpha[3]*fUpwind[15]+0.25*alpha[6]*fUpwind[14]+0.25*alpha[7]*fUpwind[13]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.25*fUpwind[10]*alpha[11]+0.25*alpha[1]*fUpwind[9]+0.25*fUpwind[1]*alpha[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[2]*alpha[8]+0.25*alpha[4]*fUpwind[5]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[13] = 0.25*alpha[2]*fUpwind[15]+0.25*alpha[5]*fUpwind[14]+0.25*alpha[0]*fUpwind[13]+0.25*alpha[7]*fUpwind[12]+0.25*fUpwind[7]*alpha[12]+0.25*alpha[9]*fUpwind[11]+0.25*fUpwind[9]*alpha[11]+0.25*alpha[1]*fUpwind[10]+0.25*alpha[3]*fUpwind[8]+0.25*fUpwind[3]*alpha[8]+0.25*alpha[4]*fUpwind[6]+0.25*fUpwind[4]*alpha[6]; 
  Ghat[14] = 0.25*alpha[1]*fUpwind[15]+0.25*alpha[0]*fUpwind[14]+0.25*alpha[5]*fUpwind[13]+0.25*alpha[6]*fUpwind[12]+0.25*fUpwind[6]*alpha[12]+0.25*alpha[8]*fUpwind[11]+0.25*fUpwind[8]*alpha[11]+0.25*alpha[2]*fUpwind[10]+0.25*alpha[3]*fUpwind[9]+0.25*fUpwind[3]*alpha[9]+0.25*alpha[4]*fUpwind[7]+0.25*fUpwind[4]*alpha[7]; 
  Ghat[15] = 0.25*alpha[0]*fUpwind[15]+0.25*alpha[1]*fUpwind[14]+0.25*alpha[2]*fUpwind[13]+0.25*alpha[3]*fUpwind[12]+0.25*fUpwind[3]*alpha[12]+0.25*alpha[4]*fUpwind[11]+0.25*fUpwind[4]*alpha[11]+0.25*alpha[5]*fUpwind[10]+0.25*alpha[6]*fUpwind[9]+0.25*fUpwind[6]*alpha[9]+0.25*alpha[7]*fUpwind[8]+0.25*fUpwind[7]*alpha[8]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += -0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += -0.7071067811865475*Ghat[8]*dv10; 
  out[13] += -0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += -0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -0.7071067811865475*Ghat[11]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 
  out[20] += -0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[8]*dv10; 
  out[22] += -1.224744871391589*Ghat[9]*dv10; 
  out[23] += -0.7071067811865475*Ghat[13]*dv10; 
  out[24] += -0.7071067811865475*Ghat[14]*dv10; 
  out[25] += -1.224744871391589*Ghat[10]*dv10; 
  out[26] += -1.224744871391589*Ghat[11]*dv10; 
  out[27] += -1.224744871391589*Ghat[12]*dv10; 
  out[28] += -0.7071067811865475*Ghat[15]*dv10; 
  out[29] += -1.224744871391589*Ghat[13]*dv10; 
  out[30] += -1.224744871391589*Ghat[14]*dv10; 
  out[31] += -1.224744871391589*Ghat[15]*dv10; 

  } else { 

  if ((-0.25*(alpha[12]+alpha[11]))+0.25*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.25*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = tensor_5x_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = tensor_5x_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (0.25*alpha[12]-0.25*(alpha[11]+alpha[9]+alpha[8])+0.25*(alpha[7]+alpha[6]+alpha[5]+alpha[4])-0.25*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = tensor_5x_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = tensor_5x_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if ((-0.25*alpha[12])+0.25*(alpha[11]+alpha[9]+alpha[8])-0.25*(alpha[7]+alpha[6])+0.25*alpha[5]-0.25*alpha[4]+0.25*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = tensor_5x_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = tensor_5x_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.25*(alpha[12]+alpha[11])-0.25*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = tensor_5x_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = tensor_5x_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if (0.25*(alpha[12]+alpha[11])-0.25*alpha[9]+0.25*alpha[8]-0.25*alpha[7]+0.25*alpha[6]-0.25*(alpha[5]+alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = tensor_5x_p1_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = tensor_5x_p1_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if ((-0.25*alpha[12])+0.25*(alpha[11]+alpha[9])-0.25*(alpha[8]+alpha[7])+0.25*alpha[6]-0.25*alpha[5]+0.25*alpha[4]-0.25*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = tensor_5x_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = tensor_5x_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (0.25*alpha[12]-0.25*(alpha[11]+alpha[9])+0.25*(alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5]+alpha[4])+0.25*(alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = tensor_5x_p1_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = tensor_5x_p1_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]))+0.25*alpha[9]-0.25*alpha[8]+0.25*alpha[7]-0.25*(alpha[6]+alpha[5])+0.25*(alpha[4]+alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = tensor_5x_p1_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = tensor_5x_p1_surfx3_eval_quad_node_7_l(fSkin); 
  } 
  if (0.25*(alpha[12]+alpha[11]+alpha[9])-0.25*alpha[8]+0.25*alpha[7]-0.25*(alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[8] = tensor_5x_p1_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = tensor_5x_p1_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if ((-0.25*alpha[12])+0.25*alpha[11]-0.25*alpha[9]+0.25*(alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5])+0.25*alpha[4]-0.25*(alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[9] = tensor_5x_p1_surfx3_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = tensor_5x_p1_surfx3_eval_quad_node_9_l(fSkin); 
  } 
  if (0.25*alpha[12]-0.25*alpha[11]+0.25*alpha[9]-0.25*(alpha[8]+alpha[7])+0.25*alpha[6]-0.25*(alpha[5]+alpha[4])+0.25*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[10] = tensor_5x_p1_surfx3_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = tensor_5x_p1_surfx3_eval_quad_node_10_l(fSkin); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]+alpha[9]))+0.25*alpha[8]-0.25*alpha[7]+0.25*alpha[6]-0.25*alpha[5]+0.25*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[11] = tensor_5x_p1_surfx3_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = tensor_5x_p1_surfx3_eval_quad_node_11_l(fSkin); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.25*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[12] = tensor_5x_p1_surfx3_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = tensor_5x_p1_surfx3_eval_quad_node_12_l(fSkin); 
  } 
  if (0.25*alpha[12]-0.25*alpha[11]+0.25*(alpha[9]+alpha[8])-0.25*(alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4])-0.25*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[13] = tensor_5x_p1_surfx3_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = tensor_5x_p1_surfx3_eval_quad_node_13_l(fSkin); 
  } 
  if ((-0.25*alpha[12])+0.25*alpha[11]-0.25*(alpha[9]+alpha[8])+0.25*(alpha[7]+alpha[6]+alpha[5])-0.25*alpha[4]+0.25*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[14] = tensor_5x_p1_surfx3_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = tensor_5x_p1_surfx3_eval_quad_node_14_l(fSkin); 
  } 
  if (0.25*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[15] = tensor_5x_p1_surfx3_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = tensor_5x_p1_surfx3_eval_quad_node_15_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_5x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[12]*fUpwind[12]+0.25*alpha[11]*fUpwind[11]+0.25*alpha[9]*fUpwind[9]+0.25*alpha[8]*fUpwind[8]+0.25*alpha[7]*fUpwind[7]+0.25*alpha[6]*fUpwind[6]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[4]*fUpwind[4]+0.25*alpha[3]*fUpwind[3]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.25*alpha[9]*fUpwind[12]+0.25*fUpwind[9]*alpha[12]+0.25*alpha[7]*fUpwind[11]+0.25*fUpwind[7]*alpha[11]+0.25*alpha[4]*fUpwind[8]+0.25*fUpwind[4]*alpha[8]+0.25*alpha[3]*fUpwind[6]+0.25*fUpwind[3]*alpha[6]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.25*alpha[8]*fUpwind[12]+0.25*fUpwind[8]*alpha[12]+0.25*alpha[6]*fUpwind[11]+0.25*fUpwind[6]*alpha[11]+0.25*alpha[4]*fUpwind[9]+0.25*fUpwind[4]*alpha[9]+0.25*alpha[3]*fUpwind[7]+0.25*fUpwind[3]*alpha[7]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.25*alpha[12]*fUpwind[15]+0.25*alpha[9]*fUpwind[14]+0.25*alpha[8]*fUpwind[13]+0.25*alpha[5]*fUpwind[11]+0.25*fUpwind[5]*alpha[11]+0.25*alpha[4]*fUpwind[10]+0.25*alpha[2]*fUpwind[7]+0.25*fUpwind[2]*alpha[7]+0.25*alpha[1]*fUpwind[6]+0.25*fUpwind[1]*alpha[6]+0.25*alpha[0]*fUpwind[3]+0.25*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.25*alpha[11]*fUpwind[15]+0.25*alpha[7]*fUpwind[14]+0.25*alpha[6]*fUpwind[13]+0.25*alpha[5]*fUpwind[12]+0.25*fUpwind[5]*alpha[12]+0.25*alpha[3]*fUpwind[10]+0.25*alpha[2]*fUpwind[9]+0.25*fUpwind[2]*alpha[9]+0.25*alpha[1]*fUpwind[8]+0.25*fUpwind[1]*alpha[8]+0.25*alpha[0]*fUpwind[4]+0.25*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.25*alpha[4]*fUpwind[12]+0.25*fUpwind[4]*alpha[12]+0.25*alpha[3]*fUpwind[11]+0.25*fUpwind[3]*alpha[11]+0.25*alpha[8]*fUpwind[9]+0.25*fUpwind[8]*alpha[9]+0.25*alpha[6]*fUpwind[7]+0.25*fUpwind[6]*alpha[7]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.25*alpha[9]*fUpwind[15]+0.25*alpha[12]*fUpwind[14]+0.25*alpha[4]*fUpwind[13]+0.25*alpha[2]*fUpwind[11]+0.25*fUpwind[2]*alpha[11]+0.25*alpha[8]*fUpwind[10]+0.25*alpha[5]*fUpwind[7]+0.25*fUpwind[5]*alpha[7]+0.25*alpha[0]*fUpwind[6]+0.25*fUpwind[0]*alpha[6]+0.25*alpha[1]*fUpwind[3]+0.25*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.25*alpha[8]*fUpwind[15]+0.25*alpha[4]*fUpwind[14]+0.25*alpha[12]*fUpwind[13]+0.25*alpha[1]*fUpwind[11]+0.25*fUpwind[1]*alpha[11]+0.25*alpha[9]*fUpwind[10]+0.25*alpha[0]*fUpwind[7]+0.25*fUpwind[0]*alpha[7]+0.25*alpha[5]*fUpwind[6]+0.25*fUpwind[5]*alpha[6]+0.25*alpha[2]*fUpwind[3]+0.25*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.25*alpha[7]*fUpwind[15]+0.25*alpha[11]*fUpwind[14]+0.25*alpha[3]*fUpwind[13]+0.25*alpha[2]*fUpwind[12]+0.25*fUpwind[2]*alpha[12]+0.25*alpha[6]*fUpwind[10]+0.25*alpha[5]*fUpwind[9]+0.25*fUpwind[5]*alpha[9]+0.25*alpha[0]*fUpwind[8]+0.25*fUpwind[0]*alpha[8]+0.25*alpha[1]*fUpwind[4]+0.25*fUpwind[1]*alpha[4]; 
  Ghat[9] = 0.25*alpha[6]*fUpwind[15]+0.25*alpha[3]*fUpwind[14]+0.25*alpha[11]*fUpwind[13]+0.25*alpha[1]*fUpwind[12]+0.25*fUpwind[1]*alpha[12]+0.25*alpha[7]*fUpwind[10]+0.25*alpha[0]*fUpwind[9]+0.25*fUpwind[0]*alpha[9]+0.25*alpha[5]*fUpwind[8]+0.25*fUpwind[5]*alpha[8]+0.25*alpha[2]*fUpwind[4]+0.25*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.25*alpha[5]*fUpwind[15]+0.25*alpha[2]*fUpwind[14]+0.25*alpha[1]*fUpwind[13]+0.25*alpha[11]*fUpwind[12]+0.25*fUpwind[11]*alpha[12]+0.25*alpha[0]*fUpwind[10]+0.25*alpha[7]*fUpwind[9]+0.25*fUpwind[7]*alpha[9]+0.25*alpha[6]*fUpwind[8]+0.25*fUpwind[6]*alpha[8]+0.25*alpha[3]*fUpwind[4]+0.25*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.25*alpha[4]*fUpwind[15]+0.25*alpha[8]*fUpwind[14]+0.25*alpha[9]*fUpwind[13]+0.25*fUpwind[10]*alpha[12]+0.25*alpha[0]*fUpwind[11]+0.25*fUpwind[0]*alpha[11]+0.25*alpha[1]*fUpwind[7]+0.25*fUpwind[1]*alpha[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[2]*alpha[6]+0.25*alpha[3]*fUpwind[5]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[12] = 0.25*alpha[3]*fUpwind[15]+0.25*alpha[6]*fUpwind[14]+0.25*alpha[7]*fUpwind[13]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.25*fUpwind[10]*alpha[11]+0.25*alpha[1]*fUpwind[9]+0.25*fUpwind[1]*alpha[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[2]*alpha[8]+0.25*alpha[4]*fUpwind[5]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[13] = 0.25*alpha[2]*fUpwind[15]+0.25*alpha[5]*fUpwind[14]+0.25*alpha[0]*fUpwind[13]+0.25*alpha[7]*fUpwind[12]+0.25*fUpwind[7]*alpha[12]+0.25*alpha[9]*fUpwind[11]+0.25*fUpwind[9]*alpha[11]+0.25*alpha[1]*fUpwind[10]+0.25*alpha[3]*fUpwind[8]+0.25*fUpwind[3]*alpha[8]+0.25*alpha[4]*fUpwind[6]+0.25*fUpwind[4]*alpha[6]; 
  Ghat[14] = 0.25*alpha[1]*fUpwind[15]+0.25*alpha[0]*fUpwind[14]+0.25*alpha[5]*fUpwind[13]+0.25*alpha[6]*fUpwind[12]+0.25*fUpwind[6]*alpha[12]+0.25*alpha[8]*fUpwind[11]+0.25*fUpwind[8]*alpha[11]+0.25*alpha[2]*fUpwind[10]+0.25*alpha[3]*fUpwind[9]+0.25*fUpwind[3]*alpha[9]+0.25*alpha[4]*fUpwind[7]+0.25*fUpwind[4]*alpha[7]; 
  Ghat[15] = 0.25*alpha[0]*fUpwind[15]+0.25*alpha[1]*fUpwind[14]+0.25*alpha[2]*fUpwind[13]+0.25*alpha[3]*fUpwind[12]+0.25*fUpwind[3]*alpha[12]+0.25*alpha[4]*fUpwind[11]+0.25*fUpwind[4]*alpha[11]+0.25*alpha[5]*fUpwind[10]+0.25*alpha[6]*fUpwind[9]+0.25*fUpwind[6]*alpha[9]+0.25*alpha[7]*fUpwind[8]+0.25*fUpwind[7]*alpha[8]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += 0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += 0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += 0.7071067811865475*Ghat[11]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 
  out[20] += 0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[8]*dv10; 
  out[22] += -1.224744871391589*Ghat[9]*dv10; 
  out[23] += 0.7071067811865475*Ghat[13]*dv10; 
  out[24] += 0.7071067811865475*Ghat[14]*dv10; 
  out[25] += -1.224744871391589*Ghat[10]*dv10; 
  out[26] += -1.224744871391589*Ghat[11]*dv10; 
  out[27] += -1.224744871391589*Ghat[12]*dv10; 
  out[28] += 0.7071067811865475*Ghat[15]*dv10; 
  out[29] += -1.224744871391589*Ghat[13]*dv10; 
  out[30] += -1.224744871391589*Ghat[14]*dv10; 
  out[31] += -1.224744871391589*Ghat[15]*dv10; 

  } 
  return 0.;

} 
