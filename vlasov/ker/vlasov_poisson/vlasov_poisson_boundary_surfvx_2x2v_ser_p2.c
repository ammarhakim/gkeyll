#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // pots:        potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:       external E and B fields (scaled by q/m).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 

  const double dv10 = 2/dxv[2]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *phi = &pots[0]; 

  double alpha[20] = {0.0}; 

  alpha[0] = -(2.4494897427831783*phi[1]*dx10); 
  alpha[1] = -(5.477225575051662*phi[4]*dx10); 
  alpha[2] = -(2.4494897427831783*phi[3]*dx10); 
  alpha[4] = -(5.477225575051662*phi[6]*dx10); 
  alpha[8] = -(2.4494897427831783*phi[7]*dx10); 

  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};
  double Ghat[20] = {0.0}; 

  if (edge == -1) { 

  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if (-(0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fEdge); 
  } 
  if (-(0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fEdge); 
  } 
  if (-(0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fSkin); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fEdge); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[4]*fUpwind[4]+0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[8]*fUpwind[12]+0.31622776601683794*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[1]*fUpwind[7]+0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.31622776601683794*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alpha[8]+0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3535533905932737*alpha[8]*fUpwind[14]+0.3535533905932737*alpha[4]*fUpwind[10]+0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.31622776601683794*alpha[2]*fUpwind[12]+0.31622776601683794*alpha[1]*fUpwind[11]+0.3162277660168379*alpha[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[7]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3535533905932737*alpha[8]*fUpwind[18]+0.3162277660168379*alpha[4]*fUpwind[17]+0.31622776601683794*alpha[1]*fUpwind[13]+0.3535533905932737*alpha[2]*fUpwind[10]+0.3535533905932737*alpha[4]*fUpwind[6]+0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*alpha[1]*fUpwind[3]; 
  Ghat[6] = 0.3162277660168379*alpha[4]*fUpwind[18]+0.31622776601683794*alpha[2]*fUpwind[14]+0.3535533905932737*alpha[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*alpha[4]*fUpwind[5]+0.3535533905932737*alpha[2]*fUpwind[3]; 
  Ghat[7] = 0.3535533905932737*alpha[2]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[7]+0.3162277660168379*alpha[4]*fUpwind[4]+0.3162277660168379*alpha[1]*fUpwind[1]; 
  Ghat[8] = 0.3535533905932737*alpha[1]*fUpwind[12]+0.22587697572631277*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[4]+0.3162277660168379*alpha[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alpha[4]*fUpwind[19]+0.3535533905932737*alpha[2]*fUpwind[16]+0.3535533905932737*alpha[1]*fUpwind[15]+0.3535533905932737*alpha[0]*fUpwind[9]; 
  Ghat[10] = 0.3162277660168379*alpha[2]*fUpwind[18]+0.3162277660168379*alpha[1]*fUpwind[17]+0.31622776601683794*alpha[4]*fUpwind[14]+0.31622776601683794*alpha[4]*fUpwind[13]+0.3162277660168379*alpha[8]*fUpwind[10]+0.3535533905932737*alpha[0]*fUpwind[10]+0.3535533905932737*alpha[1]*fUpwind[6]+0.3535533905932737*alpha[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.28284271247461906*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[8]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[11]+0.3535533905932737*alpha[2]*fUpwind[7]+0.31622776601683794*alpha[1]*fUpwind[4]+0.31622776601683794*fUpwind[1]*alpha[4]; 
  Ghat[12] = 0.22587697572631277*alpha[8]*fUpwind[12]+0.3535533905932737*alpha[0]*fUpwind[12]+0.28284271247461906*alpha[4]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind[4]+0.31622776601683794*fUpwind[2]*alpha[4]; 
  Ghat[13] = 0.3535533905932737*alpha[2]*fUpwind[17]+0.3535533905932737*alpha[0]*fUpwind[13]+0.31622776601683794*alpha[4]*fUpwind[10]+0.31622776601683794*alpha[1]*fUpwind[5]; 
  Ghat[14] = 0.3535533905932737*alpha[1]*fUpwind[18]+0.22587697572631277*alpha[8]*fUpwind[14]+0.3535533905932737*alpha[0]*fUpwind[14]+0.31622776601683794*alpha[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind[6]; 
  Ghat[15] = 0.3535533905932737*alpha[2]*fUpwind[19]+0.3535533905932737*alpha[4]*fUpwind[16]+0.3535533905932737*alpha[0]*fUpwind[15]+0.3535533905932737*alpha[1]*fUpwind[9]; 
  Ghat[16] = 0.3535533905932737*alpha[1]*fUpwind[19]+0.3162277660168379*alpha[8]*fUpwind[16]+0.3535533905932737*alpha[0]*fUpwind[16]+0.3535533905932737*alpha[4]*fUpwind[15]+0.3535533905932737*alpha[2]*fUpwind[9]; 
  Ghat[17] = 0.28284271247461906*alpha[4]*fUpwind[18]+0.3162277660168379*alpha[8]*fUpwind[17]+0.3535533905932737*alpha[0]*fUpwind[17]+0.3535533905932737*alpha[2]*fUpwind[13]+0.3162277660168379*alpha[1]*fUpwind[10]+0.3162277660168379*alpha[4]*fUpwind[5]; 
  Ghat[18] = 0.22587697572631277*alpha[8]*fUpwind[18]+0.3535533905932737*alpha[0]*fUpwind[18]+0.28284271247461906*alpha[4]*fUpwind[17]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3162277660168379*alpha[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[6]; 
  Ghat[19] = 0.3162277660168379*alpha[8]*fUpwind[19]+0.3535533905932737*alpha[0]*fUpwind[19]+0.3535533905932737*alpha[1]*fUpwind[16]+0.3535533905932737*alpha[2]*fUpwind[15]+0.3535533905932737*alpha[4]*fUpwind[9]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv10); 
  out[3] += -(1.224744871391589*Ghat[0]*dv10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv10); 
  out[6] += -(1.224744871391589*Ghat[1]*dv10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10); 
  out[8] += -(0.7071067811865475*Ghat[5]*dv10); 
  out[9] += -(0.7071067811865475*Ghat[6]*dv10); 
  out[10] += -(1.224744871391589*Ghat[3]*dv10); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv10); 
  out[12] += -(0.7071067811865475*Ghat[8]*dv10); 
  out[13] += -(1.5811388300841895*Ghat[0]*dv10); 
  out[14] += -(0.7071067811865475*Ghat[9]*dv10); 
  out[15] += -(1.224744871391589*Ghat[4]*dv10); 
  out[16] += -(0.7071067811865475*Ghat[10]*dv10); 
  out[17] += -(1.224744871391589*Ghat[5]*dv10); 
  out[18] += -(1.224744871391589*Ghat[6]*dv10); 
  out[19] += -(0.7071067811865475*Ghat[11]*dv10); 
  out[20] += -(0.7071067811865475*Ghat[12]*dv10); 
  out[21] += -(1.224744871391589*Ghat[7]*dv10); 
  out[22] += -(1.224744871391589*Ghat[8]*dv10); 
  out[23] += -(1.5811388300841898*Ghat[1]*dv10); 
  out[24] += -(1.5811388300841898*Ghat[2]*dv10); 
  out[25] += -(0.7071067811865475*Ghat[13]*dv10); 
  out[26] += -(0.7071067811865475*Ghat[14]*dv10); 
  out[27] += -(1.5811388300841898*Ghat[3]*dv10); 
  out[28] += -(0.7071067811865475*Ghat[15]*dv10); 
  out[29] += -(0.7071067811865475*Ghat[16]*dv10); 
  out[30] += -(1.224744871391589*Ghat[9]*dv10); 
  out[31] += -(1.224744871391589*Ghat[10]*dv10); 
  out[32] += -(1.224744871391589*Ghat[11]*dv10); 
  out[33] += -(1.224744871391589*Ghat[12]*dv10); 
  out[34] += -(1.5811388300841895*Ghat[4]*dv10); 
  out[35] += -(0.7071067811865475*Ghat[17]*dv10); 
  out[36] += -(0.7071067811865475*Ghat[18]*dv10); 
  out[37] += -(1.224744871391589*Ghat[13]*dv10); 
  out[38] += -(1.224744871391589*Ghat[14]*dv10); 
  out[39] += -(1.5811388300841895*Ghat[5]*dv10); 
  out[40] += -(1.5811388300841895*Ghat[6]*dv10); 
  out[41] += -(0.7071067811865475*Ghat[19]*dv10); 
  out[42] += -(1.224744871391589*Ghat[15]*dv10); 
  out[43] += -(1.224744871391589*Ghat[16]*dv10); 
  out[44] += -(1.224744871391589*Ghat[17]*dv10); 
  out[45] += -(1.224744871391589*Ghat[18]*dv10); 
  out[46] += -(1.5811388300841898*Ghat[10]*dv10); 
  out[47] += -(1.224744871391589*Ghat[19]*dv10); 

  } else { 

  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if (-(0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fSkin); 
  } 
  if (-(0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fSkin); 
  } 
  if (-(0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fEdge); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fSkin); 
    fUpwindQuad[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[4]*fUpwind[4]+0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[8]*fUpwind[12]+0.31622776601683794*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[1]*fUpwind[7]+0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.31622776601683794*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alpha[8]+0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3535533905932737*alpha[8]*fUpwind[14]+0.3535533905932737*alpha[4]*fUpwind[10]+0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.31622776601683794*alpha[2]*fUpwind[12]+0.31622776601683794*alpha[1]*fUpwind[11]+0.3162277660168379*alpha[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[7]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3535533905932737*alpha[8]*fUpwind[18]+0.3162277660168379*alpha[4]*fUpwind[17]+0.31622776601683794*alpha[1]*fUpwind[13]+0.3535533905932737*alpha[2]*fUpwind[10]+0.3535533905932737*alpha[4]*fUpwind[6]+0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*alpha[1]*fUpwind[3]; 
  Ghat[6] = 0.3162277660168379*alpha[4]*fUpwind[18]+0.31622776601683794*alpha[2]*fUpwind[14]+0.3535533905932737*alpha[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*alpha[4]*fUpwind[5]+0.3535533905932737*alpha[2]*fUpwind[3]; 
  Ghat[7] = 0.3535533905932737*alpha[2]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[7]+0.3162277660168379*alpha[4]*fUpwind[4]+0.3162277660168379*alpha[1]*fUpwind[1]; 
  Ghat[8] = 0.3535533905932737*alpha[1]*fUpwind[12]+0.22587697572631277*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[4]+0.3162277660168379*alpha[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alpha[4]*fUpwind[19]+0.3535533905932737*alpha[2]*fUpwind[16]+0.3535533905932737*alpha[1]*fUpwind[15]+0.3535533905932737*alpha[0]*fUpwind[9]; 
  Ghat[10] = 0.3162277660168379*alpha[2]*fUpwind[18]+0.3162277660168379*alpha[1]*fUpwind[17]+0.31622776601683794*alpha[4]*fUpwind[14]+0.31622776601683794*alpha[4]*fUpwind[13]+0.3162277660168379*alpha[8]*fUpwind[10]+0.3535533905932737*alpha[0]*fUpwind[10]+0.3535533905932737*alpha[1]*fUpwind[6]+0.3535533905932737*alpha[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.28284271247461906*alpha[4]*fUpwind[12]+0.3162277660168379*alpha[8]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[11]+0.3535533905932737*alpha[2]*fUpwind[7]+0.31622776601683794*alpha[1]*fUpwind[4]+0.31622776601683794*fUpwind[1]*alpha[4]; 
  Ghat[12] = 0.22587697572631277*alpha[8]*fUpwind[12]+0.3535533905932737*alpha[0]*fUpwind[12]+0.28284271247461906*alpha[4]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind[4]+0.31622776601683794*fUpwind[2]*alpha[4]; 
  Ghat[13] = 0.3535533905932737*alpha[2]*fUpwind[17]+0.3535533905932737*alpha[0]*fUpwind[13]+0.31622776601683794*alpha[4]*fUpwind[10]+0.31622776601683794*alpha[1]*fUpwind[5]; 
  Ghat[14] = 0.3535533905932737*alpha[1]*fUpwind[18]+0.22587697572631277*alpha[8]*fUpwind[14]+0.3535533905932737*alpha[0]*fUpwind[14]+0.31622776601683794*alpha[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind[6]; 
  Ghat[15] = 0.3535533905932737*alpha[2]*fUpwind[19]+0.3535533905932737*alpha[4]*fUpwind[16]+0.3535533905932737*alpha[0]*fUpwind[15]+0.3535533905932737*alpha[1]*fUpwind[9]; 
  Ghat[16] = 0.3535533905932737*alpha[1]*fUpwind[19]+0.3162277660168379*alpha[8]*fUpwind[16]+0.3535533905932737*alpha[0]*fUpwind[16]+0.3535533905932737*alpha[4]*fUpwind[15]+0.3535533905932737*alpha[2]*fUpwind[9]; 
  Ghat[17] = 0.28284271247461906*alpha[4]*fUpwind[18]+0.3162277660168379*alpha[8]*fUpwind[17]+0.3535533905932737*alpha[0]*fUpwind[17]+0.3535533905932737*alpha[2]*fUpwind[13]+0.3162277660168379*alpha[1]*fUpwind[10]+0.3162277660168379*alpha[4]*fUpwind[5]; 
  Ghat[18] = 0.22587697572631277*alpha[8]*fUpwind[18]+0.3535533905932737*alpha[0]*fUpwind[18]+0.28284271247461906*alpha[4]*fUpwind[17]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3162277660168379*alpha[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[6]; 
  Ghat[19] = 0.3162277660168379*alpha[8]*fUpwind[19]+0.3535533905932737*alpha[0]*fUpwind[19]+0.3535533905932737*alpha[1]*fUpwind[16]+0.3535533905932737*alpha[2]*fUpwind[15]+0.3535533905932737*alpha[4]*fUpwind[9]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -(1.224744871391589*Ghat[0]*dv10); 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -(1.224744871391589*Ghat[1]*dv10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10); 
  out[8] += 0.7071067811865475*Ghat[5]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -(1.224744871391589*Ghat[3]*dv10); 
  out[11] += 0.7071067811865475*Ghat[7]*dv10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 1.5811388300841895*Ghat[0]*dv10; 
  out[14] += 0.7071067811865475*Ghat[9]*dv10; 
  out[15] += -(1.224744871391589*Ghat[4]*dv10); 
  out[16] += 0.7071067811865475*Ghat[10]*dv10; 
  out[17] += -(1.224744871391589*Ghat[5]*dv10); 
  out[18] += -(1.224744871391589*Ghat[6]*dv10); 
  out[19] += 0.7071067811865475*Ghat[11]*dv10; 
  out[20] += 0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -(1.224744871391589*Ghat[7]*dv10); 
  out[22] += -(1.224744871391589*Ghat[8]*dv10); 
  out[23] += 1.5811388300841898*Ghat[1]*dv10; 
  out[24] += 1.5811388300841898*Ghat[2]*dv10; 
  out[25] += 0.7071067811865475*Ghat[13]*dv10; 
  out[26] += 0.7071067811865475*Ghat[14]*dv10; 
  out[27] += 1.5811388300841898*Ghat[3]*dv10; 
  out[28] += 0.7071067811865475*Ghat[15]*dv10; 
  out[29] += 0.7071067811865475*Ghat[16]*dv10; 
  out[30] += -(1.224744871391589*Ghat[9]*dv10); 
  out[31] += -(1.224744871391589*Ghat[10]*dv10); 
  out[32] += -(1.224744871391589*Ghat[11]*dv10); 
  out[33] += -(1.224744871391589*Ghat[12]*dv10); 
  out[34] += 1.5811388300841895*Ghat[4]*dv10; 
  out[35] += 0.7071067811865475*Ghat[17]*dv10; 
  out[36] += 0.7071067811865475*Ghat[18]*dv10; 
  out[37] += -(1.224744871391589*Ghat[13]*dv10); 
  out[38] += -(1.224744871391589*Ghat[14]*dv10); 
  out[39] += 1.5811388300841895*Ghat[5]*dv10; 
  out[40] += 1.5811388300841895*Ghat[6]*dv10; 
  out[41] += 0.7071067811865475*Ghat[19]*dv10; 
  out[42] += -(1.224744871391589*Ghat[15]*dv10); 
  out[43] += -(1.224744871391589*Ghat[16]*dv10); 
  out[44] += -(1.224744871391589*Ghat[17]*dv10); 
  out[45] += -(1.224744871391589*Ghat[18]*dv10); 
  out[46] += 1.5811388300841898*Ghat[10]*dv10; 
  out[47] += -(1.224744871391589*Ghat[19]*dv10); 

  } 
  return 0.;

} 
