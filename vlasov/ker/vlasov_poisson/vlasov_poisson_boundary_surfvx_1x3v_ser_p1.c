#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // pots:        potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:       external E and B fields (scaled by q/m).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 

  const double dv10 = 2/dxv[1]; 
  const double dx10 = 2/dxv[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *phi = &pots[0]; 

  double alpha[16] = {0.0}; 

  alpha[0] = -(3.4641016151377544*phi[1]*dx10); 

  double fUpwindQuad[18] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fSkin); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fEdge); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fSkin); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fSkin); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fEdge); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fEdge); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.3535533905932737*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.3535533905932737*alpha[0]*fUpwind[4]; 
  Ghat[5] = 0.3535533905932737*alpha[0]*fUpwind[5]; 
  Ghat[6] = 0.3535533905932737*alpha[0]*fUpwind[6]; 
  Ghat[7] = 0.3535533905932737*alpha[0]*fUpwind[7]; 
  Ghat[8] = 0.3535533905932737*alpha[0]*fUpwind[8]; 
  Ghat[9] = 0.3535533905932737*alpha[0]*fUpwind[9]; 
  Ghat[10] = 0.3535533905932737*alpha[0]*fUpwind[10]; 
  Ghat[11] = 0.3535533905932737*alpha[0]*fUpwind[11]; 
  Ghat[12] = 0.3535533905932737*alpha[0]*fUpwind[12]; 
  Ghat[13] = 0.3535533905932737*alpha[0]*fUpwind[13]; 
  Ghat[14] = 0.3535533905932737*alpha[0]*fUpwind[14]; 
  Ghat[15] = 0.3535533905932737*alpha[0]*fUpwind[15]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10); 
  out[5] += -(1.224744871391589*Ghat[1]*dv10); 
  out[6] += -(0.7071067811865475*Ghat[4]*dv10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10); 
  out[8] += -(0.7071067811865475*Ghat[5]*dv10); 
  out[9] += -(1.224744871391589*Ghat[3]*dv10); 
  out[10] += -(0.7071067811865475*Ghat[6]*dv10); 
  out[11] += -(1.224744871391589*Ghat[4]*dv10); 
  out[12] += -(1.224744871391589*Ghat[5]*dv10); 
  out[13] += -(0.7071067811865475*Ghat[7]*dv10); 
  out[14] += -(1.224744871391589*Ghat[6]*dv10); 
  out[15] += -(1.224744871391589*Ghat[7]*dv10); 
  out[16] += -(1.5811388300841895*Ghat[0]*dv10); 
  out[17] += -(1.5811388300841898*Ghat[1]*dv10); 
  out[18] += -(1.5811388300841898*Ghat[2]*dv10); 
  out[19] += -(1.5811388300841898*Ghat[3]*dv10); 
  out[20] += -(1.5811388300841895*Ghat[4]*dv10); 
  out[21] += -(1.5811388300841895*Ghat[5]*dv10); 
  out[22] += -(1.5811388300841895*Ghat[6]*dv10); 
  out[23] += -(1.5811388300841898*Ghat[7]*dv10); 
  out[24] += -(0.7071067811865475*Ghat[8]*dv10); 
  out[25] += -(0.7071067811865475*Ghat[9]*dv10); 
  out[26] += -(1.224744871391589*Ghat[8]*dv10); 
  out[27] += -(0.7071067811865475*Ghat[10]*dv10); 
  out[28] += -(1.224744871391589*Ghat[9]*dv10); 
  out[29] += -(0.7071067811865475*Ghat[11]*dv10); 
  out[30] += -(1.224744871391589*Ghat[10]*dv10); 
  out[31] += -(1.224744871391589*Ghat[11]*dv10); 
  out[32] += -(0.7071067811865475*Ghat[12]*dv10); 
  out[33] += -(0.7071067811865475*Ghat[13]*dv10); 
  out[34] += -(1.224744871391589*Ghat[12]*dv10); 
  out[35] += -(0.7071067811865475*Ghat[14]*dv10); 
  out[36] += -(1.224744871391589*Ghat[13]*dv10); 
  out[37] += -(0.7071067811865475*Ghat[15]*dv10); 
  out[38] += -(1.224744871391589*Ghat[14]*dv10); 
  out[39] += -(1.224744871391589*Ghat[15]*dv10); 

  } else { 

  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fEdge); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fSkin); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fEdge); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fEdge); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fSkin); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fSkin); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.3535533905932737*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.3535533905932737*alpha[0]*fUpwind[4]; 
  Ghat[5] = 0.3535533905932737*alpha[0]*fUpwind[5]; 
  Ghat[6] = 0.3535533905932737*alpha[0]*fUpwind[6]; 
  Ghat[7] = 0.3535533905932737*alpha[0]*fUpwind[7]; 
  Ghat[8] = 0.3535533905932737*alpha[0]*fUpwind[8]; 
  Ghat[9] = 0.3535533905932737*alpha[0]*fUpwind[9]; 
  Ghat[10] = 0.3535533905932737*alpha[0]*fUpwind[10]; 
  Ghat[11] = 0.3535533905932737*alpha[0]*fUpwind[11]; 
  Ghat[12] = 0.3535533905932737*alpha[0]*fUpwind[12]; 
  Ghat[13] = 0.3535533905932737*alpha[0]*fUpwind[13]; 
  Ghat[14] = 0.3535533905932737*alpha[0]*fUpwind[14]; 
  Ghat[15] = 0.3535533905932737*alpha[0]*fUpwind[15]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -(1.224744871391589*Ghat[1]*dv10); 
  out[6] += 0.7071067811865475*Ghat[4]*dv10; 
  out[7] += -(1.224744871391589*Ghat[2]*dv10); 
  out[8] += 0.7071067811865475*Ghat[5]*dv10; 
  out[9] += -(1.224744871391589*Ghat[3]*dv10); 
  out[10] += 0.7071067811865475*Ghat[6]*dv10; 
  out[11] += -(1.224744871391589*Ghat[4]*dv10); 
  out[12] += -(1.224744871391589*Ghat[5]*dv10); 
  out[13] += 0.7071067811865475*Ghat[7]*dv10; 
  out[14] += -(1.224744871391589*Ghat[6]*dv10); 
  out[15] += -(1.224744871391589*Ghat[7]*dv10); 
  out[16] += 1.5811388300841895*Ghat[0]*dv10; 
  out[17] += 1.5811388300841898*Ghat[1]*dv10; 
  out[18] += 1.5811388300841898*Ghat[2]*dv10; 
  out[19] += 1.5811388300841898*Ghat[3]*dv10; 
  out[20] += 1.5811388300841895*Ghat[4]*dv10; 
  out[21] += 1.5811388300841895*Ghat[5]*dv10; 
  out[22] += 1.5811388300841895*Ghat[6]*dv10; 
  out[23] += 1.5811388300841898*Ghat[7]*dv10; 
  out[24] += 0.7071067811865475*Ghat[8]*dv10; 
  out[25] += 0.7071067811865475*Ghat[9]*dv10; 
  out[26] += -(1.224744871391589*Ghat[8]*dv10); 
  out[27] += 0.7071067811865475*Ghat[10]*dv10; 
  out[28] += -(1.224744871391589*Ghat[9]*dv10); 
  out[29] += 0.7071067811865475*Ghat[11]*dv10; 
  out[30] += -(1.224744871391589*Ghat[10]*dv10); 
  out[31] += -(1.224744871391589*Ghat[11]*dv10); 
  out[32] += 0.7071067811865475*Ghat[12]*dv10; 
  out[33] += 0.7071067811865475*Ghat[13]*dv10; 
  out[34] += -(1.224744871391589*Ghat[12]*dv10); 
  out[35] += 0.7071067811865475*Ghat[14]*dv10; 
  out[36] += -(1.224744871391589*Ghat[13]*dv10); 
  out[37] += 0.7071067811865475*Ghat[15]*dv10; 
  out[38] += -(1.224744871391589*Ghat[14]*dv10); 
  out[39] += -(1.224744871391589*Ghat[15]*dv10); 

  } 
  return 0.;

} 
