#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
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

  const double *phi = &pots[0]; 

  double alpha[8] = {0.0}; 

  alpha[0] = -(2.4494897427831783*phi[1]*dx10); 
  alpha[1] = -(5.477225575051662*phi[2]*dx10); 

  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fEdge); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fEdge); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.4472135954999579*alpha[1]*fUpwind[4]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.5*alpha[1]*fUpwind[3]+0.5*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.44721359549995804*alpha[1]*fUpwind[6]+0.5*alpha[0]*fUpwind[3]+0.5*alpha[1]*fUpwind[2]; 
  Ghat[4] = 0.5*alpha[0]*fUpwind[4]+0.4472135954999579*alpha[1]*fUpwind[1]; 
  Ghat[5] = 0.5000000000000001*alpha[1]*fUpwind[7]+0.5*alpha[0]*fUpwind[5]; 
  Ghat[6] = 0.5*alpha[0]*fUpwind[6]+0.44721359549995804*alpha[1]*fUpwind[3]; 
  Ghat[7] = 0.5*alpha[0]*fUpwind[7]+0.5000000000000001*alpha[1]*fUpwind[5]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10); 
  out[4] += -(1.224744871391589*Ghat[1]*dv10); 
  out[5] += -(0.7071067811865475*Ghat[3]*dv10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += -(0.7071067811865475*Ghat[4]*dv10); 
  out[8] += -(1.5811388300841895*Ghat[0]*dv10); 
  out[9] += -(0.7071067811865475*Ghat[5]*dv10); 
  out[10] += -(1.224744871391589*Ghat[3]*dv10); 
  out[11] += -(1.224744871391589*Ghat[4]*dv10); 
  out[12] += -(1.5811388300841898*Ghat[1]*dv10); 
  out[13] += -(0.7071067811865475*Ghat[6]*dv10); 
  out[14] += -(1.5811388300841898*Ghat[2]*dv10); 
  out[15] += -(0.7071067811865475*Ghat[7]*dv10); 
  out[16] += -(1.224744871391589*Ghat[5]*dv10); 
  out[17] += -(1.224744871391589*Ghat[6]*dv10); 
  out[18] += -(1.5811388300841895*Ghat[3]*dv10); 
  out[19] += -(1.224744871391589*Ghat[7]*dv10); 

  } else { 

  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fSkin); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fSkin); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.4472135954999579*alpha[1]*fUpwind[4]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.5*alpha[1]*fUpwind[3]+0.5*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.44721359549995804*alpha[1]*fUpwind[6]+0.5*alpha[0]*fUpwind[3]+0.5*alpha[1]*fUpwind[2]; 
  Ghat[4] = 0.5*alpha[0]*fUpwind[4]+0.4472135954999579*alpha[1]*fUpwind[1]; 
  Ghat[5] = 0.5000000000000001*alpha[1]*fUpwind[7]+0.5*alpha[0]*fUpwind[5]; 
  Ghat[6] = 0.5*alpha[0]*fUpwind[6]+0.44721359549995804*alpha[1]*fUpwind[3]; 
  Ghat[7] = 0.5*alpha[0]*fUpwind[7]+0.5000000000000001*alpha[1]*fUpwind[5]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -(1.224744871391589*Ghat[1]*dv10); 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += 0.7071067811865475*Ghat[4]*dv10; 
  out[8] += 1.5811388300841895*Ghat[0]*dv10; 
  out[9] += 0.7071067811865475*Ghat[5]*dv10; 
  out[10] += -(1.224744871391589*Ghat[3]*dv10); 
  out[11] += -(1.224744871391589*Ghat[4]*dv10); 
  out[12] += 1.5811388300841898*Ghat[1]*dv10; 
  out[13] += 0.7071067811865475*Ghat[6]*dv10; 
  out[14] += 1.5811388300841898*Ghat[2]*dv10; 
  out[15] += 0.7071067811865475*Ghat[7]*dv10; 
  out[16] += -(1.224744871391589*Ghat[5]*dv10); 
  out[17] += -(1.224744871391589*Ghat[6]*dv10); 
  out[18] += 1.5811388300841895*Ghat[3]*dv10; 
  out[19] += -(1.224744871391589*Ghat[7]*dv10); 

  } 
  return 0.;

} 
