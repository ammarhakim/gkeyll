#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
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

  double alpha[3] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ex = &EBext[0]; 
  alpha[0] = Ex[0]-1.7320508075688772*phi[1]*dx10; 
  alpha[1] = Ex[1]-3.872983346207417*phi[2]*dx10; 
  alpha[2] = Ex[2]; 

  double fUpwindQuad[3] = {0.0};
  double fUpwind[3] = {0.0};
  double Ghat[3] = {0.0}; 

  if (edge == -1) { 

  if (0.6324555320336759*alpha[2]-0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2] > 0) { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (0.6324555320336759*alpha[2]+0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alpha[2]*fUpwind[2]+0.7071067811865475*alpha[1]*fUpwind[1]+0.7071067811865475*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alpha[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alpha[2]+0.7071067811865475*alpha[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.45175395145262565*alpha[2]*fUpwind[2]+0.7071067811865475*alpha[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alpha[2]+0.6324555320336759*alpha[1]*fUpwind[1]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10); 
  out[4] += -(0.7071067811865475*Ghat[2]*dv10); 
  out[5] += -(1.5811388300841895*Ghat[0]*dv10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += -(1.5811388300841898*Ghat[1]*dv10); 

  } else { 

  if (0.6324555320336759*alpha[2]-0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2] > 0) { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (0.6324555320336759*alpha[2]+0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alpha[2]*fUpwind[2]+0.7071067811865475*alpha[1]*fUpwind[1]+0.7071067811865475*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alpha[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alpha[2]+0.7071067811865475*alpha[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.45175395145262565*alpha[2]*fUpwind[2]+0.7071067811865475*alpha[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alpha[2]+0.6324555320336759*alpha[1]*fUpwind[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10); 
  out[4] += 0.7071067811865475*Ghat[2]*dv10; 
  out[5] += 1.5811388300841895*Ghat[0]*dv10; 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += 1.5811388300841898*Ghat[1]*dv10; 

  } 
  return 0.;

} 
