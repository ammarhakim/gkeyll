#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_tensor_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_boundary_surfvpar_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/fEdge:   Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in skin cell/last edge cell 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[9]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[9]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[9]; 
  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf[3] = {0.0}; 
  double F_0_UpwindQuad[3] = {0.0};
  double F_0_Upwind[3] = {0.0};;
  double Ghat_F_0[3] = {0.0}; 
  double G_1_UpwindQuad[3] = {0.0};
  double G_1_Upwind[3] = {0.0};;
  double Ghat_G_1[3] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar+0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 

  if (0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(F_0Skin); 
    G_1_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(F_0Edge); 
    G_1_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(G_1Edge); 
  } 
  if (0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2] < 0) { 
    F_0_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(F_0Skin); 
    G_1_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(F_0Edge); 
    G_1_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(G_1Edge); 
  } 
  if (0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(F_0Skin); 
    G_1_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(F_0Edge); 
    G_1_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(G_1Edge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  tensor_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  Ghat_F_0[0] = 0.7071067811865475*F_0_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*F_0_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.6324555320336759*F_0_Upwind[1]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*F_0_Upwind[2]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.4517539514526256*F_0_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*F_0_Upwind[2]+0.6324555320336759*F_0_Upwind[1]*alphaDrSurf[1]; 
  Ghat_G_1[0] = 0.7071067811865475*G_1_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*G_1_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.6324555320336759*G_1_Upwind[1]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*G_1_Upwind[2]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.4517539514526256*G_1_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*G_1_Upwind[2]+0.6324555320336759*G_1_Upwind[1]*alphaDrSurf[1]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[4] += 0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[5] += 1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_F_0[8] += 1.58113883008419*Ghat_F_0[2]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[4] += 0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[5] += 1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 1.58113883008419*Ghat_G_1[1]*dv1par; 
  out_G_1[8] += 1.58113883008419*Ghat_G_1[2]*dv1par; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar-0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 

  if (0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(F_0Edge); 
    G_1_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(F_0Skin); 
    G_1_UpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(G_1Skin); 
  } 
  if (0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2] < 0) { 
    F_0_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(F_0Edge); 
    G_1_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(F_0Skin); 
    G_1_UpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(G_1Skin); 
  } 
  if (0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(F_0Edge); 
    G_1_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(F_0Skin); 
    G_1_UpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(G_1Skin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  tensor_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  Ghat_F_0[0] = 0.7071067811865475*F_0_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*F_0_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.6324555320336759*F_0_Upwind[1]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*F_0_Upwind[2]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.4517539514526256*F_0_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*F_0_Upwind[2]+0.6324555320336759*F_0_Upwind[1]*alphaDrSurf[1]; 
  Ghat_G_1[0] = 0.7071067811865475*G_1_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*G_1_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.6324555320336759*G_1_Upwind[1]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*G_1_Upwind[2]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.4517539514526256*G_1_Upwind[2]*alphaDrSurf[2]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*G_1_Upwind[2]+0.6324555320336759*G_1_Upwind[1]*alphaDrSurf[1]; 

  out_F_0[0] += -0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += -0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[4] += -0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[5] += -1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += -1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_F_0[8] += -1.58113883008419*Ghat_F_0[2]*dv1par; 
  out_G_1[0] += -0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += -0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[4] += -0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[5] += -1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += -1.58113883008419*Ghat_G_1[1]*dv1par; 
  out_G_1[8] += -1.58113883008419*Ghat_G_1[2]*dv1par; 

  }

  return 0.;

} 
