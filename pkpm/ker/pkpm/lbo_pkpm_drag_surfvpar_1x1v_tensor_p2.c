#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_tensor_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_surfvpar_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr:      Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)] in left/center/right cells. 
  // out:           Incremented output distribution functions in center cell. 

  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[9]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[9]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[9]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[9]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[3] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_l[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf_l[2] = nuSum[2]*wvpar-0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 

  double alphaDrSurf_r[3] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_r[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf_r[2] = nuSum[2]*wvpar+0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 

  double F_0_UpwindQuad_l[3] = {0.0};
  double F_0_UpwindQuad_r[3] = {0.0};
  double F_0_Upwind_l[3] = {0.0};
  double F_0_Upwind_r[3] = {0.0};
  double Ghat_F_0_l[3] = {0.0}; 
  double Ghat_F_0_r[3] = {0.0}; 
  double G_1_UpwindQuad_l[3] = {0.0};
  double G_1_UpwindQuad_r[3] = {0.0};
  double G_1_Upwind_l[3] = {0.0};
  double G_1_Upwind_r[3] = {0.0};
  double Ghat_G_1_l[3] = {0.0}; 
  double Ghat_G_1_r[3] = {0.0}; 

  if (0.6324555320336759*alphaDrSurf_l[2]-0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(G_1c); 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]-0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(G_1r); 
  } 
  if (0.7071067811865475*alphaDrSurf_l[0]-0.7905694150420947*alphaDrSurf_l[2] < 0) { 
    F_0_UpwindQuad_l[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(G_1c); 
  } 
  if (0.7071067811865475*alphaDrSurf_r[0]-0.7905694150420947*alphaDrSurf_r[2] < 0) { 
    F_0_UpwindQuad_r[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(G_1r); 
  } 
  if (0.6324555320336759*alphaDrSurf_l[2]+0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(G_1c); 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]+0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  tensor_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  tensor_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  tensor_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  Ghat_F_0_l[0] = 0.7071067811865475*(F_0_Upwind_l[2]*alphaDrSurf_l[2]+F_0_Upwind_l[1]*alphaDrSurf_l[1]+F_0_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 0.6324555320336759*(F_0_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*F_0_Upwind_l[2])+0.7071067811865475*(F_0_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0_Upwind_l[1]); 
  Ghat_F_0_l[2] = 0.4517539514526256*F_0_Upwind_l[2]*alphaDrSurf_l[2]+0.7071067811865475*(F_0_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*F_0_Upwind_l[2])+0.6324555320336759*F_0_Upwind_l[1]*alphaDrSurf_l[1]; 
  Ghat_G_1_l[0] = 0.7071067811865475*(G_1_Upwind_l[2]*alphaDrSurf_l[2]+G_1_Upwind_l[1]*alphaDrSurf_l[1]+G_1_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 0.6324555320336759*(G_1_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*G_1_Upwind_l[2])+0.7071067811865475*(G_1_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1_Upwind_l[1]); 
  Ghat_G_1_l[2] = 0.4517539514526256*G_1_Upwind_l[2]*alphaDrSurf_l[2]+0.7071067811865475*(G_1_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*G_1_Upwind_l[2])+0.6324555320336759*G_1_Upwind_l[1]*alphaDrSurf_l[1]; 

  Ghat_F_0_r[0] = 0.7071067811865475*(F_0_Upwind_r[2]*alphaDrSurf_r[2]+F_0_Upwind_r[1]*alphaDrSurf_r[1]+F_0_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 0.6324555320336759*(F_0_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*F_0_Upwind_r[2])+0.7071067811865475*(F_0_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0_Upwind_r[1]); 
  Ghat_F_0_r[2] = 0.4517539514526256*F_0_Upwind_r[2]*alphaDrSurf_r[2]+0.7071067811865475*(F_0_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*F_0_Upwind_r[2])+0.6324555320336759*F_0_Upwind_r[1]*alphaDrSurf_r[1]; 
  Ghat_G_1_r[0] = 0.7071067811865475*(G_1_Upwind_r[2]*alphaDrSurf_r[2]+G_1_Upwind_r[1]*alphaDrSurf_r[1]+G_1_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 0.6324555320336759*(G_1_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*G_1_Upwind_r[2])+0.7071067811865475*(G_1_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1_Upwind_r[1]); 
  Ghat_G_1_r[2] = 0.4517539514526256*G_1_Upwind_r[2]*alphaDrSurf_r[2]+0.7071067811865475*(G_1_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*G_1_Upwind_r[2])+0.6324555320336759*G_1_Upwind_r[1]*alphaDrSurf_r[1]; 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_r[0]-0.7071067811865475*Ghat_F_0_l[0])*dv1par; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_r[1]-0.7071067811865475*Ghat_F_0_l[1])*dv1par; 
  out_F_0[2] += 1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dv1par; 
  out_F_0[3] += 1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dv1par; 
  out_F_0[4] += (0.7071067811865475*Ghat_F_0_r[2]-0.7071067811865475*Ghat_F_0_l[2])*dv1par; 
  out_F_0[5] += (1.58113883008419*Ghat_F_0_r[0]-1.58113883008419*Ghat_F_0_l[0])*dv1par; 
  out_F_0[6] += 1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dv1par; 
  out_F_0[7] += (1.58113883008419*Ghat_F_0_r[1]-1.58113883008419*Ghat_F_0_l[1])*dv1par; 
  out_F_0[8] += (1.58113883008419*Ghat_F_0_r[2]-1.58113883008419*Ghat_F_0_l[2])*dv1par; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_r[0]-0.7071067811865475*Ghat_G_1_l[0])*dv1par; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_r[1]-0.7071067811865475*Ghat_G_1_l[1])*dv1par; 
  out_G_1[2] += 1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dv1par; 
  out_G_1[3] += 1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dv1par; 
  out_G_1[4] += (0.7071067811865475*Ghat_G_1_r[2]-0.7071067811865475*Ghat_G_1_l[2])*dv1par; 
  out_G_1[5] += (1.58113883008419*Ghat_G_1_r[0]-1.58113883008419*Ghat_G_1_l[0])*dv1par; 
  out_G_1[6] += 1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dv1par; 
  out_G_1[7] += (1.58113883008419*Ghat_G_1_r[1]-1.58113883008419*Ghat_G_1_l[1])*dv1par; 
  out_G_1[8] += (1.58113883008419*Ghat_G_1_r[2]-1.58113883008419*Ghat_G_1_l[2])*dv1par; 

  return 0.;

} 
