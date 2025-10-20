#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_2x1v_tensor_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                Cell-center coordinates. 
  // dxv[NDIM]:              Cell spacing. 
  // div_b:                  Input volume expansion of div(b). 
  // pkpm_accel_vars:        Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // g_dist_sourceEdge/Skin: Input [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
  //                         (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]. 
  //                         in skin cell/last edge cell. First input is mirror force source, second input is vperp characteristics source. 
  // edge:                   Determines if the update is for the left edge (-1) or right edge (+1). 
  // fSkin/fEdge:            Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)] in skin cell/last edge cell. 
  // out:                    Incremented output distribution functions in center cell. 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[27]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[27]; 
  const double *F_0_sourceSkin = &fSkin[27]; 
  const double *G_1_sourceSkin = &g_dist_sourceSkin[0]; 
  const double *F_0_sourceEdge = &fEdge[27]; 
  const double *G_1_sourceEdge = &g_dist_sourceEdge[0]; 
  const double *bb_grad_u = &pkpm_accel_vars[9]; 
  const double *p_force = &pkpm_accel_vars[18]; 

  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[27]; 
  double div_b_Surf[9] = {0.0}; 
  div_b_Surf[0] = div_b[0]; 
  div_b_Surf[1] = div_b[1]; 
  div_b_Surf[2] = div_b[2]; 
  div_b_Surf[3] = div_b[3]; 
  div_b_Surf[4] = div_b[4]; 
  div_b_Surf[5] = div_b[5]; 
  div_b_Surf[6] = div_b[6]; 
  div_b_Surf[7] = div_b[7]; 
  div_b_Surf[8] = div_b[8]; 

  double alphaSurf[9] = {0.0}; 
  double F_0_UpwindQuad[9] = {0.0};
  double F_0_Upwind[9] = {0.0};;
  double Ghat_F_0[9] = {0.0}; 
  double G_1_UpwindQuad[9] = {0.0};
  double G_1_Upwind[9] = {0.0};;
  double Ghat_G_1[9] = {0.0}; 
  double F_0_div_b_UpwindQuad[9] = {0.0};
  double F_0_div_b_Upwind[9] = {0.0};;
  double Ghat_F_0_div_b[9] = {0.0}; 
  double G_1_div_b_UpwindQuad[9] = {0.0};
  double G_1_div_b_Upwind[9] = {0.0};;
  double Ghat_G_1_div_b[9] = {0.0}; 

  // get stable timestep of alpha_v = 1/rho (div(p_par b) - p_perp div(b)) - v_par bb : grad(u) 
  // from the quadrature point evaluation needed to compute upwinded distribution functions 
  double cflFreq = 0.0;
  double alphaOrd = 0.0;

  if (edge == -1) { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf[4] = (-1.0*bb_grad_u[4]*wvpar)-0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf[5] = (-1.0*bb_grad_u[5]*wvpar)-0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf[6] = (-1.0*bb_grad_u[6]*wvpar)-0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf[7] = (-1.0*bb_grad_u[7]*wvpar)-0.5*bb_grad_u[7]*dvpar+p_force[7]; 
  alphaSurf[8] = (-1.0*bb_grad_u[8]*wvpar)-0.5*bb_grad_u[8]*dvpar+p_force[8]; 

  alphaOrd = 0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]-0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]-0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0Skin); 
    G_1_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0Edge); 
    G_1_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1Edge); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])+0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])+0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0Skin); 
    G_1_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0Edge); 
    G_1_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1Edge); 
  } 
  alphaOrd = 0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]+0.6708203932499369*alphaSurf[2]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]+0.6708203932499369*alphaSurf[2]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0Skin); 
    G_1_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0Edge); 
    G_1_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1Edge); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])+0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]-0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])+0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]-0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0Skin); 
    G_1_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0Edge); 
    G_1_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1Edge); 
  } 
  alphaOrd = 0.625*alphaSurf[8]-0.5590169943749475*(alphaSurf[5]+alphaSurf[4])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.625*alphaSurf[8]-0.5590169943749475*(alphaSurf[5]+alphaSurf[4])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0Skin); 
    G_1_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0Edge); 
    G_1_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1Edge); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])-0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]+0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])-0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]+0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0Skin); 
    G_1_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0Edge); 
    G_1_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1Edge); 
  } 
  alphaOrd = 0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]-0.6708203932499369*alphaSurf[2]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]-0.6708203932499369*alphaSurf[2]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0Skin); 
    G_1_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0Edge); 
    G_1_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1Edge); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])-0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])-0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0Skin); 
    G_1_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0Edge); 
    G_1_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1Edge); 
  } 
  alphaOrd = 0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]+0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]+0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0Skin); 
    G_1_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0Edge); 
    G_1_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1Edge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  if (0.4*div_b_Surf[8]-0.5999999999999995*div_b_Surf[7]-0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])+0.9*div_b_Surf[3]-0.6708203932499369*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1_sourceEdge); 
  } 
  if ((-0.5*div_b_Surf[8])+0.75*div_b_Surf[7]-0.5590169943749475*div_b_Surf[5]+0.4472135954999579*div_b_Surf[4]-0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1_sourceEdge); 
  } 
  if (0.4*div_b_Surf[8]-0.5999999999999995*div_b_Surf[7]+0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])-0.9*div_b_Surf[3]+0.6708203932499369*div_b_Surf[2]-0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1_sourceEdge); 
  } 
  if ((-0.5*div_b_Surf[8])+0.75*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]-0.5590169943749475*div_b_Surf[4]-0.6708203932499369*div_b_Surf[2]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1_sourceEdge); 
  } 
  if (0.625*div_b_Surf[8]-0.5590169943749475*(div_b_Surf[5]+div_b_Surf[4])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1_sourceEdge); 
  } 
  if ((-0.5*div_b_Surf[8])-0.75*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]-0.5590169943749475*div_b_Surf[4]+0.6708203932499369*div_b_Surf[2]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1_sourceEdge); 
  } 
  if (0.4*div_b_Surf[8]+0.5999999999999995*div_b_Surf[7]-0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])-0.9*div_b_Surf[3]-0.6708203932499369*div_b_Surf[2]+0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1_sourceEdge); 
  } 
  if ((-0.5*div_b_Surf[8])-0.75*div_b_Surf[7]-0.5590169943749475*div_b_Surf[5]+0.4472135954999579*div_b_Surf[4]+0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1_sourceEdge); 
  } 
  if (0.4*div_b_Surf[8]+0.5999999999999995*div_b_Surf[7]+0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])+0.9*div_b_Surf[3]+0.6708203932499369*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1_sourceEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(F_0_div_b_UpwindQuad, F_0_div_b_Upwind); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_div_b_UpwindQuad, G_1_div_b_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[8]*alphaSurf[8]+0.5*F_0_Upwind[7]*alphaSurf[7]+0.5*F_0_Upwind[6]*alphaSurf[6]+0.5*F_0_Upwind[5]*alphaSurf[5]+0.5*F_0_Upwind[4]*alphaSurf[4]+0.5*F_0_Upwind[3]*alphaSurf[3]+0.5*F_0_Upwind[2]*alphaSurf[2]+0.5*F_0_Upwind[1]*alphaSurf[1]+0.5*F_0_Upwind[0]*alphaSurf[0]; 
  Ghat_F_0[1] = 0.447213595499958*F_0_Upwind[7]*alphaSurf[8]+0.447213595499958*alphaSurf[7]*F_0_Upwind[8]+0.5000000000000001*F_0_Upwind[5]*alphaSurf[7]+0.5000000000000001*alphaSurf[5]*F_0_Upwind[7]+0.447213595499958*F_0_Upwind[3]*alphaSurf[6]+0.447213595499958*alphaSurf[3]*F_0_Upwind[6]+0.4472135954999579*F_0_Upwind[1]*alphaSurf[4]+0.4472135954999579*alphaSurf[1]*F_0_Upwind[4]+0.5*F_0_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.447213595499958*F_0_Upwind[6]*alphaSurf[8]+0.447213595499958*alphaSurf[6]*F_0_Upwind[8]+0.447213595499958*F_0_Upwind[3]*alphaSurf[7]+0.447213595499958*alphaSurf[3]*F_0_Upwind[7]+0.5000000000000001*F_0_Upwind[4]*alphaSurf[6]+0.5000000000000001*alphaSurf[4]*F_0_Upwind[6]+0.4472135954999579*F_0_Upwind[2]*alphaSurf[5]+0.4472135954999579*alphaSurf[2]*F_0_Upwind[5]+0.5*F_0_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.4*F_0_Upwind[3]*alphaSurf[8]+0.4*alphaSurf[3]*F_0_Upwind[8]+0.4*F_0_Upwind[6]*alphaSurf[7]+0.447213595499958*F_0_Upwind[2]*alphaSurf[7]+0.4*alphaSurf[6]*F_0_Upwind[7]+0.447213595499958*alphaSurf[2]*F_0_Upwind[7]+0.447213595499958*F_0_Upwind[1]*alphaSurf[6]+0.447213595499958*alphaSurf[1]*F_0_Upwind[6]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[5]+0.4472135954999579*alphaSurf[3]*F_0_Upwind[5]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[4]+0.4472135954999579*alphaSurf[3]*F_0_Upwind[4]+0.5*F_0_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*F_0_Upwind[2]; 
  Ghat_F_0[4] = 0.31943828249997*F_0_Upwind[8]*alphaSurf[8]+0.5*F_0_Upwind[5]*alphaSurf[8]+0.5*alphaSurf[5]*F_0_Upwind[8]+0.4472135954999579*F_0_Upwind[7]*alphaSurf[7]+0.31943828249997*F_0_Upwind[6]*alphaSurf[6]+0.5000000000000001*F_0_Upwind[2]*alphaSurf[6]+0.5000000000000001*alphaSurf[2]*F_0_Upwind[6]+0.31943828249997*F_0_Upwind[4]*alphaSurf[4]+0.5*F_0_Upwind[0]*alphaSurf[4]+0.5*alphaSurf[0]*F_0_Upwind[4]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[3]+0.4472135954999579*F_0_Upwind[1]*alphaSurf[1]; 
  Ghat_F_0[5] = 0.31943828249997*F_0_Upwind[8]*alphaSurf[8]+0.5*F_0_Upwind[4]*alphaSurf[8]+0.5*alphaSurf[4]*F_0_Upwind[8]+0.31943828249997*F_0_Upwind[7]*alphaSurf[7]+0.5000000000000001*F_0_Upwind[1]*alphaSurf[7]+0.5000000000000001*alphaSurf[1]*F_0_Upwind[7]+0.4472135954999579*F_0_Upwind[6]*alphaSurf[6]+0.31943828249997*F_0_Upwind[5]*alphaSurf[5]+0.5*F_0_Upwind[0]*alphaSurf[5]+0.5*alphaSurf[0]*F_0_Upwind[5]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[3]+0.4472135954999579*F_0_Upwind[2]*alphaSurf[2]; 
  Ghat_F_0[6] = 0.2857142857142857*F_0_Upwind[6]*alphaSurf[8]+0.447213595499958*F_0_Upwind[2]*alphaSurf[8]+0.2857142857142857*alphaSurf[6]*F_0_Upwind[8]+0.447213595499958*alphaSurf[2]*F_0_Upwind[8]+0.4*F_0_Upwind[3]*alphaSurf[7]+0.4*alphaSurf[3]*F_0_Upwind[7]+0.4472135954999579*F_0_Upwind[5]*alphaSurf[6]+0.31943828249997*F_0_Upwind[4]*alphaSurf[6]+0.5*F_0_Upwind[0]*alphaSurf[6]+0.4472135954999579*alphaSurf[5]*F_0_Upwind[6]+0.31943828249997*alphaSurf[4]*F_0_Upwind[6]+0.5*alphaSurf[0]*F_0_Upwind[6]+0.5000000000000001*F_0_Upwind[2]*alphaSurf[4]+0.5000000000000001*alphaSurf[2]*F_0_Upwind[4]+0.447213595499958*F_0_Upwind[1]*alphaSurf[3]+0.447213595499958*alphaSurf[1]*F_0_Upwind[3]; 
  Ghat_F_0[7] = 0.2857142857142857*F_0_Upwind[7]*alphaSurf[8]+0.447213595499958*F_0_Upwind[1]*alphaSurf[8]+0.2857142857142857*alphaSurf[7]*F_0_Upwind[8]+0.447213595499958*alphaSurf[1]*F_0_Upwind[8]+0.31943828249997*F_0_Upwind[5]*alphaSurf[7]+0.4472135954999579*F_0_Upwind[4]*alphaSurf[7]+0.5*F_0_Upwind[0]*alphaSurf[7]+0.31943828249997*alphaSurf[5]*F_0_Upwind[7]+0.4472135954999579*alphaSurf[4]*F_0_Upwind[7]+0.5*alphaSurf[0]*F_0_Upwind[7]+0.4*F_0_Upwind[3]*alphaSurf[6]+0.4*alphaSurf[3]*F_0_Upwind[6]+0.5000000000000001*F_0_Upwind[1]*alphaSurf[5]+0.5000000000000001*alphaSurf[1]*F_0_Upwind[5]+0.447213595499958*F_0_Upwind[2]*alphaSurf[3]+0.447213595499958*alphaSurf[2]*F_0_Upwind[3]; 
  Ghat_F_0[8] = 0.2040816326530612*F_0_Upwind[8]*alphaSurf[8]+0.31943828249997*F_0_Upwind[5]*alphaSurf[8]+0.31943828249997*F_0_Upwind[4]*alphaSurf[8]+0.5*F_0_Upwind[0]*alphaSurf[8]+0.31943828249997*alphaSurf[5]*F_0_Upwind[8]+0.31943828249997*alphaSurf[4]*F_0_Upwind[8]+0.5*alphaSurf[0]*F_0_Upwind[8]+0.2857142857142857*F_0_Upwind[7]*alphaSurf[7]+0.447213595499958*F_0_Upwind[1]*alphaSurf[7]+0.447213595499958*alphaSurf[1]*F_0_Upwind[7]+0.2857142857142857*F_0_Upwind[6]*alphaSurf[6]+0.447213595499958*F_0_Upwind[2]*alphaSurf[6]+0.447213595499958*alphaSurf[2]*F_0_Upwind[6]+0.5*F_0_Upwind[4]*alphaSurf[5]+0.5*alphaSurf[4]*F_0_Upwind[5]+0.4*F_0_Upwind[3]*alphaSurf[3]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[8]*alphaSurf[8]+0.5*G_1_Upwind[7]*alphaSurf[7]+0.5*G_1_Upwind[6]*alphaSurf[6]+0.5*G_1_Upwind[5]*alphaSurf[5]+0.5*G_1_Upwind[4]*alphaSurf[4]+0.5*G_1_Upwind[3]*alphaSurf[3]+0.5*G_1_Upwind[2]*alphaSurf[2]+0.5*G_1_Upwind[1]*alphaSurf[1]+0.5*G_1_Upwind[0]*alphaSurf[0]; 
  Ghat_G_1[1] = 0.447213595499958*G_1_Upwind[7]*alphaSurf[8]+0.447213595499958*alphaSurf[7]*G_1_Upwind[8]+0.5000000000000001*G_1_Upwind[5]*alphaSurf[7]+0.5000000000000001*alphaSurf[5]*G_1_Upwind[7]+0.447213595499958*G_1_Upwind[3]*alphaSurf[6]+0.447213595499958*alphaSurf[3]*G_1_Upwind[6]+0.4472135954999579*G_1_Upwind[1]*alphaSurf[4]+0.4472135954999579*alphaSurf[1]*G_1_Upwind[4]+0.5*G_1_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.447213595499958*G_1_Upwind[6]*alphaSurf[8]+0.447213595499958*alphaSurf[6]*G_1_Upwind[8]+0.447213595499958*G_1_Upwind[3]*alphaSurf[7]+0.447213595499958*alphaSurf[3]*G_1_Upwind[7]+0.5000000000000001*G_1_Upwind[4]*alphaSurf[6]+0.5000000000000001*alphaSurf[4]*G_1_Upwind[6]+0.4472135954999579*G_1_Upwind[2]*alphaSurf[5]+0.4472135954999579*alphaSurf[2]*G_1_Upwind[5]+0.5*G_1_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.4*G_1_Upwind[3]*alphaSurf[8]+0.4*alphaSurf[3]*G_1_Upwind[8]+0.4*G_1_Upwind[6]*alphaSurf[7]+0.447213595499958*G_1_Upwind[2]*alphaSurf[7]+0.4*alphaSurf[6]*G_1_Upwind[7]+0.447213595499958*alphaSurf[2]*G_1_Upwind[7]+0.447213595499958*G_1_Upwind[1]*alphaSurf[6]+0.447213595499958*alphaSurf[1]*G_1_Upwind[6]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[5]+0.4472135954999579*alphaSurf[3]*G_1_Upwind[5]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[4]+0.4472135954999579*alphaSurf[3]*G_1_Upwind[4]+0.5*G_1_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*G_1_Upwind[2]; 
  Ghat_G_1[4] = 0.31943828249997*G_1_Upwind[8]*alphaSurf[8]+0.5*G_1_Upwind[5]*alphaSurf[8]+0.5*alphaSurf[5]*G_1_Upwind[8]+0.4472135954999579*G_1_Upwind[7]*alphaSurf[7]+0.31943828249997*G_1_Upwind[6]*alphaSurf[6]+0.5000000000000001*G_1_Upwind[2]*alphaSurf[6]+0.5000000000000001*alphaSurf[2]*G_1_Upwind[6]+0.31943828249997*G_1_Upwind[4]*alphaSurf[4]+0.5*G_1_Upwind[0]*alphaSurf[4]+0.5*alphaSurf[0]*G_1_Upwind[4]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[3]+0.4472135954999579*G_1_Upwind[1]*alphaSurf[1]; 
  Ghat_G_1[5] = 0.31943828249997*G_1_Upwind[8]*alphaSurf[8]+0.5*G_1_Upwind[4]*alphaSurf[8]+0.5*alphaSurf[4]*G_1_Upwind[8]+0.31943828249997*G_1_Upwind[7]*alphaSurf[7]+0.5000000000000001*G_1_Upwind[1]*alphaSurf[7]+0.5000000000000001*alphaSurf[1]*G_1_Upwind[7]+0.4472135954999579*G_1_Upwind[6]*alphaSurf[6]+0.31943828249997*G_1_Upwind[5]*alphaSurf[5]+0.5*G_1_Upwind[0]*alphaSurf[5]+0.5*alphaSurf[0]*G_1_Upwind[5]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[3]+0.4472135954999579*G_1_Upwind[2]*alphaSurf[2]; 
  Ghat_G_1[6] = 0.2857142857142857*G_1_Upwind[6]*alphaSurf[8]+0.447213595499958*G_1_Upwind[2]*alphaSurf[8]+0.2857142857142857*alphaSurf[6]*G_1_Upwind[8]+0.447213595499958*alphaSurf[2]*G_1_Upwind[8]+0.4*G_1_Upwind[3]*alphaSurf[7]+0.4*alphaSurf[3]*G_1_Upwind[7]+0.4472135954999579*G_1_Upwind[5]*alphaSurf[6]+0.31943828249997*G_1_Upwind[4]*alphaSurf[6]+0.5*G_1_Upwind[0]*alphaSurf[6]+0.4472135954999579*alphaSurf[5]*G_1_Upwind[6]+0.31943828249997*alphaSurf[4]*G_1_Upwind[6]+0.5*alphaSurf[0]*G_1_Upwind[6]+0.5000000000000001*G_1_Upwind[2]*alphaSurf[4]+0.5000000000000001*alphaSurf[2]*G_1_Upwind[4]+0.447213595499958*G_1_Upwind[1]*alphaSurf[3]+0.447213595499958*alphaSurf[1]*G_1_Upwind[3]; 
  Ghat_G_1[7] = 0.2857142857142857*G_1_Upwind[7]*alphaSurf[8]+0.447213595499958*G_1_Upwind[1]*alphaSurf[8]+0.2857142857142857*alphaSurf[7]*G_1_Upwind[8]+0.447213595499958*alphaSurf[1]*G_1_Upwind[8]+0.31943828249997*G_1_Upwind[5]*alphaSurf[7]+0.4472135954999579*G_1_Upwind[4]*alphaSurf[7]+0.5*G_1_Upwind[0]*alphaSurf[7]+0.31943828249997*alphaSurf[5]*G_1_Upwind[7]+0.4472135954999579*alphaSurf[4]*G_1_Upwind[7]+0.5*alphaSurf[0]*G_1_Upwind[7]+0.4*G_1_Upwind[3]*alphaSurf[6]+0.4*alphaSurf[3]*G_1_Upwind[6]+0.5000000000000001*G_1_Upwind[1]*alphaSurf[5]+0.5000000000000001*alphaSurf[1]*G_1_Upwind[5]+0.447213595499958*G_1_Upwind[2]*alphaSurf[3]+0.447213595499958*alphaSurf[2]*G_1_Upwind[3]; 
  Ghat_G_1[8] = 0.2040816326530612*G_1_Upwind[8]*alphaSurf[8]+0.31943828249997*G_1_Upwind[5]*alphaSurf[8]+0.31943828249997*G_1_Upwind[4]*alphaSurf[8]+0.5*G_1_Upwind[0]*alphaSurf[8]+0.31943828249997*alphaSurf[5]*G_1_Upwind[8]+0.31943828249997*alphaSurf[4]*G_1_Upwind[8]+0.5*alphaSurf[0]*G_1_Upwind[8]+0.2857142857142857*G_1_Upwind[7]*alphaSurf[7]+0.447213595499958*G_1_Upwind[1]*alphaSurf[7]+0.447213595499958*alphaSurf[1]*G_1_Upwind[7]+0.2857142857142857*G_1_Upwind[6]*alphaSurf[6]+0.447213595499958*G_1_Upwind[2]*alphaSurf[6]+0.447213595499958*alphaSurf[2]*G_1_Upwind[6]+0.5*G_1_Upwind[4]*alphaSurf[5]+0.5*alphaSurf[4]*G_1_Upwind[5]+0.4*G_1_Upwind[3]*alphaSurf[3]; 
  Ghat_F_0_div_b[0] = 0.5*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.5*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.5*F_0_div_b_Upwind[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b[1] = 0.447213595499958*F_0_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*F_0_div_b_Upwind[8]+0.5000000000000001*F_0_div_b_Upwind[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*F_0_div_b_Upwind[7]+0.447213595499958*F_0_div_b_Upwind[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind[6]+0.4472135954999579*F_0_div_b_Upwind[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*F_0_div_b_Upwind[4]+0.5*F_0_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[1]; 
  Ghat_F_0_div_b[2] = 0.447213595499958*F_0_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*F_0_div_b_Upwind[8]+0.447213595499958*F_0_div_b_Upwind[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind[7]+0.5000000000000001*F_0_div_b_Upwind[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*F_0_div_b_Upwind[6]+0.4472135954999579*F_0_div_b_Upwind[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*F_0_div_b_Upwind[5]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[2]; 
  Ghat_F_0_div_b[3] = 0.4*F_0_div_b_Upwind[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*F_0_div_b_Upwind[8]+0.4*F_0_div_b_Upwind[6]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*F_0_div_b_Upwind[7]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[7]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[6]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind[5]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind[4]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[2]; 
  Ghat_F_0_div_b[4] = 0.31943828249997*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*F_0_div_b_Upwind[8]+0.4472135954999579*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.31943828249997*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.5000000000000001*F_0_div_b_Upwind[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind[6]+0.31943828249997*F_0_div_b_Upwind[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[4]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind[1]*div_b_Surf[1]; 
  Ghat_F_0_div_b[5] = 0.31943828249997*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*F_0_div_b_Upwind[8]+0.31943828249997*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.5000000000000001*F_0_div_b_Upwind[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind[7]+0.4472135954999579*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[5]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind[2]*div_b_Surf[2]; 
  Ghat_F_0_div_b[6] = 0.2857142857142857*F_0_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*F_0_div_b_Upwind[8]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[8]+0.4*F_0_div_b_Upwind[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*F_0_div_b_Upwind[7]+0.4472135954999579*F_0_div_b_Upwind[5]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind[4]*div_b_Surf[6]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*F_0_div_b_Upwind[6]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind[6]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[6]+0.5000000000000001*F_0_div_b_Upwind[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind[4]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[3]; 
  Ghat_F_0_div_b[7] = 0.2857142857142857*F_0_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*F_0_div_b_Upwind[8]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[8]+0.31943828249997*F_0_div_b_Upwind[5]*div_b_Surf[7]+0.4472135954999579*F_0_div_b_Upwind[4]*div_b_Surf[7]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind[7]+0.4472135954999579*div_b_Surf[4]*F_0_div_b_Upwind[7]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[7]+0.4*F_0_div_b_Upwind[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*F_0_div_b_Upwind[6]+0.5000000000000001*F_0_div_b_Upwind[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind[5]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[3]; 
  Ghat_F_0_div_b[8] = 0.2040816326530612*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind[5]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind[4]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind[8]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind[8]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[8]+0.2857142857142857*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[7]+0.2857142857142857*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[6]+0.5*F_0_div_b_Upwind[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*F_0_div_b_Upwind[5]+0.4*F_0_div_b_Upwind[3]*div_b_Surf[3]; 
  Ghat_G_1_div_b[0] = 0.5*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.5*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.5*G_1_div_b_Upwind[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b[1] = 0.447213595499958*G_1_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*G_1_div_b_Upwind[8]+0.5000000000000001*G_1_div_b_Upwind[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*G_1_div_b_Upwind[7]+0.447213595499958*G_1_div_b_Upwind[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind[6]+0.4472135954999579*G_1_div_b_Upwind[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*G_1_div_b_Upwind[4]+0.5*G_1_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[1]; 
  Ghat_G_1_div_b[2] = 0.447213595499958*G_1_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*G_1_div_b_Upwind[8]+0.447213595499958*G_1_div_b_Upwind[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind[7]+0.5000000000000001*G_1_div_b_Upwind[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*G_1_div_b_Upwind[6]+0.4472135954999579*G_1_div_b_Upwind[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*G_1_div_b_Upwind[5]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[2]; 
  Ghat_G_1_div_b[3] = 0.4*G_1_div_b_Upwind[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*G_1_div_b_Upwind[8]+0.4*G_1_div_b_Upwind[6]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*G_1_div_b_Upwind[7]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[7]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[6]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind[5]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind[4]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[2]; 
  Ghat_G_1_div_b[4] = 0.31943828249997*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*G_1_div_b_Upwind[8]+0.4472135954999579*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.31943828249997*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.5000000000000001*G_1_div_b_Upwind[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind[6]+0.31943828249997*G_1_div_b_Upwind[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[4]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind[1]*div_b_Surf[1]; 
  Ghat_G_1_div_b[5] = 0.31943828249997*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*G_1_div_b_Upwind[8]+0.31943828249997*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.5000000000000001*G_1_div_b_Upwind[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind[7]+0.4472135954999579*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[5]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind[2]*div_b_Surf[2]; 
  Ghat_G_1_div_b[6] = 0.2857142857142857*G_1_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*G_1_div_b_Upwind[8]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[8]+0.4*G_1_div_b_Upwind[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*G_1_div_b_Upwind[7]+0.4472135954999579*G_1_div_b_Upwind[5]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind[4]*div_b_Surf[6]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*G_1_div_b_Upwind[6]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind[6]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[6]+0.5000000000000001*G_1_div_b_Upwind[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind[4]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[3]; 
  Ghat_G_1_div_b[7] = 0.2857142857142857*G_1_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*G_1_div_b_Upwind[8]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[8]+0.31943828249997*G_1_div_b_Upwind[5]*div_b_Surf[7]+0.4472135954999579*G_1_div_b_Upwind[4]*div_b_Surf[7]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind[7]+0.4472135954999579*div_b_Surf[4]*G_1_div_b_Upwind[7]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[7]+0.4*G_1_div_b_Upwind[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*G_1_div_b_Upwind[6]+0.5000000000000001*G_1_div_b_Upwind[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind[5]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[3]; 
  Ghat_G_1_div_b[8] = 0.2040816326530612*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind[5]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind[4]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind[8]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind[8]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[8]+0.2857142857142857*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[7]+0.2857142857142857*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[6]+0.5*G_1_div_b_Upwind[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*G_1_div_b_Upwind[5]+0.4*G_1_div_b_Upwind[3]*div_b_Surf[3]; 

  out_F_0[0] += (-0.7071067811865475*Ghat_F_0_div_b[0]*dv1par)-0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += (-0.7071067811865475*Ghat_F_0_div_b[1]*dv1par)-0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += (-0.7071067811865475*Ghat_F_0_div_b[2]*dv1par)-0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += (-1.224744871391589*Ghat_F_0_div_b[0]*dv1par)-1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += (-0.7071067811865475*Ghat_F_0_div_b[3]*dv1par)-0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += (-1.224744871391589*Ghat_F_0_div_b[1]*dv1par)-1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += (-1.224744871391589*Ghat_F_0_div_b[2]*dv1par)-1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += (-0.7071067811865475*Ghat_F_0_div_b[4]*dv1par)-0.7071067811865475*Ghat_F_0[4]*dv1par; 
  out_F_0[8] += (-0.7071067811865475*Ghat_F_0_div_b[5]*dv1par)-0.7071067811865475*Ghat_F_0[5]*dv1par; 
  out_F_0[9] += (-1.58113883008419*Ghat_F_0_div_b[0]*dv1par)-1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[10] += (-1.224744871391589*Ghat_F_0_div_b[3]*dv1par)-1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[11] += (-0.7071067811865475*Ghat_F_0_div_b[6]*dv1par)-0.7071067811865475*Ghat_F_0[6]*dv1par; 
  out_F_0[12] += (-0.7071067811865475*Ghat_F_0_div_b[7]*dv1par)-0.7071067811865475*Ghat_F_0[7]*dv1par; 
  out_F_0[13] += (-1.224744871391589*Ghat_F_0_div_b[4]*dv1par)-1.224744871391589*Ghat_F_0[4]*dv1par; 
  out_F_0[14] += (-1.224744871391589*Ghat_F_0_div_b[5]*dv1par)-1.224744871391589*Ghat_F_0[5]*dv1par; 
  out_F_0[15] += (-1.58113883008419*Ghat_F_0_div_b[1]*dv1par)-1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_F_0[16] += (-1.58113883008419*Ghat_F_0_div_b[2]*dv1par)-1.58113883008419*Ghat_F_0[2]*dv1par; 
  out_F_0[17] += (-1.224744871391589*Ghat_F_0_div_b[6]*dv1par)-1.224744871391589*Ghat_F_0[6]*dv1par; 
  out_F_0[18] += (-1.224744871391589*Ghat_F_0_div_b[7]*dv1par)-1.224744871391589*Ghat_F_0[7]*dv1par; 
  out_F_0[19] += (-1.58113883008419*Ghat_F_0_div_b[3]*dv1par)-1.58113883008419*Ghat_F_0[3]*dv1par; 
  out_F_0[20] += (-0.7071067811865475*Ghat_F_0_div_b[8]*dv1par)-0.7071067811865475*Ghat_F_0[8]*dv1par; 
  out_F_0[21] += (-1.58113883008419*Ghat_F_0_div_b[4]*dv1par)-1.58113883008419*Ghat_F_0[4]*dv1par; 
  out_F_0[22] += (-1.58113883008419*Ghat_F_0_div_b[5]*dv1par)-1.58113883008419*Ghat_F_0[5]*dv1par; 
  out_F_0[23] += (-1.224744871391589*Ghat_F_0_div_b[8]*dv1par)-1.224744871391589*Ghat_F_0[8]*dv1par; 
  out_F_0[24] += (-1.58113883008419*Ghat_F_0_div_b[6]*dv1par)-1.58113883008419*Ghat_F_0[6]*dv1par; 
  out_F_0[25] += (-1.58113883008419*Ghat_F_0_div_b[7]*dv1par)-1.58113883008419*Ghat_F_0[7]*dv1par; 
  out_F_0[26] += (-1.58113883008419*Ghat_F_0_div_b[8]*dv1par)-1.58113883008419*Ghat_F_0[8]*dv1par; 
  out_G_1[0] += (-0.7071067811865475*Ghat_G_1_div_b[0]*dv1par)-0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += (-0.7071067811865475*Ghat_G_1_div_b[1]*dv1par)-0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += (-0.7071067811865475*Ghat_G_1_div_b[2]*dv1par)-0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += (-1.224744871391589*Ghat_G_1_div_b[0]*dv1par)-1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += (-0.7071067811865475*Ghat_G_1_div_b[3]*dv1par)-0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += (-1.224744871391589*Ghat_G_1_div_b[1]*dv1par)-1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += (-1.224744871391589*Ghat_G_1_div_b[2]*dv1par)-1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += (-0.7071067811865475*Ghat_G_1_div_b[4]*dv1par)-0.7071067811865475*Ghat_G_1[4]*dv1par; 
  out_G_1[8] += (-0.7071067811865475*Ghat_G_1_div_b[5]*dv1par)-0.7071067811865475*Ghat_G_1[5]*dv1par; 
  out_G_1[9] += (-1.58113883008419*Ghat_G_1_div_b[0]*dv1par)-1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[10] += (-1.224744871391589*Ghat_G_1_div_b[3]*dv1par)-1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[11] += (-0.7071067811865475*Ghat_G_1_div_b[6]*dv1par)-0.7071067811865475*Ghat_G_1[6]*dv1par; 
  out_G_1[12] += (-0.7071067811865475*Ghat_G_1_div_b[7]*dv1par)-0.7071067811865475*Ghat_G_1[7]*dv1par; 
  out_G_1[13] += (-1.224744871391589*Ghat_G_1_div_b[4]*dv1par)-1.224744871391589*Ghat_G_1[4]*dv1par; 
  out_G_1[14] += (-1.224744871391589*Ghat_G_1_div_b[5]*dv1par)-1.224744871391589*Ghat_G_1[5]*dv1par; 
  out_G_1[15] += (-1.58113883008419*Ghat_G_1_div_b[1]*dv1par)-1.58113883008419*Ghat_G_1[1]*dv1par; 
  out_G_1[16] += (-1.58113883008419*Ghat_G_1_div_b[2]*dv1par)-1.58113883008419*Ghat_G_1[2]*dv1par; 
  out_G_1[17] += (-1.224744871391589*Ghat_G_1_div_b[6]*dv1par)-1.224744871391589*Ghat_G_1[6]*dv1par; 
  out_G_1[18] += (-1.224744871391589*Ghat_G_1_div_b[7]*dv1par)-1.224744871391589*Ghat_G_1[7]*dv1par; 
  out_G_1[19] += (-1.58113883008419*Ghat_G_1_div_b[3]*dv1par)-1.58113883008419*Ghat_G_1[3]*dv1par; 
  out_G_1[20] += (-0.7071067811865475*Ghat_G_1_div_b[8]*dv1par)-0.7071067811865475*Ghat_G_1[8]*dv1par; 
  out_G_1[21] += (-1.58113883008419*Ghat_G_1_div_b[4]*dv1par)-1.58113883008419*Ghat_G_1[4]*dv1par; 
  out_G_1[22] += (-1.58113883008419*Ghat_G_1_div_b[5]*dv1par)-1.58113883008419*Ghat_G_1[5]*dv1par; 
  out_G_1[23] += (-1.224744871391589*Ghat_G_1_div_b[8]*dv1par)-1.224744871391589*Ghat_G_1[8]*dv1par; 
  out_G_1[24] += (-1.58113883008419*Ghat_G_1_div_b[6]*dv1par)-1.58113883008419*Ghat_G_1[6]*dv1par; 
  out_G_1[25] += (-1.58113883008419*Ghat_G_1_div_b[7]*dv1par)-1.58113883008419*Ghat_G_1[7]*dv1par; 
  out_G_1[26] += (-1.58113883008419*Ghat_G_1_div_b[8]*dv1par)-1.58113883008419*Ghat_G_1[8]*dv1par; 

  } else { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf[4] = (-1.0*bb_grad_u[4]*wvpar)+0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf[5] = (-1.0*bb_grad_u[5]*wvpar)+0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf[6] = (-1.0*bb_grad_u[6]*wvpar)+0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf[7] = (-1.0*bb_grad_u[7]*wvpar)+0.5*bb_grad_u[7]*dvpar+p_force[7]; 
  alphaSurf[8] = (-1.0*bb_grad_u[8]*wvpar)+0.5*bb_grad_u[8]*dvpar+p_force[8]; 

  alphaOrd = 0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]-0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]-0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0Edge); 
    G_1_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0Skin); 
    G_1_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1Skin); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])+0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])+0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0Edge); 
    G_1_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0Skin); 
    G_1_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1Skin); 
  } 
  alphaOrd = 0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]+0.6708203932499369*alphaSurf[2]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]-0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]+0.6708203932499369*alphaSurf[2]-0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0Edge); 
    G_1_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0Skin); 
    G_1_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1Skin); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])+0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]-0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])+0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]-0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0Edge); 
    G_1_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0Skin); 
    G_1_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1Skin); 
  } 
  alphaOrd = 0.625*alphaSurf[8]-0.5590169943749475*(alphaSurf[5]+alphaSurf[4])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.625*alphaSurf[8]-0.5590169943749475*(alphaSurf[5]+alphaSurf[4])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0Edge); 
    G_1_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0Skin); 
    G_1_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1Skin); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])-0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]+0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])-0.75*alphaSurf[6]+0.4472135954999579*alphaSurf[5]-0.5590169943749475*alphaSurf[4]+0.6708203932499369*alphaSurf[2]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0Edge); 
    G_1_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0Skin); 
    G_1_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1Skin); 
  } 
  alphaOrd = 0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]-0.6708203932499369*alphaSurf[2]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]-0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])-0.9*alphaSurf[3]-0.6708203932499369*alphaSurf[2]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0Edge); 
    G_1_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0Skin); 
    G_1_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1Skin); 
  } 
  alphaOrd = (-0.5*alphaSurf[8])-0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if ((-0.5*alphaSurf[8])-0.75*alphaSurf[7]-0.5590169943749475*alphaSurf[5]+0.4472135954999579*alphaSurf[4]+0.6708203932499369*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0Edge); 
    G_1_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0Skin); 
    G_1_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1Skin); 
  } 
  alphaOrd = 0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]+0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.4*alphaSurf[8]+0.5999999999999995*alphaSurf[7]+0.5999999999999999*alphaSurf[6]+0.4472135954999579*(alphaSurf[5]+alphaSurf[4])+0.9*alphaSurf[3]+0.6708203932499369*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0Edge); 
    G_1_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0Skin); 
    G_1_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1Skin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  if (0.4*div_b_Surf[8]-0.5999999999999995*div_b_Surf[7]-0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])+0.9*div_b_Surf[3]-0.6708203932499369*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1_sourceSkin); 
  } 
  if ((-0.5*div_b_Surf[8])+0.75*div_b_Surf[7]-0.5590169943749475*div_b_Surf[5]+0.4472135954999579*div_b_Surf[4]-0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1_sourceSkin); 
  } 
  if (0.4*div_b_Surf[8]-0.5999999999999995*div_b_Surf[7]+0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])-0.9*div_b_Surf[3]+0.6708203932499369*div_b_Surf[2]-0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1_sourceSkin); 
  } 
  if ((-0.5*div_b_Surf[8])+0.75*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]-0.5590169943749475*div_b_Surf[4]-0.6708203932499369*div_b_Surf[2]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1_sourceSkin); 
  } 
  if (0.625*div_b_Surf[8]-0.5590169943749475*(div_b_Surf[5]+div_b_Surf[4])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1_sourceSkin); 
  } 
  if ((-0.5*div_b_Surf[8])-0.75*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]-0.5590169943749475*div_b_Surf[4]+0.6708203932499369*div_b_Surf[2]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1_sourceSkin); 
  } 
  if (0.4*div_b_Surf[8]+0.5999999999999995*div_b_Surf[7]-0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])-0.9*div_b_Surf[3]-0.6708203932499369*div_b_Surf[2]+0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1_sourceSkin); 
  } 
  if ((-0.5*div_b_Surf[8])-0.75*div_b_Surf[7]-0.5590169943749475*div_b_Surf[5]+0.4472135954999579*div_b_Surf[4]+0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1_sourceSkin); 
  } 
  if (0.4*div_b_Surf[8]+0.5999999999999995*div_b_Surf[7]+0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])+0.9*div_b_Surf[3]+0.6708203932499369*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1_sourceSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(F_0_div_b_UpwindQuad, F_0_div_b_Upwind); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_div_b_UpwindQuad, G_1_div_b_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[8]*alphaSurf[8]+0.5*F_0_Upwind[7]*alphaSurf[7]+0.5*F_0_Upwind[6]*alphaSurf[6]+0.5*F_0_Upwind[5]*alphaSurf[5]+0.5*F_0_Upwind[4]*alphaSurf[4]+0.5*F_0_Upwind[3]*alphaSurf[3]+0.5*F_0_Upwind[2]*alphaSurf[2]+0.5*F_0_Upwind[1]*alphaSurf[1]+0.5*F_0_Upwind[0]*alphaSurf[0]; 
  Ghat_F_0[1] = 0.447213595499958*F_0_Upwind[7]*alphaSurf[8]+0.447213595499958*alphaSurf[7]*F_0_Upwind[8]+0.5000000000000001*F_0_Upwind[5]*alphaSurf[7]+0.5000000000000001*alphaSurf[5]*F_0_Upwind[7]+0.447213595499958*F_0_Upwind[3]*alphaSurf[6]+0.447213595499958*alphaSurf[3]*F_0_Upwind[6]+0.4472135954999579*F_0_Upwind[1]*alphaSurf[4]+0.4472135954999579*alphaSurf[1]*F_0_Upwind[4]+0.5*F_0_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.447213595499958*F_0_Upwind[6]*alphaSurf[8]+0.447213595499958*alphaSurf[6]*F_0_Upwind[8]+0.447213595499958*F_0_Upwind[3]*alphaSurf[7]+0.447213595499958*alphaSurf[3]*F_0_Upwind[7]+0.5000000000000001*F_0_Upwind[4]*alphaSurf[6]+0.5000000000000001*alphaSurf[4]*F_0_Upwind[6]+0.4472135954999579*F_0_Upwind[2]*alphaSurf[5]+0.4472135954999579*alphaSurf[2]*F_0_Upwind[5]+0.5*F_0_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.4*F_0_Upwind[3]*alphaSurf[8]+0.4*alphaSurf[3]*F_0_Upwind[8]+0.4*F_0_Upwind[6]*alphaSurf[7]+0.447213595499958*F_0_Upwind[2]*alphaSurf[7]+0.4*alphaSurf[6]*F_0_Upwind[7]+0.447213595499958*alphaSurf[2]*F_0_Upwind[7]+0.447213595499958*F_0_Upwind[1]*alphaSurf[6]+0.447213595499958*alphaSurf[1]*F_0_Upwind[6]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[5]+0.4472135954999579*alphaSurf[3]*F_0_Upwind[5]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[4]+0.4472135954999579*alphaSurf[3]*F_0_Upwind[4]+0.5*F_0_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*F_0_Upwind[2]; 
  Ghat_F_0[4] = 0.31943828249997*F_0_Upwind[8]*alphaSurf[8]+0.5*F_0_Upwind[5]*alphaSurf[8]+0.5*alphaSurf[5]*F_0_Upwind[8]+0.4472135954999579*F_0_Upwind[7]*alphaSurf[7]+0.31943828249997*F_0_Upwind[6]*alphaSurf[6]+0.5000000000000001*F_0_Upwind[2]*alphaSurf[6]+0.5000000000000001*alphaSurf[2]*F_0_Upwind[6]+0.31943828249997*F_0_Upwind[4]*alphaSurf[4]+0.5*F_0_Upwind[0]*alphaSurf[4]+0.5*alphaSurf[0]*F_0_Upwind[4]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[3]+0.4472135954999579*F_0_Upwind[1]*alphaSurf[1]; 
  Ghat_F_0[5] = 0.31943828249997*F_0_Upwind[8]*alphaSurf[8]+0.5*F_0_Upwind[4]*alphaSurf[8]+0.5*alphaSurf[4]*F_0_Upwind[8]+0.31943828249997*F_0_Upwind[7]*alphaSurf[7]+0.5000000000000001*F_0_Upwind[1]*alphaSurf[7]+0.5000000000000001*alphaSurf[1]*F_0_Upwind[7]+0.4472135954999579*F_0_Upwind[6]*alphaSurf[6]+0.31943828249997*F_0_Upwind[5]*alphaSurf[5]+0.5*F_0_Upwind[0]*alphaSurf[5]+0.5*alphaSurf[0]*F_0_Upwind[5]+0.4472135954999579*F_0_Upwind[3]*alphaSurf[3]+0.4472135954999579*F_0_Upwind[2]*alphaSurf[2]; 
  Ghat_F_0[6] = 0.2857142857142857*F_0_Upwind[6]*alphaSurf[8]+0.447213595499958*F_0_Upwind[2]*alphaSurf[8]+0.2857142857142857*alphaSurf[6]*F_0_Upwind[8]+0.447213595499958*alphaSurf[2]*F_0_Upwind[8]+0.4*F_0_Upwind[3]*alphaSurf[7]+0.4*alphaSurf[3]*F_0_Upwind[7]+0.4472135954999579*F_0_Upwind[5]*alphaSurf[6]+0.31943828249997*F_0_Upwind[4]*alphaSurf[6]+0.5*F_0_Upwind[0]*alphaSurf[6]+0.4472135954999579*alphaSurf[5]*F_0_Upwind[6]+0.31943828249997*alphaSurf[4]*F_0_Upwind[6]+0.5*alphaSurf[0]*F_0_Upwind[6]+0.5000000000000001*F_0_Upwind[2]*alphaSurf[4]+0.5000000000000001*alphaSurf[2]*F_0_Upwind[4]+0.447213595499958*F_0_Upwind[1]*alphaSurf[3]+0.447213595499958*alphaSurf[1]*F_0_Upwind[3]; 
  Ghat_F_0[7] = 0.2857142857142857*F_0_Upwind[7]*alphaSurf[8]+0.447213595499958*F_0_Upwind[1]*alphaSurf[8]+0.2857142857142857*alphaSurf[7]*F_0_Upwind[8]+0.447213595499958*alphaSurf[1]*F_0_Upwind[8]+0.31943828249997*F_0_Upwind[5]*alphaSurf[7]+0.4472135954999579*F_0_Upwind[4]*alphaSurf[7]+0.5*F_0_Upwind[0]*alphaSurf[7]+0.31943828249997*alphaSurf[5]*F_0_Upwind[7]+0.4472135954999579*alphaSurf[4]*F_0_Upwind[7]+0.5*alphaSurf[0]*F_0_Upwind[7]+0.4*F_0_Upwind[3]*alphaSurf[6]+0.4*alphaSurf[3]*F_0_Upwind[6]+0.5000000000000001*F_0_Upwind[1]*alphaSurf[5]+0.5000000000000001*alphaSurf[1]*F_0_Upwind[5]+0.447213595499958*F_0_Upwind[2]*alphaSurf[3]+0.447213595499958*alphaSurf[2]*F_0_Upwind[3]; 
  Ghat_F_0[8] = 0.2040816326530612*F_0_Upwind[8]*alphaSurf[8]+0.31943828249997*F_0_Upwind[5]*alphaSurf[8]+0.31943828249997*F_0_Upwind[4]*alphaSurf[8]+0.5*F_0_Upwind[0]*alphaSurf[8]+0.31943828249997*alphaSurf[5]*F_0_Upwind[8]+0.31943828249997*alphaSurf[4]*F_0_Upwind[8]+0.5*alphaSurf[0]*F_0_Upwind[8]+0.2857142857142857*F_0_Upwind[7]*alphaSurf[7]+0.447213595499958*F_0_Upwind[1]*alphaSurf[7]+0.447213595499958*alphaSurf[1]*F_0_Upwind[7]+0.2857142857142857*F_0_Upwind[6]*alphaSurf[6]+0.447213595499958*F_0_Upwind[2]*alphaSurf[6]+0.447213595499958*alphaSurf[2]*F_0_Upwind[6]+0.5*F_0_Upwind[4]*alphaSurf[5]+0.5*alphaSurf[4]*F_0_Upwind[5]+0.4*F_0_Upwind[3]*alphaSurf[3]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[8]*alphaSurf[8]+0.5*G_1_Upwind[7]*alphaSurf[7]+0.5*G_1_Upwind[6]*alphaSurf[6]+0.5*G_1_Upwind[5]*alphaSurf[5]+0.5*G_1_Upwind[4]*alphaSurf[4]+0.5*G_1_Upwind[3]*alphaSurf[3]+0.5*G_1_Upwind[2]*alphaSurf[2]+0.5*G_1_Upwind[1]*alphaSurf[1]+0.5*G_1_Upwind[0]*alphaSurf[0]; 
  Ghat_G_1[1] = 0.447213595499958*G_1_Upwind[7]*alphaSurf[8]+0.447213595499958*alphaSurf[7]*G_1_Upwind[8]+0.5000000000000001*G_1_Upwind[5]*alphaSurf[7]+0.5000000000000001*alphaSurf[5]*G_1_Upwind[7]+0.447213595499958*G_1_Upwind[3]*alphaSurf[6]+0.447213595499958*alphaSurf[3]*G_1_Upwind[6]+0.4472135954999579*G_1_Upwind[1]*alphaSurf[4]+0.4472135954999579*alphaSurf[1]*G_1_Upwind[4]+0.5*G_1_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.447213595499958*G_1_Upwind[6]*alphaSurf[8]+0.447213595499958*alphaSurf[6]*G_1_Upwind[8]+0.447213595499958*G_1_Upwind[3]*alphaSurf[7]+0.447213595499958*alphaSurf[3]*G_1_Upwind[7]+0.5000000000000001*G_1_Upwind[4]*alphaSurf[6]+0.5000000000000001*alphaSurf[4]*G_1_Upwind[6]+0.4472135954999579*G_1_Upwind[2]*alphaSurf[5]+0.4472135954999579*alphaSurf[2]*G_1_Upwind[5]+0.5*G_1_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.4*G_1_Upwind[3]*alphaSurf[8]+0.4*alphaSurf[3]*G_1_Upwind[8]+0.4*G_1_Upwind[6]*alphaSurf[7]+0.447213595499958*G_1_Upwind[2]*alphaSurf[7]+0.4*alphaSurf[6]*G_1_Upwind[7]+0.447213595499958*alphaSurf[2]*G_1_Upwind[7]+0.447213595499958*G_1_Upwind[1]*alphaSurf[6]+0.447213595499958*alphaSurf[1]*G_1_Upwind[6]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[5]+0.4472135954999579*alphaSurf[3]*G_1_Upwind[5]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[4]+0.4472135954999579*alphaSurf[3]*G_1_Upwind[4]+0.5*G_1_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*G_1_Upwind[2]; 
  Ghat_G_1[4] = 0.31943828249997*G_1_Upwind[8]*alphaSurf[8]+0.5*G_1_Upwind[5]*alphaSurf[8]+0.5*alphaSurf[5]*G_1_Upwind[8]+0.4472135954999579*G_1_Upwind[7]*alphaSurf[7]+0.31943828249997*G_1_Upwind[6]*alphaSurf[6]+0.5000000000000001*G_1_Upwind[2]*alphaSurf[6]+0.5000000000000001*alphaSurf[2]*G_1_Upwind[6]+0.31943828249997*G_1_Upwind[4]*alphaSurf[4]+0.5*G_1_Upwind[0]*alphaSurf[4]+0.5*alphaSurf[0]*G_1_Upwind[4]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[3]+0.4472135954999579*G_1_Upwind[1]*alphaSurf[1]; 
  Ghat_G_1[5] = 0.31943828249997*G_1_Upwind[8]*alphaSurf[8]+0.5*G_1_Upwind[4]*alphaSurf[8]+0.5*alphaSurf[4]*G_1_Upwind[8]+0.31943828249997*G_1_Upwind[7]*alphaSurf[7]+0.5000000000000001*G_1_Upwind[1]*alphaSurf[7]+0.5000000000000001*alphaSurf[1]*G_1_Upwind[7]+0.4472135954999579*G_1_Upwind[6]*alphaSurf[6]+0.31943828249997*G_1_Upwind[5]*alphaSurf[5]+0.5*G_1_Upwind[0]*alphaSurf[5]+0.5*alphaSurf[0]*G_1_Upwind[5]+0.4472135954999579*G_1_Upwind[3]*alphaSurf[3]+0.4472135954999579*G_1_Upwind[2]*alphaSurf[2]; 
  Ghat_G_1[6] = 0.2857142857142857*G_1_Upwind[6]*alphaSurf[8]+0.447213595499958*G_1_Upwind[2]*alphaSurf[8]+0.2857142857142857*alphaSurf[6]*G_1_Upwind[8]+0.447213595499958*alphaSurf[2]*G_1_Upwind[8]+0.4*G_1_Upwind[3]*alphaSurf[7]+0.4*alphaSurf[3]*G_1_Upwind[7]+0.4472135954999579*G_1_Upwind[5]*alphaSurf[6]+0.31943828249997*G_1_Upwind[4]*alphaSurf[6]+0.5*G_1_Upwind[0]*alphaSurf[6]+0.4472135954999579*alphaSurf[5]*G_1_Upwind[6]+0.31943828249997*alphaSurf[4]*G_1_Upwind[6]+0.5*alphaSurf[0]*G_1_Upwind[6]+0.5000000000000001*G_1_Upwind[2]*alphaSurf[4]+0.5000000000000001*alphaSurf[2]*G_1_Upwind[4]+0.447213595499958*G_1_Upwind[1]*alphaSurf[3]+0.447213595499958*alphaSurf[1]*G_1_Upwind[3]; 
  Ghat_G_1[7] = 0.2857142857142857*G_1_Upwind[7]*alphaSurf[8]+0.447213595499958*G_1_Upwind[1]*alphaSurf[8]+0.2857142857142857*alphaSurf[7]*G_1_Upwind[8]+0.447213595499958*alphaSurf[1]*G_1_Upwind[8]+0.31943828249997*G_1_Upwind[5]*alphaSurf[7]+0.4472135954999579*G_1_Upwind[4]*alphaSurf[7]+0.5*G_1_Upwind[0]*alphaSurf[7]+0.31943828249997*alphaSurf[5]*G_1_Upwind[7]+0.4472135954999579*alphaSurf[4]*G_1_Upwind[7]+0.5*alphaSurf[0]*G_1_Upwind[7]+0.4*G_1_Upwind[3]*alphaSurf[6]+0.4*alphaSurf[3]*G_1_Upwind[6]+0.5000000000000001*G_1_Upwind[1]*alphaSurf[5]+0.5000000000000001*alphaSurf[1]*G_1_Upwind[5]+0.447213595499958*G_1_Upwind[2]*alphaSurf[3]+0.447213595499958*alphaSurf[2]*G_1_Upwind[3]; 
  Ghat_G_1[8] = 0.2040816326530612*G_1_Upwind[8]*alphaSurf[8]+0.31943828249997*G_1_Upwind[5]*alphaSurf[8]+0.31943828249997*G_1_Upwind[4]*alphaSurf[8]+0.5*G_1_Upwind[0]*alphaSurf[8]+0.31943828249997*alphaSurf[5]*G_1_Upwind[8]+0.31943828249997*alphaSurf[4]*G_1_Upwind[8]+0.5*alphaSurf[0]*G_1_Upwind[8]+0.2857142857142857*G_1_Upwind[7]*alphaSurf[7]+0.447213595499958*G_1_Upwind[1]*alphaSurf[7]+0.447213595499958*alphaSurf[1]*G_1_Upwind[7]+0.2857142857142857*G_1_Upwind[6]*alphaSurf[6]+0.447213595499958*G_1_Upwind[2]*alphaSurf[6]+0.447213595499958*alphaSurf[2]*G_1_Upwind[6]+0.5*G_1_Upwind[4]*alphaSurf[5]+0.5*alphaSurf[4]*G_1_Upwind[5]+0.4*G_1_Upwind[3]*alphaSurf[3]; 
  Ghat_F_0_div_b[0] = 0.5*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.5*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.5*F_0_div_b_Upwind[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b[1] = 0.447213595499958*F_0_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*F_0_div_b_Upwind[8]+0.5000000000000001*F_0_div_b_Upwind[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*F_0_div_b_Upwind[7]+0.447213595499958*F_0_div_b_Upwind[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind[6]+0.4472135954999579*F_0_div_b_Upwind[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*F_0_div_b_Upwind[4]+0.5*F_0_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[1]; 
  Ghat_F_0_div_b[2] = 0.447213595499958*F_0_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*F_0_div_b_Upwind[8]+0.447213595499958*F_0_div_b_Upwind[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind[7]+0.5000000000000001*F_0_div_b_Upwind[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*F_0_div_b_Upwind[6]+0.4472135954999579*F_0_div_b_Upwind[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*F_0_div_b_Upwind[5]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[2]; 
  Ghat_F_0_div_b[3] = 0.4*F_0_div_b_Upwind[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*F_0_div_b_Upwind[8]+0.4*F_0_div_b_Upwind[6]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*F_0_div_b_Upwind[7]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[7]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[6]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind[5]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind[4]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[2]; 
  Ghat_F_0_div_b[4] = 0.31943828249997*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*F_0_div_b_Upwind[8]+0.4472135954999579*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.31943828249997*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.5000000000000001*F_0_div_b_Upwind[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind[6]+0.31943828249997*F_0_div_b_Upwind[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[4]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind[1]*div_b_Surf[1]; 
  Ghat_F_0_div_b[5] = 0.31943828249997*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*F_0_div_b_Upwind[8]+0.31943828249997*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.5000000000000001*F_0_div_b_Upwind[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind[7]+0.4472135954999579*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[5]+0.4472135954999579*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind[2]*div_b_Surf[2]; 
  Ghat_F_0_div_b[6] = 0.2857142857142857*F_0_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*F_0_div_b_Upwind[8]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[8]+0.4*F_0_div_b_Upwind[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*F_0_div_b_Upwind[7]+0.4472135954999579*F_0_div_b_Upwind[5]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind[4]*div_b_Surf[6]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*F_0_div_b_Upwind[6]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind[6]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[6]+0.5000000000000001*F_0_div_b_Upwind[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind[4]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[3]; 
  Ghat_F_0_div_b[7] = 0.2857142857142857*F_0_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*F_0_div_b_Upwind[8]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[8]+0.31943828249997*F_0_div_b_Upwind[5]*div_b_Surf[7]+0.4472135954999579*F_0_div_b_Upwind[4]*div_b_Surf[7]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind[7]+0.4472135954999579*div_b_Surf[4]*F_0_div_b_Upwind[7]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[7]+0.4*F_0_div_b_Upwind[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*F_0_div_b_Upwind[6]+0.5000000000000001*F_0_div_b_Upwind[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind[5]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[3]; 
  Ghat_F_0_div_b[8] = 0.2040816326530612*F_0_div_b_Upwind[8]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind[5]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind[4]*div_b_Surf[8]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind[8]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind[8]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[8]+0.2857142857142857*F_0_div_b_Upwind[7]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind[7]+0.2857142857142857*F_0_div_b_Upwind[6]*div_b_Surf[6]+0.447213595499958*F_0_div_b_Upwind[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind[6]+0.5*F_0_div_b_Upwind[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*F_0_div_b_Upwind[5]+0.4*F_0_div_b_Upwind[3]*div_b_Surf[3]; 
  Ghat_G_1_div_b[0] = 0.5*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.5*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.5*G_1_div_b_Upwind[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b[1] = 0.447213595499958*G_1_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*G_1_div_b_Upwind[8]+0.5000000000000001*G_1_div_b_Upwind[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*G_1_div_b_Upwind[7]+0.447213595499958*G_1_div_b_Upwind[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind[6]+0.4472135954999579*G_1_div_b_Upwind[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*G_1_div_b_Upwind[4]+0.5*G_1_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[1]; 
  Ghat_G_1_div_b[2] = 0.447213595499958*G_1_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*G_1_div_b_Upwind[8]+0.447213595499958*G_1_div_b_Upwind[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind[7]+0.5000000000000001*G_1_div_b_Upwind[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*G_1_div_b_Upwind[6]+0.4472135954999579*G_1_div_b_Upwind[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*G_1_div_b_Upwind[5]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[2]; 
  Ghat_G_1_div_b[3] = 0.4*G_1_div_b_Upwind[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*G_1_div_b_Upwind[8]+0.4*G_1_div_b_Upwind[6]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*G_1_div_b_Upwind[7]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[7]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[6]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind[5]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind[4]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[2]; 
  Ghat_G_1_div_b[4] = 0.31943828249997*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*G_1_div_b_Upwind[8]+0.4472135954999579*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.31943828249997*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.5000000000000001*G_1_div_b_Upwind[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind[6]+0.31943828249997*G_1_div_b_Upwind[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[4]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind[1]*div_b_Surf[1]; 
  Ghat_G_1_div_b[5] = 0.31943828249997*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*G_1_div_b_Upwind[8]+0.31943828249997*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.5000000000000001*G_1_div_b_Upwind[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind[7]+0.4472135954999579*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[5]+0.4472135954999579*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind[2]*div_b_Surf[2]; 
  Ghat_G_1_div_b[6] = 0.2857142857142857*G_1_div_b_Upwind[6]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*G_1_div_b_Upwind[8]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[8]+0.4*G_1_div_b_Upwind[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*G_1_div_b_Upwind[7]+0.4472135954999579*G_1_div_b_Upwind[5]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind[4]*div_b_Surf[6]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*G_1_div_b_Upwind[6]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind[6]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[6]+0.5000000000000001*G_1_div_b_Upwind[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind[4]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[3]; 
  Ghat_G_1_div_b[7] = 0.2857142857142857*G_1_div_b_Upwind[7]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*G_1_div_b_Upwind[8]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[8]+0.31943828249997*G_1_div_b_Upwind[5]*div_b_Surf[7]+0.4472135954999579*G_1_div_b_Upwind[4]*div_b_Surf[7]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind[7]+0.4472135954999579*div_b_Surf[4]*G_1_div_b_Upwind[7]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[7]+0.4*G_1_div_b_Upwind[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*G_1_div_b_Upwind[6]+0.5000000000000001*G_1_div_b_Upwind[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind[5]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[3]; 
  Ghat_G_1_div_b[8] = 0.2040816326530612*G_1_div_b_Upwind[8]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind[5]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind[4]*div_b_Surf[8]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind[8]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind[8]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[8]+0.2857142857142857*G_1_div_b_Upwind[7]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind[7]+0.2857142857142857*G_1_div_b_Upwind[6]*div_b_Surf[6]+0.447213595499958*G_1_div_b_Upwind[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind[6]+0.5*G_1_div_b_Upwind[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*G_1_div_b_Upwind[5]+0.4*G_1_div_b_Upwind[3]*div_b_Surf[3]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0_div_b[0]*dv1par+0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0_div_b[1]*dv1par+0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 0.7071067811865475*Ghat_F_0_div_b[2]*dv1par+0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += (-1.224744871391589*Ghat_F_0_div_b[0]*dv1par)-1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += 0.7071067811865475*Ghat_F_0_div_b[3]*dv1par+0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += (-1.224744871391589*Ghat_F_0_div_b[1]*dv1par)-1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += (-1.224744871391589*Ghat_F_0_div_b[2]*dv1par)-1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 0.7071067811865475*Ghat_F_0_div_b[4]*dv1par+0.7071067811865475*Ghat_F_0[4]*dv1par; 
  out_F_0[8] += 0.7071067811865475*Ghat_F_0_div_b[5]*dv1par+0.7071067811865475*Ghat_F_0[5]*dv1par; 
  out_F_0[9] += 1.58113883008419*Ghat_F_0_div_b[0]*dv1par+1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[10] += (-1.224744871391589*Ghat_F_0_div_b[3]*dv1par)-1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[11] += 0.7071067811865475*Ghat_F_0_div_b[6]*dv1par+0.7071067811865475*Ghat_F_0[6]*dv1par; 
  out_F_0[12] += 0.7071067811865475*Ghat_F_0_div_b[7]*dv1par+0.7071067811865475*Ghat_F_0[7]*dv1par; 
  out_F_0[13] += (-1.224744871391589*Ghat_F_0_div_b[4]*dv1par)-1.224744871391589*Ghat_F_0[4]*dv1par; 
  out_F_0[14] += (-1.224744871391589*Ghat_F_0_div_b[5]*dv1par)-1.224744871391589*Ghat_F_0[5]*dv1par; 
  out_F_0[15] += 1.58113883008419*Ghat_F_0_div_b[1]*dv1par+1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_F_0[16] += 1.58113883008419*Ghat_F_0_div_b[2]*dv1par+1.58113883008419*Ghat_F_0[2]*dv1par; 
  out_F_0[17] += (-1.224744871391589*Ghat_F_0_div_b[6]*dv1par)-1.224744871391589*Ghat_F_0[6]*dv1par; 
  out_F_0[18] += (-1.224744871391589*Ghat_F_0_div_b[7]*dv1par)-1.224744871391589*Ghat_F_0[7]*dv1par; 
  out_F_0[19] += 1.58113883008419*Ghat_F_0_div_b[3]*dv1par+1.58113883008419*Ghat_F_0[3]*dv1par; 
  out_F_0[20] += 0.7071067811865475*Ghat_F_0_div_b[8]*dv1par+0.7071067811865475*Ghat_F_0[8]*dv1par; 
  out_F_0[21] += 1.58113883008419*Ghat_F_0_div_b[4]*dv1par+1.58113883008419*Ghat_F_0[4]*dv1par; 
  out_F_0[22] += 1.58113883008419*Ghat_F_0_div_b[5]*dv1par+1.58113883008419*Ghat_F_0[5]*dv1par; 
  out_F_0[23] += (-1.224744871391589*Ghat_F_0_div_b[8]*dv1par)-1.224744871391589*Ghat_F_0[8]*dv1par; 
  out_F_0[24] += 1.58113883008419*Ghat_F_0_div_b[6]*dv1par+1.58113883008419*Ghat_F_0[6]*dv1par; 
  out_F_0[25] += 1.58113883008419*Ghat_F_0_div_b[7]*dv1par+1.58113883008419*Ghat_F_0[7]*dv1par; 
  out_F_0[26] += 1.58113883008419*Ghat_F_0_div_b[8]*dv1par+1.58113883008419*Ghat_F_0[8]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1_div_b[0]*dv1par+0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1_div_b[1]*dv1par+0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 0.7071067811865475*Ghat_G_1_div_b[2]*dv1par+0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += (-1.224744871391589*Ghat_G_1_div_b[0]*dv1par)-1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += 0.7071067811865475*Ghat_G_1_div_b[3]*dv1par+0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += (-1.224744871391589*Ghat_G_1_div_b[1]*dv1par)-1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += (-1.224744871391589*Ghat_G_1_div_b[2]*dv1par)-1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 0.7071067811865475*Ghat_G_1_div_b[4]*dv1par+0.7071067811865475*Ghat_G_1[4]*dv1par; 
  out_G_1[8] += 0.7071067811865475*Ghat_G_1_div_b[5]*dv1par+0.7071067811865475*Ghat_G_1[5]*dv1par; 
  out_G_1[9] += 1.58113883008419*Ghat_G_1_div_b[0]*dv1par+1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[10] += (-1.224744871391589*Ghat_G_1_div_b[3]*dv1par)-1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[11] += 0.7071067811865475*Ghat_G_1_div_b[6]*dv1par+0.7071067811865475*Ghat_G_1[6]*dv1par; 
  out_G_1[12] += 0.7071067811865475*Ghat_G_1_div_b[7]*dv1par+0.7071067811865475*Ghat_G_1[7]*dv1par; 
  out_G_1[13] += (-1.224744871391589*Ghat_G_1_div_b[4]*dv1par)-1.224744871391589*Ghat_G_1[4]*dv1par; 
  out_G_1[14] += (-1.224744871391589*Ghat_G_1_div_b[5]*dv1par)-1.224744871391589*Ghat_G_1[5]*dv1par; 
  out_G_1[15] += 1.58113883008419*Ghat_G_1_div_b[1]*dv1par+1.58113883008419*Ghat_G_1[1]*dv1par; 
  out_G_1[16] += 1.58113883008419*Ghat_G_1_div_b[2]*dv1par+1.58113883008419*Ghat_G_1[2]*dv1par; 
  out_G_1[17] += (-1.224744871391589*Ghat_G_1_div_b[6]*dv1par)-1.224744871391589*Ghat_G_1[6]*dv1par; 
  out_G_1[18] += (-1.224744871391589*Ghat_G_1_div_b[7]*dv1par)-1.224744871391589*Ghat_G_1[7]*dv1par; 
  out_G_1[19] += 1.58113883008419*Ghat_G_1_div_b[3]*dv1par+1.58113883008419*Ghat_G_1[3]*dv1par; 
  out_G_1[20] += 0.7071067811865475*Ghat_G_1_div_b[8]*dv1par+0.7071067811865475*Ghat_G_1[8]*dv1par; 
  out_G_1[21] += 1.58113883008419*Ghat_G_1_div_b[4]*dv1par+1.58113883008419*Ghat_G_1[4]*dv1par; 
  out_G_1[22] += 1.58113883008419*Ghat_G_1_div_b[5]*dv1par+1.58113883008419*Ghat_G_1[5]*dv1par; 
  out_G_1[23] += (-1.224744871391589*Ghat_G_1_div_b[8]*dv1par)-1.224744871391589*Ghat_G_1[8]*dv1par; 
  out_G_1[24] += 1.58113883008419*Ghat_G_1_div_b[6]*dv1par+1.58113883008419*Ghat_G_1[6]*dv1par; 
  out_G_1[25] += 1.58113883008419*Ghat_G_1_div_b[7]*dv1par+1.58113883008419*Ghat_G_1[7]*dv1par; 
  out_G_1[26] += 1.58113883008419*Ghat_G_1_div_b[8]*dv1par+1.58113883008419*Ghat_G_1[8]*dv1par; 

  } 

  return 2.5*dv1par*cflFreq;

} 
