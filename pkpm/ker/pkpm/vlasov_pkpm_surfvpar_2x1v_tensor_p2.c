#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_pkpm_surfvpar_2x1v_tensor_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:            Cell-center coordinates. 
  // dxv[NDIM]:          Cell spacing. 
  // div_b:              Input volume expansion of div(b). 
  // pkpm_accel_vars:    Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // g_dist_sourcel/c/r: Input [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
  //                     (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]. 
  //                     in left/center/right cells. First input is mirror force source, second input is vperp characteristics source. 
  // fl/fc/fr:           Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)] in left/center/right cells. 
  // out:                Incremented output distribution functions in center cell. 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[27]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[27]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[27]; 

  const double *F_0_sourcel = &fl[27]; 
  const double *G_1_sourcel = &g_dist_sourcel[0]; 
  const double *F_0_sourcec = &fc[27]; 
  const double *G_1_sourcec = &g_dist_sourcec[0]; 
  const double *F_0_sourcer = &fr[27]; 
  const double *G_1_sourcer = &g_dist_sourcer[0]; 

  const double *bb_grad_u = &pkpm_accel_vars[9]; 
  const double *p_force = &pkpm_accel_vars[18]; 

  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[27]; 

  double alphaSurf_l[9] = {0.0}; 
  alphaSurf_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf_l[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf_l[4] = (-1.0*bb_grad_u[4]*wvpar)+0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf_l[5] = (-1.0*bb_grad_u[5]*wvpar)+0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf_l[6] = (-1.0*bb_grad_u[6]*wvpar)+0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf_l[7] = (-1.0*bb_grad_u[7]*wvpar)+0.5*bb_grad_u[7]*dvpar+p_force[7]; 
  alphaSurf_l[8] = (-1.0*bb_grad_u[8]*wvpar)+0.5*bb_grad_u[8]*dvpar+p_force[8]; 

  double alphaSurf_r[9] = {0.0}; 
  alphaSurf_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf_r[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf_r[4] = (-1.0*bb_grad_u[4]*wvpar)-0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf_r[5] = (-1.0*bb_grad_u[5]*wvpar)-0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf_r[6] = (-1.0*bb_grad_u[6]*wvpar)-0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf_r[7] = (-1.0*bb_grad_u[7]*wvpar)-0.5*bb_grad_u[7]*dvpar+p_force[7]; 
  alphaSurf_r[8] = (-1.0*bb_grad_u[8]*wvpar)-0.5*bb_grad_u[8]*dvpar+p_force[8]; 

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

  double F_0_UpwindQuad_l[9] = {0.0};
  double F_0_UpwindQuad_r[9] = {0.0};
  double F_0_Upwind_l[9] = {0.0};
  double F_0_Upwind_r[9] = {0.0};
  double Ghat_F_0_l[9] = {0.0}; 
  double Ghat_F_0_r[9] = {0.0}; 
  double G_1_UpwindQuad_l[9] = {0.0};
  double G_1_UpwindQuad_r[9] = {0.0};
  double G_1_Upwind_l[9] = {0.0};
  double G_1_Upwind_r[9] = {0.0};
  double Ghat_G_1_l[9] = {0.0}; 
  double Ghat_G_1_r[9] = {0.0}; 

  // get stable timestep of alpha_v = 1/rho (div(p_par b) - p_perp div(b)) - v_par bb : grad(u) 
  // from the quadrature point evaluation needed to compute upwinded distribution functions 
  double cflFreq = 0.0;
  double alphaOrd = 0.0;

  alphaOrd = 0.4*alphaSurf_l[8]-0.5999999999999995*alphaSurf_l[7]-0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])+0.9*alphaSurf_l[3]-0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1c); 
  } 
  alphaOrd = 0.4*alphaSurf_r[8]-0.5999999999999995*alphaSurf_r[7]-0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])+0.9*alphaSurf_r[3]-0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1r); 
  } 
  alphaOrd = (-0.5*alphaSurf_l[8])+0.75*alphaSurf_l[7]-0.5590169943749475*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1c); 
  } 
  alphaOrd = (-0.5*alphaSurf_r[8])+0.75*alphaSurf_r[7]-0.5590169943749475*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1r); 
  } 
  alphaOrd = 0.4*alphaSurf_l[8]-0.5999999999999995*alphaSurf_l[7]+0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])-0.9*alphaSurf_l[3]+0.6708203932499369*alphaSurf_l[2]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1c); 
  } 
  alphaOrd = 0.4*alphaSurf_r[8]-0.5999999999999995*alphaSurf_r[7]+0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])-0.9*alphaSurf_r[3]+0.6708203932499369*alphaSurf_r[2]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1r); 
  } 
  alphaOrd = (-0.5*alphaSurf_l[8])+0.75*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[5]-0.5590169943749475*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0l); 
    G_1_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0c); 
    G_1_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1c); 
  } 
  alphaOrd = (-0.5*alphaSurf_r[8])+0.75*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[5]-0.5590169943749475*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0c); 
    G_1_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0r); 
    G_1_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1r); 
  } 
  alphaOrd = 0.625*alphaSurf_l[8]-0.5590169943749475*(alphaSurf_l[5]+alphaSurf_l[4])+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0l); 
    G_1_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0c); 
    G_1_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1c); 
  } 
  alphaOrd = 0.625*alphaSurf_r[8]-0.5590169943749475*(alphaSurf_r[5]+alphaSurf_r[4])+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0c); 
    G_1_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0r); 
    G_1_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1r); 
  } 
  alphaOrd = (-0.5*alphaSurf_l[8])-0.75*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[5]-0.5590169943749475*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0l); 
    G_1_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0c); 
    G_1_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1c); 
  } 
  alphaOrd = (-0.5*alphaSurf_r[8])-0.75*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[5]-0.5590169943749475*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0c); 
    G_1_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0r); 
    G_1_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1r); 
  } 
  alphaOrd = 0.4*alphaSurf_l[8]+0.5999999999999995*alphaSurf_l[7]-0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])-0.9*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0l); 
    G_1_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0c); 
    G_1_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1c); 
  } 
  alphaOrd = 0.4*alphaSurf_r[8]+0.5999999999999995*alphaSurf_r[7]-0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])-0.9*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0c); 
    G_1_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0r); 
    G_1_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1r); 
  } 
  alphaOrd = (-0.5*alphaSurf_l[8])-0.75*alphaSurf_l[7]-0.5590169943749475*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0l); 
    G_1_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0c); 
    G_1_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1c); 
  } 
  alphaOrd = (-0.5*alphaSurf_r[8])-0.75*alphaSurf_r[7]-0.5590169943749475*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0c); 
    G_1_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0r); 
    G_1_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1r); 
  } 
  alphaOrd = 0.4*alphaSurf_l[8]+0.5999999999999995*alphaSurf_l[7]+0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])+0.9*alphaSurf_l[3]+0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0l); 
    G_1_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0c); 
    G_1_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1c); 
  } 
  alphaOrd = 0.4*alphaSurf_r[8]+0.5999999999999995*alphaSurf_r[7]+0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])+0.9*alphaSurf_r[3]+0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (alphaOrd > 0) { 
    F_0_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0c); 
    G_1_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0r); 
    G_1_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  double F_0_div_b_UpwindQuad_l[9] = {0.0};
  double F_0_div_b_UpwindQuad_r[9] = {0.0};
  double F_0_div_b_Upwind_l[9] = {0.0};
  double F_0_div_b_Upwind_r[9] = {0.0};
  double Ghat_F_0_div_b_l[9] = {0.0}; 
  double Ghat_F_0_div_b_r[9] = {0.0}; 

  double G_1_div_b_UpwindQuad_l[9] = {0.0};
  double G_1_div_b_UpwindQuad_r[9] = {0.0};
  double G_1_div_b_Upwind_l[9] = {0.0};
  double G_1_div_b_Upwind_r[9] = {0.0};
  double Ghat_G_1_div_b_l[9] = {0.0}; 
  double Ghat_G_1_div_b_r[9] = {0.0}; 

  if (0.4*div_b_Surf[8]-0.5999999999999995*div_b_Surf[7]-0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])+0.9*div_b_Surf[3]-0.6708203932499369*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(G_1_sourcer); 
  } 
  if ((-0.5*div_b_Surf[8])+0.75*div_b_Surf[7]-0.5590169943749475*div_b_Surf[5]+0.4472135954999579*div_b_Surf[4]-0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(G_1_sourcer); 
  } 
  if (0.4*div_b_Surf[8]-0.5999999999999995*div_b_Surf[7]+0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])-0.9*div_b_Surf[3]+0.6708203932499369*div_b_Surf[2]-0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(G_1_sourcer); 
  } 
  if ((-0.5*div_b_Surf[8])+0.75*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]-0.5590169943749475*div_b_Surf[4]-0.6708203932499369*div_b_Surf[2]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(G_1_sourcer); 
  } 
  if (0.625*div_b_Surf[8]-0.5590169943749475*(div_b_Surf[5]+div_b_Surf[4])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(G_1_sourcer); 
  } 
  if ((-0.5*div_b_Surf[8])-0.75*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]-0.5590169943749475*div_b_Surf[4]+0.6708203932499369*div_b_Surf[2]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(G_1_sourcer); 
  } 
  if (0.4*div_b_Surf[8]+0.5999999999999995*div_b_Surf[7]-0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])-0.9*div_b_Surf[3]-0.6708203932499369*div_b_Surf[2]+0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(G_1_sourcer); 
  } 
  if ((-0.5*div_b_Surf[8])-0.75*div_b_Surf[7]-0.5590169943749475*div_b_Surf[5]+0.4472135954999579*div_b_Surf[4]+0.6708203932499369*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(G_1_sourcer); 
  } 
  if (0.4*div_b_Surf[8]+0.5999999999999995*div_b_Surf[7]+0.5999999999999999*div_b_Surf[6]+0.4472135954999579*(div_b_Surf[5]+div_b_Surf[4])+0.9*div_b_Surf[3]+0.6708203932499369*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(G_1_sourcer); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(F_0_div_b_UpwindQuad_l, F_0_div_b_Upwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(F_0_div_b_UpwindQuad_r, F_0_div_b_Upwind_r); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_div_b_UpwindQuad_l, G_1_div_b_Upwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_div_b_UpwindQuad_r, G_1_div_b_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*F_0_Upwind_l[8]*alphaSurf_l[8]+0.5*F_0_Upwind_l[7]*alphaSurf_l[7]+0.5*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5*F_0_Upwind_l[5]*alphaSurf_l[5]+0.5*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[3]*alphaSurf_l[3]+0.5*F_0_Upwind_l[2]*alphaSurf_l[2]+0.5*F_0_Upwind_l[1]*alphaSurf_l[1]+0.5*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.447213595499958*F_0_Upwind_l[7]*alphaSurf_l[8]+0.447213595499958*alphaSurf_l[7]*F_0_Upwind_l[8]+0.5000000000000001*F_0_Upwind_l[5]*alphaSurf_l[7]+0.5000000000000001*alphaSurf_l[5]*F_0_Upwind_l[7]+0.447213595499958*F_0_Upwind_l[3]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[6]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[1]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.447213595499958*F_0_Upwind_l[6]*alphaSurf_l[8]+0.447213595499958*alphaSurf_l[6]*F_0_Upwind_l[8]+0.447213595499958*F_0_Upwind_l[3]*alphaSurf_l[7]+0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[7]+0.5000000000000001*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[4]*F_0_Upwind_l[6]+0.4472135954999579*F_0_Upwind_l[2]*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[2]*F_0_Upwind_l[5]+0.5*F_0_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.4*F_0_Upwind_l[3]*alphaSurf_l[8]+0.4*alphaSurf_l[3]*F_0_Upwind_l[8]+0.4*F_0_Upwind_l[6]*alphaSurf_l[7]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[7]+0.4*alphaSurf_l[6]*F_0_Upwind_l[7]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[7]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[6]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[4] = 0.31943828249997*F_0_Upwind_l[8]*alphaSurf_l[8]+0.5*F_0_Upwind_l[5]*alphaSurf_l[8]+0.5*alphaSurf_l[5]*F_0_Upwind_l[8]+0.4472135954999579*F_0_Upwind_l[7]*alphaSurf_l[7]+0.31943828249997*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*F_0_Upwind_l[4]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_F_0_l[5] = 0.31943828249997*F_0_Upwind_l[8]*alphaSurf_l[8]+0.5*F_0_Upwind_l[4]*alphaSurf_l[8]+0.5*alphaSurf_l[4]*F_0_Upwind_l[8]+0.31943828249997*F_0_Upwind_l[7]*alphaSurf_l[7]+0.5000000000000001*F_0_Upwind_l[1]*alphaSurf_l[7]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[6]*alphaSurf_l[6]+0.31943828249997*F_0_Upwind_l[5]*alphaSurf_l[5]+0.5*F_0_Upwind_l[0]*alphaSurf_l[5]+0.5*alphaSurf_l[0]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_F_0_l[6] = 0.2857142857142857*F_0_Upwind_l[6]*alphaSurf_l[8]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[8]+0.2857142857142857*alphaSurf_l[6]*F_0_Upwind_l[8]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[8]+0.4*F_0_Upwind_l[3]*alphaSurf_l[7]+0.4*alphaSurf_l[3]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[5]*alphaSurf_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5*F_0_Upwind_l[0]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[5]*F_0_Upwind_l[6]+0.31943828249997*alphaSurf_l[4]*F_0_Upwind_l[6]+0.5*alphaSurf_l[0]*F_0_Upwind_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[4]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[4]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[7] = 0.2857142857142857*F_0_Upwind_l[7]*alphaSurf_l[8]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[8]+0.2857142857142857*alphaSurf_l[7]*F_0_Upwind_l[8]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[8]+0.31943828249997*F_0_Upwind_l[5]*alphaSurf_l[7]+0.4472135954999579*F_0_Upwind_l[4]*alphaSurf_l[7]+0.5*F_0_Upwind_l[0]*alphaSurf_l[7]+0.31943828249997*alphaSurf_l[5]*F_0_Upwind_l[7]+0.4472135954999579*alphaSurf_l[4]*F_0_Upwind_l[7]+0.5*alphaSurf_l[0]*F_0_Upwind_l[7]+0.4*F_0_Upwind_l[3]*alphaSurf_l[6]+0.4*alphaSurf_l[3]*F_0_Upwind_l[6]+0.5000000000000001*F_0_Upwind_l[1]*alphaSurf_l[5]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[5]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[8] = 0.2040816326530612*F_0_Upwind_l[8]*alphaSurf_l[8]+0.31943828249997*F_0_Upwind_l[5]*alphaSurf_l[8]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[8]+0.5*F_0_Upwind_l[0]*alphaSurf_l[8]+0.31943828249997*alphaSurf_l[5]*F_0_Upwind_l[8]+0.31943828249997*alphaSurf_l[4]*F_0_Upwind_l[8]+0.5*alphaSurf_l[0]*F_0_Upwind_l[8]+0.2857142857142857*F_0_Upwind_l[7]*alphaSurf_l[7]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[7]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[7]+0.2857142857142857*F_0_Upwind_l[6]*alphaSurf_l[6]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[6]+0.5*F_0_Upwind_l[4]*alphaSurf_l[5]+0.5*alphaSurf_l[4]*F_0_Upwind_l[5]+0.4*F_0_Upwind_l[3]*alphaSurf_l[3]; 
  Ghat_G_1_l[0] = 0.5*G_1_Upwind_l[8]*alphaSurf_l[8]+0.5*G_1_Upwind_l[7]*alphaSurf_l[7]+0.5*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5*G_1_Upwind_l[5]*alphaSurf_l[5]+0.5*G_1_Upwind_l[4]*alphaSurf_l[4]+0.5*G_1_Upwind_l[3]*alphaSurf_l[3]+0.5*G_1_Upwind_l[2]*alphaSurf_l[2]+0.5*G_1_Upwind_l[1]*alphaSurf_l[1]+0.5*G_1_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_G_1_l[1] = 0.447213595499958*G_1_Upwind_l[7]*alphaSurf_l[8]+0.447213595499958*alphaSurf_l[7]*G_1_Upwind_l[8]+0.5000000000000001*G_1_Upwind_l[5]*alphaSurf_l[7]+0.5000000000000001*alphaSurf_l[5]*G_1_Upwind_l[7]+0.447213595499958*G_1_Upwind_l[3]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[3]*G_1_Upwind_l[6]+0.4472135954999579*G_1_Upwind_l[1]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[1]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.447213595499958*G_1_Upwind_l[6]*alphaSurf_l[8]+0.447213595499958*alphaSurf_l[6]*G_1_Upwind_l[8]+0.447213595499958*G_1_Upwind_l[3]*alphaSurf_l[7]+0.447213595499958*alphaSurf_l[3]*G_1_Upwind_l[7]+0.5000000000000001*G_1_Upwind_l[4]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[4]*G_1_Upwind_l[6]+0.4472135954999579*G_1_Upwind_l[2]*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[2]*G_1_Upwind_l[5]+0.5*G_1_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[3] = 0.4*G_1_Upwind_l[3]*alphaSurf_l[8]+0.4*alphaSurf_l[3]*G_1_Upwind_l[8]+0.4*G_1_Upwind_l[6]*alphaSurf_l[7]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[7]+0.4*alphaSurf_l[6]*G_1_Upwind_l[7]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[7]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[6]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[3]*G_1_Upwind_l[5]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[4] = 0.31943828249997*G_1_Upwind_l[8]*alphaSurf_l[8]+0.5*G_1_Upwind_l[5]*alphaSurf_l[8]+0.5*alphaSurf_l[5]*G_1_Upwind_l[8]+0.4472135954999579*G_1_Upwind_l[7]*alphaSurf_l[7]+0.31943828249997*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5000000000000001*G_1_Upwind_l[2]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*G_1_Upwind_l[6]+0.31943828249997*G_1_Upwind_l[4]*alphaSurf_l[4]+0.5*G_1_Upwind_l[0]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*G_1_Upwind_l[4]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*G_1_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_G_1_l[5] = 0.31943828249997*G_1_Upwind_l[8]*alphaSurf_l[8]+0.5*G_1_Upwind_l[4]*alphaSurf_l[8]+0.5*alphaSurf_l[4]*G_1_Upwind_l[8]+0.31943828249997*G_1_Upwind_l[7]*alphaSurf_l[7]+0.5000000000000001*G_1_Upwind_l[1]*alphaSurf_l[7]+0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[7]+0.4472135954999579*G_1_Upwind_l[6]*alphaSurf_l[6]+0.31943828249997*G_1_Upwind_l[5]*alphaSurf_l[5]+0.5*G_1_Upwind_l[0]*alphaSurf_l[5]+0.5*alphaSurf_l[0]*G_1_Upwind_l[5]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*G_1_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_G_1_l[6] = 0.2857142857142857*G_1_Upwind_l[6]*alphaSurf_l[8]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[8]+0.2857142857142857*alphaSurf_l[6]*G_1_Upwind_l[8]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[8]+0.4*G_1_Upwind_l[3]*alphaSurf_l[7]+0.4*alphaSurf_l[3]*G_1_Upwind_l[7]+0.4472135954999579*G_1_Upwind_l[5]*alphaSurf_l[6]+0.31943828249997*G_1_Upwind_l[4]*alphaSurf_l[6]+0.5*G_1_Upwind_l[0]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[5]*G_1_Upwind_l[6]+0.31943828249997*alphaSurf_l[4]*G_1_Upwind_l[6]+0.5*alphaSurf_l[0]*G_1_Upwind_l[6]+0.5000000000000001*G_1_Upwind_l[2]*alphaSurf_l[4]+0.5000000000000001*alphaSurf_l[2]*G_1_Upwind_l[4]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[7] = 0.2857142857142857*G_1_Upwind_l[7]*alphaSurf_l[8]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[8]+0.2857142857142857*alphaSurf_l[7]*G_1_Upwind_l[8]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[8]+0.31943828249997*G_1_Upwind_l[5]*alphaSurf_l[7]+0.4472135954999579*G_1_Upwind_l[4]*alphaSurf_l[7]+0.5*G_1_Upwind_l[0]*alphaSurf_l[7]+0.31943828249997*alphaSurf_l[5]*G_1_Upwind_l[7]+0.4472135954999579*alphaSurf_l[4]*G_1_Upwind_l[7]+0.5*alphaSurf_l[0]*G_1_Upwind_l[7]+0.4*G_1_Upwind_l[3]*alphaSurf_l[6]+0.4*alphaSurf_l[3]*G_1_Upwind_l[6]+0.5000000000000001*G_1_Upwind_l[1]*alphaSurf_l[5]+0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[5]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[8] = 0.2040816326530612*G_1_Upwind_l[8]*alphaSurf_l[8]+0.31943828249997*G_1_Upwind_l[5]*alphaSurf_l[8]+0.31943828249997*G_1_Upwind_l[4]*alphaSurf_l[8]+0.5*G_1_Upwind_l[0]*alphaSurf_l[8]+0.31943828249997*alphaSurf_l[5]*G_1_Upwind_l[8]+0.31943828249997*alphaSurf_l[4]*G_1_Upwind_l[8]+0.5*alphaSurf_l[0]*G_1_Upwind_l[8]+0.2857142857142857*G_1_Upwind_l[7]*alphaSurf_l[7]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[7]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[7]+0.2857142857142857*G_1_Upwind_l[6]*alphaSurf_l[6]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[6]+0.5*G_1_Upwind_l[4]*alphaSurf_l[5]+0.5*alphaSurf_l[4]*G_1_Upwind_l[5]+0.4*G_1_Upwind_l[3]*alphaSurf_l[3]; 
  Ghat_F_0_div_b_l[0] = 0.5*F_0_div_b_Upwind_l[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_l[7]*div_b_Surf[7]+0.5*F_0_div_b_Upwind_l[6]*div_b_Surf[6]+0.5*F_0_div_b_Upwind_l[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind_l[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind_l[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind_l[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind_l[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_l[1] = 0.447213595499958*F_0_div_b_Upwind_l[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*F_0_div_b_Upwind_l[8]+0.5000000000000001*F_0_div_b_Upwind_l[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*F_0_div_b_Upwind_l[7]+0.447213595499958*F_0_div_b_Upwind_l[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind_l[6]+0.4472135954999579*F_0_div_b_Upwind_l[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*F_0_div_b_Upwind_l[4]+0.5*F_0_div_b_Upwind_l[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind_l[3]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[1]; 
  Ghat_F_0_div_b_l[2] = 0.447213595499958*F_0_div_b_Upwind_l[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*F_0_div_b_Upwind_l[8]+0.447213595499958*F_0_div_b_Upwind_l[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind_l[7]+0.5000000000000001*F_0_div_b_Upwind_l[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*F_0_div_b_Upwind_l[6]+0.4472135954999579*F_0_div_b_Upwind_l[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*F_0_div_b_Upwind_l[5]+0.5*F_0_div_b_Upwind_l[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_l[3]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[2]; 
  Ghat_F_0_div_b_l[3] = 0.4*F_0_div_b_Upwind_l[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*F_0_div_b_Upwind_l[8]+0.4*F_0_div_b_Upwind_l[6]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind_l[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*F_0_div_b_Upwind_l[7]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_l[7]+0.447213595499958*F_0_div_b_Upwind_l[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_l[6]+0.4472135954999579*F_0_div_b_Upwind_l[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind_l[5]+0.4472135954999579*F_0_div_b_Upwind_l[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind_l[4]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[3]+0.5*F_0_div_b_Upwind_l[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_l[2]; 
  Ghat_F_0_div_b_l[4] = 0.31943828249997*F_0_div_b_Upwind_l[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_l[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*F_0_div_b_Upwind_l[8]+0.4472135954999579*F_0_div_b_Upwind_l[7]*div_b_Surf[7]+0.31943828249997*F_0_div_b_Upwind_l[6]*div_b_Surf[6]+0.5000000000000001*F_0_div_b_Upwind_l[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind_l[6]+0.31943828249997*F_0_div_b_Upwind_l[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[4]+0.4472135954999579*F_0_div_b_Upwind_l[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind_l[1]*div_b_Surf[1]; 
  Ghat_F_0_div_b_l[5] = 0.31943828249997*F_0_div_b_Upwind_l[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_l[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*F_0_div_b_Upwind_l[8]+0.31943828249997*F_0_div_b_Upwind_l[7]*div_b_Surf[7]+0.5000000000000001*F_0_div_b_Upwind_l[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind_l[7]+0.4472135954999579*F_0_div_b_Upwind_l[6]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind_l[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[5]+0.4472135954999579*F_0_div_b_Upwind_l[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind_l[2]*div_b_Surf[2]; 
  Ghat_F_0_div_b_l[6] = 0.2857142857142857*F_0_div_b_Upwind_l[6]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind_l[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*F_0_div_b_Upwind_l[8]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_l[8]+0.4*F_0_div_b_Upwind_l[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*F_0_div_b_Upwind_l[7]+0.4472135954999579*F_0_div_b_Upwind_l[5]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind_l[4]*div_b_Surf[6]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*F_0_div_b_Upwind_l[6]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind_l[6]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[6]+0.5000000000000001*F_0_div_b_Upwind_l[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind_l[4]+0.447213595499958*F_0_div_b_Upwind_l[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_l[3]; 
  Ghat_F_0_div_b_l[7] = 0.2857142857142857*F_0_div_b_Upwind_l[7]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind_l[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*F_0_div_b_Upwind_l[8]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_l[8]+0.31943828249997*F_0_div_b_Upwind_l[5]*div_b_Surf[7]+0.4472135954999579*F_0_div_b_Upwind_l[4]*div_b_Surf[7]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind_l[7]+0.4472135954999579*div_b_Surf[4]*F_0_div_b_Upwind_l[7]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[7]+0.4*F_0_div_b_Upwind_l[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*F_0_div_b_Upwind_l[6]+0.5000000000000001*F_0_div_b_Upwind_l[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind_l[5]+0.447213595499958*F_0_div_b_Upwind_l[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_l[3]; 
  Ghat_F_0_div_b_l[8] = 0.2040816326530612*F_0_div_b_Upwind_l[8]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind_l[5]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind_l[4]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind_l[8]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind_l[8]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[8]+0.2857142857142857*F_0_div_b_Upwind_l[7]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind_l[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_l[7]+0.2857142857142857*F_0_div_b_Upwind_l[6]*div_b_Surf[6]+0.447213595499958*F_0_div_b_Upwind_l[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_l[6]+0.5*F_0_div_b_Upwind_l[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*F_0_div_b_Upwind_l[5]+0.4*F_0_div_b_Upwind_l[3]*div_b_Surf[3]; 
  Ghat_G_1_div_b_l[0] = 0.5*G_1_div_b_Upwind_l[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_l[7]*div_b_Surf[7]+0.5*G_1_div_b_Upwind_l[6]*div_b_Surf[6]+0.5*G_1_div_b_Upwind_l[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind_l[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind_l[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind_l[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind_l[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_l[1] = 0.447213595499958*G_1_div_b_Upwind_l[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*G_1_div_b_Upwind_l[8]+0.5000000000000001*G_1_div_b_Upwind_l[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*G_1_div_b_Upwind_l[7]+0.447213595499958*G_1_div_b_Upwind_l[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind_l[6]+0.4472135954999579*G_1_div_b_Upwind_l[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*G_1_div_b_Upwind_l[4]+0.5*G_1_div_b_Upwind_l[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind_l[3]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[1]; 
  Ghat_G_1_div_b_l[2] = 0.447213595499958*G_1_div_b_Upwind_l[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*G_1_div_b_Upwind_l[8]+0.447213595499958*G_1_div_b_Upwind_l[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind_l[7]+0.5000000000000001*G_1_div_b_Upwind_l[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*G_1_div_b_Upwind_l[6]+0.4472135954999579*G_1_div_b_Upwind_l[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*G_1_div_b_Upwind_l[5]+0.5*G_1_div_b_Upwind_l[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_l[3]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[2]; 
  Ghat_G_1_div_b_l[3] = 0.4*G_1_div_b_Upwind_l[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*G_1_div_b_Upwind_l[8]+0.4*G_1_div_b_Upwind_l[6]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind_l[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*G_1_div_b_Upwind_l[7]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_l[7]+0.447213595499958*G_1_div_b_Upwind_l[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_l[6]+0.4472135954999579*G_1_div_b_Upwind_l[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind_l[5]+0.4472135954999579*G_1_div_b_Upwind_l[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind_l[4]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[3]+0.5*G_1_div_b_Upwind_l[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_l[2]; 
  Ghat_G_1_div_b_l[4] = 0.31943828249997*G_1_div_b_Upwind_l[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_l[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*G_1_div_b_Upwind_l[8]+0.4472135954999579*G_1_div_b_Upwind_l[7]*div_b_Surf[7]+0.31943828249997*G_1_div_b_Upwind_l[6]*div_b_Surf[6]+0.5000000000000001*G_1_div_b_Upwind_l[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind_l[6]+0.31943828249997*G_1_div_b_Upwind_l[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[4]+0.4472135954999579*G_1_div_b_Upwind_l[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind_l[1]*div_b_Surf[1]; 
  Ghat_G_1_div_b_l[5] = 0.31943828249997*G_1_div_b_Upwind_l[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_l[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*G_1_div_b_Upwind_l[8]+0.31943828249997*G_1_div_b_Upwind_l[7]*div_b_Surf[7]+0.5000000000000001*G_1_div_b_Upwind_l[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind_l[7]+0.4472135954999579*G_1_div_b_Upwind_l[6]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind_l[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[5]+0.4472135954999579*G_1_div_b_Upwind_l[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind_l[2]*div_b_Surf[2]; 
  Ghat_G_1_div_b_l[6] = 0.2857142857142857*G_1_div_b_Upwind_l[6]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind_l[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*G_1_div_b_Upwind_l[8]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_l[8]+0.4*G_1_div_b_Upwind_l[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*G_1_div_b_Upwind_l[7]+0.4472135954999579*G_1_div_b_Upwind_l[5]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind_l[4]*div_b_Surf[6]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*G_1_div_b_Upwind_l[6]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind_l[6]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[6]+0.5000000000000001*G_1_div_b_Upwind_l[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind_l[4]+0.447213595499958*G_1_div_b_Upwind_l[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_l[3]; 
  Ghat_G_1_div_b_l[7] = 0.2857142857142857*G_1_div_b_Upwind_l[7]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind_l[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*G_1_div_b_Upwind_l[8]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_l[8]+0.31943828249997*G_1_div_b_Upwind_l[5]*div_b_Surf[7]+0.4472135954999579*G_1_div_b_Upwind_l[4]*div_b_Surf[7]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind_l[7]+0.4472135954999579*div_b_Surf[4]*G_1_div_b_Upwind_l[7]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[7]+0.4*G_1_div_b_Upwind_l[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*G_1_div_b_Upwind_l[6]+0.5000000000000001*G_1_div_b_Upwind_l[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind_l[5]+0.447213595499958*G_1_div_b_Upwind_l[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_l[3]; 
  Ghat_G_1_div_b_l[8] = 0.2040816326530612*G_1_div_b_Upwind_l[8]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind_l[5]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind_l[4]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind_l[8]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind_l[8]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[8]+0.2857142857142857*G_1_div_b_Upwind_l[7]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind_l[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_l[7]+0.2857142857142857*G_1_div_b_Upwind_l[6]*div_b_Surf[6]+0.447213595499958*G_1_div_b_Upwind_l[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_l[6]+0.5*G_1_div_b_Upwind_l[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*G_1_div_b_Upwind_l[5]+0.4*G_1_div_b_Upwind_l[3]*div_b_Surf[3]; 

  Ghat_F_0_r[0] = 0.5*F_0_Upwind_r[8]*alphaSurf_r[8]+0.5*F_0_Upwind_r[7]*alphaSurf_r[7]+0.5*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5*F_0_Upwind_r[5]*alphaSurf_r[5]+0.5*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[3]*alphaSurf_r[3]+0.5*F_0_Upwind_r[2]*alphaSurf_r[2]+0.5*F_0_Upwind_r[1]*alphaSurf_r[1]+0.5*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.447213595499958*F_0_Upwind_r[7]*alphaSurf_r[8]+0.447213595499958*alphaSurf_r[7]*F_0_Upwind_r[8]+0.5000000000000001*F_0_Upwind_r[5]*alphaSurf_r[7]+0.5000000000000001*alphaSurf_r[5]*F_0_Upwind_r[7]+0.447213595499958*F_0_Upwind_r[3]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[6]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[1]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.447213595499958*F_0_Upwind_r[6]*alphaSurf_r[8]+0.447213595499958*alphaSurf_r[6]*F_0_Upwind_r[8]+0.447213595499958*F_0_Upwind_r[3]*alphaSurf_r[7]+0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[7]+0.5000000000000001*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[4]*F_0_Upwind_r[6]+0.4472135954999579*F_0_Upwind_r[2]*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[2]*F_0_Upwind_r[5]+0.5*F_0_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.4*F_0_Upwind_r[3]*alphaSurf_r[8]+0.4*alphaSurf_r[3]*F_0_Upwind_r[8]+0.4*F_0_Upwind_r[6]*alphaSurf_r[7]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[7]+0.4*alphaSurf_r[6]*F_0_Upwind_r[7]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[7]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[6]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[4] = 0.31943828249997*F_0_Upwind_r[8]*alphaSurf_r[8]+0.5*F_0_Upwind_r[5]*alphaSurf_r[8]+0.5*alphaSurf_r[5]*F_0_Upwind_r[8]+0.4472135954999579*F_0_Upwind_r[7]*alphaSurf_r[7]+0.31943828249997*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*F_0_Upwind_r[4]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_F_0_r[5] = 0.31943828249997*F_0_Upwind_r[8]*alphaSurf_r[8]+0.5*F_0_Upwind_r[4]*alphaSurf_r[8]+0.5*alphaSurf_r[4]*F_0_Upwind_r[8]+0.31943828249997*F_0_Upwind_r[7]*alphaSurf_r[7]+0.5000000000000001*F_0_Upwind_r[1]*alphaSurf_r[7]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[6]*alphaSurf_r[6]+0.31943828249997*F_0_Upwind_r[5]*alphaSurf_r[5]+0.5*F_0_Upwind_r[0]*alphaSurf_r[5]+0.5*alphaSurf_r[0]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_F_0_r[6] = 0.2857142857142857*F_0_Upwind_r[6]*alphaSurf_r[8]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[8]+0.2857142857142857*alphaSurf_r[6]*F_0_Upwind_r[8]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[8]+0.4*F_0_Upwind_r[3]*alphaSurf_r[7]+0.4*alphaSurf_r[3]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[5]*alphaSurf_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5*F_0_Upwind_r[0]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[5]*F_0_Upwind_r[6]+0.31943828249997*alphaSurf_r[4]*F_0_Upwind_r[6]+0.5*alphaSurf_r[0]*F_0_Upwind_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[4]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[4]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[7] = 0.2857142857142857*F_0_Upwind_r[7]*alphaSurf_r[8]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[8]+0.2857142857142857*alphaSurf_r[7]*F_0_Upwind_r[8]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[8]+0.31943828249997*F_0_Upwind_r[5]*alphaSurf_r[7]+0.4472135954999579*F_0_Upwind_r[4]*alphaSurf_r[7]+0.5*F_0_Upwind_r[0]*alphaSurf_r[7]+0.31943828249997*alphaSurf_r[5]*F_0_Upwind_r[7]+0.4472135954999579*alphaSurf_r[4]*F_0_Upwind_r[7]+0.5*alphaSurf_r[0]*F_0_Upwind_r[7]+0.4*F_0_Upwind_r[3]*alphaSurf_r[6]+0.4*alphaSurf_r[3]*F_0_Upwind_r[6]+0.5000000000000001*F_0_Upwind_r[1]*alphaSurf_r[5]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[5]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[8] = 0.2040816326530612*F_0_Upwind_r[8]*alphaSurf_r[8]+0.31943828249997*F_0_Upwind_r[5]*alphaSurf_r[8]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[8]+0.5*F_0_Upwind_r[0]*alphaSurf_r[8]+0.31943828249997*alphaSurf_r[5]*F_0_Upwind_r[8]+0.31943828249997*alphaSurf_r[4]*F_0_Upwind_r[8]+0.5*alphaSurf_r[0]*F_0_Upwind_r[8]+0.2857142857142857*F_0_Upwind_r[7]*alphaSurf_r[7]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[7]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[7]+0.2857142857142857*F_0_Upwind_r[6]*alphaSurf_r[6]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[6]+0.5*F_0_Upwind_r[4]*alphaSurf_r[5]+0.5*alphaSurf_r[4]*F_0_Upwind_r[5]+0.4*F_0_Upwind_r[3]*alphaSurf_r[3]; 
  Ghat_G_1_r[0] = 0.5*G_1_Upwind_r[8]*alphaSurf_r[8]+0.5*G_1_Upwind_r[7]*alphaSurf_r[7]+0.5*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5*G_1_Upwind_r[5]*alphaSurf_r[5]+0.5*G_1_Upwind_r[4]*alphaSurf_r[4]+0.5*G_1_Upwind_r[3]*alphaSurf_r[3]+0.5*G_1_Upwind_r[2]*alphaSurf_r[2]+0.5*G_1_Upwind_r[1]*alphaSurf_r[1]+0.5*G_1_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_G_1_r[1] = 0.447213595499958*G_1_Upwind_r[7]*alphaSurf_r[8]+0.447213595499958*alphaSurf_r[7]*G_1_Upwind_r[8]+0.5000000000000001*G_1_Upwind_r[5]*alphaSurf_r[7]+0.5000000000000001*alphaSurf_r[5]*G_1_Upwind_r[7]+0.447213595499958*G_1_Upwind_r[3]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[3]*G_1_Upwind_r[6]+0.4472135954999579*G_1_Upwind_r[1]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[1]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.447213595499958*G_1_Upwind_r[6]*alphaSurf_r[8]+0.447213595499958*alphaSurf_r[6]*G_1_Upwind_r[8]+0.447213595499958*G_1_Upwind_r[3]*alphaSurf_r[7]+0.447213595499958*alphaSurf_r[3]*G_1_Upwind_r[7]+0.5000000000000001*G_1_Upwind_r[4]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[4]*G_1_Upwind_r[6]+0.4472135954999579*G_1_Upwind_r[2]*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[2]*G_1_Upwind_r[5]+0.5*G_1_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[3] = 0.4*G_1_Upwind_r[3]*alphaSurf_r[8]+0.4*alphaSurf_r[3]*G_1_Upwind_r[8]+0.4*G_1_Upwind_r[6]*alphaSurf_r[7]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[7]+0.4*alphaSurf_r[6]*G_1_Upwind_r[7]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[7]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[6]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[3]*G_1_Upwind_r[5]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[4] = 0.31943828249997*G_1_Upwind_r[8]*alphaSurf_r[8]+0.5*G_1_Upwind_r[5]*alphaSurf_r[8]+0.5*alphaSurf_r[5]*G_1_Upwind_r[8]+0.4472135954999579*G_1_Upwind_r[7]*alphaSurf_r[7]+0.31943828249997*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5000000000000001*G_1_Upwind_r[2]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*G_1_Upwind_r[6]+0.31943828249997*G_1_Upwind_r[4]*alphaSurf_r[4]+0.5*G_1_Upwind_r[0]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*G_1_Upwind_r[4]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*G_1_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_G_1_r[5] = 0.31943828249997*G_1_Upwind_r[8]*alphaSurf_r[8]+0.5*G_1_Upwind_r[4]*alphaSurf_r[8]+0.5*alphaSurf_r[4]*G_1_Upwind_r[8]+0.31943828249997*G_1_Upwind_r[7]*alphaSurf_r[7]+0.5000000000000001*G_1_Upwind_r[1]*alphaSurf_r[7]+0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[7]+0.4472135954999579*G_1_Upwind_r[6]*alphaSurf_r[6]+0.31943828249997*G_1_Upwind_r[5]*alphaSurf_r[5]+0.5*G_1_Upwind_r[0]*alphaSurf_r[5]+0.5*alphaSurf_r[0]*G_1_Upwind_r[5]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*G_1_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_G_1_r[6] = 0.2857142857142857*G_1_Upwind_r[6]*alphaSurf_r[8]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[8]+0.2857142857142857*alphaSurf_r[6]*G_1_Upwind_r[8]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[8]+0.4*G_1_Upwind_r[3]*alphaSurf_r[7]+0.4*alphaSurf_r[3]*G_1_Upwind_r[7]+0.4472135954999579*G_1_Upwind_r[5]*alphaSurf_r[6]+0.31943828249997*G_1_Upwind_r[4]*alphaSurf_r[6]+0.5*G_1_Upwind_r[0]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[5]*G_1_Upwind_r[6]+0.31943828249997*alphaSurf_r[4]*G_1_Upwind_r[6]+0.5*alphaSurf_r[0]*G_1_Upwind_r[6]+0.5000000000000001*G_1_Upwind_r[2]*alphaSurf_r[4]+0.5000000000000001*alphaSurf_r[2]*G_1_Upwind_r[4]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[7] = 0.2857142857142857*G_1_Upwind_r[7]*alphaSurf_r[8]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[8]+0.2857142857142857*alphaSurf_r[7]*G_1_Upwind_r[8]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[8]+0.31943828249997*G_1_Upwind_r[5]*alphaSurf_r[7]+0.4472135954999579*G_1_Upwind_r[4]*alphaSurf_r[7]+0.5*G_1_Upwind_r[0]*alphaSurf_r[7]+0.31943828249997*alphaSurf_r[5]*G_1_Upwind_r[7]+0.4472135954999579*alphaSurf_r[4]*G_1_Upwind_r[7]+0.5*alphaSurf_r[0]*G_1_Upwind_r[7]+0.4*G_1_Upwind_r[3]*alphaSurf_r[6]+0.4*alphaSurf_r[3]*G_1_Upwind_r[6]+0.5000000000000001*G_1_Upwind_r[1]*alphaSurf_r[5]+0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[5]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[8] = 0.2040816326530612*G_1_Upwind_r[8]*alphaSurf_r[8]+0.31943828249997*G_1_Upwind_r[5]*alphaSurf_r[8]+0.31943828249997*G_1_Upwind_r[4]*alphaSurf_r[8]+0.5*G_1_Upwind_r[0]*alphaSurf_r[8]+0.31943828249997*alphaSurf_r[5]*G_1_Upwind_r[8]+0.31943828249997*alphaSurf_r[4]*G_1_Upwind_r[8]+0.5*alphaSurf_r[0]*G_1_Upwind_r[8]+0.2857142857142857*G_1_Upwind_r[7]*alphaSurf_r[7]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[7]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[7]+0.2857142857142857*G_1_Upwind_r[6]*alphaSurf_r[6]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[6]+0.5*G_1_Upwind_r[4]*alphaSurf_r[5]+0.5*alphaSurf_r[4]*G_1_Upwind_r[5]+0.4*G_1_Upwind_r[3]*alphaSurf_r[3]; 
  Ghat_F_0_div_b_r[0] = 0.5*F_0_div_b_Upwind_r[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_r[7]*div_b_Surf[7]+0.5*F_0_div_b_Upwind_r[6]*div_b_Surf[6]+0.5*F_0_div_b_Upwind_r[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind_r[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind_r[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind_r[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind_r[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_r[1] = 0.447213595499958*F_0_div_b_Upwind_r[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*F_0_div_b_Upwind_r[8]+0.5000000000000001*F_0_div_b_Upwind_r[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*F_0_div_b_Upwind_r[7]+0.447213595499958*F_0_div_b_Upwind_r[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind_r[6]+0.4472135954999579*F_0_div_b_Upwind_r[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*F_0_div_b_Upwind_r[4]+0.5*F_0_div_b_Upwind_r[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind_r[3]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[1]; 
  Ghat_F_0_div_b_r[2] = 0.447213595499958*F_0_div_b_Upwind_r[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*F_0_div_b_Upwind_r[8]+0.447213595499958*F_0_div_b_Upwind_r[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*F_0_div_b_Upwind_r[7]+0.5000000000000001*F_0_div_b_Upwind_r[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*F_0_div_b_Upwind_r[6]+0.4472135954999579*F_0_div_b_Upwind_r[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*F_0_div_b_Upwind_r[5]+0.5*F_0_div_b_Upwind_r[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_r[3]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[2]; 
  Ghat_F_0_div_b_r[3] = 0.4*F_0_div_b_Upwind_r[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*F_0_div_b_Upwind_r[8]+0.4*F_0_div_b_Upwind_r[6]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind_r[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*F_0_div_b_Upwind_r[7]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_r[7]+0.447213595499958*F_0_div_b_Upwind_r[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_r[6]+0.4472135954999579*F_0_div_b_Upwind_r[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind_r[5]+0.4472135954999579*F_0_div_b_Upwind_r[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*F_0_div_b_Upwind_r[4]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[3]+0.5*F_0_div_b_Upwind_r[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_r[2]; 
  Ghat_F_0_div_b_r[4] = 0.31943828249997*F_0_div_b_Upwind_r[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_r[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*F_0_div_b_Upwind_r[8]+0.4472135954999579*F_0_div_b_Upwind_r[7]*div_b_Surf[7]+0.31943828249997*F_0_div_b_Upwind_r[6]*div_b_Surf[6]+0.5000000000000001*F_0_div_b_Upwind_r[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind_r[6]+0.31943828249997*F_0_div_b_Upwind_r[4]*div_b_Surf[4]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[4]+0.4472135954999579*F_0_div_b_Upwind_r[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind_r[1]*div_b_Surf[1]; 
  Ghat_F_0_div_b_r[5] = 0.31943828249997*F_0_div_b_Upwind_r[8]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_r[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*F_0_div_b_Upwind_r[8]+0.31943828249997*F_0_div_b_Upwind_r[7]*div_b_Surf[7]+0.5000000000000001*F_0_div_b_Upwind_r[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind_r[7]+0.4472135954999579*F_0_div_b_Upwind_r[6]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind_r[5]*div_b_Surf[5]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[5]+0.4472135954999579*F_0_div_b_Upwind_r[3]*div_b_Surf[3]+0.4472135954999579*F_0_div_b_Upwind_r[2]*div_b_Surf[2]; 
  Ghat_F_0_div_b_r[6] = 0.2857142857142857*F_0_div_b_Upwind_r[6]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind_r[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*F_0_div_b_Upwind_r[8]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_r[8]+0.4*F_0_div_b_Upwind_r[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*F_0_div_b_Upwind_r[7]+0.4472135954999579*F_0_div_b_Upwind_r[5]*div_b_Surf[6]+0.31943828249997*F_0_div_b_Upwind_r[4]*div_b_Surf[6]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*F_0_div_b_Upwind_r[6]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind_r[6]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[6]+0.5000000000000001*F_0_div_b_Upwind_r[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*F_0_div_b_Upwind_r[4]+0.447213595499958*F_0_div_b_Upwind_r[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_r[3]; 
  Ghat_F_0_div_b_r[7] = 0.2857142857142857*F_0_div_b_Upwind_r[7]*div_b_Surf[8]+0.447213595499958*F_0_div_b_Upwind_r[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*F_0_div_b_Upwind_r[8]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_r[8]+0.31943828249997*F_0_div_b_Upwind_r[5]*div_b_Surf[7]+0.4472135954999579*F_0_div_b_Upwind_r[4]*div_b_Surf[7]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind_r[7]+0.4472135954999579*div_b_Surf[4]*F_0_div_b_Upwind_r[7]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[7]+0.4*F_0_div_b_Upwind_r[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*F_0_div_b_Upwind_r[6]+0.5000000000000001*F_0_div_b_Upwind_r[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*F_0_div_b_Upwind_r[5]+0.447213595499958*F_0_div_b_Upwind_r[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_r[3]; 
  Ghat_F_0_div_b_r[8] = 0.2040816326530612*F_0_div_b_Upwind_r[8]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind_r[5]*div_b_Surf[8]+0.31943828249997*F_0_div_b_Upwind_r[4]*div_b_Surf[8]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*F_0_div_b_Upwind_r[8]+0.31943828249997*div_b_Surf[4]*F_0_div_b_Upwind_r[8]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[8]+0.2857142857142857*F_0_div_b_Upwind_r[7]*div_b_Surf[7]+0.447213595499958*F_0_div_b_Upwind_r[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*F_0_div_b_Upwind_r[7]+0.2857142857142857*F_0_div_b_Upwind_r[6]*div_b_Surf[6]+0.447213595499958*F_0_div_b_Upwind_r[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*F_0_div_b_Upwind_r[6]+0.5*F_0_div_b_Upwind_r[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*F_0_div_b_Upwind_r[5]+0.4*F_0_div_b_Upwind_r[3]*div_b_Surf[3]; 
  Ghat_G_1_div_b_r[0] = 0.5*G_1_div_b_Upwind_r[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_r[7]*div_b_Surf[7]+0.5*G_1_div_b_Upwind_r[6]*div_b_Surf[6]+0.5*G_1_div_b_Upwind_r[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind_r[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind_r[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind_r[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind_r[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_r[1] = 0.447213595499958*G_1_div_b_Upwind_r[7]*div_b_Surf[8]+0.447213595499958*div_b_Surf[7]*G_1_div_b_Upwind_r[8]+0.5000000000000001*G_1_div_b_Upwind_r[5]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[5]*G_1_div_b_Upwind_r[7]+0.447213595499958*G_1_div_b_Upwind_r[3]*div_b_Surf[6]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind_r[6]+0.4472135954999579*G_1_div_b_Upwind_r[1]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[1]*G_1_div_b_Upwind_r[4]+0.5*G_1_div_b_Upwind_r[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind_r[3]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[1]; 
  Ghat_G_1_div_b_r[2] = 0.447213595499958*G_1_div_b_Upwind_r[6]*div_b_Surf[8]+0.447213595499958*div_b_Surf[6]*G_1_div_b_Upwind_r[8]+0.447213595499958*G_1_div_b_Upwind_r[3]*div_b_Surf[7]+0.447213595499958*div_b_Surf[3]*G_1_div_b_Upwind_r[7]+0.5000000000000001*G_1_div_b_Upwind_r[4]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[4]*G_1_div_b_Upwind_r[6]+0.4472135954999579*G_1_div_b_Upwind_r[2]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[2]*G_1_div_b_Upwind_r[5]+0.5*G_1_div_b_Upwind_r[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_r[3]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[2]; 
  Ghat_G_1_div_b_r[3] = 0.4*G_1_div_b_Upwind_r[3]*div_b_Surf[8]+0.4*div_b_Surf[3]*G_1_div_b_Upwind_r[8]+0.4*G_1_div_b_Upwind_r[6]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind_r[2]*div_b_Surf[7]+0.4*div_b_Surf[6]*G_1_div_b_Upwind_r[7]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_r[7]+0.447213595499958*G_1_div_b_Upwind_r[1]*div_b_Surf[6]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_r[6]+0.4472135954999579*G_1_div_b_Upwind_r[3]*div_b_Surf[5]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind_r[5]+0.4472135954999579*G_1_div_b_Upwind_r[3]*div_b_Surf[4]+0.4472135954999579*div_b_Surf[3]*G_1_div_b_Upwind_r[4]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[3]+0.5*G_1_div_b_Upwind_r[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_r[2]; 
  Ghat_G_1_div_b_r[4] = 0.31943828249997*G_1_div_b_Upwind_r[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_r[5]*div_b_Surf[8]+0.5*div_b_Surf[5]*G_1_div_b_Upwind_r[8]+0.4472135954999579*G_1_div_b_Upwind_r[7]*div_b_Surf[7]+0.31943828249997*G_1_div_b_Upwind_r[6]*div_b_Surf[6]+0.5000000000000001*G_1_div_b_Upwind_r[2]*div_b_Surf[6]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind_r[6]+0.31943828249997*G_1_div_b_Upwind_r[4]*div_b_Surf[4]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[4]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[4]+0.4472135954999579*G_1_div_b_Upwind_r[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind_r[1]*div_b_Surf[1]; 
  Ghat_G_1_div_b_r[5] = 0.31943828249997*G_1_div_b_Upwind_r[8]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_r[4]*div_b_Surf[8]+0.5*div_b_Surf[4]*G_1_div_b_Upwind_r[8]+0.31943828249997*G_1_div_b_Upwind_r[7]*div_b_Surf[7]+0.5000000000000001*G_1_div_b_Upwind_r[1]*div_b_Surf[7]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind_r[7]+0.4472135954999579*G_1_div_b_Upwind_r[6]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind_r[5]*div_b_Surf[5]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[5]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[5]+0.4472135954999579*G_1_div_b_Upwind_r[3]*div_b_Surf[3]+0.4472135954999579*G_1_div_b_Upwind_r[2]*div_b_Surf[2]; 
  Ghat_G_1_div_b_r[6] = 0.2857142857142857*G_1_div_b_Upwind_r[6]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind_r[2]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[6]*G_1_div_b_Upwind_r[8]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_r[8]+0.4*G_1_div_b_Upwind_r[3]*div_b_Surf[7]+0.4*div_b_Surf[3]*G_1_div_b_Upwind_r[7]+0.4472135954999579*G_1_div_b_Upwind_r[5]*div_b_Surf[6]+0.31943828249997*G_1_div_b_Upwind_r[4]*div_b_Surf[6]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[6]+0.4472135954999579*div_b_Surf[5]*G_1_div_b_Upwind_r[6]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind_r[6]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[6]+0.5000000000000001*G_1_div_b_Upwind_r[2]*div_b_Surf[4]+0.5000000000000001*div_b_Surf[2]*G_1_div_b_Upwind_r[4]+0.447213595499958*G_1_div_b_Upwind_r[1]*div_b_Surf[3]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_r[3]; 
  Ghat_G_1_div_b_r[7] = 0.2857142857142857*G_1_div_b_Upwind_r[7]*div_b_Surf[8]+0.447213595499958*G_1_div_b_Upwind_r[1]*div_b_Surf[8]+0.2857142857142857*div_b_Surf[7]*G_1_div_b_Upwind_r[8]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_r[8]+0.31943828249997*G_1_div_b_Upwind_r[5]*div_b_Surf[7]+0.4472135954999579*G_1_div_b_Upwind_r[4]*div_b_Surf[7]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[7]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind_r[7]+0.4472135954999579*div_b_Surf[4]*G_1_div_b_Upwind_r[7]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[7]+0.4*G_1_div_b_Upwind_r[3]*div_b_Surf[6]+0.4*div_b_Surf[3]*G_1_div_b_Upwind_r[6]+0.5000000000000001*G_1_div_b_Upwind_r[1]*div_b_Surf[5]+0.5000000000000001*div_b_Surf[1]*G_1_div_b_Upwind_r[5]+0.447213595499958*G_1_div_b_Upwind_r[2]*div_b_Surf[3]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_r[3]; 
  Ghat_G_1_div_b_r[8] = 0.2040816326530612*G_1_div_b_Upwind_r[8]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind_r[5]*div_b_Surf[8]+0.31943828249997*G_1_div_b_Upwind_r[4]*div_b_Surf[8]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[8]+0.31943828249997*div_b_Surf[5]*G_1_div_b_Upwind_r[8]+0.31943828249997*div_b_Surf[4]*G_1_div_b_Upwind_r[8]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[8]+0.2857142857142857*G_1_div_b_Upwind_r[7]*div_b_Surf[7]+0.447213595499958*G_1_div_b_Upwind_r[1]*div_b_Surf[7]+0.447213595499958*div_b_Surf[1]*G_1_div_b_Upwind_r[7]+0.2857142857142857*G_1_div_b_Upwind_r[6]*div_b_Surf[6]+0.447213595499958*G_1_div_b_Upwind_r[2]*div_b_Surf[6]+0.447213595499958*div_b_Surf[2]*G_1_div_b_Upwind_r[6]+0.5*G_1_div_b_Upwind_r[4]*div_b_Surf[5]+0.5*div_b_Surf[4]*G_1_div_b_Upwind_r[5]+0.4*G_1_div_b_Upwind_r[3]*div_b_Surf[3]; 

  out_F_0[0] += ((-0.7071067811865475*Ghat_F_0_r[0])+0.7071067811865475*Ghat_F_0_l[0]-0.7071067811865475*Ghat_F_0_div_b_r[0]+0.7071067811865475*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[1] += ((-0.7071067811865475*Ghat_F_0_r[1])+0.7071067811865475*Ghat_F_0_l[1]-0.7071067811865475*Ghat_F_0_div_b_r[1]+0.7071067811865475*Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[2] += ((-0.7071067811865475*Ghat_F_0_r[2])+0.7071067811865475*Ghat_F_0_l[2]-0.7071067811865475*Ghat_F_0_div_b_r[2]+0.7071067811865475*Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[3] += -1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0]+Ghat_F_0_div_b_r[0]+Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[4] += ((-0.7071067811865475*Ghat_F_0_r[3])+0.7071067811865475*Ghat_F_0_l[3]-0.7071067811865475*Ghat_F_0_div_b_r[3]+0.7071067811865475*Ghat_F_0_div_b_l[3])*dv1par; 
  out_F_0[5] += -1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1]+Ghat_F_0_div_b_r[1]+Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[6] += -1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2]+Ghat_F_0_div_b_r[2]+Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[7] += ((-0.7071067811865475*Ghat_F_0_r[4])+0.7071067811865475*Ghat_F_0_l[4]-0.7071067811865475*Ghat_F_0_div_b_r[4]+0.7071067811865475*Ghat_F_0_div_b_l[4])*dv1par; 
  out_F_0[8] += ((-0.7071067811865475*Ghat_F_0_r[5])+0.7071067811865475*Ghat_F_0_l[5]-0.7071067811865475*Ghat_F_0_div_b_r[5]+0.7071067811865475*Ghat_F_0_div_b_l[5])*dv1par; 
  out_F_0[9] += ((-1.58113883008419*Ghat_F_0_r[0])+1.58113883008419*Ghat_F_0_l[0]-1.58113883008419*Ghat_F_0_div_b_r[0]+1.58113883008419*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[10] += -1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3]+Ghat_F_0_div_b_r[3]+Ghat_F_0_div_b_l[3])*dv1par; 
  out_F_0[11] += ((-0.7071067811865475*Ghat_F_0_r[6])+0.7071067811865475*Ghat_F_0_l[6]-0.7071067811865475*Ghat_F_0_div_b_r[6]+0.7071067811865475*Ghat_F_0_div_b_l[6])*dv1par; 
  out_F_0[12] += ((-0.7071067811865475*Ghat_F_0_r[7])+0.7071067811865475*Ghat_F_0_l[7]-0.7071067811865475*Ghat_F_0_div_b_r[7]+0.7071067811865475*Ghat_F_0_div_b_l[7])*dv1par; 
  out_F_0[13] += -1.224744871391589*(Ghat_F_0_r[4]+Ghat_F_0_l[4]+Ghat_F_0_div_b_r[4]+Ghat_F_0_div_b_l[4])*dv1par; 
  out_F_0[14] += -1.224744871391589*(Ghat_F_0_r[5]+Ghat_F_0_l[5]+Ghat_F_0_div_b_r[5]+Ghat_F_0_div_b_l[5])*dv1par; 
  out_F_0[15] += ((-1.58113883008419*Ghat_F_0_r[1])+1.58113883008419*Ghat_F_0_l[1]-1.58113883008419*Ghat_F_0_div_b_r[1]+1.58113883008419*Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[16] += ((-1.58113883008419*Ghat_F_0_r[2])+1.58113883008419*Ghat_F_0_l[2]-1.58113883008419*Ghat_F_0_div_b_r[2]+1.58113883008419*Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[17] += -1.224744871391589*(Ghat_F_0_r[6]+Ghat_F_0_l[6]+Ghat_F_0_div_b_r[6]+Ghat_F_0_div_b_l[6])*dv1par; 
  out_F_0[18] += -1.224744871391589*(Ghat_F_0_r[7]+Ghat_F_0_l[7]+Ghat_F_0_div_b_r[7]+Ghat_F_0_div_b_l[7])*dv1par; 
  out_F_0[19] += ((-1.58113883008419*Ghat_F_0_r[3])+1.58113883008419*Ghat_F_0_l[3]-1.58113883008419*Ghat_F_0_div_b_r[3]+1.58113883008419*Ghat_F_0_div_b_l[3])*dv1par; 
  out_F_0[20] += ((-0.7071067811865475*Ghat_F_0_r[8])+0.7071067811865475*Ghat_F_0_l[8]-0.7071067811865475*Ghat_F_0_div_b_r[8]+0.7071067811865475*Ghat_F_0_div_b_l[8])*dv1par; 
  out_F_0[21] += ((-1.58113883008419*Ghat_F_0_r[4])+1.58113883008419*Ghat_F_0_l[4]-1.58113883008419*Ghat_F_0_div_b_r[4]+1.58113883008419*Ghat_F_0_div_b_l[4])*dv1par; 
  out_F_0[22] += ((-1.58113883008419*Ghat_F_0_r[5])+1.58113883008419*Ghat_F_0_l[5]-1.58113883008419*Ghat_F_0_div_b_r[5]+1.58113883008419*Ghat_F_0_div_b_l[5])*dv1par; 
  out_F_0[23] += -1.224744871391589*(Ghat_F_0_r[8]+Ghat_F_0_l[8]+Ghat_F_0_div_b_r[8]+Ghat_F_0_div_b_l[8])*dv1par; 
  out_F_0[24] += ((-1.58113883008419*Ghat_F_0_r[6])+1.58113883008419*Ghat_F_0_l[6]-1.58113883008419*Ghat_F_0_div_b_r[6]+1.58113883008419*Ghat_F_0_div_b_l[6])*dv1par; 
  out_F_0[25] += ((-1.58113883008419*Ghat_F_0_r[7])+1.58113883008419*Ghat_F_0_l[7]-1.58113883008419*Ghat_F_0_div_b_r[7]+1.58113883008419*Ghat_F_0_div_b_l[7])*dv1par; 
  out_F_0[26] += ((-1.58113883008419*Ghat_F_0_r[8])+1.58113883008419*Ghat_F_0_l[8]-1.58113883008419*Ghat_F_0_div_b_r[8]+1.58113883008419*Ghat_F_0_div_b_l[8])*dv1par; 
  out_G_1[0] += ((-0.7071067811865475*Ghat_G_1_r[0])+0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_div_b_r[0]+0.7071067811865475*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[1] += ((-0.7071067811865475*Ghat_G_1_r[1])+0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_div_b_r[1]+0.7071067811865475*Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[2] += ((-0.7071067811865475*Ghat_G_1_r[2])+0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_div_b_r[2]+0.7071067811865475*Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[3] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0]+Ghat_G_1_div_b_r[0]+Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[4] += ((-0.7071067811865475*Ghat_G_1_r[3])+0.7071067811865475*Ghat_G_1_l[3]-0.7071067811865475*Ghat_G_1_div_b_r[3]+0.7071067811865475*Ghat_G_1_div_b_l[3])*dv1par; 
  out_G_1[5] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1]+Ghat_G_1_div_b_r[1]+Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[6] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2]+Ghat_G_1_div_b_r[2]+Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[7] += ((-0.7071067811865475*Ghat_G_1_r[4])+0.7071067811865475*Ghat_G_1_l[4]-0.7071067811865475*Ghat_G_1_div_b_r[4]+0.7071067811865475*Ghat_G_1_div_b_l[4])*dv1par; 
  out_G_1[8] += ((-0.7071067811865475*Ghat_G_1_r[5])+0.7071067811865475*Ghat_G_1_l[5]-0.7071067811865475*Ghat_G_1_div_b_r[5]+0.7071067811865475*Ghat_G_1_div_b_l[5])*dv1par; 
  out_G_1[9] += ((-1.58113883008419*Ghat_G_1_r[0])+1.58113883008419*Ghat_G_1_l[0]-1.58113883008419*Ghat_G_1_div_b_r[0]+1.58113883008419*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[10] += -1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3]+Ghat_G_1_div_b_r[3]+Ghat_G_1_div_b_l[3])*dv1par; 
  out_G_1[11] += ((-0.7071067811865475*Ghat_G_1_r[6])+0.7071067811865475*Ghat_G_1_l[6]-0.7071067811865475*Ghat_G_1_div_b_r[6]+0.7071067811865475*Ghat_G_1_div_b_l[6])*dv1par; 
  out_G_1[12] += ((-0.7071067811865475*Ghat_G_1_r[7])+0.7071067811865475*Ghat_G_1_l[7]-0.7071067811865475*Ghat_G_1_div_b_r[7]+0.7071067811865475*Ghat_G_1_div_b_l[7])*dv1par; 
  out_G_1[13] += -1.224744871391589*(Ghat_G_1_r[4]+Ghat_G_1_l[4]+Ghat_G_1_div_b_r[4]+Ghat_G_1_div_b_l[4])*dv1par; 
  out_G_1[14] += -1.224744871391589*(Ghat_G_1_r[5]+Ghat_G_1_l[5]+Ghat_G_1_div_b_r[5]+Ghat_G_1_div_b_l[5])*dv1par; 
  out_G_1[15] += ((-1.58113883008419*Ghat_G_1_r[1])+1.58113883008419*Ghat_G_1_l[1]-1.58113883008419*Ghat_G_1_div_b_r[1]+1.58113883008419*Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[16] += ((-1.58113883008419*Ghat_G_1_r[2])+1.58113883008419*Ghat_G_1_l[2]-1.58113883008419*Ghat_G_1_div_b_r[2]+1.58113883008419*Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[17] += -1.224744871391589*(Ghat_G_1_r[6]+Ghat_G_1_l[6]+Ghat_G_1_div_b_r[6]+Ghat_G_1_div_b_l[6])*dv1par; 
  out_G_1[18] += -1.224744871391589*(Ghat_G_1_r[7]+Ghat_G_1_l[7]+Ghat_G_1_div_b_r[7]+Ghat_G_1_div_b_l[7])*dv1par; 
  out_G_1[19] += ((-1.58113883008419*Ghat_G_1_r[3])+1.58113883008419*Ghat_G_1_l[3]-1.58113883008419*Ghat_G_1_div_b_r[3]+1.58113883008419*Ghat_G_1_div_b_l[3])*dv1par; 
  out_G_1[20] += ((-0.7071067811865475*Ghat_G_1_r[8])+0.7071067811865475*Ghat_G_1_l[8]-0.7071067811865475*Ghat_G_1_div_b_r[8]+0.7071067811865475*Ghat_G_1_div_b_l[8])*dv1par; 
  out_G_1[21] += ((-1.58113883008419*Ghat_G_1_r[4])+1.58113883008419*Ghat_G_1_l[4]-1.58113883008419*Ghat_G_1_div_b_r[4]+1.58113883008419*Ghat_G_1_div_b_l[4])*dv1par; 
  out_G_1[22] += ((-1.58113883008419*Ghat_G_1_r[5])+1.58113883008419*Ghat_G_1_l[5]-1.58113883008419*Ghat_G_1_div_b_r[5]+1.58113883008419*Ghat_G_1_div_b_l[5])*dv1par; 
  out_G_1[23] += -1.224744871391589*(Ghat_G_1_r[8]+Ghat_G_1_l[8]+Ghat_G_1_div_b_r[8]+Ghat_G_1_div_b_l[8])*dv1par; 
  out_G_1[24] += ((-1.58113883008419*Ghat_G_1_r[6])+1.58113883008419*Ghat_G_1_l[6]-1.58113883008419*Ghat_G_1_div_b_r[6]+1.58113883008419*Ghat_G_1_div_b_l[6])*dv1par; 
  out_G_1[25] += ((-1.58113883008419*Ghat_G_1_r[7])+1.58113883008419*Ghat_G_1_l[7]-1.58113883008419*Ghat_G_1_div_b_r[7]+1.58113883008419*Ghat_G_1_div_b_l[7])*dv1par; 
  out_G_1[26] += ((-1.58113883008419*Ghat_G_1_r[8])+1.58113883008419*Ghat_G_1_l[8]-1.58113883008419*Ghat_G_1_div_b_r[8]+1.58113883008419*Ghat_G_1_div_b_l[8])*dv1par; 

  return 2.5*dv1par*cflFreq;

} 
