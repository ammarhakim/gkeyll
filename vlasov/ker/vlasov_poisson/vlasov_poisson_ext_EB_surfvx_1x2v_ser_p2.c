#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 

  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dx10 = 2/dxv[0]; 

  double alpha[8] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ex = &EBext[0]; 
  const double *Ey = &EBext[3]; 
  const double *Bx = &EBext[9]; 
  const double *By = &EBext[12]; 
  const double *Bz = &EBext[15]; 

  alpha[0] = 1.4142135623730951*Bz[0]*wv2-2.4494897427831783*phi[1]*dx10+1.4142135623730951*Ex[0]; 
  alpha[1] = 1.4142135623730951*Bz[1]*wv2-5.477225575051662*phi[2]*dx10+1.4142135623730951*Ex[1]; 
  alpha[2] = 0.408248290463863*Bz[0]*dv2; 
  alpha[3] = 0.408248290463863*Bz[1]*dv2; 
  alpha[4] = 1.4142135623730951*Bz[2]*wv2+1.4142135623730951*Ex[2]; 
  alpha[6] = 0.40824829046386296*Bz[2]*dv2; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[8] = {0.0};;
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if (-(0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fr); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fr); 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fr); 
  } 
  if (0.5*alpha[0]-0.5590169943749475*alpha[4] > 0) { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fr); 
  } 
  if (-(0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fr); 
  } 
  if (-(0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fr); 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fr); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alpha[6]*fUpwind_l[6]+0.5*alpha[4]*fUpwind_l[4]+0.5*alpha[3]*fUpwind_l[3]+0.5*alpha[2]*fUpwind_l[2]+0.5*alpha[1]*fUpwind_l[1]+0.5*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.44721359549995804*alpha[3]*fUpwind_l[6]+0.44721359549995804*fUpwind_l[3]*alpha[6]+0.4472135954999579*alpha[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alpha[4]+0.5*alpha[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alpha[3]+0.5*alpha[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.44721359549995804*alpha[3]*fUpwind_l[7]+0.5000000000000001*alpha[4]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[4]*alpha[6]+0.4472135954999579*alpha[2]*fUpwind_l[5]+0.5*alpha[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alpha[3]+0.5*alpha[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.4*alpha[6]*fUpwind_l[7]+0.44721359549995804*alpha[2]*fUpwind_l[7]+0.44721359549995804*alpha[1]*fUpwind_l[6]+0.44721359549995804*fUpwind_l[1]*alpha[6]+0.4472135954999579*alpha[3]*fUpwind_l[5]+0.4472135954999579*alpha[3]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[3]*alpha[4]+0.5*alpha[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alpha[3]+0.5*alpha[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alpha[2]; 
  Ghat_l[4] = 0.31943828249996997*alpha[6]*fUpwind_l[6]+0.5000000000000001*alpha[2]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alpha[6]+0.31943828249996997*alpha[4]*fUpwind_l[4]+0.5*alpha[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alpha[4]+0.4472135954999579*alpha[3]*fUpwind_l[3]+0.4472135954999579*alpha[1]*fUpwind_l[1]; 
  Ghat_l[5] = 0.5000000000000001*alpha[1]*fUpwind_l[7]+0.4472135954999579*alpha[6]*fUpwind_l[6]+0.5*alpha[0]*fUpwind_l[5]+0.4472135954999579*alpha[3]*fUpwind_l[3]+0.4472135954999579*alpha[2]*fUpwind_l[2]; 
  Ghat_l[6] = 0.4*alpha[3]*fUpwind_l[7]+0.31943828249996997*alpha[4]*fUpwind_l[6]+0.5*alpha[0]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[5]*alpha[6]+0.31943828249996997*fUpwind_l[4]*alpha[6]+0.5*fUpwind_l[0]*alpha[6]+0.5000000000000001*alpha[2]*fUpwind_l[4]+0.5000000000000001*fUpwind_l[2]*alpha[4]+0.44721359549995804*alpha[1]*fUpwind_l[3]+0.44721359549995804*fUpwind_l[1]*alpha[3]; 
  Ghat_l[7] = 0.4472135954999579*alpha[4]*fUpwind_l[7]+0.5*alpha[0]*fUpwind_l[7]+0.4*alpha[3]*fUpwind_l[6]+0.4*fUpwind_l[3]*alpha[6]+0.5000000000000001*alpha[1]*fUpwind_l[5]+0.44721359549995804*alpha[2]*fUpwind_l[3]+0.44721359549995804*fUpwind_l[2]*alpha[3]; 

  Ghat_r[0] = 0.5*alpha[6]*fUpwind_r[6]+0.5*alpha[4]*fUpwind_r[4]+0.5*alpha[3]*fUpwind_r[3]+0.5*alpha[2]*fUpwind_r[2]+0.5*alpha[1]*fUpwind_r[1]+0.5*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.44721359549995804*alpha[3]*fUpwind_r[6]+0.44721359549995804*fUpwind_r[3]*alpha[6]+0.4472135954999579*alpha[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alpha[4]+0.5*alpha[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alpha[3]+0.5*alpha[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.44721359549995804*alpha[3]*fUpwind_r[7]+0.5000000000000001*alpha[4]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[4]*alpha[6]+0.4472135954999579*alpha[2]*fUpwind_r[5]+0.5*alpha[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alpha[3]+0.5*alpha[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.4*alpha[6]*fUpwind_r[7]+0.44721359549995804*alpha[2]*fUpwind_r[7]+0.44721359549995804*alpha[1]*fUpwind_r[6]+0.44721359549995804*fUpwind_r[1]*alpha[6]+0.4472135954999579*alpha[3]*fUpwind_r[5]+0.4472135954999579*alpha[3]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[3]*alpha[4]+0.5*alpha[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alpha[3]+0.5*alpha[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alpha[2]; 
  Ghat_r[4] = 0.31943828249996997*alpha[6]*fUpwind_r[6]+0.5000000000000001*alpha[2]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alpha[6]+0.31943828249996997*alpha[4]*fUpwind_r[4]+0.5*alpha[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alpha[4]+0.4472135954999579*alpha[3]*fUpwind_r[3]+0.4472135954999579*alpha[1]*fUpwind_r[1]; 
  Ghat_r[5] = 0.5000000000000001*alpha[1]*fUpwind_r[7]+0.4472135954999579*alpha[6]*fUpwind_r[6]+0.5*alpha[0]*fUpwind_r[5]+0.4472135954999579*alpha[3]*fUpwind_r[3]+0.4472135954999579*alpha[2]*fUpwind_r[2]; 
  Ghat_r[6] = 0.4*alpha[3]*fUpwind_r[7]+0.31943828249996997*alpha[4]*fUpwind_r[6]+0.5*alpha[0]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[5]*alpha[6]+0.31943828249996997*fUpwind_r[4]*alpha[6]+0.5*fUpwind_r[0]*alpha[6]+0.5000000000000001*alpha[2]*fUpwind_r[4]+0.5000000000000001*fUpwind_r[2]*alpha[4]+0.44721359549995804*alpha[1]*fUpwind_r[3]+0.44721359549995804*fUpwind_r[1]*alpha[3]; 
  Ghat_r[7] = 0.4472135954999579*alpha[4]*fUpwind_r[7]+0.5*alpha[0]*fUpwind_r[7]+0.4*alpha[3]*fUpwind_r[6]+0.4*fUpwind_r[3]*alpha[6]+0.5000000000000001*alpha[1]*fUpwind_r[5]+0.44721359549995804*alpha[2]*fUpwind_r[3]+0.44721359549995804*fUpwind_r[2]*alpha[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10); 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10); 
  out[5] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[6] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10); 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[8] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[10] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10); 
  out[11] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10); 
  out[12] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[14] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv10; 
  out[15] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[16] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10); 
  out[17] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10); 
  out[18] += (1.5811388300841895*Ghat_l[3]-1.5811388300841895*Ghat_r[3])*dv10; 
  out[19] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10); 

  return 0.;

} 
