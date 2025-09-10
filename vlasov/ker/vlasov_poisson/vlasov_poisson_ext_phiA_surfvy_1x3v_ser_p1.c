#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 

  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double dx10 = 2/dxv[0]; 

  double alpha[16] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ax = &pots[2]; 
  const double *Ay = &pots[4]; 
  const double *Az = &pots[6]; 

  alpha[0] = -(3.4641016151377544*Ay[1]*dx10*wv1); 
  alpha[2] = -(1.0*Ay[1]*dv1*dx10); 

  double fUpwindQuad_l[18] = {0.0};
  double fUpwindQuad_r[18] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad_l[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad_l[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad_l[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_l(fr); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_l(fr); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad_l[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad_l[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad_l[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_l(fr); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_l(fr); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_l(fr); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_l(fr); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_l(fr); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]; 
  Ghat_l[2] = 0.3162277660168379*alpha[2]*fUpwind_l[8]+0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*alpha[0]*fUpwind_l[3]; 
  Ghat_l[4] = 0.31622776601683794*alpha[2]*fUpwind_l[9]+0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[5]; 
  Ghat_l[6] = 0.31622776601683794*alpha[2]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*alpha[2]*fUpwind_l[3]; 
  Ghat_l[7] = 0.3162277660168379*alpha[2]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*alpha[2]*fUpwind_l[5]; 
  Ghat_l[8] = 0.3535533905932737*alpha[0]*fUpwind_l[8]+0.3162277660168379*alpha[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha[0]*fUpwind_l[9]+0.31622776601683794*alpha[2]*fUpwind_l[4]; 
  Ghat_l[10] = 0.3535533905932737*alpha[0]*fUpwind_l[10]+0.31622776601683794*alpha[2]*fUpwind_l[6]; 
  Ghat_l[11] = 0.3535533905932737*alpha[0]*fUpwind_l[11]+0.3162277660168379*alpha[2]*fUpwind_l[7]; 
  Ghat_l[12] = 0.3535533905932737*alpha[2]*fUpwind_l[14]+0.3535533905932737*alpha[0]*fUpwind_l[12]; 
  Ghat_l[13] = 0.3535533905932737*alpha[2]*fUpwind_l[15]+0.3535533905932737*alpha[0]*fUpwind_l[13]; 
  Ghat_l[14] = 0.3535533905932737*alpha[0]*fUpwind_l[14]+0.3535533905932737*alpha[2]*fUpwind_l[12]; 
  Ghat_l[15] = 0.3535533905932737*alpha[0]*fUpwind_l[15]+0.3535533905932737*alpha[2]*fUpwind_l[13]; 

  Ghat_r[0] = 0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]; 
  Ghat_r[2] = 0.3162277660168379*alpha[2]*fUpwind_r[8]+0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*alpha[0]*fUpwind_r[3]; 
  Ghat_r[4] = 0.31622776601683794*alpha[2]*fUpwind_r[9]+0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[5]; 
  Ghat_r[6] = 0.31622776601683794*alpha[2]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*alpha[2]*fUpwind_r[3]; 
  Ghat_r[7] = 0.3162277660168379*alpha[2]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*alpha[2]*fUpwind_r[5]; 
  Ghat_r[8] = 0.3535533905932737*alpha[0]*fUpwind_r[8]+0.3162277660168379*alpha[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha[0]*fUpwind_r[9]+0.31622776601683794*alpha[2]*fUpwind_r[4]; 
  Ghat_r[10] = 0.3535533905932737*alpha[0]*fUpwind_r[10]+0.31622776601683794*alpha[2]*fUpwind_r[6]; 
  Ghat_r[11] = 0.3535533905932737*alpha[0]*fUpwind_r[11]+0.3162277660168379*alpha[2]*fUpwind_r[7]; 
  Ghat_r[12] = 0.3535533905932737*alpha[2]*fUpwind_r[14]+0.3535533905932737*alpha[0]*fUpwind_r[12]; 
  Ghat_r[13] = 0.3535533905932737*alpha[2]*fUpwind_r[15]+0.3535533905932737*alpha[0]*fUpwind_r[13]; 
  Ghat_r[14] = 0.3535533905932737*alpha[0]*fUpwind_r[14]+0.3535533905932737*alpha[2]*fUpwind_r[12]; 
  Ghat_r[15] = 0.3535533905932737*alpha[0]*fUpwind_r[15]+0.3535533905932737*alpha[2]*fUpwind_r[13]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11); 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11); 
  out[7] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11); 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[10] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11); 
  out[11] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11); 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[13] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11); 
  out[14] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11); 
  out[15] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11); 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11; 
  out[17] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11; 
  out[18] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11); 
  out[19] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11; 
  out[20] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11); 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11; 
  out[22] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11); 
  out[23] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11); 
  out[24] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv11; 
  out[25] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv11; 
  out[26] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv11; 
  out[27] += (1.5811388300841898*Ghat_l[3]-1.5811388300841898*Ghat_r[3])*dv11; 
  out[28] += (1.5811388300841895*Ghat_l[4]-1.5811388300841895*Ghat_r[4])*dv11; 
  out[29] += (1.5811388300841895*Ghat_l[5]-1.5811388300841895*Ghat_r[5])*dv11; 
  out[30] += (1.5811388300841895*Ghat_l[6]-1.5811388300841895*Ghat_r[6])*dv11; 
  out[31] += (1.5811388300841898*Ghat_l[7]-1.5811388300841898*Ghat_r[7])*dv11; 
  out[32] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11; 
  out[33] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11; 
  out[34] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11; 
  out[35] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11); 
  out[36] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11; 
  out[37] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11); 
  out[38] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11); 
  out[39] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11); 

  return 0.;

} 
