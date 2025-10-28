#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double dx10 = 2/dxv[0]; 

  double alpha[20] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ex = &EBext[0]; 
  const double *Ey = &EBext[3]; 
  const double *Ez = &EBext[6]; 
  const double *Bx = &EBext[9]; 
  const double *By = &EBext[12]; 
  const double *Bz = &EBext[15]; 

  alpha[0] = -(2.0*By[0]*wv3)+2.0*Bz[0]*wv2-3.4641016151377544*phi[1]*dx10+2.0*Ex[0]; 
  alpha[1] = -(2.0*By[1]*wv3)+2.0*Bz[1]*wv2-7.745966692414834*phi[2]*dx10+2.0*Ex[1]; 
  alpha[2] = 0.5773502691896258*Bz[0]*dv2; 
  alpha[3] = -(0.5773502691896258*By[0]*dv3); 
  alpha[4] = 0.5773502691896258*Bz[1]*dv2; 
  alpha[5] = -(0.5773502691896258*By[1]*dv3); 
  alpha[7] = -(2.0*By[2]*wv3)+2.0*Bz[2]*wv2+2.0*Ex[2]; 
  alpha[11] = 0.5773502691896257*Bz[2]*dv2; 
  alpha[13] = -(0.5773502691896257*By[2]*dv3); 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (-(0.42426406871192845*(alpha[13]+alpha[11]))+0.3162277660168379*alpha[7]+0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx2_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx2_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[11])+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx2_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx2_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx2_eval_quad_node_1_l(fr); 
  } 
  if (0.42426406871192845*alpha[13]-0.42426406871192845*alpha[11]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx2_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx2_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx2_eval_quad_node_2_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[13])+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx2_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx2_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx2_eval_quad_node_3_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx2_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx2_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx2_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx2_eval_quad_node_4_l(fr); 
  } 
  if (0.42426406871192845*alpha[13]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx2_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx2_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx2_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx2_eval_quad_node_5_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[13])+0.42426406871192845*alpha[11]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx2_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx2_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx2_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx2_eval_quad_node_6_l(fr); 
  } 
  if (0.42426406871192845*alpha[11]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx2_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx2_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx2_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx2_eval_quad_node_7_l(fr); 
  } 
  if (0.42426406871192845*(alpha[13]+alpha[11])+0.3162277660168379*alpha[7]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx2_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx2_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx2_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx2_eval_quad_node_8_l(fr); 
  } 
  if (0.5303300858899104*(alpha[13]+alpha[11])-0.3952847075210473*alpha[7]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx2_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx2_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx2_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx2_eval_quad_node_9_l(fr); 
  } 
  if (0.5303300858899104*alpha[11]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx2_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx2_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx2_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx2_eval_quad_node_10_l(fr); 
  } 
  if (-(0.5303300858899104*alpha[13])+0.5303300858899104*alpha[11]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx2_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx2_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx2_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx2_eval_quad_node_11_l(fr); 
  } 
  if (0.5303300858899104*alpha[13]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx2_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx2_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx2_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx2_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7] > 0) { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx2_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx2_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx2_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx2_eval_quad_node_13_l(fr); 
  } 
  if (-(0.5303300858899104*alpha[13])-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx2_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx2_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx2_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx2_eval_quad_node_14_l(fr); 
  } 
  if (0.5303300858899104*alpha[13]-0.5303300858899104*alpha[11]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx2_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx2_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx2_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx2_eval_quad_node_15_l(fr); 
  } 
  if (-(0.5303300858899104*alpha[11])-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx2_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx2_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx2_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx2_eval_quad_node_16_l(fr); 
  } 
  if (-(0.5303300858899104*(alpha[13]+alpha[11]))-0.3952847075210473*alpha[7]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx2_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx2_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx2_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx2_eval_quad_node_17_l(fr); 
  } 
  if (-(0.42426406871192845*(alpha[13]+alpha[11]))+0.3162277660168379*alpha[7]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx2_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx2_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx2_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx2_eval_quad_node_18_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[11])+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx2_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx2_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx2_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx2_eval_quad_node_19_l(fr); 
  } 
  if (0.42426406871192845*alpha[13]-0.42426406871192845*alpha[11]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx2_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx2_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx2_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx2_eval_quad_node_20_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[13])+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx2_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx2_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx2_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx2_eval_quad_node_21_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx2_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx2_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx2_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx2_eval_quad_node_22_l(fr); 
  } 
  if (0.42426406871192845*alpha[13]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx2_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx2_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx2_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx2_eval_quad_node_23_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[13])+0.42426406871192845*alpha[11]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx2_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx2_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx2_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx2_eval_quad_node_24_l(fr); 
  } 
  if (0.42426406871192845*alpha[11]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx2_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx2_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx2_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx2_eval_quad_node_25_l(fr); 
  } 
  if (0.42426406871192845*(alpha[13]+alpha[11])+0.3162277660168379*alpha[7]+0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx2_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx2_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx2_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx2_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[5]*fUpwind_l[5]+0.3535533905932737*alpha[4]*fUpwind_l[4]+0.3535533905932737*alpha[3]*fUpwind_l[3]+0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[1]*fUpwind_l[1]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.31622776601683794*alpha[5]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[5]*alpha[13]+0.31622776601683794*alpha[4]*fUpwind_l[11]+0.31622776601683794*fUpwind_l[4]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.3535533905932737*alpha[13]*fUpwind_l[17]+0.31622776601683794*alpha[4]*fUpwind_l[12]+0.3535533905932737*alpha[7]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[7]*alpha[11]+0.3535533905932737*alpha[5]*fUpwind_l[10]+0.3162277660168379*alpha[2]*fUpwind_l[8]+0.3535533905932737*alpha[3]*fUpwind_l[6]+0.3535533905932737*alpha[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3535533905932737*alpha[11]*fUpwind_l[17]+0.31622776601683794*alpha[5]*fUpwind_l[15]+0.3535533905932737*alpha[7]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_l[10]+0.3162277660168379*alpha[3]*fUpwind_l[9]+0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*alpha[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.3162277660168379*alpha[5]*fUpwind_l[17]+0.31622776601683794*fUpwind_l[10]*alpha[13]+0.28284271247461906*alpha[11]*fUpwind_l[12]+0.31622776601683794*alpha[2]*fUpwind_l[12]+0.31622776601683794*alpha[1]*fUpwind_l[11]+0.31622776601683794*fUpwind_l[1]*alpha[11]+0.3535533905932737*alpha[3]*fUpwind_l[10]+0.3162277660168379*alpha[4]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_l[6]+0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.3162277660168379*alpha[4]*fUpwind_l[17]+0.28284271247461906*alpha[13]*fUpwind_l[15]+0.31622776601683794*alpha[3]*fUpwind_l[15]+0.31622776601683794*alpha[1]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[1]*alpha[13]+0.31622776601683794*fUpwind_l[10]*alpha[11]+0.3535533905932737*alpha[2]*fUpwind_l[10]+0.3162277660168379*alpha[5]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_l[6]+0.3535533905932737*alpha[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alpha[3]; 
  Ghat_l[6] = 0.3162277660168379*alpha[5]*fUpwind_l[19]+0.3162277660168379*alpha[4]*fUpwind_l[18]+0.3535533905932737*alpha[7]*fUpwind_l[17]+0.31622776601683794*alpha[3]*fUpwind_l[16]+0.31622776601683794*alpha[2]*fUpwind_l[14]+0.3535533905932737*alpha[11]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[11]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*alpha[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alpha[3]; 
  Ghat_l[7] = 0.22587697572631277*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[3]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alpha[13]+0.22587697572631277*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[2]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2]*alpha[11]+0.22587697572631277*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[1]*fUpwind_l[1]; 
  Ghat_l[8] = 0.3535533905932737*alpha[5]*fUpwind_l[18]+0.3535533905932737*alpha[3]*fUpwind_l[14]+0.3535533905932737*alpha[1]*fUpwind_l[12]+0.3162277660168379*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha[4]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[16]+0.3535533905932737*alpha[1]*fUpwind_l[15]+0.3162277660168379*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[0]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[3]*fUpwind_l[3]; 
  Ghat_l[10] = 0.282842712474619*alpha[13]*fUpwind_l[19]+0.3162277660168379*alpha[3]*fUpwind_l[19]+0.282842712474619*alpha[11]*fUpwind_l[18]+0.3162277660168379*alpha[2]*fUpwind_l[18]+0.3162277660168379*alpha[1]*fUpwind_l[17]+0.31622776601683794*alpha[5]*fUpwind_l[16]+0.31622776601683794*alpha[4]*fUpwind_l[14]+0.31622776601683794*alpha[4]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[4]*alpha[13]+0.31622776601683794*alpha[5]*fUpwind_l[11]+0.31622776601683794*fUpwind_l[5]*alpha[11]+0.3162277660168379*alpha[7]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[10]+0.3535533905932737*alpha[1]*fUpwind_l[6]+0.3535533905932737*alpha[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.22587697572631277*alpha[13]*fUpwind_l[17]+0.3535533905932737*alpha[3]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[6]*alpha[13]+0.28284271247461906*alpha[4]*fUpwind_l[12]+0.22587697572631277*alpha[7]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[8]*alpha[11]+0.22587697572631277*fUpwind_l[7]*alpha[11]+0.3535533905932737*fUpwind_l[0]*alpha[11]+0.31622776601683794*alpha[5]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_l[4]+0.31622776601683794*fUpwind_l[1]*alpha[4]; 
  Ghat_l[12] = 0.3162277660168379*alpha[13]*fUpwind_l[18]+0.3535533905932737*alpha[3]*fUpwind_l[18]+0.3535533905932737*alpha[5]*fUpwind_l[14]+0.3162277660168379*alpha[7]*fUpwind_l[12]+0.3535533905932737*alpha[0]*fUpwind_l[12]+0.28284271247461906*alpha[4]*fUpwind_l[11]+0.28284271247461906*fUpwind_l[4]*alpha[11]+0.3535533905932737*alpha[1]*fUpwind_l[8]+0.31622776601683794*alpha[2]*fUpwind_l[4]+0.31622776601683794*fUpwind_l[2]*alpha[4]; 
  Ghat_l[13] = 0.22587697572631277*alpha[11]*fUpwind_l[17]+0.3535533905932737*alpha[2]*fUpwind_l[17]+0.28284271247461906*alpha[5]*fUpwind_l[15]+0.22587697572631277*alpha[7]*fUpwind_l[13]+0.3535533905932737*alpha[0]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[9]*alpha[13]+0.22587697572631277*fUpwind_l[7]*alpha[13]+0.3535533905932737*fUpwind_l[0]*alpha[13]+0.3535533905932737*fUpwind_l[6]*alpha[11]+0.31622776601683794*alpha[4]*fUpwind_l[10]+0.3535533905932737*alpha[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_l[5]+0.31622776601683794*fUpwind_l[1]*alpha[5]; 
  Ghat_l[14] = 0.3535533905932737*alpha[1]*fUpwind_l[18]+0.3162277660168379*alpha[11]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[14]+0.3535533905932737*alpha[5]*fUpwind_l[12]+0.31622776601683794*alpha[4]*fUpwind_l[10]+0.3535533905932737*alpha[3]*fUpwind_l[8]+0.31622776601683794*alpha[2]*fUpwind_l[6]; 
  Ghat_l[15] = 0.3162277660168379*alpha[11]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[19]+0.3535533905932737*alpha[4]*fUpwind_l[16]+0.3162277660168379*alpha[7]*fUpwind_l[15]+0.3535533905932737*alpha[0]*fUpwind_l[15]+0.28284271247461906*alpha[5]*fUpwind_l[13]+0.28284271247461906*fUpwind_l[5]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_l[9]+0.31622776601683794*alpha[3]*fUpwind_l[5]+0.31622776601683794*fUpwind_l[3]*alpha[5]; 
  Ghat_l[16] = 0.3535533905932737*alpha[1]*fUpwind_l[19]+0.3162277660168379*alpha[13]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[16]+0.3535533905932737*alpha[4]*fUpwind_l[15]+0.31622776601683794*alpha[5]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[9]+0.31622776601683794*alpha[3]*fUpwind_l[6]; 
  Ghat_l[17] = 0.28284271247461906*alpha[5]*fUpwind_l[19]+0.28284271247461906*alpha[4]*fUpwind_l[18]+0.22587697572631277*alpha[7]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[17]+0.3162277660168379*alpha[13]*fUpwind_l[16]+0.3162277660168379*alpha[11]*fUpwind_l[14]+0.22587697572631277*alpha[11]*fUpwind_l[13]+0.3535533905932737*alpha[2]*fUpwind_l[13]+0.22587697572631277*fUpwind_l[11]*alpha[13]+0.3535533905932737*fUpwind_l[2]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[3]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[4]*alpha[5]; 
  Ghat_l[18] = 0.3162277660168379*alpha[7]*fUpwind_l[18]+0.3535533905932737*alpha[0]*fUpwind_l[18]+0.28284271247461906*alpha[4]*fUpwind_l[17]+0.3535533905932737*alpha[1]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_l[12]+0.282842712474619*fUpwind_l[10]*alpha[11]+0.3162277660168379*alpha[2]*fUpwind_l[10]+0.3535533905932737*alpha[5]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[6]; 
  Ghat_l[19] = 0.3162277660168379*alpha[7]*fUpwind_l[19]+0.3535533905932737*alpha[0]*fUpwind_l[19]+0.28284271247461906*alpha[5]*fUpwind_l[17]+0.3535533905932737*alpha[1]*fUpwind_l[16]+0.3162277660168379*alpha[11]*fUpwind_l[15]+0.3535533905932737*alpha[2]*fUpwind_l[15]+0.282842712474619*fUpwind_l[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_l[10]+0.3535533905932737*alpha[4]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[6]; 

  Ghat_r[0] = 0.3535533905932737*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[5]*fUpwind_r[5]+0.3535533905932737*alpha[4]*fUpwind_r[4]+0.3535533905932737*alpha[3]*fUpwind_r[3]+0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[1]*fUpwind_r[1]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.31622776601683794*alpha[5]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[5]*alpha[13]+0.31622776601683794*alpha[4]*fUpwind_r[11]+0.31622776601683794*fUpwind_r[4]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.3535533905932737*alpha[13]*fUpwind_r[17]+0.31622776601683794*alpha[4]*fUpwind_r[12]+0.3535533905932737*alpha[7]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[7]*alpha[11]+0.3535533905932737*alpha[5]*fUpwind_r[10]+0.3162277660168379*alpha[2]*fUpwind_r[8]+0.3535533905932737*alpha[3]*fUpwind_r[6]+0.3535533905932737*alpha[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3535533905932737*alpha[11]*fUpwind_r[17]+0.31622776601683794*alpha[5]*fUpwind_r[15]+0.3535533905932737*alpha[7]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_r[10]+0.3162277660168379*alpha[3]*fUpwind_r[9]+0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*alpha[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.3162277660168379*alpha[5]*fUpwind_r[17]+0.31622776601683794*fUpwind_r[10]*alpha[13]+0.28284271247461906*alpha[11]*fUpwind_r[12]+0.31622776601683794*alpha[2]*fUpwind_r[12]+0.31622776601683794*alpha[1]*fUpwind_r[11]+0.31622776601683794*fUpwind_r[1]*alpha[11]+0.3535533905932737*alpha[3]*fUpwind_r[10]+0.3162277660168379*alpha[4]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_r[6]+0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.3162277660168379*alpha[4]*fUpwind_r[17]+0.28284271247461906*alpha[13]*fUpwind_r[15]+0.31622776601683794*alpha[3]*fUpwind_r[15]+0.31622776601683794*alpha[1]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[1]*alpha[13]+0.31622776601683794*fUpwind_r[10]*alpha[11]+0.3535533905932737*alpha[2]*fUpwind_r[10]+0.3162277660168379*alpha[5]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_r[6]+0.3535533905932737*alpha[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alpha[3]; 
  Ghat_r[6] = 0.3162277660168379*alpha[5]*fUpwind_r[19]+0.3162277660168379*alpha[4]*fUpwind_r[18]+0.3535533905932737*alpha[7]*fUpwind_r[17]+0.31622776601683794*alpha[3]*fUpwind_r[16]+0.31622776601683794*alpha[2]*fUpwind_r[14]+0.3535533905932737*alpha[11]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[11]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*alpha[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alpha[3]; 
  Ghat_r[7] = 0.22587697572631277*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[3]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alpha[13]+0.22587697572631277*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[2]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2]*alpha[11]+0.22587697572631277*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[1]*fUpwind_r[1]; 
  Ghat_r[8] = 0.3535533905932737*alpha[5]*fUpwind_r[18]+0.3535533905932737*alpha[3]*fUpwind_r[14]+0.3535533905932737*alpha[1]*fUpwind_r[12]+0.3162277660168379*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha[4]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[16]+0.3535533905932737*alpha[1]*fUpwind_r[15]+0.3162277660168379*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[0]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[3]*fUpwind_r[3]; 
  Ghat_r[10] = 0.282842712474619*alpha[13]*fUpwind_r[19]+0.3162277660168379*alpha[3]*fUpwind_r[19]+0.282842712474619*alpha[11]*fUpwind_r[18]+0.3162277660168379*alpha[2]*fUpwind_r[18]+0.3162277660168379*alpha[1]*fUpwind_r[17]+0.31622776601683794*alpha[5]*fUpwind_r[16]+0.31622776601683794*alpha[4]*fUpwind_r[14]+0.31622776601683794*alpha[4]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[4]*alpha[13]+0.31622776601683794*alpha[5]*fUpwind_r[11]+0.31622776601683794*fUpwind_r[5]*alpha[11]+0.3162277660168379*alpha[7]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[10]+0.3535533905932737*alpha[1]*fUpwind_r[6]+0.3535533905932737*alpha[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.22587697572631277*alpha[13]*fUpwind_r[17]+0.3535533905932737*alpha[3]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[6]*alpha[13]+0.28284271247461906*alpha[4]*fUpwind_r[12]+0.22587697572631277*alpha[7]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[8]*alpha[11]+0.22587697572631277*fUpwind_r[7]*alpha[11]+0.3535533905932737*fUpwind_r[0]*alpha[11]+0.31622776601683794*alpha[5]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_r[4]+0.31622776601683794*fUpwind_r[1]*alpha[4]; 
  Ghat_r[12] = 0.3162277660168379*alpha[13]*fUpwind_r[18]+0.3535533905932737*alpha[3]*fUpwind_r[18]+0.3535533905932737*alpha[5]*fUpwind_r[14]+0.3162277660168379*alpha[7]*fUpwind_r[12]+0.3535533905932737*alpha[0]*fUpwind_r[12]+0.28284271247461906*alpha[4]*fUpwind_r[11]+0.28284271247461906*fUpwind_r[4]*alpha[11]+0.3535533905932737*alpha[1]*fUpwind_r[8]+0.31622776601683794*alpha[2]*fUpwind_r[4]+0.31622776601683794*fUpwind_r[2]*alpha[4]; 
  Ghat_r[13] = 0.22587697572631277*alpha[11]*fUpwind_r[17]+0.3535533905932737*alpha[2]*fUpwind_r[17]+0.28284271247461906*alpha[5]*fUpwind_r[15]+0.22587697572631277*alpha[7]*fUpwind_r[13]+0.3535533905932737*alpha[0]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[9]*alpha[13]+0.22587697572631277*fUpwind_r[7]*alpha[13]+0.3535533905932737*fUpwind_r[0]*alpha[13]+0.3535533905932737*fUpwind_r[6]*alpha[11]+0.31622776601683794*alpha[4]*fUpwind_r[10]+0.3535533905932737*alpha[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_r[5]+0.31622776601683794*fUpwind_r[1]*alpha[5]; 
  Ghat_r[14] = 0.3535533905932737*alpha[1]*fUpwind_r[18]+0.3162277660168379*alpha[11]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[14]+0.3535533905932737*alpha[5]*fUpwind_r[12]+0.31622776601683794*alpha[4]*fUpwind_r[10]+0.3535533905932737*alpha[3]*fUpwind_r[8]+0.31622776601683794*alpha[2]*fUpwind_r[6]; 
  Ghat_r[15] = 0.3162277660168379*alpha[11]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[19]+0.3535533905932737*alpha[4]*fUpwind_r[16]+0.3162277660168379*alpha[7]*fUpwind_r[15]+0.3535533905932737*alpha[0]*fUpwind_r[15]+0.28284271247461906*alpha[5]*fUpwind_r[13]+0.28284271247461906*fUpwind_r[5]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_r[9]+0.31622776601683794*alpha[3]*fUpwind_r[5]+0.31622776601683794*fUpwind_r[3]*alpha[5]; 
  Ghat_r[16] = 0.3535533905932737*alpha[1]*fUpwind_r[19]+0.3162277660168379*alpha[13]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[16]+0.3535533905932737*alpha[4]*fUpwind_r[15]+0.31622776601683794*alpha[5]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[9]+0.31622776601683794*alpha[3]*fUpwind_r[6]; 
  Ghat_r[17] = 0.28284271247461906*alpha[5]*fUpwind_r[19]+0.28284271247461906*alpha[4]*fUpwind_r[18]+0.22587697572631277*alpha[7]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[17]+0.3162277660168379*alpha[13]*fUpwind_r[16]+0.3162277660168379*alpha[11]*fUpwind_r[14]+0.22587697572631277*alpha[11]*fUpwind_r[13]+0.3535533905932737*alpha[2]*fUpwind_r[13]+0.22587697572631277*fUpwind_r[11]*alpha[13]+0.3535533905932737*fUpwind_r[2]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[3]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[4]*alpha[5]; 
  Ghat_r[18] = 0.3162277660168379*alpha[7]*fUpwind_r[18]+0.3535533905932737*alpha[0]*fUpwind_r[18]+0.28284271247461906*alpha[4]*fUpwind_r[17]+0.3535533905932737*alpha[1]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_r[12]+0.282842712474619*fUpwind_r[10]*alpha[11]+0.3162277660168379*alpha[2]*fUpwind_r[10]+0.3535533905932737*alpha[5]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[6]; 
  Ghat_r[19] = 0.3162277660168379*alpha[7]*fUpwind_r[19]+0.3535533905932737*alpha[0]*fUpwind_r[19]+0.28284271247461906*alpha[5]*fUpwind_r[17]+0.3535533905932737*alpha[1]*fUpwind_r[16]+0.3162277660168379*alpha[11]*fUpwind_r[15]+0.3535533905932737*alpha[2]*fUpwind_r[15]+0.282842712474619*fUpwind_r[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_r[10]+0.3535533905932737*alpha[4]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[6]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10); 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10); 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[7] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10); 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[9] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10); 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[12] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[15] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10); 
  out[16] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10); 
  out[17] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[18] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10); 
  out[19] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10); 
  out[20] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv10; 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[22] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv10; 
  out[23] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[24] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10); 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[26] += (1.5811388300841898*Ghat_l[3]-1.5811388300841898*Ghat_r[3])*dv10; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[29] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10); 
  out[30] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv10; 
  out[31] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10); 
  out[32] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10); 
  out[33] += (1.5811388300841895*Ghat_l[4]-1.5811388300841895*Ghat_r[4])*dv10; 
  out[34] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10); 
  out[35] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10); 
  out[36] += (1.5811388300841895*Ghat_l[5]-1.5811388300841895*Ghat_r[5])*dv10; 
  out[37] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv10; 
  out[38] += (1.5811388300841895*Ghat_l[6]-1.5811388300841895*Ghat_r[6])*dv10; 
  out[39] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv10; 
  out[40] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10); 
  out[41] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10); 
  out[42] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv10; 
  out[43] += -(1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv10); 
  out[44] += -(1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv10); 
  out[45] += (1.5811388300841898*Ghat_l[10]-1.5811388300841898*Ghat_r[10])*dv10; 
  out[46] += -(1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv10); 
  out[47] += -(1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv10); 

  return 0.;

} 
