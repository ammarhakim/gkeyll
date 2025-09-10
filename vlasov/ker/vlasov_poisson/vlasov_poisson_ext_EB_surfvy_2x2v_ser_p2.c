#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_EB_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 

  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 

  double alpha[20] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ex = &EBext[0]; 
  const double *Ey = &EBext[8]; 
  const double *Bx = &EBext[24]; 
  const double *By = &EBext[32]; 
  const double *Bz = &EBext[40]; 

  alpha[0] = -(1.4142135623730951*Bz[0]*wv1)-2.4494897427831783*phi[2]*dx11+1.4142135623730951*Ey[0]; 
  alpha[1] = -(1.4142135623730951*Bz[1]*wv1)-2.4494897427831783*phi[3]*dx11+1.4142135623730951*Ey[1]; 
  alpha[2] = -(1.4142135623730951*Bz[2]*wv1)-5.477225575051662*phi[5]*dx11+1.4142135623730951*Ey[2]; 
  alpha[3] = -(0.408248290463863*Bz[0]*dv1); 
  alpha[4] = -(1.4142135623730951*Bz[3]*wv1)-5.477225575051662*phi[7]*dx11+1.4142135623730951*Ey[3]; 
  alpha[5] = -(0.408248290463863*Bz[1]*dv1); 
  alpha[6] = -(0.408248290463863*Bz[2]*dv1); 
  alpha[7] = -(1.4142135623730951*Bz[4]*wv1)-2.4494897427831783*phi[6]*dx11+1.4142135623730951*Ey[4]; 
  alpha[8] = 1.4142135623730951*Ey[5]-1.4142135623730951*Bz[5]*wv1; 
  alpha[10] = -(0.408248290463863*Bz[3]*dv1); 
  alpha[11] = 1.4142135623730951*Ey[6]-1.4142135623730951*Bz[6]*wv1; 
  alpha[12] = 1.4142135623730951*Ey[7]-1.4142135623730951*Bz[7]*wv1; 
  alpha[13] = -(0.40824829046386296*Bz[4]*dv1); 
  alpha[14] = -(0.40824829046386296*Bz[5]*dv1); 
  alpha[17] = -(0.40824829046386296*Bz[6]*dv1); 
  alpha[18] = -(0.40824829046386296*Bz[7]*dv1); 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (0.5692099788303082*(alpha[18]+alpha[17])-0.42426406871192845*(alpha[14]+alpha[13])-0.42426406871192807*alpha[12]-0.42426406871192845*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fr); 
  } 
  if (-(0.42426406871192807*alpha[12])-0.42426406871192845*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fr); 
  } 
  if (-(0.5692099788303082*(alpha[18]+alpha[17]))+0.42426406871192845*(alpha[14]+alpha[13])-0.42426406871192807*alpha[12]-0.42426406871192845*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fr); 
  } 
  if (-(0.711512473537885*alpha[18])+0.5303300858899104*alpha[14]-0.42426406871192845*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fr); 
  } 
  if (0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fr); 
  } 
  if (0.711512473537885*alpha[18]-0.5303300858899104*alpha[14]+0.42426406871192845*alpha[13]+0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fr); 
  } 
  if (0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]-0.42426406871192845*(alpha[14]+alpha[13])-0.42426406871192807*alpha[12]+0.42426406871192845*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fr); 
  } 
  if (-(0.42426406871192807*alpha[12])+0.42426406871192845*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fr); 
  } 
  if (-(0.5692099788303082*alpha[18])+0.5692099788303082*alpha[17]+0.42426406871192845*(alpha[14]+alpha[13])-0.42426406871192807*alpha[12]+0.42426406871192845*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fr); 
  } 
  if (-(0.711512473537885*alpha[17])-0.42426406871192845*alpha[14]+0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fr); 
  } 
  if (0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fr); 
  } 
  if (0.711512473537885*alpha[17]+0.42426406871192845*alpha[14]-0.5303300858899104*alpha[13]+0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fr); 
  } 
  if (0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*(alpha[8]+alpha[7]) > 0) { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fr); 
  } 
  if (-(0.5303300858899104*(alpha[14]+alpha[13]))-0.3952847075210473*(alpha[8]+alpha[7])+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fr); 
  } 
  if (0.711512473537885*alpha[17]-0.42426406871192845*alpha[14]+0.5303300858899104*alpha[13]-0.5303300858899104*alpha[11]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fr); 
  } 
  if (-(0.5303300858899104*alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fr); 
  } 
  if (-(0.711512473537885*alpha[17])+0.42426406871192845*alpha[14]-0.5303300858899104*(alpha[13]+alpha[11])+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fr); 
  } 
  if (-(0.5692099788303082*alpha[18])+0.5692099788303082*alpha[17]-0.42426406871192845*(alpha[14]+alpha[13])+0.42426406871192807*alpha[12]-0.42426406871192845*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fr); 
  } 
  if (0.42426406871192807*alpha[12]-0.42426406871192845*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fr); 
  } 
  if (0.5692099788303082*alpha[18]-0.5692099788303082*alpha[17]+0.42426406871192845*(alpha[14]+alpha[13])+0.42426406871192807*alpha[12]-0.42426406871192845*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fr); 
  } 
  if (0.711512473537885*alpha[18]+0.5303300858899104*alpha[14]-0.42426406871192845*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fr); 
  } 
  if (-(0.5303300858899104*alpha[12])-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fr); 
  } 
  if (-(0.711512473537885*alpha[18])-0.5303300858899104*alpha[14]+0.42426406871192845*alpha[13]-0.5303300858899104*alpha[12]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fr); 
  } 
  if (-(0.5692099788303082*(alpha[18]+alpha[17]))-0.42426406871192845*(alpha[14]+alpha[13])+0.42426406871192807*alpha[12]+0.42426406871192845*alpha[11]-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fr); 
  } 
  if (0.42426406871192807*alpha[12]+0.42426406871192845*alpha[11]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fr); 
  } 
  if (0.5692099788303082*(alpha[18]+alpha[17])+0.42426406871192845*(alpha[14]+alpha[13])+0.42426406871192807*alpha[12]+0.42426406871192845*alpha[11]+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[18]*fUpwind_l[18]+0.3535533905932737*alpha[17]*fUpwind_l[17]+0.3535533905932737*alpha[14]*fUpwind_l[14]+0.3535533905932737*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[12]*fUpwind_l[12]+0.3535533905932737*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[10]*fUpwind_l[10]+0.3535533905932737*alpha[8]*fUpwind_l[8]+0.3535533905932737*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[6]*fUpwind_l[6]+0.3535533905932737*alpha[5]*fUpwind_l[5]+0.3535533905932737*alpha[4]*fUpwind_l[4]+0.3535533905932737*alpha[3]*fUpwind_l[3]+0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[1]*fUpwind_l[1]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alpha[14]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[14]*alpha[18]+0.3162277660168379*alpha[10]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[10]*alpha[17]+0.31622776601683794*alpha[5]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[5]*alpha[13]+0.3535533905932737*alpha[8]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[8]*alpha[12]+0.31622776601683794*alpha[4]*fUpwind_l[11]+0.31622776601683794*fUpwind_l[4]*alpha[11]+0.3535533905932737*alpha[6]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alpha[10]+0.3162277660168379*alpha[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.3162277660168379*alpha[10]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[10]*alpha[18]+0.3535533905932737*alpha[13]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[13]*alpha[17]+0.31622776601683794*alpha[6]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[6]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_l[12]+0.31622776601683794*fUpwind_l[4]*alpha[12]+0.3535533905932737*alpha[7]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[7]*alpha[11]+0.3535533905932737*alpha[5]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[5]*alpha[10]+0.3162277660168379*alpha[2]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[2]*alpha[8]+0.3535533905932737*alpha[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3162277660168379*alpha[10]*fUpwind_l[19]+0.3535533905932737*alpha[12]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[12]*alpha[18]+0.3535533905932737*alpha[11]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[11]*alpha[17]+0.31622776601683794*alpha[6]*fUpwind_l[16]+0.31622776601683794*alpha[5]*fUpwind_l[15]+0.3535533905932737*alpha[8]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[8]*alpha[14]+0.3535533905932737*alpha[7]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[4]*alpha[10]+0.3162277660168379*alpha[3]*fUpwind_l[9]+0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.28284271247461906*alpha[17]*fUpwind_l[18]+0.3162277660168379*alpha[6]*fUpwind_l[18]+0.28284271247461906*fUpwind_l[17]*alpha[18]+0.3162277660168379*fUpwind_l[6]*alpha[18]+0.3162277660168379*alpha[5]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[5]*alpha[17]+0.31622776601683794*alpha[10]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[10]*alpha[14]+0.31622776601683794*alpha[10]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[10]*alpha[13]+0.28284271247461906*alpha[11]*fUpwind_l[12]+0.31622776601683794*alpha[2]*fUpwind_l[12]+0.28284271247461906*fUpwind_l[11]*alpha[12]+0.31622776601683794*fUpwind_l[2]*alpha[12]+0.31622776601683794*alpha[1]*fUpwind_l[11]+0.31622776601683794*fUpwind_l[1]*alpha[11]+0.3535533905932737*alpha[3]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alpha[10]+0.3162277660168379*alpha[4]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[4]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.28284271247461906*alpha[17]*fUpwind_l[19]+0.3162277660168379*alpha[6]*fUpwind_l[19]+0.3535533905932737*alpha[8]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[8]*alpha[18]+0.3162277660168379*alpha[4]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[4]*alpha[17]+0.31622776601683794*alpha[10]*fUpwind_l[16]+0.28284271247461906*alpha[13]*fUpwind_l[15]+0.31622776601683794*alpha[3]*fUpwind_l[15]+0.3535533905932737*alpha[12]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[12]*alpha[14]+0.31622776601683794*alpha[1]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[1]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_l[11]+0.31622776601683794*fUpwind_l[10]*alpha[11]+0.3535533905932737*alpha[2]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[2]*alpha[10]+0.3162277660168379*alpha[5]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alpha[3]; 
  Ghat_l[6] = 0.28284271247461906*alpha[18]*fUpwind_l[19]+0.3162277660168379*alpha[5]*fUpwind_l[19]+0.3162277660168379*alpha[4]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[4]*alpha[18]+0.3535533905932737*alpha[7]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[7]*alpha[17]+0.28284271247461906*alpha[14]*fUpwind_l[16]+0.31622776601683794*alpha[3]*fUpwind_l[16]+0.31622776601683794*alpha[10]*fUpwind_l[15]+0.31622776601683794*alpha[2]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[2]*alpha[14]+0.3535533905932737*alpha[11]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[11]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_l[12]+0.31622776601683794*fUpwind_l[10]*alpha[12]+0.3535533905932737*alpha[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[1]*alpha[10]+0.3162277660168379*alpha[6]*fUpwind_l[9]+0.3162277660168379*alpha[6]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[6]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alpha[3]; 
  Ghat_l[7] = 0.3162277660168379*alpha[18]*fUpwind_l[18]+0.22587697572631277*alpha[17]*fUpwind_l[17]+0.3535533905932737*alpha[6]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[6]*alpha[17]+0.22587697572631277*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[3]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alpha[13]+0.3162277660168379*alpha[12]*fUpwind_l[12]+0.22587697572631277*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[2]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2]*alpha[11]+0.3162277660168379*alpha[10]*fUpwind_l[10]+0.22587697572631277*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[1]*fUpwind_l[1]; 
  Ghat_l[8] = 0.22587697572631277*alpha[18]*fUpwind_l[18]+0.3535533905932737*alpha[5]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[5]*alpha[18]+0.3162277660168379*alpha[17]*fUpwind_l[17]+0.22587697572631277*alpha[14]*fUpwind_l[14]+0.3535533905932737*alpha[3]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[3]*alpha[14]+0.22587697572631277*alpha[12]*fUpwind_l[12]+0.3535533905932737*alpha[1]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[1]*alpha[12]+0.3162277660168379*alpha[11]*fUpwind_l[11]+0.3162277660168379*alpha[10]*fUpwind_l[10]+0.22587697572631277*alpha[8]*fUpwind_l[8]+0.3535533905932737*alpha[0]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[0]*alpha[8]+0.3162277660168379*alpha[6]*fUpwind_l[6]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha[4]*fUpwind_l[19]+0.3162277660168379*alpha[18]*fUpwind_l[18]+0.3162277660168379*alpha[17]*fUpwind_l[17]+0.3535533905932737*alpha[2]*fUpwind_l[16]+0.3535533905932737*alpha[1]*fUpwind_l[15]+0.3162277660168379*alpha[14]*fUpwind_l[14]+0.3162277660168379*alpha[13]*fUpwind_l[13]+0.3162277660168379*alpha[10]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[9]+0.3162277660168379*alpha[6]*fUpwind_l[6]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[3]*fUpwind_l[3]; 
  Ghat_l[10] = 0.282842712474619*alpha[14]*fUpwind_l[19]+0.282842712474619*alpha[13]*fUpwind_l[19]+0.3162277660168379*alpha[3]*fUpwind_l[19]+0.282842712474619*alpha[11]*fUpwind_l[18]+0.3162277660168379*alpha[2]*fUpwind_l[18]+0.282842712474619*fUpwind_l[16]*alpha[18]+0.282842712474619*fUpwind_l[11]*alpha[18]+0.3162277660168379*fUpwind_l[2]*alpha[18]+0.282842712474619*alpha[12]*fUpwind_l[17]+0.3162277660168379*alpha[1]*fUpwind_l[17]+0.282842712474619*fUpwind_l[15]*alpha[17]+0.282842712474619*fUpwind_l[12]*alpha[17]+0.3162277660168379*fUpwind_l[1]*alpha[17]+0.31622776601683794*alpha[5]*fUpwind_l[16]+0.31622776601683794*alpha[6]*fUpwind_l[15]+0.31622776601683794*alpha[4]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[4]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[4]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_l[12]+0.31622776601683794*fUpwind_l[6]*alpha[12]+0.31622776601683794*alpha[5]*fUpwind_l[11]+0.31622776601683794*fUpwind_l[5]*alpha[11]+0.3162277660168379*alpha[8]*fUpwind_l[10]+0.3162277660168379*alpha[7]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[9]*alpha[10]+0.3162277660168379*fUpwind_l[8]*alpha[10]+0.3162277660168379*fUpwind_l[7]*alpha[10]+0.3535533905932737*fUpwind_l[0]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.282842712474619*alpha[10]*fUpwind_l[18]+0.282842712474619*fUpwind_l[10]*alpha[18]+0.3162277660168379*alpha[14]*fUpwind_l[17]+0.22587697572631277*alpha[13]*fUpwind_l[17]+0.3535533905932737*alpha[3]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[14]*alpha[17]+0.22587697572631277*fUpwind_l[13]*alpha[17]+0.3535533905932737*fUpwind_l[3]*alpha[17]+0.3535533905932737*alpha[6]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[6]*alpha[13]+0.28284271247461906*alpha[4]*fUpwind_l[12]+0.28284271247461906*fUpwind_l[4]*alpha[12]+0.3162277660168379*alpha[8]*fUpwind_l[11]+0.22587697572631277*alpha[7]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[8]*alpha[11]+0.22587697572631277*fUpwind_l[7]*alpha[11]+0.3535533905932737*fUpwind_l[0]*alpha[11]+0.31622776601683794*alpha[5]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_l[4]+0.31622776601683794*fUpwind_l[1]*alpha[4]; 
  Ghat_l[12] = 0.22587697572631277*alpha[14]*fUpwind_l[18]+0.3162277660168379*alpha[13]*fUpwind_l[18]+0.3535533905932737*alpha[3]*fUpwind_l[18]+0.22587697572631277*fUpwind_l[14]*alpha[18]+0.3162277660168379*fUpwind_l[13]*alpha[18]+0.3535533905932737*fUpwind_l[3]*alpha[18]+0.282842712474619*alpha[10]*fUpwind_l[17]+0.282842712474619*fUpwind_l[10]*alpha[17]+0.3535533905932737*alpha[5]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[5]*alpha[14]+0.22587697572631277*alpha[8]*fUpwind_l[12]+0.3162277660168379*alpha[7]*fUpwind_l[12]+0.3535533905932737*alpha[0]*fUpwind_l[12]+0.22587697572631277*fUpwind_l[8]*alpha[12]+0.3162277660168379*fUpwind_l[7]*alpha[12]+0.3535533905932737*fUpwind_l[0]*alpha[12]+0.28284271247461906*alpha[4]*fUpwind_l[11]+0.28284271247461906*fUpwind_l[4]*alpha[11]+0.31622776601683794*alpha[6]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[1]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_l[4]+0.31622776601683794*fUpwind_l[2]*alpha[4]; 
  Ghat_l[13] = 0.282842712474619*alpha[10]*fUpwind_l[19]+0.3162277660168379*alpha[12]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[12]*alpha[18]+0.22587697572631277*alpha[11]*fUpwind_l[17]+0.3535533905932737*alpha[2]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[16]*alpha[17]+0.22587697572631277*fUpwind_l[11]*alpha[17]+0.3535533905932737*fUpwind_l[2]*alpha[17]+0.28284271247461906*alpha[5]*fUpwind_l[15]+0.22587697572631277*alpha[7]*fUpwind_l[13]+0.3535533905932737*alpha[0]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[9]*alpha[13]+0.22587697572631277*fUpwind_l[7]*alpha[13]+0.3535533905932737*fUpwind_l[0]*alpha[13]+0.3535533905932737*alpha[6]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[6]*alpha[11]+0.31622776601683794*alpha[4]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_l[5]+0.31622776601683794*fUpwind_l[1]*alpha[5]; 
  Ghat_l[14] = 0.282842712474619*alpha[10]*fUpwind_l[19]+0.22587697572631277*alpha[12]*fUpwind_l[18]+0.3535533905932737*alpha[1]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[15]*alpha[18]+0.22587697572631277*fUpwind_l[12]*alpha[18]+0.3535533905932737*fUpwind_l[1]*alpha[18]+0.3162277660168379*alpha[11]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[11]*alpha[17]+0.28284271247461906*alpha[6]*fUpwind_l[16]+0.22587697572631277*alpha[8]*fUpwind_l[14]+0.3535533905932737*alpha[0]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[9]*alpha[14]+0.22587697572631277*fUpwind_l[8]*alpha[14]+0.3535533905932737*fUpwind_l[0]*alpha[14]+0.3535533905932737*alpha[5]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[5]*alpha[12]+0.31622776601683794*alpha[4]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[3]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_l[6]+0.31622776601683794*fUpwind_l[2]*alpha[6]; 
  Ghat_l[15] = 0.3162277660168379*alpha[11]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[19]+0.3162277660168379*alpha[14]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[14]*alpha[18]+0.282842712474619*alpha[10]*fUpwind_l[17]+0.282842712474619*fUpwind_l[10]*alpha[17]+0.3535533905932737*alpha[4]*fUpwind_l[16]+0.3162277660168379*alpha[7]*fUpwind_l[15]+0.3535533905932737*alpha[0]*fUpwind_l[15]+0.28284271247461906*alpha[5]*fUpwind_l[13]+0.28284271247461906*fUpwind_l[5]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_l[9]+0.31622776601683794*alpha[3]*fUpwind_l[5]+0.31622776601683794*fUpwind_l[3]*alpha[5]; 
  Ghat_l[16] = 0.3162277660168379*alpha[12]*fUpwind_l[19]+0.3535533905932737*alpha[1]*fUpwind_l[19]+0.282842712474619*alpha[10]*fUpwind_l[18]+0.282842712474619*fUpwind_l[10]*alpha[18]+0.3162277660168379*alpha[13]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[13]*alpha[17]+0.3162277660168379*alpha[8]*fUpwind_l[16]+0.3535533905932737*alpha[0]*fUpwind_l[16]+0.3535533905932737*alpha[4]*fUpwind_l[15]+0.28284271247461906*alpha[6]*fUpwind_l[14]+0.28284271247461906*fUpwind_l[6]*alpha[14]+0.31622776601683794*alpha[5]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_l[9]+0.31622776601683794*alpha[3]*fUpwind_l[6]+0.31622776601683794*fUpwind_l[3]*alpha[6]; 
  Ghat_l[17] = 0.2529822128134704*alpha[18]*fUpwind_l[19]+0.28284271247461906*alpha[5]*fUpwind_l[19]+0.28284271247461906*alpha[4]*fUpwind_l[18]+0.28284271247461906*fUpwind_l[4]*alpha[18]+0.3162277660168379*alpha[8]*fUpwind_l[17]+0.22587697572631277*alpha[7]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[9]*alpha[17]+0.3162277660168379*fUpwind_l[8]*alpha[17]+0.22587697572631277*fUpwind_l[7]*alpha[17]+0.3535533905932737*fUpwind_l[0]*alpha[17]+0.3162277660168379*alpha[13]*fUpwind_l[16]+0.282842712474619*alpha[10]*fUpwind_l[15]+0.3162277660168379*alpha[11]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[11]*alpha[14]+0.22587697572631277*alpha[11]*fUpwind_l[13]+0.3535533905932737*alpha[2]*fUpwind_l[13]+0.22587697572631277*fUpwind_l[11]*alpha[13]+0.3535533905932737*fUpwind_l[2]*alpha[13]+0.282842712474619*alpha[10]*fUpwind_l[12]+0.282842712474619*fUpwind_l[10]*alpha[12]+0.3535533905932737*alpha[3]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[3]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[1]*alpha[10]+0.3535533905932737*alpha[6]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[4]*alpha[5]; 
  Ghat_l[18] = 0.2529822128134704*alpha[17]*fUpwind_l[19]+0.28284271247461906*alpha[6]*fUpwind_l[19]+0.22587697572631277*alpha[8]*fUpwind_l[18]+0.3162277660168379*alpha[7]*fUpwind_l[18]+0.3535533905932737*alpha[0]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[9]*alpha[18]+0.22587697572631277*fUpwind_l[8]*alpha[18]+0.3162277660168379*fUpwind_l[7]*alpha[18]+0.3535533905932737*fUpwind_l[0]*alpha[18]+0.28284271247461906*alpha[4]*fUpwind_l[17]+0.28284271247461906*fUpwind_l[4]*alpha[17]+0.282842712474619*alpha[10]*fUpwind_l[16]+0.3162277660168379*alpha[14]*fUpwind_l[15]+0.22587697572631277*alpha[12]*fUpwind_l[14]+0.3535533905932737*alpha[1]*fUpwind_l[14]+0.22587697572631277*fUpwind_l[12]*alpha[14]+0.3535533905932737*fUpwind_l[1]*alpha[14]+0.3162277660168379*alpha[12]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[3]*alpha[12]+0.282842712474619*alpha[10]*fUpwind_l[11]+0.282842712474619*fUpwind_l[10]*alpha[11]+0.3162277660168379*alpha[2]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[2]*alpha[10]+0.3535533905932737*alpha[5]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[4]*alpha[6]; 
  Ghat_l[19] = 0.3162277660168379*alpha[8]*fUpwind_l[19]+0.3162277660168379*alpha[7]*fUpwind_l[19]+0.3535533905932737*alpha[0]*fUpwind_l[19]+0.2529822128134704*alpha[17]*fUpwind_l[18]+0.28284271247461906*alpha[6]*fUpwind_l[18]+0.2529822128134704*fUpwind_l[17]*alpha[18]+0.28284271247461906*fUpwind_l[6]*alpha[18]+0.28284271247461906*alpha[5]*fUpwind_l[17]+0.28284271247461906*fUpwind_l[5]*alpha[17]+0.3162277660168379*alpha[12]*fUpwind_l[16]+0.3535533905932737*alpha[1]*fUpwind_l[16]+0.3162277660168379*alpha[11]*fUpwind_l[15]+0.3535533905932737*alpha[2]*fUpwind_l[15]+0.282842712474619*alpha[10]*fUpwind_l[14]+0.282842712474619*fUpwind_l[10]*alpha[14]+0.282842712474619*alpha[10]*fUpwind_l[13]+0.282842712474619*fUpwind_l[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[3]*alpha[10]+0.3535533905932737*alpha[4]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[5]*alpha[6]; 

  Ghat_r[0] = 0.3535533905932737*alpha[18]*fUpwind_r[18]+0.3535533905932737*alpha[17]*fUpwind_r[17]+0.3535533905932737*alpha[14]*fUpwind_r[14]+0.3535533905932737*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[12]*fUpwind_r[12]+0.3535533905932737*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[10]*fUpwind_r[10]+0.3535533905932737*alpha[8]*fUpwind_r[8]+0.3535533905932737*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[6]*fUpwind_r[6]+0.3535533905932737*alpha[5]*fUpwind_r[5]+0.3535533905932737*alpha[4]*fUpwind_r[4]+0.3535533905932737*alpha[3]*fUpwind_r[3]+0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[1]*fUpwind_r[1]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alpha[14]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[14]*alpha[18]+0.3162277660168379*alpha[10]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[10]*alpha[17]+0.31622776601683794*alpha[5]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[5]*alpha[13]+0.3535533905932737*alpha[8]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[8]*alpha[12]+0.31622776601683794*alpha[4]*fUpwind_r[11]+0.31622776601683794*fUpwind_r[4]*alpha[11]+0.3535533905932737*alpha[6]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alpha[10]+0.3162277660168379*alpha[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.3162277660168379*alpha[10]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[10]*alpha[18]+0.3535533905932737*alpha[13]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[13]*alpha[17]+0.31622776601683794*alpha[6]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[6]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_r[12]+0.31622776601683794*fUpwind_r[4]*alpha[12]+0.3535533905932737*alpha[7]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[7]*alpha[11]+0.3535533905932737*alpha[5]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[5]*alpha[10]+0.3162277660168379*alpha[2]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[2]*alpha[8]+0.3535533905932737*alpha[3]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3162277660168379*alpha[10]*fUpwind_r[19]+0.3535533905932737*alpha[12]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[12]*alpha[18]+0.3535533905932737*alpha[11]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[11]*alpha[17]+0.31622776601683794*alpha[6]*fUpwind_r[16]+0.31622776601683794*alpha[5]*fUpwind_r[15]+0.3535533905932737*alpha[8]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[8]*alpha[14]+0.3535533905932737*alpha[7]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[4]*alpha[10]+0.3162277660168379*alpha[3]*fUpwind_r[9]+0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.28284271247461906*alpha[17]*fUpwind_r[18]+0.3162277660168379*alpha[6]*fUpwind_r[18]+0.28284271247461906*fUpwind_r[17]*alpha[18]+0.3162277660168379*fUpwind_r[6]*alpha[18]+0.3162277660168379*alpha[5]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[5]*alpha[17]+0.31622776601683794*alpha[10]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[10]*alpha[14]+0.31622776601683794*alpha[10]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[10]*alpha[13]+0.28284271247461906*alpha[11]*fUpwind_r[12]+0.31622776601683794*alpha[2]*fUpwind_r[12]+0.28284271247461906*fUpwind_r[11]*alpha[12]+0.31622776601683794*fUpwind_r[2]*alpha[12]+0.31622776601683794*alpha[1]*fUpwind_r[11]+0.31622776601683794*fUpwind_r[1]*alpha[11]+0.3535533905932737*alpha[3]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alpha[10]+0.3162277660168379*alpha[4]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[4]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.28284271247461906*alpha[17]*fUpwind_r[19]+0.3162277660168379*alpha[6]*fUpwind_r[19]+0.3535533905932737*alpha[8]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[8]*alpha[18]+0.3162277660168379*alpha[4]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[4]*alpha[17]+0.31622776601683794*alpha[10]*fUpwind_r[16]+0.28284271247461906*alpha[13]*fUpwind_r[15]+0.31622776601683794*alpha[3]*fUpwind_r[15]+0.3535533905932737*alpha[12]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[12]*alpha[14]+0.31622776601683794*alpha[1]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[1]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_r[11]+0.31622776601683794*fUpwind_r[10]*alpha[11]+0.3535533905932737*alpha[2]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[2]*alpha[10]+0.3162277660168379*alpha[5]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alpha[3]; 
  Ghat_r[6] = 0.28284271247461906*alpha[18]*fUpwind_r[19]+0.3162277660168379*alpha[5]*fUpwind_r[19]+0.3162277660168379*alpha[4]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[4]*alpha[18]+0.3535533905932737*alpha[7]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[7]*alpha[17]+0.28284271247461906*alpha[14]*fUpwind_r[16]+0.31622776601683794*alpha[3]*fUpwind_r[16]+0.31622776601683794*alpha[10]*fUpwind_r[15]+0.31622776601683794*alpha[2]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[2]*alpha[14]+0.3535533905932737*alpha[11]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[11]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_r[12]+0.31622776601683794*fUpwind_r[10]*alpha[12]+0.3535533905932737*alpha[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[1]*alpha[10]+0.3162277660168379*alpha[6]*fUpwind_r[9]+0.3162277660168379*alpha[6]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[6]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alpha[3]; 
  Ghat_r[7] = 0.3162277660168379*alpha[18]*fUpwind_r[18]+0.22587697572631277*alpha[17]*fUpwind_r[17]+0.3535533905932737*alpha[6]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[6]*alpha[17]+0.22587697572631277*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[3]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alpha[13]+0.3162277660168379*alpha[12]*fUpwind_r[12]+0.22587697572631277*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[2]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2]*alpha[11]+0.3162277660168379*alpha[10]*fUpwind_r[10]+0.22587697572631277*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[1]*fUpwind_r[1]; 
  Ghat_r[8] = 0.22587697572631277*alpha[18]*fUpwind_r[18]+0.3535533905932737*alpha[5]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[5]*alpha[18]+0.3162277660168379*alpha[17]*fUpwind_r[17]+0.22587697572631277*alpha[14]*fUpwind_r[14]+0.3535533905932737*alpha[3]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[3]*alpha[14]+0.22587697572631277*alpha[12]*fUpwind_r[12]+0.3535533905932737*alpha[1]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[1]*alpha[12]+0.3162277660168379*alpha[11]*fUpwind_r[11]+0.3162277660168379*alpha[10]*fUpwind_r[10]+0.22587697572631277*alpha[8]*fUpwind_r[8]+0.3535533905932737*alpha[0]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[0]*alpha[8]+0.3162277660168379*alpha[6]*fUpwind_r[6]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha[4]*fUpwind_r[19]+0.3162277660168379*alpha[18]*fUpwind_r[18]+0.3162277660168379*alpha[17]*fUpwind_r[17]+0.3535533905932737*alpha[2]*fUpwind_r[16]+0.3535533905932737*alpha[1]*fUpwind_r[15]+0.3162277660168379*alpha[14]*fUpwind_r[14]+0.3162277660168379*alpha[13]*fUpwind_r[13]+0.3162277660168379*alpha[10]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[9]+0.3162277660168379*alpha[6]*fUpwind_r[6]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[3]*fUpwind_r[3]; 
  Ghat_r[10] = 0.282842712474619*alpha[14]*fUpwind_r[19]+0.282842712474619*alpha[13]*fUpwind_r[19]+0.3162277660168379*alpha[3]*fUpwind_r[19]+0.282842712474619*alpha[11]*fUpwind_r[18]+0.3162277660168379*alpha[2]*fUpwind_r[18]+0.282842712474619*fUpwind_r[16]*alpha[18]+0.282842712474619*fUpwind_r[11]*alpha[18]+0.3162277660168379*fUpwind_r[2]*alpha[18]+0.282842712474619*alpha[12]*fUpwind_r[17]+0.3162277660168379*alpha[1]*fUpwind_r[17]+0.282842712474619*fUpwind_r[15]*alpha[17]+0.282842712474619*fUpwind_r[12]*alpha[17]+0.3162277660168379*fUpwind_r[1]*alpha[17]+0.31622776601683794*alpha[5]*fUpwind_r[16]+0.31622776601683794*alpha[6]*fUpwind_r[15]+0.31622776601683794*alpha[4]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[4]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[4]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_r[12]+0.31622776601683794*fUpwind_r[6]*alpha[12]+0.31622776601683794*alpha[5]*fUpwind_r[11]+0.31622776601683794*fUpwind_r[5]*alpha[11]+0.3162277660168379*alpha[8]*fUpwind_r[10]+0.3162277660168379*alpha[7]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[9]*alpha[10]+0.3162277660168379*fUpwind_r[8]*alpha[10]+0.3162277660168379*fUpwind_r[7]*alpha[10]+0.3535533905932737*fUpwind_r[0]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.282842712474619*alpha[10]*fUpwind_r[18]+0.282842712474619*fUpwind_r[10]*alpha[18]+0.3162277660168379*alpha[14]*fUpwind_r[17]+0.22587697572631277*alpha[13]*fUpwind_r[17]+0.3535533905932737*alpha[3]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[14]*alpha[17]+0.22587697572631277*fUpwind_r[13]*alpha[17]+0.3535533905932737*fUpwind_r[3]*alpha[17]+0.3535533905932737*alpha[6]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[6]*alpha[13]+0.28284271247461906*alpha[4]*fUpwind_r[12]+0.28284271247461906*fUpwind_r[4]*alpha[12]+0.3162277660168379*alpha[8]*fUpwind_r[11]+0.22587697572631277*alpha[7]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[8]*alpha[11]+0.22587697572631277*fUpwind_r[7]*alpha[11]+0.3535533905932737*fUpwind_r[0]*alpha[11]+0.31622776601683794*alpha[5]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_r[4]+0.31622776601683794*fUpwind_r[1]*alpha[4]; 
  Ghat_r[12] = 0.22587697572631277*alpha[14]*fUpwind_r[18]+0.3162277660168379*alpha[13]*fUpwind_r[18]+0.3535533905932737*alpha[3]*fUpwind_r[18]+0.22587697572631277*fUpwind_r[14]*alpha[18]+0.3162277660168379*fUpwind_r[13]*alpha[18]+0.3535533905932737*fUpwind_r[3]*alpha[18]+0.282842712474619*alpha[10]*fUpwind_r[17]+0.282842712474619*fUpwind_r[10]*alpha[17]+0.3535533905932737*alpha[5]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[5]*alpha[14]+0.22587697572631277*alpha[8]*fUpwind_r[12]+0.3162277660168379*alpha[7]*fUpwind_r[12]+0.3535533905932737*alpha[0]*fUpwind_r[12]+0.22587697572631277*fUpwind_r[8]*alpha[12]+0.3162277660168379*fUpwind_r[7]*alpha[12]+0.3535533905932737*fUpwind_r[0]*alpha[12]+0.28284271247461906*alpha[4]*fUpwind_r[11]+0.28284271247461906*fUpwind_r[4]*alpha[11]+0.31622776601683794*alpha[6]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[1]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_r[4]+0.31622776601683794*fUpwind_r[2]*alpha[4]; 
  Ghat_r[13] = 0.282842712474619*alpha[10]*fUpwind_r[19]+0.3162277660168379*alpha[12]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[12]*alpha[18]+0.22587697572631277*alpha[11]*fUpwind_r[17]+0.3535533905932737*alpha[2]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[16]*alpha[17]+0.22587697572631277*fUpwind_r[11]*alpha[17]+0.3535533905932737*fUpwind_r[2]*alpha[17]+0.28284271247461906*alpha[5]*fUpwind_r[15]+0.22587697572631277*alpha[7]*fUpwind_r[13]+0.3535533905932737*alpha[0]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[9]*alpha[13]+0.22587697572631277*fUpwind_r[7]*alpha[13]+0.3535533905932737*fUpwind_r[0]*alpha[13]+0.3535533905932737*alpha[6]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[6]*alpha[11]+0.31622776601683794*alpha[4]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_r[5]+0.31622776601683794*fUpwind_r[1]*alpha[5]; 
  Ghat_r[14] = 0.282842712474619*alpha[10]*fUpwind_r[19]+0.22587697572631277*alpha[12]*fUpwind_r[18]+0.3535533905932737*alpha[1]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[15]*alpha[18]+0.22587697572631277*fUpwind_r[12]*alpha[18]+0.3535533905932737*fUpwind_r[1]*alpha[18]+0.3162277660168379*alpha[11]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[11]*alpha[17]+0.28284271247461906*alpha[6]*fUpwind_r[16]+0.22587697572631277*alpha[8]*fUpwind_r[14]+0.3535533905932737*alpha[0]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[9]*alpha[14]+0.22587697572631277*fUpwind_r[8]*alpha[14]+0.3535533905932737*fUpwind_r[0]*alpha[14]+0.3535533905932737*alpha[5]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[5]*alpha[12]+0.31622776601683794*alpha[4]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[3]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_r[6]+0.31622776601683794*fUpwind_r[2]*alpha[6]; 
  Ghat_r[15] = 0.3162277660168379*alpha[11]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[19]+0.3162277660168379*alpha[14]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[14]*alpha[18]+0.282842712474619*alpha[10]*fUpwind_r[17]+0.282842712474619*fUpwind_r[10]*alpha[17]+0.3535533905932737*alpha[4]*fUpwind_r[16]+0.3162277660168379*alpha[7]*fUpwind_r[15]+0.3535533905932737*alpha[0]*fUpwind_r[15]+0.28284271247461906*alpha[5]*fUpwind_r[13]+0.28284271247461906*fUpwind_r[5]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_r[9]+0.31622776601683794*alpha[3]*fUpwind_r[5]+0.31622776601683794*fUpwind_r[3]*alpha[5]; 
  Ghat_r[16] = 0.3162277660168379*alpha[12]*fUpwind_r[19]+0.3535533905932737*alpha[1]*fUpwind_r[19]+0.282842712474619*alpha[10]*fUpwind_r[18]+0.282842712474619*fUpwind_r[10]*alpha[18]+0.3162277660168379*alpha[13]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[13]*alpha[17]+0.3162277660168379*alpha[8]*fUpwind_r[16]+0.3535533905932737*alpha[0]*fUpwind_r[16]+0.3535533905932737*alpha[4]*fUpwind_r[15]+0.28284271247461906*alpha[6]*fUpwind_r[14]+0.28284271247461906*fUpwind_r[6]*alpha[14]+0.31622776601683794*alpha[5]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_r[9]+0.31622776601683794*alpha[3]*fUpwind_r[6]+0.31622776601683794*fUpwind_r[3]*alpha[6]; 
  Ghat_r[17] = 0.2529822128134704*alpha[18]*fUpwind_r[19]+0.28284271247461906*alpha[5]*fUpwind_r[19]+0.28284271247461906*alpha[4]*fUpwind_r[18]+0.28284271247461906*fUpwind_r[4]*alpha[18]+0.3162277660168379*alpha[8]*fUpwind_r[17]+0.22587697572631277*alpha[7]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[9]*alpha[17]+0.3162277660168379*fUpwind_r[8]*alpha[17]+0.22587697572631277*fUpwind_r[7]*alpha[17]+0.3535533905932737*fUpwind_r[0]*alpha[17]+0.3162277660168379*alpha[13]*fUpwind_r[16]+0.282842712474619*alpha[10]*fUpwind_r[15]+0.3162277660168379*alpha[11]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[11]*alpha[14]+0.22587697572631277*alpha[11]*fUpwind_r[13]+0.3535533905932737*alpha[2]*fUpwind_r[13]+0.22587697572631277*fUpwind_r[11]*alpha[13]+0.3535533905932737*fUpwind_r[2]*alpha[13]+0.282842712474619*alpha[10]*fUpwind_r[12]+0.282842712474619*fUpwind_r[10]*alpha[12]+0.3535533905932737*alpha[3]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[3]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[1]*alpha[10]+0.3535533905932737*alpha[6]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[4]*alpha[5]; 
  Ghat_r[18] = 0.2529822128134704*alpha[17]*fUpwind_r[19]+0.28284271247461906*alpha[6]*fUpwind_r[19]+0.22587697572631277*alpha[8]*fUpwind_r[18]+0.3162277660168379*alpha[7]*fUpwind_r[18]+0.3535533905932737*alpha[0]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[9]*alpha[18]+0.22587697572631277*fUpwind_r[8]*alpha[18]+0.3162277660168379*fUpwind_r[7]*alpha[18]+0.3535533905932737*fUpwind_r[0]*alpha[18]+0.28284271247461906*alpha[4]*fUpwind_r[17]+0.28284271247461906*fUpwind_r[4]*alpha[17]+0.282842712474619*alpha[10]*fUpwind_r[16]+0.3162277660168379*alpha[14]*fUpwind_r[15]+0.22587697572631277*alpha[12]*fUpwind_r[14]+0.3535533905932737*alpha[1]*fUpwind_r[14]+0.22587697572631277*fUpwind_r[12]*alpha[14]+0.3535533905932737*fUpwind_r[1]*alpha[14]+0.3162277660168379*alpha[12]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[3]*alpha[12]+0.282842712474619*alpha[10]*fUpwind_r[11]+0.282842712474619*fUpwind_r[10]*alpha[11]+0.3162277660168379*alpha[2]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[2]*alpha[10]+0.3535533905932737*alpha[5]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[4]*alpha[6]; 
  Ghat_r[19] = 0.3162277660168379*alpha[8]*fUpwind_r[19]+0.3162277660168379*alpha[7]*fUpwind_r[19]+0.3535533905932737*alpha[0]*fUpwind_r[19]+0.2529822128134704*alpha[17]*fUpwind_r[18]+0.28284271247461906*alpha[6]*fUpwind_r[18]+0.2529822128134704*fUpwind_r[17]*alpha[18]+0.28284271247461906*fUpwind_r[6]*alpha[18]+0.28284271247461906*alpha[5]*fUpwind_r[17]+0.28284271247461906*fUpwind_r[5]*alpha[17]+0.3162277660168379*alpha[12]*fUpwind_r[16]+0.3535533905932737*alpha[1]*fUpwind_r[16]+0.3162277660168379*alpha[11]*fUpwind_r[15]+0.3535533905932737*alpha[2]*fUpwind_r[15]+0.282842712474619*alpha[10]*fUpwind_r[14]+0.282842712474619*fUpwind_r[10]*alpha[14]+0.282842712474619*alpha[10]*fUpwind_r[13]+0.282842712474619*fUpwind_r[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[3]*alpha[10]+0.3535533905932737*alpha[4]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[5]*alpha[6]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[4] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11); 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[8] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11); 
  out[9] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11); 
  out[10] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11); 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11; 
  out[14] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv11; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11; 
  out[16] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11); 
  out[17] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11); 
  out[18] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11); 
  out[19] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11; 
  out[21] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11; 
  out[22] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11; 
  out[23] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11; 
  out[24] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv11; 
  out[25] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11); 
  out[26] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11); 
  out[27] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11); 
  out[28] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv11; 
  out[29] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv11; 
  out[30] += (1.5811388300841898*Ghat_l[3]-1.5811388300841898*Ghat_r[3])*dv11; 
  out[31] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11); 
  out[32] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv11; 
  out[33] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv11; 
  out[34] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv11; 
  out[35] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11); 
  out[36] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11); 
  out[37] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11); 
  out[38] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11); 
  out[39] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11); 
  out[40] += -(1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv11); 
  out[41] += (1.5811388300841895*Ghat_l[4]-1.5811388300841895*Ghat_r[4])*dv11; 
  out[42] += (1.5811388300841895*Ghat_l[5]-1.5811388300841895*Ghat_r[5])*dv11; 
  out[43] += (1.5811388300841895*Ghat_l[6]-1.5811388300841895*Ghat_r[6])*dv11; 
  out[44] += -(1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv11); 
  out[45] += -(1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv11); 
  out[46] += -(1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv11); 
  out[47] += (1.5811388300841898*Ghat_l[10]-1.5811388300841898*Ghat_r[10])*dv11; 

  return 0.;

} 
