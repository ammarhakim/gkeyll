#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 

  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 

  double alpha[20] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ax = &pots[8]; 
  const double *Ay = &pots[16]; 

  alpha[0] = -(2.4494897427831783*Ax[2]*dx11*wv2)+2.4494897427831783*Ay[1]*dx10*wv2-2.4494897427831783*phi[1]*dx10; 
  alpha[1] = -(2.4494897427831783*Ax[3]*dx11*wv2)+5.477225575051662*Ay[4]*dx10*wv2-5.477225575051662*phi[4]*dx10; 
  alpha[2] = -(5.477225575051662*Ax[5]*dx11*wv2)+2.4494897427831783*Ay[3]*dx10*wv2-2.4494897427831783*phi[3]*dx10; 
  alpha[3] = 0.7071067811865475*Ay[1]*dv2*dx10-0.7071067811865475*Ax[2]*dv2*dx11; 
  alpha[4] = -(5.477225575051662*Ax[7]*dx11*wv2)+5.477225575051662*Ay[6]*dx10*wv2-5.477225575051662*phi[6]*dx10; 
  alpha[5] = 1.5811388300841895*Ay[4]*dv2*dx10-0.7071067811865475*Ax[3]*dv2*dx11; 
  alpha[6] = 0.7071067811865475*Ay[3]*dv2*dx10-1.5811388300841895*Ax[5]*dv2*dx11; 
  alpha[7] = -(2.4494897427831783*Ax[6]*dx11*wv2); 
  alpha[8] = 2.4494897427831783*Ay[7]*dx10*wv2-2.4494897427831783*phi[7]*dx10; 
  alpha[10] = 1.5811388300841898*Ay[6]*dv2*dx10-1.5811388300841898*Ax[7]*dv2*dx11; 
  alpha[13] = -(0.7071067811865475*Ax[6]*dv2*dx11); 
  alpha[14] = 0.7071067811865475*Ay[7]*dv2*dx10; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (-(0.42426406871192845*(alpha[14]+alpha[13]))-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fr); 
  } 
  if (0.42426406871192845*(alpha[14]+alpha[13])+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.5303300858899104*alpha[14]-0.42426406871192845*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (-(0.3952847075210473*alpha[8])+0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (-(0.5303300858899104*alpha[14])+0.42426406871192845*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (-(0.42426406871192845*(alpha[14]+alpha[13]))+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fr); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fr); 
  } 
  if (0.42426406871192845*(alpha[14]+alpha[13])-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[14])+0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fr); 
  } 
  if (0.42426406871192845*alpha[14]-0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fr); 
  } 
  if (0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*(alpha[8]+alpha[7]) > 0) { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fr); 
  } 
  if (-(0.5303300858899104*(alpha[14]+alpha[13]))-0.3952847075210473*(alpha[8]+alpha[7])+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fr); 
  } 
  if (-(0.42426406871192845*alpha[14])+0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fr); 
  } 
  if (0.42426406871192845*alpha[14]-0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fr); 
  } 
  if (-(0.42426406871192845*(alpha[14]+alpha[13]))+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fr); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fr); 
  } 
  if (0.42426406871192845*(alpha[14]+alpha[13])-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fr); 
  } 
  if (0.5303300858899104*alpha[14]-0.42426406871192845*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fr); 
  } 
  if (-(0.3952847075210473*alpha[8])+0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fr); 
  } 
  if (-(0.5303300858899104*alpha[14])+0.42426406871192845*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fr); 
  } 
  if (-(0.42426406871192845*(alpha[14]+alpha[13]))-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fr); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fr); 
  } 
  if (0.42426406871192845*(alpha[14]+alpha[13])+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[14]*fUpwind_l[14]+0.3535533905932737*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[10]*fUpwind_l[10]+0.3535533905932737*alpha[8]*fUpwind_l[8]+0.3535533905932737*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[6]*fUpwind_l[6]+0.3535533905932737*alpha[5]*fUpwind_l[5]+0.3535533905932737*alpha[4]*fUpwind_l[4]+0.3535533905932737*alpha[3]*fUpwind_l[3]+0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[1]*fUpwind_l[1]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alpha[14]*fUpwind_l[18]+0.3162277660168379*alpha[10]*fUpwind_l[17]+0.31622776601683794*alpha[5]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[5]*alpha[13]+0.3535533905932737*alpha[8]*fUpwind_l[12]+0.31622776601683794*alpha[4]*fUpwind_l[11]+0.3535533905932737*alpha[6]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alpha[10]+0.3162277660168379*alpha[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.3162277660168379*alpha[10]*fUpwind_l[18]+0.3535533905932737*alpha[13]*fUpwind_l[17]+0.31622776601683794*alpha[6]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[6]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_l[12]+0.3535533905932737*alpha[7]*fUpwind_l[11]+0.3535533905932737*alpha[5]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[5]*alpha[10]+0.3162277660168379*alpha[2]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[2]*alpha[8]+0.3535533905932737*alpha[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3162277660168379*alpha[10]*fUpwind_l[19]+0.31622776601683794*alpha[6]*fUpwind_l[16]+0.31622776601683794*alpha[5]*fUpwind_l[15]+0.3535533905932737*alpha[8]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[8]*alpha[14]+0.3535533905932737*alpha[7]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[4]*alpha[10]+0.3162277660168379*alpha[3]*fUpwind_l[9]+0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.3162277660168379*alpha[6]*fUpwind_l[18]+0.3162277660168379*alpha[5]*fUpwind_l[17]+0.31622776601683794*alpha[10]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[10]*alpha[14]+0.31622776601683794*alpha[10]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[10]*alpha[13]+0.31622776601683794*alpha[2]*fUpwind_l[12]+0.31622776601683794*alpha[1]*fUpwind_l[11]+0.3535533905932737*alpha[3]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alpha[10]+0.3162277660168379*alpha[4]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[4]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.3162277660168379*alpha[6]*fUpwind_l[19]+0.3535533905932737*alpha[8]*fUpwind_l[18]+0.3162277660168379*alpha[4]*fUpwind_l[17]+0.31622776601683794*alpha[10]*fUpwind_l[16]+0.28284271247461906*alpha[13]*fUpwind_l[15]+0.31622776601683794*alpha[3]*fUpwind_l[15]+0.3535533905932737*fUpwind_l[12]*alpha[14]+0.31622776601683794*alpha[1]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[1]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_l[11]+0.3535533905932737*alpha[2]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[2]*alpha[10]+0.3162277660168379*alpha[5]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alpha[3]; 
  Ghat_l[6] = 0.3162277660168379*alpha[5]*fUpwind_l[19]+0.3162277660168379*alpha[4]*fUpwind_l[18]+0.3535533905932737*alpha[7]*fUpwind_l[17]+0.28284271247461906*alpha[14]*fUpwind_l[16]+0.31622776601683794*alpha[3]*fUpwind_l[16]+0.31622776601683794*alpha[10]*fUpwind_l[15]+0.31622776601683794*alpha[2]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[2]*alpha[14]+0.3535533905932737*fUpwind_l[11]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_l[12]+0.3535533905932737*alpha[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[1]*alpha[10]+0.3162277660168379*alpha[6]*fUpwind_l[9]+0.3162277660168379*alpha[6]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[6]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alpha[3]; 
  Ghat_l[7] = 0.3535533905932737*alpha[6]*fUpwind_l[17]+0.22587697572631277*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[3]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alpha[13]+0.3535533905932737*alpha[2]*fUpwind_l[11]+0.3162277660168379*alpha[10]*fUpwind_l[10]+0.22587697572631277*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[1]*fUpwind_l[1]; 
  Ghat_l[8] = 0.3535533905932737*alpha[5]*fUpwind_l[18]+0.22587697572631277*alpha[14]*fUpwind_l[14]+0.3535533905932737*alpha[3]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[3]*alpha[14]+0.3535533905932737*alpha[1]*fUpwind_l[12]+0.3162277660168379*alpha[10]*fUpwind_l[10]+0.22587697572631277*alpha[8]*fUpwind_l[8]+0.3535533905932737*alpha[0]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[0]*alpha[8]+0.3162277660168379*alpha[6]*fUpwind_l[6]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha[4]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[16]+0.3535533905932737*alpha[1]*fUpwind_l[15]+0.3162277660168379*alpha[14]*fUpwind_l[14]+0.3162277660168379*alpha[13]*fUpwind_l[13]+0.3162277660168379*alpha[10]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[9]+0.3162277660168379*alpha[6]*fUpwind_l[6]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[3]*fUpwind_l[3]; 
  Ghat_l[10] = 0.282842712474619*alpha[14]*fUpwind_l[19]+0.282842712474619*alpha[13]*fUpwind_l[19]+0.3162277660168379*alpha[3]*fUpwind_l[19]+0.3162277660168379*alpha[2]*fUpwind_l[18]+0.3162277660168379*alpha[1]*fUpwind_l[17]+0.31622776601683794*alpha[5]*fUpwind_l[16]+0.31622776601683794*alpha[6]*fUpwind_l[15]+0.31622776601683794*alpha[4]*fUpwind_l[14]+0.31622776601683794*fUpwind_l[4]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_l[13]+0.31622776601683794*fUpwind_l[4]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_l[12]+0.31622776601683794*alpha[5]*fUpwind_l[11]+0.3162277660168379*alpha[8]*fUpwind_l[10]+0.3162277660168379*alpha[7]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[9]*alpha[10]+0.3162277660168379*fUpwind_l[8]*alpha[10]+0.3162277660168379*fUpwind_l[7]*alpha[10]+0.3535533905932737*fUpwind_l[0]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.282842712474619*alpha[10]*fUpwind_l[18]+0.3162277660168379*alpha[14]*fUpwind_l[17]+0.22587697572631277*alpha[13]*fUpwind_l[17]+0.3535533905932737*alpha[3]*fUpwind_l[17]+0.3535533905932737*alpha[6]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[6]*alpha[13]+0.28284271247461906*alpha[4]*fUpwind_l[12]+0.3162277660168379*alpha[8]*fUpwind_l[11]+0.22587697572631277*alpha[7]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[11]+0.31622776601683794*alpha[5]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_l[4]+0.31622776601683794*fUpwind_l[1]*alpha[4]; 
  Ghat_l[12] = 0.22587697572631277*alpha[14]*fUpwind_l[18]+0.3162277660168379*alpha[13]*fUpwind_l[18]+0.3535533905932737*alpha[3]*fUpwind_l[18]+0.282842712474619*alpha[10]*fUpwind_l[17]+0.3535533905932737*alpha[5]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[5]*alpha[14]+0.22587697572631277*alpha[8]*fUpwind_l[12]+0.3162277660168379*alpha[7]*fUpwind_l[12]+0.3535533905932737*alpha[0]*fUpwind_l[12]+0.28284271247461906*alpha[4]*fUpwind_l[11]+0.31622776601683794*alpha[6]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[1]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_l[4]+0.31622776601683794*fUpwind_l[2]*alpha[4]; 
  Ghat_l[13] = 0.282842712474619*alpha[10]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[17]+0.28284271247461906*alpha[5]*fUpwind_l[15]+0.22587697572631277*alpha[7]*fUpwind_l[13]+0.3535533905932737*alpha[0]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[9]*alpha[13]+0.22587697572631277*fUpwind_l[7]*alpha[13]+0.3535533905932737*fUpwind_l[0]*alpha[13]+0.3535533905932737*alpha[6]*fUpwind_l[11]+0.31622776601683794*alpha[4]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_l[5]+0.31622776601683794*fUpwind_l[1]*alpha[5]; 
  Ghat_l[14] = 0.282842712474619*alpha[10]*fUpwind_l[19]+0.3535533905932737*alpha[1]*fUpwind_l[18]+0.28284271247461906*alpha[6]*fUpwind_l[16]+0.22587697572631277*alpha[8]*fUpwind_l[14]+0.3535533905932737*alpha[0]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[9]*alpha[14]+0.22587697572631277*fUpwind_l[8]*alpha[14]+0.3535533905932737*fUpwind_l[0]*alpha[14]+0.3535533905932737*alpha[5]*fUpwind_l[12]+0.31622776601683794*alpha[4]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[3]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_l[6]+0.31622776601683794*fUpwind_l[2]*alpha[6]; 
  Ghat_l[15] = 0.3535533905932737*alpha[2]*fUpwind_l[19]+0.3162277660168379*alpha[14]*fUpwind_l[18]+0.282842712474619*alpha[10]*fUpwind_l[17]+0.3535533905932737*alpha[4]*fUpwind_l[16]+0.3162277660168379*alpha[7]*fUpwind_l[15]+0.3535533905932737*alpha[0]*fUpwind_l[15]+0.28284271247461906*alpha[5]*fUpwind_l[13]+0.28284271247461906*fUpwind_l[5]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_l[9]+0.31622776601683794*alpha[3]*fUpwind_l[5]+0.31622776601683794*fUpwind_l[3]*alpha[5]; 
  Ghat_l[16] = 0.3535533905932737*alpha[1]*fUpwind_l[19]+0.282842712474619*alpha[10]*fUpwind_l[18]+0.3162277660168379*alpha[13]*fUpwind_l[17]+0.3162277660168379*alpha[8]*fUpwind_l[16]+0.3535533905932737*alpha[0]*fUpwind_l[16]+0.3535533905932737*alpha[4]*fUpwind_l[15]+0.28284271247461906*alpha[6]*fUpwind_l[14]+0.28284271247461906*fUpwind_l[6]*alpha[14]+0.31622776601683794*alpha[5]*fUpwind_l[10]+0.31622776601683794*fUpwind_l[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_l[9]+0.31622776601683794*alpha[3]*fUpwind_l[6]+0.31622776601683794*fUpwind_l[3]*alpha[6]; 
  Ghat_l[17] = 0.28284271247461906*alpha[5]*fUpwind_l[19]+0.28284271247461906*alpha[4]*fUpwind_l[18]+0.3162277660168379*alpha[8]*fUpwind_l[17]+0.22587697572631277*alpha[7]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[17]+0.3162277660168379*alpha[13]*fUpwind_l[16]+0.282842712474619*alpha[10]*fUpwind_l[15]+0.3162277660168379*fUpwind_l[11]*alpha[14]+0.3535533905932737*alpha[2]*fUpwind_l[13]+0.22587697572631277*fUpwind_l[11]*alpha[13]+0.3535533905932737*fUpwind_l[2]*alpha[13]+0.282842712474619*alpha[10]*fUpwind_l[12]+0.3535533905932737*alpha[3]*fUpwind_l[11]+0.3162277660168379*alpha[1]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[1]*alpha[10]+0.3535533905932737*alpha[6]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[4]*alpha[5]; 
  Ghat_l[18] = 0.28284271247461906*alpha[6]*fUpwind_l[19]+0.22587697572631277*alpha[8]*fUpwind_l[18]+0.3162277660168379*alpha[7]*fUpwind_l[18]+0.3535533905932737*alpha[0]*fUpwind_l[18]+0.28284271247461906*alpha[4]*fUpwind_l[17]+0.282842712474619*alpha[10]*fUpwind_l[16]+0.3162277660168379*alpha[14]*fUpwind_l[15]+0.3535533905932737*alpha[1]*fUpwind_l[14]+0.22587697572631277*fUpwind_l[12]*alpha[14]+0.3535533905932737*fUpwind_l[1]*alpha[14]+0.3162277660168379*fUpwind_l[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_l[12]+0.282842712474619*alpha[10]*fUpwind_l[11]+0.3162277660168379*alpha[2]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[2]*alpha[10]+0.3535533905932737*alpha[5]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[4]*alpha[6]; 
  Ghat_l[19] = 0.3162277660168379*alpha[8]*fUpwind_l[19]+0.3162277660168379*alpha[7]*fUpwind_l[19]+0.3535533905932737*alpha[0]*fUpwind_l[19]+0.28284271247461906*alpha[6]*fUpwind_l[18]+0.28284271247461906*alpha[5]*fUpwind_l[17]+0.3535533905932737*alpha[1]*fUpwind_l[16]+0.3535533905932737*alpha[2]*fUpwind_l[15]+0.282842712474619*alpha[10]*fUpwind_l[14]+0.282842712474619*fUpwind_l[10]*alpha[14]+0.282842712474619*alpha[10]*fUpwind_l[13]+0.282842712474619*fUpwind_l[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[3]*alpha[10]+0.3535533905932737*alpha[4]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[5]*alpha[6]; 

  Ghat_r[0] = 0.3535533905932737*alpha[14]*fUpwind_r[14]+0.3535533905932737*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[10]*fUpwind_r[10]+0.3535533905932737*alpha[8]*fUpwind_r[8]+0.3535533905932737*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[6]*fUpwind_r[6]+0.3535533905932737*alpha[5]*fUpwind_r[5]+0.3535533905932737*alpha[4]*fUpwind_r[4]+0.3535533905932737*alpha[3]*fUpwind_r[3]+0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[1]*fUpwind_r[1]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alpha[14]*fUpwind_r[18]+0.3162277660168379*alpha[10]*fUpwind_r[17]+0.31622776601683794*alpha[5]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[5]*alpha[13]+0.3535533905932737*alpha[8]*fUpwind_r[12]+0.31622776601683794*alpha[4]*fUpwind_r[11]+0.3535533905932737*alpha[6]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alpha[10]+0.3162277660168379*alpha[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.3162277660168379*alpha[10]*fUpwind_r[18]+0.3535533905932737*alpha[13]*fUpwind_r[17]+0.31622776601683794*alpha[6]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[6]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_r[12]+0.3535533905932737*alpha[7]*fUpwind_r[11]+0.3535533905932737*alpha[5]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[5]*alpha[10]+0.3162277660168379*alpha[2]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[2]*alpha[8]+0.3535533905932737*alpha[3]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3162277660168379*alpha[10]*fUpwind_r[19]+0.31622776601683794*alpha[6]*fUpwind_r[16]+0.31622776601683794*alpha[5]*fUpwind_r[15]+0.3535533905932737*alpha[8]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[8]*alpha[14]+0.3535533905932737*alpha[7]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[4]*alpha[10]+0.3162277660168379*alpha[3]*fUpwind_r[9]+0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.3162277660168379*alpha[6]*fUpwind_r[18]+0.3162277660168379*alpha[5]*fUpwind_r[17]+0.31622776601683794*alpha[10]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[10]*alpha[14]+0.31622776601683794*alpha[10]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[10]*alpha[13]+0.31622776601683794*alpha[2]*fUpwind_r[12]+0.31622776601683794*alpha[1]*fUpwind_r[11]+0.3535533905932737*alpha[3]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alpha[10]+0.3162277660168379*alpha[4]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[4]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.3162277660168379*alpha[6]*fUpwind_r[19]+0.3535533905932737*alpha[8]*fUpwind_r[18]+0.3162277660168379*alpha[4]*fUpwind_r[17]+0.31622776601683794*alpha[10]*fUpwind_r[16]+0.28284271247461906*alpha[13]*fUpwind_r[15]+0.31622776601683794*alpha[3]*fUpwind_r[15]+0.3535533905932737*fUpwind_r[12]*alpha[14]+0.31622776601683794*alpha[1]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[1]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_r[11]+0.3535533905932737*alpha[2]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[2]*alpha[10]+0.3162277660168379*alpha[5]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alpha[3]; 
  Ghat_r[6] = 0.3162277660168379*alpha[5]*fUpwind_r[19]+0.3162277660168379*alpha[4]*fUpwind_r[18]+0.3535533905932737*alpha[7]*fUpwind_r[17]+0.28284271247461906*alpha[14]*fUpwind_r[16]+0.31622776601683794*alpha[3]*fUpwind_r[16]+0.31622776601683794*alpha[10]*fUpwind_r[15]+0.31622776601683794*alpha[2]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[2]*alpha[14]+0.3535533905932737*fUpwind_r[11]*alpha[13]+0.31622776601683794*alpha[10]*fUpwind_r[12]+0.3535533905932737*alpha[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[1]*alpha[10]+0.3162277660168379*alpha[6]*fUpwind_r[9]+0.3162277660168379*alpha[6]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[6]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alpha[3]; 
  Ghat_r[7] = 0.3535533905932737*alpha[6]*fUpwind_r[17]+0.22587697572631277*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[3]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alpha[13]+0.3535533905932737*alpha[2]*fUpwind_r[11]+0.3162277660168379*alpha[10]*fUpwind_r[10]+0.22587697572631277*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[1]*fUpwind_r[1]; 
  Ghat_r[8] = 0.3535533905932737*alpha[5]*fUpwind_r[18]+0.22587697572631277*alpha[14]*fUpwind_r[14]+0.3535533905932737*alpha[3]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[3]*alpha[14]+0.3535533905932737*alpha[1]*fUpwind_r[12]+0.3162277660168379*alpha[10]*fUpwind_r[10]+0.22587697572631277*alpha[8]*fUpwind_r[8]+0.3535533905932737*alpha[0]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[0]*alpha[8]+0.3162277660168379*alpha[6]*fUpwind_r[6]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha[4]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[16]+0.3535533905932737*alpha[1]*fUpwind_r[15]+0.3162277660168379*alpha[14]*fUpwind_r[14]+0.3162277660168379*alpha[13]*fUpwind_r[13]+0.3162277660168379*alpha[10]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[9]+0.3162277660168379*alpha[6]*fUpwind_r[6]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[3]*fUpwind_r[3]; 
  Ghat_r[10] = 0.282842712474619*alpha[14]*fUpwind_r[19]+0.282842712474619*alpha[13]*fUpwind_r[19]+0.3162277660168379*alpha[3]*fUpwind_r[19]+0.3162277660168379*alpha[2]*fUpwind_r[18]+0.3162277660168379*alpha[1]*fUpwind_r[17]+0.31622776601683794*alpha[5]*fUpwind_r[16]+0.31622776601683794*alpha[6]*fUpwind_r[15]+0.31622776601683794*alpha[4]*fUpwind_r[14]+0.31622776601683794*fUpwind_r[4]*alpha[14]+0.31622776601683794*alpha[4]*fUpwind_r[13]+0.31622776601683794*fUpwind_r[4]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_r[12]+0.31622776601683794*alpha[5]*fUpwind_r[11]+0.3162277660168379*alpha[8]*fUpwind_r[10]+0.3162277660168379*alpha[7]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[9]*alpha[10]+0.3162277660168379*fUpwind_r[8]*alpha[10]+0.3162277660168379*fUpwind_r[7]*alpha[10]+0.3535533905932737*fUpwind_r[0]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.282842712474619*alpha[10]*fUpwind_r[18]+0.3162277660168379*alpha[14]*fUpwind_r[17]+0.22587697572631277*alpha[13]*fUpwind_r[17]+0.3535533905932737*alpha[3]*fUpwind_r[17]+0.3535533905932737*alpha[6]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[6]*alpha[13]+0.28284271247461906*alpha[4]*fUpwind_r[12]+0.3162277660168379*alpha[8]*fUpwind_r[11]+0.22587697572631277*alpha[7]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[11]+0.31622776601683794*alpha[5]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_r[4]+0.31622776601683794*fUpwind_r[1]*alpha[4]; 
  Ghat_r[12] = 0.22587697572631277*alpha[14]*fUpwind_r[18]+0.3162277660168379*alpha[13]*fUpwind_r[18]+0.3535533905932737*alpha[3]*fUpwind_r[18]+0.282842712474619*alpha[10]*fUpwind_r[17]+0.3535533905932737*alpha[5]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[5]*alpha[14]+0.22587697572631277*alpha[8]*fUpwind_r[12]+0.3162277660168379*alpha[7]*fUpwind_r[12]+0.3535533905932737*alpha[0]*fUpwind_r[12]+0.28284271247461906*alpha[4]*fUpwind_r[11]+0.31622776601683794*alpha[6]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[1]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_r[4]+0.31622776601683794*fUpwind_r[2]*alpha[4]; 
  Ghat_r[13] = 0.282842712474619*alpha[10]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[17]+0.28284271247461906*alpha[5]*fUpwind_r[15]+0.22587697572631277*alpha[7]*fUpwind_r[13]+0.3535533905932737*alpha[0]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[9]*alpha[13]+0.22587697572631277*fUpwind_r[7]*alpha[13]+0.3535533905932737*fUpwind_r[0]*alpha[13]+0.3535533905932737*alpha[6]*fUpwind_r[11]+0.31622776601683794*alpha[4]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alpha[7]+0.31622776601683794*alpha[1]*fUpwind_r[5]+0.31622776601683794*fUpwind_r[1]*alpha[5]; 
  Ghat_r[14] = 0.282842712474619*alpha[10]*fUpwind_r[19]+0.3535533905932737*alpha[1]*fUpwind_r[18]+0.28284271247461906*alpha[6]*fUpwind_r[16]+0.22587697572631277*alpha[8]*fUpwind_r[14]+0.3535533905932737*alpha[0]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[9]*alpha[14]+0.22587697572631277*fUpwind_r[8]*alpha[14]+0.3535533905932737*fUpwind_r[0]*alpha[14]+0.3535533905932737*alpha[5]*fUpwind_r[12]+0.31622776601683794*alpha[4]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[4]*alpha[10]+0.3535533905932737*alpha[3]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[3]*alpha[8]+0.31622776601683794*alpha[2]*fUpwind_r[6]+0.31622776601683794*fUpwind_r[2]*alpha[6]; 
  Ghat_r[15] = 0.3535533905932737*alpha[2]*fUpwind_r[19]+0.3162277660168379*alpha[14]*fUpwind_r[18]+0.282842712474619*alpha[10]*fUpwind_r[17]+0.3535533905932737*alpha[4]*fUpwind_r[16]+0.3162277660168379*alpha[7]*fUpwind_r[15]+0.3535533905932737*alpha[0]*fUpwind_r[15]+0.28284271247461906*alpha[5]*fUpwind_r[13]+0.28284271247461906*fUpwind_r[5]*alpha[13]+0.31622776601683794*alpha[6]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[6]*alpha[10]+0.3535533905932737*alpha[1]*fUpwind_r[9]+0.31622776601683794*alpha[3]*fUpwind_r[5]+0.31622776601683794*fUpwind_r[3]*alpha[5]; 
  Ghat_r[16] = 0.3535533905932737*alpha[1]*fUpwind_r[19]+0.282842712474619*alpha[10]*fUpwind_r[18]+0.3162277660168379*alpha[13]*fUpwind_r[17]+0.3162277660168379*alpha[8]*fUpwind_r[16]+0.3535533905932737*alpha[0]*fUpwind_r[16]+0.3535533905932737*alpha[4]*fUpwind_r[15]+0.28284271247461906*alpha[6]*fUpwind_r[14]+0.28284271247461906*fUpwind_r[6]*alpha[14]+0.31622776601683794*alpha[5]*fUpwind_r[10]+0.31622776601683794*fUpwind_r[5]*alpha[10]+0.3535533905932737*alpha[2]*fUpwind_r[9]+0.31622776601683794*alpha[3]*fUpwind_r[6]+0.31622776601683794*fUpwind_r[3]*alpha[6]; 
  Ghat_r[17] = 0.28284271247461906*alpha[5]*fUpwind_r[19]+0.28284271247461906*alpha[4]*fUpwind_r[18]+0.3162277660168379*alpha[8]*fUpwind_r[17]+0.22587697572631277*alpha[7]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[17]+0.3162277660168379*alpha[13]*fUpwind_r[16]+0.282842712474619*alpha[10]*fUpwind_r[15]+0.3162277660168379*fUpwind_r[11]*alpha[14]+0.3535533905932737*alpha[2]*fUpwind_r[13]+0.22587697572631277*fUpwind_r[11]*alpha[13]+0.3535533905932737*fUpwind_r[2]*alpha[13]+0.282842712474619*alpha[10]*fUpwind_r[12]+0.3535533905932737*alpha[3]*fUpwind_r[11]+0.3162277660168379*alpha[1]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[1]*alpha[10]+0.3535533905932737*alpha[6]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[4]*alpha[5]; 
  Ghat_r[18] = 0.28284271247461906*alpha[6]*fUpwind_r[19]+0.22587697572631277*alpha[8]*fUpwind_r[18]+0.3162277660168379*alpha[7]*fUpwind_r[18]+0.3535533905932737*alpha[0]*fUpwind_r[18]+0.28284271247461906*alpha[4]*fUpwind_r[17]+0.282842712474619*alpha[10]*fUpwind_r[16]+0.3162277660168379*alpha[14]*fUpwind_r[15]+0.3535533905932737*alpha[1]*fUpwind_r[14]+0.22587697572631277*fUpwind_r[12]*alpha[14]+0.3535533905932737*fUpwind_r[1]*alpha[14]+0.3162277660168379*fUpwind_r[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_r[12]+0.282842712474619*alpha[10]*fUpwind_r[11]+0.3162277660168379*alpha[2]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[2]*alpha[10]+0.3535533905932737*alpha[5]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[4]*alpha[6]; 
  Ghat_r[19] = 0.3162277660168379*alpha[8]*fUpwind_r[19]+0.3162277660168379*alpha[7]*fUpwind_r[19]+0.3535533905932737*alpha[0]*fUpwind_r[19]+0.28284271247461906*alpha[6]*fUpwind_r[18]+0.28284271247461906*alpha[5]*fUpwind_r[17]+0.3535533905932737*alpha[1]*fUpwind_r[16]+0.3535533905932737*alpha[2]*fUpwind_r[15]+0.282842712474619*alpha[10]*fUpwind_r[14]+0.282842712474619*fUpwind_r[10]*alpha[14]+0.282842712474619*alpha[10]*fUpwind_r[13]+0.282842712474619*fUpwind_r[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[3]*alpha[10]+0.3535533905932737*alpha[4]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[5]*alpha[6]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10); 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10); 
  out[7] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10); 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[10] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10); 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[13] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[15] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10); 
  out[16] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[17] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10); 
  out[18] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10); 
  out[19] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[21] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10); 
  out[22] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10); 
  out[23] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv10; 
  out[24] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv10; 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[26] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[27] += (1.5811388300841898*Ghat_l[3]-1.5811388300841898*Ghat_r[3])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[29] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv10; 
  out[30] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10); 
  out[31] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10); 
  out[32] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10); 
  out[33] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10); 
  out[34] += (1.5811388300841895*Ghat_l[4]-1.5811388300841895*Ghat_r[4])*dv10; 
  out[35] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv10; 
  out[36] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv10; 
  out[37] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10); 
  out[38] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10); 
  out[39] += (1.5811388300841895*Ghat_l[5]-1.5811388300841895*Ghat_r[5])*dv10; 
  out[40] += (1.5811388300841895*Ghat_l[6]-1.5811388300841895*Ghat_r[6])*dv10; 
  out[41] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv10; 
  out[42] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10); 
  out[43] += -(1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv10); 
  out[44] += -(1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv10); 
  out[45] += -(1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv10); 
  out[46] += (1.5811388300841898*Ghat_l[10]-1.5811388300841898*Ghat_r[10])*dv10; 
  out[47] += -(1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv10); 

  return 0.;

} 
