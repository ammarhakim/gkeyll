#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_tensor_4x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E0 = &field[0]; 
  const double *B0 = &field[9]; 
  const double *B1 = &field[12]; 
  const double *B2 = &field[15]; 

  double alpha[27] = {0.0}; 

  alpha[0] = (-2.0*B1[0]*wv3)+2.0*B2[0]*wv2+2.0*E0[0]; 
  alpha[1] = (-2.0*B1[1]*wv3)+2.0*B2[1]*wv2+2.0*E0[1]; 
  alpha[2] = 0.5773502691896258*B2[0]*dv2; 
  alpha[3] = -0.5773502691896258*B1[0]*dv3; 
  alpha[4] = 0.5773502691896258*B2[1]*dv2; 
  alpha[5] = -0.5773502691896258*B1[1]*dv3; 
  alpha[7] = (-2.0*B1[2]*wv3)+2.0*B2[2]*wv2+2.0*E0[2]; 
  alpha[11] = 0.5773502691896257*B2[2]*dv2; 
  alpha[13] = -0.5773502691896257*B1[2]*dv3; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[27] = {0.0};;
  double fUpwind_r[27] = {0.0};
  double Ghat_l[27] = {0.0}; 
  double Ghat_r[27] = {0.0}; 

  if ((-0.4242640687119285*(alpha[13]+alpha[11]))+0.3162277660168379*alpha[7]+0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[0] = tensor_4x_p2_surfx2_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = tensor_4x_p2_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = tensor_4x_p2_surfx2_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = tensor_4x_p2_surfx2_eval_quad_node_0_l(fr); 
  } 
  if ((-0.4242640687119285*alpha[11])+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[1] = tensor_4x_p2_surfx2_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = tensor_4x_p2_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = tensor_4x_p2_surfx2_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = tensor_4x_p2_surfx2_eval_quad_node_1_l(fr); 
  } 
  if (0.4242640687119285*alpha[13]-0.4242640687119285*alpha[11]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[2] = tensor_4x_p2_surfx2_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = tensor_4x_p2_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = tensor_4x_p2_surfx2_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = tensor_4x_p2_surfx2_eval_quad_node_2_l(fr); 
  } 
  if ((-0.4242640687119285*alpha[13])+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[3] = tensor_4x_p2_surfx2_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = tensor_4x_p2_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = tensor_4x_p2_surfx2_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = tensor_4x_p2_surfx2_eval_quad_node_3_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[4] = tensor_4x_p2_surfx2_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = tensor_4x_p2_surfx2_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = tensor_4x_p2_surfx2_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = tensor_4x_p2_surfx2_eval_quad_node_4_l(fr); 
  } 
  if (0.4242640687119285*alpha[13]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[5] = tensor_4x_p2_surfx2_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = tensor_4x_p2_surfx2_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = tensor_4x_p2_surfx2_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = tensor_4x_p2_surfx2_eval_quad_node_5_l(fr); 
  } 
  if ((-0.4242640687119285*alpha[13])+0.4242640687119285*alpha[11]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[6] = tensor_4x_p2_surfx2_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = tensor_4x_p2_surfx2_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = tensor_4x_p2_surfx2_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = tensor_4x_p2_surfx2_eval_quad_node_6_l(fr); 
  } 
  if (0.4242640687119285*alpha[11]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[7] = tensor_4x_p2_surfx2_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = tensor_4x_p2_surfx2_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = tensor_4x_p2_surfx2_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = tensor_4x_p2_surfx2_eval_quad_node_7_l(fr); 
  } 
  if (0.4242640687119285*(alpha[13]+alpha[11])+0.3162277660168379*alpha[7]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[8] = tensor_4x_p2_surfx2_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = tensor_4x_p2_surfx2_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = tensor_4x_p2_surfx2_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = tensor_4x_p2_surfx2_eval_quad_node_8_l(fr); 
  } 
  if (0.5303300858899104*(alpha[13]+alpha[11])-0.3952847075210473*alpha[7]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[9] = tensor_4x_p2_surfx2_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = tensor_4x_p2_surfx2_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = tensor_4x_p2_surfx2_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = tensor_4x_p2_surfx2_eval_quad_node_9_l(fr); 
  } 
  if (0.5303300858899104*alpha[11]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[10] = tensor_4x_p2_surfx2_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = tensor_4x_p2_surfx2_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = tensor_4x_p2_surfx2_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = tensor_4x_p2_surfx2_eval_quad_node_10_l(fr); 
  } 
  if ((-0.5303300858899104*alpha[13])+0.5303300858899104*alpha[11]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[11] = tensor_4x_p2_surfx2_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = tensor_4x_p2_surfx2_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = tensor_4x_p2_surfx2_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = tensor_4x_p2_surfx2_eval_quad_node_11_l(fr); 
  } 
  if (0.5303300858899104*alpha[13]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[12] = tensor_4x_p2_surfx2_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = tensor_4x_p2_surfx2_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = tensor_4x_p2_surfx2_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = tensor_4x_p2_surfx2_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7] > 0) { 
    fUpwindQuad_l[13] = tensor_4x_p2_surfx2_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = tensor_4x_p2_surfx2_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = tensor_4x_p2_surfx2_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = tensor_4x_p2_surfx2_eval_quad_node_13_l(fr); 
  } 
  if ((-0.5303300858899104*alpha[13])-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[14] = tensor_4x_p2_surfx2_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = tensor_4x_p2_surfx2_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = tensor_4x_p2_surfx2_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = tensor_4x_p2_surfx2_eval_quad_node_14_l(fr); 
  } 
  if (0.5303300858899104*alpha[13]-0.5303300858899104*alpha[11]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[15] = tensor_4x_p2_surfx2_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = tensor_4x_p2_surfx2_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = tensor_4x_p2_surfx2_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = tensor_4x_p2_surfx2_eval_quad_node_15_l(fr); 
  } 
  if ((-0.5303300858899104*alpha[11])-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[16] = tensor_4x_p2_surfx2_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = tensor_4x_p2_surfx2_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = tensor_4x_p2_surfx2_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = tensor_4x_p2_surfx2_eval_quad_node_16_l(fr); 
  } 
  if ((-0.5303300858899104*(alpha[13]+alpha[11]))-0.3952847075210473*alpha[7]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[17] = tensor_4x_p2_surfx2_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = tensor_4x_p2_surfx2_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = tensor_4x_p2_surfx2_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = tensor_4x_p2_surfx2_eval_quad_node_17_l(fr); 
  } 
  if ((-0.4242640687119285*(alpha[13]+alpha[11]))+0.3162277660168379*alpha[7]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[18] = tensor_4x_p2_surfx2_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = tensor_4x_p2_surfx2_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = tensor_4x_p2_surfx2_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = tensor_4x_p2_surfx2_eval_quad_node_18_l(fr); 
  } 
  if ((-0.4242640687119285*alpha[11])+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[19] = tensor_4x_p2_surfx2_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = tensor_4x_p2_surfx2_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = tensor_4x_p2_surfx2_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = tensor_4x_p2_surfx2_eval_quad_node_19_l(fr); 
  } 
  if (0.4242640687119285*alpha[13]-0.4242640687119285*alpha[11]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[20] = tensor_4x_p2_surfx2_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = tensor_4x_p2_surfx2_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = tensor_4x_p2_surfx2_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = tensor_4x_p2_surfx2_eval_quad_node_20_l(fr); 
  } 
  if ((-0.4242640687119285*alpha[13])+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[21] = tensor_4x_p2_surfx2_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = tensor_4x_p2_surfx2_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = tensor_4x_p2_surfx2_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = tensor_4x_p2_surfx2_eval_quad_node_21_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[22] = tensor_4x_p2_surfx2_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = tensor_4x_p2_surfx2_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = tensor_4x_p2_surfx2_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = tensor_4x_p2_surfx2_eval_quad_node_22_l(fr); 
  } 
  if (0.4242640687119285*alpha[13]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[23] = tensor_4x_p2_surfx2_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = tensor_4x_p2_surfx2_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = tensor_4x_p2_surfx2_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = tensor_4x_p2_surfx2_eval_quad_node_23_l(fr); 
  } 
  if ((-0.4242640687119285*alpha[13])+0.4242640687119285*alpha[11]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[24] = tensor_4x_p2_surfx2_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = tensor_4x_p2_surfx2_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = tensor_4x_p2_surfx2_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = tensor_4x_p2_surfx2_eval_quad_node_24_l(fr); 
  } 
  if (0.4242640687119285*alpha[11]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[25] = tensor_4x_p2_surfx2_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = tensor_4x_p2_surfx2_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = tensor_4x_p2_surfx2_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = tensor_4x_p2_surfx2_eval_quad_node_25_l(fr); 
  } 
  if (0.4242640687119285*(alpha[13]+alpha[11])+0.3162277660168379*alpha[7]+0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[26] = tensor_4x_p2_surfx2_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = tensor_4x_p2_surfx2_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = tensor_4x_p2_surfx2_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = tensor_4x_p2_surfx2_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_4x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  tensor_4x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[5]*fUpwind_l[5]+0.3535533905932737*alpha[4]*fUpwind_l[4]+0.3535533905932737*alpha[3]*fUpwind_l[3]+0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[1]*fUpwind_l[1]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3162277660168379*alpha[5]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[5]*alpha[13]+0.3162277660168379*alpha[4]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[4]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.3162277660168379*alpha[11]*fUpwind_l[20]+0.3535533905932737*alpha[13]*fUpwind_l[17]+0.3162277660168379*alpha[4]*fUpwind_l[12]+0.3535533905932737*alpha[7]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[7]*alpha[11]+0.3535533905932737*alpha[5]*fUpwind_l[10]+0.3162277660168379*alpha[2]*fUpwind_l[8]+0.3535533905932737*alpha[3]*fUpwind_l[6]+0.3535533905932737*alpha[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3162277660168379*alpha[13]*fUpwind_l[21]+0.3535533905932737*alpha[11]*fUpwind_l[17]+0.3162277660168379*alpha[5]*fUpwind_l[15]+0.3535533905932737*alpha[7]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_l[10]+0.3162277660168379*alpha[3]*fUpwind_l[9]+0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*alpha[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.2828427124746191*alpha[4]*fUpwind_l[20]+0.3162277660168379*alpha[5]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[10]*alpha[13]+0.2828427124746191*alpha[11]*fUpwind_l[12]+0.3162277660168379*alpha[2]*fUpwind_l[12]+0.3162277660168379*alpha[1]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[1]*alpha[11]+0.3535533905932737*alpha[3]*fUpwind_l[10]+0.3162277660168379*alpha[4]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_l[6]+0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.2828427124746191*alpha[5]*fUpwind_l[21]+0.3162277660168379*alpha[4]*fUpwind_l[17]+0.2828427124746191*alpha[13]*fUpwind_l[15]+0.3162277660168379*alpha[3]*fUpwind_l[15]+0.3162277660168379*alpha[1]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[1]*alpha[13]+0.3162277660168379*fUpwind_l[10]*alpha[11]+0.3535533905932737*alpha[2]*fUpwind_l[10]+0.3162277660168379*alpha[5]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_l[6]+0.3535533905932737*alpha[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alpha[3]; 
  Ghat_l[6] = 0.3162277660168379*alpha[13]*fUpwind_l[24]+0.3162277660168379*alpha[11]*fUpwind_l[23]+0.3162277660168379*alpha[5]*fUpwind_l[19]+0.3162277660168379*alpha[4]*fUpwind_l[18]+0.3535533905932737*alpha[7]*fUpwind_l[17]+0.3162277660168379*alpha[3]*fUpwind_l[16]+0.3162277660168379*alpha[2]*fUpwind_l[14]+0.3535533905932737*alpha[11]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[11]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*alpha[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alpha[3]; 
  Ghat_l[7] = 0.2258769757263128*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[3]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alpha[13]+0.2258769757263128*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[2]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2]*alpha[11]+0.2258769757263128*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[1]*fUpwind_l[1]; 
  Ghat_l[8] = 0.3535533905932737*alpha[13]*fUpwind_l[23]+0.3535533905932737*alpha[7]*fUpwind_l[20]+0.3535533905932737*alpha[5]*fUpwind_l[18]+0.3535533905932737*alpha[3]*fUpwind_l[14]+0.3535533905932737*alpha[1]*fUpwind_l[12]+0.3162277660168379*alpha[11]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha[11]*fUpwind_l[24]+0.3535533905932737*alpha[7]*fUpwind_l[21]+0.3535533905932737*alpha[4]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[16]+0.3535533905932737*alpha[1]*fUpwind_l[15]+0.3162277660168379*alpha[13]*fUpwind_l[13]+0.3535533905932737*alpha[0]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[3]*fUpwind_l[3]; 
  Ghat_l[10] = 0.2828427124746191*alpha[5]*fUpwind_l[24]+0.2828427124746191*alpha[4]*fUpwind_l[23]+0.282842712474619*alpha[13]*fUpwind_l[19]+0.3162277660168379*alpha[3]*fUpwind_l[19]+0.282842712474619*alpha[11]*fUpwind_l[18]+0.3162277660168379*alpha[2]*fUpwind_l[18]+0.3162277660168379*alpha[1]*fUpwind_l[17]+0.3162277660168379*alpha[5]*fUpwind_l[16]+0.3162277660168379*alpha[4]*fUpwind_l[14]+0.3162277660168379*alpha[4]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[4]*alpha[13]+0.3162277660168379*alpha[5]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[5]*alpha[11]+0.3162277660168379*alpha[7]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[10]+0.3535533905932737*alpha[1]*fUpwind_l[6]+0.3535533905932737*alpha[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.2020305089104422*alpha[11]*fUpwind_l[20]+0.3162277660168379*alpha[2]*fUpwind_l[20]+0.2258769757263128*alpha[13]*fUpwind_l[17]+0.3535533905932737*alpha[3]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[6]*alpha[13]+0.2828427124746191*alpha[4]*fUpwind_l[12]+0.2258769757263128*alpha[7]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[8]*alpha[11]+0.2258769757263128*fUpwind_l[7]*alpha[11]+0.3535533905932737*fUpwind_l[0]*alpha[11]+0.3162277660168379*alpha[5]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[1]*alpha[4]; 
  Ghat_l[12] = 0.3162277660168379*alpha[5]*fUpwind_l[23]+0.3162277660168379*alpha[1]*fUpwind_l[20]+0.3162277660168379*alpha[13]*fUpwind_l[18]+0.3535533905932737*alpha[3]*fUpwind_l[18]+0.3535533905932737*alpha[5]*fUpwind_l[14]+0.3162277660168379*alpha[7]*fUpwind_l[12]+0.3535533905932737*alpha[0]*fUpwind_l[12]+0.2828427124746191*alpha[4]*fUpwind_l[11]+0.2828427124746191*fUpwind_l[4]*alpha[11]+0.3535533905932737*alpha[1]*fUpwind_l[8]+0.3162277660168379*alpha[2]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[2]*alpha[4]; 
  Ghat_l[13] = 0.2020305089104422*alpha[13]*fUpwind_l[21]+0.3162277660168379*alpha[3]*fUpwind_l[21]+0.2258769757263128*alpha[11]*fUpwind_l[17]+0.3535533905932737*alpha[2]*fUpwind_l[17]+0.2828427124746191*alpha[5]*fUpwind_l[15]+0.2258769757263128*alpha[7]*fUpwind_l[13]+0.3535533905932737*alpha[0]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[9]*alpha[13]+0.2258769757263128*fUpwind_l[7]*alpha[13]+0.3535533905932737*fUpwind_l[0]*alpha[13]+0.3535533905932737*fUpwind_l[6]*alpha[11]+0.3162277660168379*alpha[4]*fUpwind_l[10]+0.3535533905932737*alpha[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[1]*alpha[5]; 
  Ghat_l[14] = 0.3162277660168379*alpha[13]*fUpwind_l[26]+0.3162277660168379*alpha[5]*fUpwind_l[25]+0.3535533905932737*alpha[7]*fUpwind_l[23]+0.3162277660168379*alpha[3]*fUpwind_l[22]+0.3535533905932737*alpha[13]*fUpwind_l[20]+0.3535533905932737*alpha[1]*fUpwind_l[18]+0.3162277660168379*alpha[11]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[14]+0.3535533905932737*alpha[5]*fUpwind_l[12]+0.3162277660168379*alpha[4]*fUpwind_l[10]+0.3535533905932737*alpha[3]*fUpwind_l[8]+0.3162277660168379*alpha[2]*fUpwind_l[6]; 
  Ghat_l[15] = 0.3162277660168379*alpha[4]*fUpwind_l[24]+0.3162277660168379*alpha[1]*fUpwind_l[21]+0.3162277660168379*alpha[11]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[19]+0.3535533905932737*alpha[4]*fUpwind_l[16]+0.3162277660168379*alpha[7]*fUpwind_l[15]+0.3535533905932737*alpha[0]*fUpwind_l[15]+0.2828427124746191*alpha[5]*fUpwind_l[13]+0.2828427124746191*fUpwind_l[5]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_l[9]+0.3162277660168379*alpha[3]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[3]*alpha[5]; 
  Ghat_l[16] = 0.3162277660168379*alpha[11]*fUpwind_l[26]+0.3162277660168379*alpha[4]*fUpwind_l[25]+0.3535533905932737*alpha[7]*fUpwind_l[24]+0.3162277660168379*alpha[2]*fUpwind_l[22]+0.3535533905932737*alpha[11]*fUpwind_l[21]+0.3535533905932737*alpha[1]*fUpwind_l[19]+0.3162277660168379*alpha[13]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[16]+0.3535533905932737*alpha[4]*fUpwind_l[15]+0.3162277660168379*alpha[5]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[9]+0.3162277660168379*alpha[3]*fUpwind_l[6]; 
  Ghat_l[17] = 0.2020305089104422*alpha[13]*fUpwind_l[24]+0.3162277660168379*alpha[3]*fUpwind_l[24]+0.2020305089104422*alpha[11]*fUpwind_l[23]+0.3162277660168379*alpha[2]*fUpwind_l[23]+0.2828427124746191*alpha[5]*fUpwind_l[19]+0.2828427124746191*alpha[4]*fUpwind_l[18]+0.2258769757263128*alpha[7]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[17]+0.3162277660168379*alpha[13]*fUpwind_l[16]+0.3162277660168379*alpha[11]*fUpwind_l[14]+0.2258769757263128*alpha[11]*fUpwind_l[13]+0.3535533905932737*alpha[2]*fUpwind_l[13]+0.2258769757263128*fUpwind_l[11]*alpha[13]+0.3535533905932737*fUpwind_l[2]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[3]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[4]*alpha[5]; 
  Ghat_l[18] = 0.2828427124746191*alpha[5]*fUpwind_l[26]+0.282842712474619*alpha[13]*fUpwind_l[25]+0.3162277660168379*alpha[3]*fUpwind_l[25]+0.3162277660168379*alpha[1]*fUpwind_l[23]+0.3162277660168379*alpha[5]*fUpwind_l[22]+0.3162277660168379*alpha[5]*fUpwind_l[20]+0.3162277660168379*alpha[7]*fUpwind_l[18]+0.3535533905932737*alpha[0]*fUpwind_l[18]+0.2828427124746191*alpha[4]*fUpwind_l[17]+0.3535533905932737*alpha[1]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_l[12]+0.282842712474619*fUpwind_l[10]*alpha[11]+0.3162277660168379*alpha[2]*fUpwind_l[10]+0.3535533905932737*alpha[5]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[6]; 
  Ghat_l[19] = 0.2828427124746191*alpha[4]*fUpwind_l[26]+0.282842712474619*alpha[11]*fUpwind_l[25]+0.3162277660168379*alpha[2]*fUpwind_l[25]+0.3162277660168379*alpha[1]*fUpwind_l[24]+0.3162277660168379*alpha[4]*fUpwind_l[22]+0.3162277660168379*alpha[4]*fUpwind_l[21]+0.3162277660168379*alpha[7]*fUpwind_l[19]+0.3535533905932737*alpha[0]*fUpwind_l[19]+0.2828427124746191*alpha[5]*fUpwind_l[17]+0.3535533905932737*alpha[1]*fUpwind_l[16]+0.3162277660168379*alpha[11]*fUpwind_l[15]+0.3535533905932737*alpha[2]*fUpwind_l[15]+0.282842712474619*fUpwind_l[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_l[10]+0.3535533905932737*alpha[4]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[6]; 
  Ghat_l[20] = 0.2258769757263128*alpha[13]*fUpwind_l[23]+0.3535533905932737*alpha[3]*fUpwind_l[23]+0.2258769757263128*alpha[7]*fUpwind_l[20]+0.3535533905932737*alpha[0]*fUpwind_l[20]+0.3162277660168379*alpha[5]*fUpwind_l[18]+0.3535533905932737*alpha[13]*fUpwind_l[14]+0.3162277660168379*alpha[1]*fUpwind_l[12]+0.2020305089104422*alpha[11]*fUpwind_l[11]+0.3162277660168379*alpha[2]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[2]*alpha[11]+0.3535533905932737*alpha[7]*fUpwind_l[8]+0.2828427124746191*alpha[4]*fUpwind_l[4]; 
  Ghat_l[21] = 0.2258769757263128*alpha[11]*fUpwind_l[24]+0.3535533905932737*alpha[2]*fUpwind_l[24]+0.2258769757263128*alpha[7]*fUpwind_l[21]+0.3535533905932737*alpha[0]*fUpwind_l[21]+0.3162277660168379*alpha[4]*fUpwind_l[19]+0.3535533905932737*alpha[11]*fUpwind_l[16]+0.3162277660168379*alpha[1]*fUpwind_l[15]+0.2020305089104422*alpha[13]*fUpwind_l[13]+0.3162277660168379*alpha[3]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[3]*alpha[13]+0.3535533905932737*alpha[7]*fUpwind_l[9]+0.2828427124746191*alpha[5]*fUpwind_l[5]; 
  Ghat_l[22] = 0.3535533905932737*alpha[7]*fUpwind_l[26]+0.3535533905932737*alpha[1]*fUpwind_l[25]+0.3162277660168379*alpha[11]*fUpwind_l[24]+0.3162277660168379*alpha[13]*fUpwind_l[23]+0.3535533905932737*alpha[0]*fUpwind_l[22]+0.3162277660168379*alpha[4]*fUpwind_l[19]+0.3162277660168379*alpha[5]*fUpwind_l[18]+0.3162277660168379*alpha[2]*fUpwind_l[16]+0.3162277660168379*alpha[3]*fUpwind_l[14]; 
  Ghat_l[23] = 0.2020305089104422*alpha[13]*fUpwind_l[26]+0.3162277660168379*alpha[3]*fUpwind_l[26]+0.2828427124746191*alpha[5]*fUpwind_l[25]+0.2258769757263128*alpha[7]*fUpwind_l[23]+0.3535533905932737*alpha[0]*fUpwind_l[23]+0.3162277660168379*alpha[13]*fUpwind_l[22]+0.2258769757263128*alpha[13]*fUpwind_l[20]+0.3535533905932737*alpha[3]*fUpwind_l[20]+0.3162277660168379*alpha[1]*fUpwind_l[18]+0.2020305089104422*alpha[11]*fUpwind_l[17]+0.3162277660168379*alpha[2]*fUpwind_l[17]+0.3535533905932737*alpha[7]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[8]*alpha[13]+0.3162277660168379*alpha[5]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[6]*alpha[11]+0.2828427124746191*alpha[4]*fUpwind_l[10]; 
  Ghat_l[24] = 0.2020305089104422*alpha[11]*fUpwind_l[26]+0.3162277660168379*alpha[2]*fUpwind_l[26]+0.2828427124746191*alpha[4]*fUpwind_l[25]+0.2258769757263128*alpha[7]*fUpwind_l[24]+0.3535533905932737*alpha[0]*fUpwind_l[24]+0.3162277660168379*alpha[11]*fUpwind_l[22]+0.2258769757263128*alpha[11]*fUpwind_l[21]+0.3535533905932737*alpha[2]*fUpwind_l[21]+0.3162277660168379*alpha[1]*fUpwind_l[19]+0.2020305089104422*alpha[13]*fUpwind_l[17]+0.3162277660168379*alpha[3]*fUpwind_l[17]+0.3535533905932737*alpha[7]*fUpwind_l[16]+0.3162277660168379*alpha[4]*fUpwind_l[15]+0.3162277660168379*fUpwind_l[6]*alpha[13]+0.3535533905932737*fUpwind_l[9]*alpha[11]+0.2828427124746191*alpha[5]*fUpwind_l[10]; 
  Ghat_l[25] = 0.3162277660168379*alpha[1]*fUpwind_l[26]+0.3162277660168379*alpha[7]*fUpwind_l[25]+0.3535533905932737*alpha[0]*fUpwind_l[25]+0.2828427124746191*alpha[4]*fUpwind_l[24]+0.2828427124746191*alpha[5]*fUpwind_l[23]+0.3535533905932737*alpha[1]*fUpwind_l[22]+0.282842712474619*alpha[11]*fUpwind_l[19]+0.3162277660168379*alpha[2]*fUpwind_l[19]+0.282842712474619*alpha[13]*fUpwind_l[18]+0.3162277660168379*alpha[3]*fUpwind_l[18]+0.3162277660168379*alpha[4]*fUpwind_l[16]+0.3162277660168379*alpha[5]*fUpwind_l[14]; 
  Ghat_l[26] = 0.2258769757263128*alpha[7]*fUpwind_l[26]+0.3535533905932737*alpha[0]*fUpwind_l[26]+0.3162277660168379*alpha[1]*fUpwind_l[25]+0.2020305089104422*alpha[11]*fUpwind_l[24]+0.3162277660168379*alpha[2]*fUpwind_l[24]+0.2020305089104422*alpha[13]*fUpwind_l[23]+0.3162277660168379*alpha[3]*fUpwind_l[23]+0.3535533905932737*alpha[7]*fUpwind_l[22]+0.2828427124746191*alpha[4]*fUpwind_l[19]+0.2828427124746191*alpha[5]*fUpwind_l[18]+0.3162277660168379*alpha[11]*fUpwind_l[16]+0.3162277660168379*alpha[13]*fUpwind_l[14]; 

  Ghat_r[0] = 0.3535533905932737*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[5]*fUpwind_r[5]+0.3535533905932737*alpha[4]*fUpwind_r[4]+0.3535533905932737*alpha[3]*fUpwind_r[3]+0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[1]*fUpwind_r[1]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3162277660168379*alpha[5]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[5]*alpha[13]+0.3162277660168379*alpha[4]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[4]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.3162277660168379*alpha[11]*fUpwind_r[20]+0.3535533905932737*alpha[13]*fUpwind_r[17]+0.3162277660168379*alpha[4]*fUpwind_r[12]+0.3535533905932737*alpha[7]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[7]*alpha[11]+0.3535533905932737*alpha[5]*fUpwind_r[10]+0.3162277660168379*alpha[2]*fUpwind_r[8]+0.3535533905932737*alpha[3]*fUpwind_r[6]+0.3535533905932737*alpha[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3162277660168379*alpha[13]*fUpwind_r[21]+0.3535533905932737*alpha[11]*fUpwind_r[17]+0.3162277660168379*alpha[5]*fUpwind_r[15]+0.3535533905932737*alpha[7]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[7]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind_r[10]+0.3162277660168379*alpha[3]*fUpwind_r[9]+0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*alpha[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.2828427124746191*alpha[4]*fUpwind_r[20]+0.3162277660168379*alpha[5]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[10]*alpha[13]+0.2828427124746191*alpha[11]*fUpwind_r[12]+0.3162277660168379*alpha[2]*fUpwind_r[12]+0.3162277660168379*alpha[1]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[1]*alpha[11]+0.3535533905932737*alpha[3]*fUpwind_r[10]+0.3162277660168379*alpha[4]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind_r[6]+0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.2828427124746191*alpha[5]*fUpwind_r[21]+0.3162277660168379*alpha[4]*fUpwind_r[17]+0.2828427124746191*alpha[13]*fUpwind_r[15]+0.3162277660168379*alpha[3]*fUpwind_r[15]+0.3162277660168379*alpha[1]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[1]*alpha[13]+0.3162277660168379*fUpwind_r[10]*alpha[11]+0.3535533905932737*alpha[2]*fUpwind_r[10]+0.3162277660168379*alpha[5]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_r[6]+0.3535533905932737*alpha[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alpha[3]; 
  Ghat_r[6] = 0.3162277660168379*alpha[13]*fUpwind_r[24]+0.3162277660168379*alpha[11]*fUpwind_r[23]+0.3162277660168379*alpha[5]*fUpwind_r[19]+0.3162277660168379*alpha[4]*fUpwind_r[18]+0.3535533905932737*alpha[7]*fUpwind_r[17]+0.3162277660168379*alpha[3]*fUpwind_r[16]+0.3162277660168379*alpha[2]*fUpwind_r[14]+0.3535533905932737*alpha[11]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[11]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*alpha[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alpha[3]; 
  Ghat_r[7] = 0.2258769757263128*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[3]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alpha[13]+0.2258769757263128*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[2]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2]*alpha[11]+0.2258769757263128*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[1]*fUpwind_r[1]; 
  Ghat_r[8] = 0.3535533905932737*alpha[13]*fUpwind_r[23]+0.3535533905932737*alpha[7]*fUpwind_r[20]+0.3535533905932737*alpha[5]*fUpwind_r[18]+0.3535533905932737*alpha[3]*fUpwind_r[14]+0.3535533905932737*alpha[1]*fUpwind_r[12]+0.3162277660168379*alpha[11]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha[11]*fUpwind_r[24]+0.3535533905932737*alpha[7]*fUpwind_r[21]+0.3535533905932737*alpha[4]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[16]+0.3535533905932737*alpha[1]*fUpwind_r[15]+0.3162277660168379*alpha[13]*fUpwind_r[13]+0.3535533905932737*alpha[0]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[3]*fUpwind_r[3]; 
  Ghat_r[10] = 0.2828427124746191*alpha[5]*fUpwind_r[24]+0.2828427124746191*alpha[4]*fUpwind_r[23]+0.282842712474619*alpha[13]*fUpwind_r[19]+0.3162277660168379*alpha[3]*fUpwind_r[19]+0.282842712474619*alpha[11]*fUpwind_r[18]+0.3162277660168379*alpha[2]*fUpwind_r[18]+0.3162277660168379*alpha[1]*fUpwind_r[17]+0.3162277660168379*alpha[5]*fUpwind_r[16]+0.3162277660168379*alpha[4]*fUpwind_r[14]+0.3162277660168379*alpha[4]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[4]*alpha[13]+0.3162277660168379*alpha[5]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[5]*alpha[11]+0.3162277660168379*alpha[7]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[10]+0.3535533905932737*alpha[1]*fUpwind_r[6]+0.3535533905932737*alpha[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.2020305089104422*alpha[11]*fUpwind_r[20]+0.3162277660168379*alpha[2]*fUpwind_r[20]+0.2258769757263128*alpha[13]*fUpwind_r[17]+0.3535533905932737*alpha[3]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[6]*alpha[13]+0.2828427124746191*alpha[4]*fUpwind_r[12]+0.2258769757263128*alpha[7]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[8]*alpha[11]+0.2258769757263128*fUpwind_r[7]*alpha[11]+0.3535533905932737*fUpwind_r[0]*alpha[11]+0.3162277660168379*alpha[5]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[1]*alpha[4]; 
  Ghat_r[12] = 0.3162277660168379*alpha[5]*fUpwind_r[23]+0.3162277660168379*alpha[1]*fUpwind_r[20]+0.3162277660168379*alpha[13]*fUpwind_r[18]+0.3535533905932737*alpha[3]*fUpwind_r[18]+0.3535533905932737*alpha[5]*fUpwind_r[14]+0.3162277660168379*alpha[7]*fUpwind_r[12]+0.3535533905932737*alpha[0]*fUpwind_r[12]+0.2828427124746191*alpha[4]*fUpwind_r[11]+0.2828427124746191*fUpwind_r[4]*alpha[11]+0.3535533905932737*alpha[1]*fUpwind_r[8]+0.3162277660168379*alpha[2]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[2]*alpha[4]; 
  Ghat_r[13] = 0.2020305089104422*alpha[13]*fUpwind_r[21]+0.3162277660168379*alpha[3]*fUpwind_r[21]+0.2258769757263128*alpha[11]*fUpwind_r[17]+0.3535533905932737*alpha[2]*fUpwind_r[17]+0.2828427124746191*alpha[5]*fUpwind_r[15]+0.2258769757263128*alpha[7]*fUpwind_r[13]+0.3535533905932737*alpha[0]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[9]*alpha[13]+0.2258769757263128*fUpwind_r[7]*alpha[13]+0.3535533905932737*fUpwind_r[0]*alpha[13]+0.3535533905932737*fUpwind_r[6]*alpha[11]+0.3162277660168379*alpha[4]*fUpwind_r[10]+0.3535533905932737*alpha[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[1]*alpha[5]; 
  Ghat_r[14] = 0.3162277660168379*alpha[13]*fUpwind_r[26]+0.3162277660168379*alpha[5]*fUpwind_r[25]+0.3535533905932737*alpha[7]*fUpwind_r[23]+0.3162277660168379*alpha[3]*fUpwind_r[22]+0.3535533905932737*alpha[13]*fUpwind_r[20]+0.3535533905932737*alpha[1]*fUpwind_r[18]+0.3162277660168379*alpha[11]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[14]+0.3535533905932737*alpha[5]*fUpwind_r[12]+0.3162277660168379*alpha[4]*fUpwind_r[10]+0.3535533905932737*alpha[3]*fUpwind_r[8]+0.3162277660168379*alpha[2]*fUpwind_r[6]; 
  Ghat_r[15] = 0.3162277660168379*alpha[4]*fUpwind_r[24]+0.3162277660168379*alpha[1]*fUpwind_r[21]+0.3162277660168379*alpha[11]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[19]+0.3535533905932737*alpha[4]*fUpwind_r[16]+0.3162277660168379*alpha[7]*fUpwind_r[15]+0.3535533905932737*alpha[0]*fUpwind_r[15]+0.2828427124746191*alpha[5]*fUpwind_r[13]+0.2828427124746191*fUpwind_r[5]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind_r[9]+0.3162277660168379*alpha[3]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[3]*alpha[5]; 
  Ghat_r[16] = 0.3162277660168379*alpha[11]*fUpwind_r[26]+0.3162277660168379*alpha[4]*fUpwind_r[25]+0.3535533905932737*alpha[7]*fUpwind_r[24]+0.3162277660168379*alpha[2]*fUpwind_r[22]+0.3535533905932737*alpha[11]*fUpwind_r[21]+0.3535533905932737*alpha[1]*fUpwind_r[19]+0.3162277660168379*alpha[13]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[16]+0.3535533905932737*alpha[4]*fUpwind_r[15]+0.3162277660168379*alpha[5]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[9]+0.3162277660168379*alpha[3]*fUpwind_r[6]; 
  Ghat_r[17] = 0.2020305089104422*alpha[13]*fUpwind_r[24]+0.3162277660168379*alpha[3]*fUpwind_r[24]+0.2020305089104422*alpha[11]*fUpwind_r[23]+0.3162277660168379*alpha[2]*fUpwind_r[23]+0.2828427124746191*alpha[5]*fUpwind_r[19]+0.2828427124746191*alpha[4]*fUpwind_r[18]+0.2258769757263128*alpha[7]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[17]+0.3162277660168379*alpha[13]*fUpwind_r[16]+0.3162277660168379*alpha[11]*fUpwind_r[14]+0.2258769757263128*alpha[11]*fUpwind_r[13]+0.3535533905932737*alpha[2]*fUpwind_r[13]+0.2258769757263128*fUpwind_r[11]*alpha[13]+0.3535533905932737*fUpwind_r[2]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[3]*alpha[11]+0.3162277660168379*alpha[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[4]*alpha[5]; 
  Ghat_r[18] = 0.2828427124746191*alpha[5]*fUpwind_r[26]+0.282842712474619*alpha[13]*fUpwind_r[25]+0.3162277660168379*alpha[3]*fUpwind_r[25]+0.3162277660168379*alpha[1]*fUpwind_r[23]+0.3162277660168379*alpha[5]*fUpwind_r[22]+0.3162277660168379*alpha[5]*fUpwind_r[20]+0.3162277660168379*alpha[7]*fUpwind_r[18]+0.3535533905932737*alpha[0]*fUpwind_r[18]+0.2828427124746191*alpha[4]*fUpwind_r[17]+0.3535533905932737*alpha[1]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[12]*alpha[13]+0.3535533905932737*alpha[3]*fUpwind_r[12]+0.282842712474619*fUpwind_r[10]*alpha[11]+0.3162277660168379*alpha[2]*fUpwind_r[10]+0.3535533905932737*alpha[5]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[6]; 
  Ghat_r[19] = 0.2828427124746191*alpha[4]*fUpwind_r[26]+0.282842712474619*alpha[11]*fUpwind_r[25]+0.3162277660168379*alpha[2]*fUpwind_r[25]+0.3162277660168379*alpha[1]*fUpwind_r[24]+0.3162277660168379*alpha[4]*fUpwind_r[22]+0.3162277660168379*alpha[4]*fUpwind_r[21]+0.3162277660168379*alpha[7]*fUpwind_r[19]+0.3535533905932737*alpha[0]*fUpwind_r[19]+0.2828427124746191*alpha[5]*fUpwind_r[17]+0.3535533905932737*alpha[1]*fUpwind_r[16]+0.3162277660168379*alpha[11]*fUpwind_r[15]+0.3535533905932737*alpha[2]*fUpwind_r[15]+0.282842712474619*fUpwind_r[10]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind_r[10]+0.3535533905932737*alpha[4]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[6]; 
  Ghat_r[20] = 0.2258769757263128*alpha[13]*fUpwind_r[23]+0.3535533905932737*alpha[3]*fUpwind_r[23]+0.2258769757263128*alpha[7]*fUpwind_r[20]+0.3535533905932737*alpha[0]*fUpwind_r[20]+0.3162277660168379*alpha[5]*fUpwind_r[18]+0.3535533905932737*alpha[13]*fUpwind_r[14]+0.3162277660168379*alpha[1]*fUpwind_r[12]+0.2020305089104422*alpha[11]*fUpwind_r[11]+0.3162277660168379*alpha[2]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[2]*alpha[11]+0.3535533905932737*alpha[7]*fUpwind_r[8]+0.2828427124746191*alpha[4]*fUpwind_r[4]; 
  Ghat_r[21] = 0.2258769757263128*alpha[11]*fUpwind_r[24]+0.3535533905932737*alpha[2]*fUpwind_r[24]+0.2258769757263128*alpha[7]*fUpwind_r[21]+0.3535533905932737*alpha[0]*fUpwind_r[21]+0.3162277660168379*alpha[4]*fUpwind_r[19]+0.3535533905932737*alpha[11]*fUpwind_r[16]+0.3162277660168379*alpha[1]*fUpwind_r[15]+0.2020305089104422*alpha[13]*fUpwind_r[13]+0.3162277660168379*alpha[3]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[3]*alpha[13]+0.3535533905932737*alpha[7]*fUpwind_r[9]+0.2828427124746191*alpha[5]*fUpwind_r[5]; 
  Ghat_r[22] = 0.3535533905932737*alpha[7]*fUpwind_r[26]+0.3535533905932737*alpha[1]*fUpwind_r[25]+0.3162277660168379*alpha[11]*fUpwind_r[24]+0.3162277660168379*alpha[13]*fUpwind_r[23]+0.3535533905932737*alpha[0]*fUpwind_r[22]+0.3162277660168379*alpha[4]*fUpwind_r[19]+0.3162277660168379*alpha[5]*fUpwind_r[18]+0.3162277660168379*alpha[2]*fUpwind_r[16]+0.3162277660168379*alpha[3]*fUpwind_r[14]; 
  Ghat_r[23] = 0.2020305089104422*alpha[13]*fUpwind_r[26]+0.3162277660168379*alpha[3]*fUpwind_r[26]+0.2828427124746191*alpha[5]*fUpwind_r[25]+0.2258769757263128*alpha[7]*fUpwind_r[23]+0.3535533905932737*alpha[0]*fUpwind_r[23]+0.3162277660168379*alpha[13]*fUpwind_r[22]+0.2258769757263128*alpha[13]*fUpwind_r[20]+0.3535533905932737*alpha[3]*fUpwind_r[20]+0.3162277660168379*alpha[1]*fUpwind_r[18]+0.2020305089104422*alpha[11]*fUpwind_r[17]+0.3162277660168379*alpha[2]*fUpwind_r[17]+0.3535533905932737*alpha[7]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[8]*alpha[13]+0.3162277660168379*alpha[5]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[6]*alpha[11]+0.2828427124746191*alpha[4]*fUpwind_r[10]; 
  Ghat_r[24] = 0.2020305089104422*alpha[11]*fUpwind_r[26]+0.3162277660168379*alpha[2]*fUpwind_r[26]+0.2828427124746191*alpha[4]*fUpwind_r[25]+0.2258769757263128*alpha[7]*fUpwind_r[24]+0.3535533905932737*alpha[0]*fUpwind_r[24]+0.3162277660168379*alpha[11]*fUpwind_r[22]+0.2258769757263128*alpha[11]*fUpwind_r[21]+0.3535533905932737*alpha[2]*fUpwind_r[21]+0.3162277660168379*alpha[1]*fUpwind_r[19]+0.2020305089104422*alpha[13]*fUpwind_r[17]+0.3162277660168379*alpha[3]*fUpwind_r[17]+0.3535533905932737*alpha[7]*fUpwind_r[16]+0.3162277660168379*alpha[4]*fUpwind_r[15]+0.3162277660168379*fUpwind_r[6]*alpha[13]+0.3535533905932737*fUpwind_r[9]*alpha[11]+0.2828427124746191*alpha[5]*fUpwind_r[10]; 
  Ghat_r[25] = 0.3162277660168379*alpha[1]*fUpwind_r[26]+0.3162277660168379*alpha[7]*fUpwind_r[25]+0.3535533905932737*alpha[0]*fUpwind_r[25]+0.2828427124746191*alpha[4]*fUpwind_r[24]+0.2828427124746191*alpha[5]*fUpwind_r[23]+0.3535533905932737*alpha[1]*fUpwind_r[22]+0.282842712474619*alpha[11]*fUpwind_r[19]+0.3162277660168379*alpha[2]*fUpwind_r[19]+0.282842712474619*alpha[13]*fUpwind_r[18]+0.3162277660168379*alpha[3]*fUpwind_r[18]+0.3162277660168379*alpha[4]*fUpwind_r[16]+0.3162277660168379*alpha[5]*fUpwind_r[14]; 
  Ghat_r[26] = 0.2258769757263128*alpha[7]*fUpwind_r[26]+0.3535533905932737*alpha[0]*fUpwind_r[26]+0.3162277660168379*alpha[1]*fUpwind_r[25]+0.2020305089104422*alpha[11]*fUpwind_r[24]+0.3162277660168379*alpha[2]*fUpwind_r[24]+0.2020305089104422*alpha[13]*fUpwind_r[23]+0.3162277660168379*alpha[3]*fUpwind_r[23]+0.3535533905932737*alpha[7]*fUpwind_r[22]+0.2828427124746191*alpha[4]*fUpwind_r[19]+0.2828427124746191*alpha[5]*fUpwind_r[18]+0.3162277660168379*alpha[11]*fUpwind_r[16]+0.3162277660168379*alpha[13]*fUpwind_r[14]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[12] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[17] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[20] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[22] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv10; 
  out[23] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[24] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[26] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv10; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[29] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[30] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv10; 
  out[31] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[32] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[33] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv10; 
  out[34] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[35] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[36] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv10; 
  out[37] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv10; 
  out[38] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv10; 
  out[39] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv10; 
  out[40] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[41] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 
  out[42] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv10; 
  out[43] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv10; 
  out[44] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv10; 
  out[45] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dv10; 
  out[46] += (1.58113883008419*Ghat_l[8]-1.58113883008419*Ghat_r[8])*dv10; 
  out[47] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dv10; 
  out[48] += (1.58113883008419*Ghat_l[9]-1.58113883008419*Ghat_r[9])*dv10; 
  out[49] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dv10; 
  out[50] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv10; 
  out[51] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dv10; 
  out[52] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv10; 
  out[53] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv10; 
  out[54] += (1.58113883008419*Ghat_l[11]-1.58113883008419*Ghat_r[11])*dv10; 
  out[55] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dv10; 
  out[56] += (1.58113883008419*Ghat_l[12]-1.58113883008419*Ghat_r[12])*dv10; 
  out[57] += (1.58113883008419*Ghat_l[13]-1.58113883008419*Ghat_r[13])*dv10; 
  out[58] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dv10; 
  out[59] += (1.58113883008419*Ghat_l[14]-1.58113883008419*Ghat_r[14])*dv10; 
  out[60] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dv10; 
  out[61] += (1.58113883008419*Ghat_l[15]-1.58113883008419*Ghat_r[15])*dv10; 
  out[62] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dv10; 
  out[63] += (1.58113883008419*Ghat_l[16]-1.58113883008419*Ghat_r[16])*dv10; 
  out[64] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dv10; 
  out[65] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dv10; 
  out[66] += (1.58113883008419*Ghat_l[17]-1.58113883008419*Ghat_r[17])*dv10; 
  out[67] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dv10; 
  out[68] += (1.58113883008419*Ghat_l[18]-1.58113883008419*Ghat_r[18])*dv10; 
  out[69] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dv10; 
  out[70] += (1.58113883008419*Ghat_l[19]-1.58113883008419*Ghat_r[19])*dv10; 
  out[71] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dv10; 
  out[72] += (1.58113883008419*Ghat_l[20]-1.58113883008419*Ghat_r[20])*dv10; 
  out[73] += (1.58113883008419*Ghat_l[21]-1.58113883008419*Ghat_r[21])*dv10; 
  out[74] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dv10; 
  out[75] += (1.58113883008419*Ghat_l[22]-1.58113883008419*Ghat_r[22])*dv10; 
  out[76] += (1.58113883008419*Ghat_l[23]-1.58113883008419*Ghat_r[23])*dv10; 
  out[77] += (1.58113883008419*Ghat_l[24]-1.58113883008419*Ghat_r[24])*dv10; 
  out[78] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dv10; 
  out[79] += (1.58113883008419*Ghat_l[25]-1.58113883008419*Ghat_r[25])*dv10; 
  out[80] += (1.58113883008419*Ghat_l[26]-1.58113883008419*Ghat_r[26])*dv10; 

  return 0.;

} 
