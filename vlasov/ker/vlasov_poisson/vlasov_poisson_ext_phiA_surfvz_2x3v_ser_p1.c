#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_phiA_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 

  const double dv12 = 2/dxv[4]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 

  double alpha[32] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ax = &pots[4]; 
  const double *Ay = &pots[8]; 
  const double *Az = &pots[12]; 

  alpha[0] = -(3.4641016151377544*Az[2]*dx11*wv2)-3.4641016151377544*Az[1]*dx10*wv1; 
  alpha[1] = -(3.4641016151377544*Az[3]*dx11*wv2); 
  alpha[2] = -(3.4641016151377544*Az[3]*dx10*wv1); 
  alpha[3] = -(1.0*Az[1]*dv1*dx10); 
  alpha[4] = -(1.0*Az[2]*dv2*dx11); 
  alpha[7] = -(1.0*Az[3]*dv1*dx10); 
  alpha[8] = -(1.0*Az[3]*dv2*dx11); 

  double fUpwindQuad_l[36] = {0.0};
  double fUpwindQuad_r[36] = {0.0};
  double fUpwind_l[32] = {0.0};;
  double fUpwind_r[32] = {0.0};
  double Ghat_l[32] = {0.0}; 
  double Ghat_r[32] = {0.0}; 

  if (0.33541019662496846*(alpha[8]+alpha[7])-0.33541019662496846*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[0] = hyb_2x3v_p1_surfx5_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_2x3v_p1_surfx5_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x3v_p1_surfx5_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_2x3v_p1_surfx5_eval_quad_node_0_l(fr); 
  } 
  if (0.33541019662496846*alpha[7]-0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[1] = hyb_2x3v_p1_surfx5_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_2x3v_p1_surfx5_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = hyb_2x3v_p1_surfx5_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_2x3v_p1_surfx5_eval_quad_node_1_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[8])+0.33541019662496846*(alpha[7]+alpha[4])-0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[2] = hyb_2x3v_p1_surfx5_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_2x3v_p1_surfx5_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = hyb_2x3v_p1_surfx5_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_2x3v_p1_surfx5_eval_quad_node_2_l(fr); 
  } 
  if (0.33541019662496846*alpha[8]-0.33541019662496846*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[3] = hyb_2x3v_p1_surfx5_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_2x3v_p1_surfx5_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_2x3v_p1_surfx5_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_2x3v_p1_surfx5_eval_quad_node_3_l(fr); 
  } 
  if (0.25*alpha[0]-0.25*(alpha[2]+alpha[1]) > 0) { 
    fUpwindQuad_l[4] = hyb_2x3v_p1_surfx5_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_2x3v_p1_surfx5_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = hyb_2x3v_p1_surfx5_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_2x3v_p1_surfx5_eval_quad_node_4_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[8])+0.33541019662496846*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[5] = hyb_2x3v_p1_surfx5_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_2x3v_p1_surfx5_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = hyb_2x3v_p1_surfx5_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_2x3v_p1_surfx5_eval_quad_node_5_l(fr); 
  } 
  if (0.33541019662496846*alpha[8]-0.33541019662496846*(alpha[7]+alpha[4])+0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[6] = hyb_2x3v_p1_surfx5_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = hyb_2x3v_p1_surfx5_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = hyb_2x3v_p1_surfx5_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = hyb_2x3v_p1_surfx5_eval_quad_node_6_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[7])+0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[7] = hyb_2x3v_p1_surfx5_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = hyb_2x3v_p1_surfx5_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = hyb_2x3v_p1_surfx5_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = hyb_2x3v_p1_surfx5_eval_quad_node_7_l(fr); 
  } 
  if (-(0.33541019662496846*(alpha[8]+alpha[7]))+0.33541019662496846*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[8] = hyb_2x3v_p1_surfx5_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = hyb_2x3v_p1_surfx5_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = hyb_2x3v_p1_surfx5_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = hyb_2x3v_p1_surfx5_eval_quad_node_8_l(fr); 
  } 
  if (0.33541019662496846*alpha[8]-0.33541019662496846*(alpha[7]+alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[9] = hyb_2x3v_p1_surfx5_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = hyb_2x3v_p1_surfx5_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = hyb_2x3v_p1_surfx5_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = hyb_2x3v_p1_surfx5_eval_quad_node_9_l(fr); 
  } 
  if (-(0.33541019662496846*(alpha[7]+alpha[3]))+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[10] = hyb_2x3v_p1_surfx5_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = hyb_2x3v_p1_surfx5_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = hyb_2x3v_p1_surfx5_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = hyb_2x3v_p1_surfx5_eval_quad_node_10_l(fr); 
  } 
  if (-(0.33541019662496846*(alpha[8]+alpha[7]))+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[11] = hyb_2x3v_p1_surfx5_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = hyb_2x3v_p1_surfx5_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = hyb_2x3v_p1_surfx5_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = hyb_2x3v_p1_surfx5_eval_quad_node_11_l(fr); 
  } 
  if (0.33541019662496846*alpha[8]-0.33541019662496846*alpha[4]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[12] = hyb_2x3v_p1_surfx5_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = hyb_2x3v_p1_surfx5_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = hyb_2x3v_p1_surfx5_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = hyb_2x3v_p1_surfx5_eval_quad_node_12_l(fr); 
  } 
  if (0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[13] = hyb_2x3v_p1_surfx5_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = hyb_2x3v_p1_surfx5_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = hyb_2x3v_p1_surfx5_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = hyb_2x3v_p1_surfx5_eval_quad_node_13_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[8])+0.33541019662496846*alpha[4]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[14] = hyb_2x3v_p1_surfx5_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = hyb_2x3v_p1_surfx5_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = hyb_2x3v_p1_surfx5_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = hyb_2x3v_p1_surfx5_eval_quad_node_14_l(fr); 
  } 
  if (0.33541019662496846*(alpha[8]+alpha[7])-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[15] = hyb_2x3v_p1_surfx5_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = hyb_2x3v_p1_surfx5_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = hyb_2x3v_p1_surfx5_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = hyb_2x3v_p1_surfx5_eval_quad_node_15_l(fr); 
  } 
  if (0.33541019662496846*(alpha[7]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[16] = hyb_2x3v_p1_surfx5_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = hyb_2x3v_p1_surfx5_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = hyb_2x3v_p1_surfx5_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = hyb_2x3v_p1_surfx5_eval_quad_node_16_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[8])+0.33541019662496846*(alpha[7]+alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[17] = hyb_2x3v_p1_surfx5_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = hyb_2x3v_p1_surfx5_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = hyb_2x3v_p1_surfx5_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = hyb_2x3v_p1_surfx5_eval_quad_node_17_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[8])+0.33541019662496846*alpha[7]-0.33541019662496846*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[18] = hyb_2x3v_p1_surfx5_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = hyb_2x3v_p1_surfx5_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = hyb_2x3v_p1_surfx5_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = hyb_2x3v_p1_surfx5_eval_quad_node_18_l(fr); 
  } 
  if (0.33541019662496846*alpha[7]-0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[19] = hyb_2x3v_p1_surfx5_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = hyb_2x3v_p1_surfx5_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = hyb_2x3v_p1_surfx5_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = hyb_2x3v_p1_surfx5_eval_quad_node_19_l(fr); 
  } 
  if (0.33541019662496846*(alpha[8]+alpha[7]+alpha[4])-0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[20] = hyb_2x3v_p1_surfx5_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = hyb_2x3v_p1_surfx5_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = hyb_2x3v_p1_surfx5_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = hyb_2x3v_p1_surfx5_eval_quad_node_20_l(fr); 
  } 
  if (-(0.33541019662496846*(alpha[8]+alpha[4]))-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[21] = hyb_2x3v_p1_surfx5_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = hyb_2x3v_p1_surfx5_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = hyb_2x3v_p1_surfx5_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = hyb_2x3v_p1_surfx5_eval_quad_node_21_l(fr); 
  } 
  if (0.25*(alpha[1]+alpha[0])-0.25*alpha[2] > 0) { 
    fUpwindQuad_l[22] = hyb_2x3v_p1_surfx5_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = hyb_2x3v_p1_surfx5_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = hyb_2x3v_p1_surfx5_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = hyb_2x3v_p1_surfx5_eval_quad_node_22_l(fr); 
  } 
  if (0.33541019662496846*(alpha[8]+alpha[4])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[23] = hyb_2x3v_p1_surfx5_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = hyb_2x3v_p1_surfx5_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = hyb_2x3v_p1_surfx5_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = hyb_2x3v_p1_surfx5_eval_quad_node_23_l(fr); 
  } 
  if (-(0.33541019662496846*(alpha[8]+alpha[7]+alpha[4]))+0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[24] = hyb_2x3v_p1_surfx5_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = hyb_2x3v_p1_surfx5_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = hyb_2x3v_p1_surfx5_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = hyb_2x3v_p1_surfx5_eval_quad_node_24_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[7])+0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[25] = hyb_2x3v_p1_surfx5_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = hyb_2x3v_p1_surfx5_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = hyb_2x3v_p1_surfx5_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = hyb_2x3v_p1_surfx5_eval_quad_node_25_l(fr); 
  } 
  if (0.33541019662496846*alpha[8]-0.33541019662496846*alpha[7]+0.33541019662496846*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[26] = hyb_2x3v_p1_surfx5_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = hyb_2x3v_p1_surfx5_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = hyb_2x3v_p1_surfx5_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = hyb_2x3v_p1_surfx5_eval_quad_node_26_l(fr); 
  } 
  if (0.25*(alpha[2]+alpha[1]+alpha[0])-0.33541019662496846*(alpha[8]+alpha[7]+alpha[4]+alpha[3]) > 0) { 
    fUpwindQuad_l[27] = hyb_2x3v_p1_surfx5_eval_quad_node_27_r(fl); 
    fUpwindQuad_r[27] = hyb_2x3v_p1_surfx5_eval_quad_node_27_r(fc); 
  } else { 
    fUpwindQuad_l[27] = hyb_2x3v_p1_surfx5_eval_quad_node_27_l(fc); 
    fUpwindQuad_r[27] = hyb_2x3v_p1_surfx5_eval_quad_node_27_l(fr); 
  } 
  if (0.25*(alpha[2]+alpha[1]+alpha[0])-0.33541019662496846*(alpha[7]+alpha[3]) > 0) { 
    fUpwindQuad_l[28] = hyb_2x3v_p1_surfx5_eval_quad_node_28_r(fl); 
    fUpwindQuad_r[28] = hyb_2x3v_p1_surfx5_eval_quad_node_28_r(fc); 
  } else { 
    fUpwindQuad_l[28] = hyb_2x3v_p1_surfx5_eval_quad_node_28_l(fc); 
    fUpwindQuad_r[28] = hyb_2x3v_p1_surfx5_eval_quad_node_28_l(fr); 
  } 
  if (0.33541019662496846*alpha[8]-0.33541019662496846*alpha[7]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[29] = hyb_2x3v_p1_surfx5_eval_quad_node_29_r(fl); 
    fUpwindQuad_r[29] = hyb_2x3v_p1_surfx5_eval_quad_node_29_r(fc); 
  } else { 
    fUpwindQuad_l[29] = hyb_2x3v_p1_surfx5_eval_quad_node_29_l(fc); 
    fUpwindQuad_r[29] = hyb_2x3v_p1_surfx5_eval_quad_node_29_l(fr); 
  } 
  if (0.25*(alpha[2]+alpha[1]+alpha[0])-0.33541019662496846*(alpha[8]+alpha[4]) > 0) { 
    fUpwindQuad_l[30] = hyb_2x3v_p1_surfx5_eval_quad_node_30_r(fl); 
    fUpwindQuad_r[30] = hyb_2x3v_p1_surfx5_eval_quad_node_30_r(fc); 
  } else { 
    fUpwindQuad_l[30] = hyb_2x3v_p1_surfx5_eval_quad_node_30_l(fc); 
    fUpwindQuad_r[30] = hyb_2x3v_p1_surfx5_eval_quad_node_30_l(fr); 
  } 
  if (0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[31] = hyb_2x3v_p1_surfx5_eval_quad_node_31_r(fl); 
    fUpwindQuad_r[31] = hyb_2x3v_p1_surfx5_eval_quad_node_31_r(fc); 
  } else { 
    fUpwindQuad_l[31] = hyb_2x3v_p1_surfx5_eval_quad_node_31_l(fc); 
    fUpwindQuad_r[31] = hyb_2x3v_p1_surfx5_eval_quad_node_31_l(fr); 
  } 
  if (0.33541019662496846*(alpha[8]+alpha[4])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[32] = hyb_2x3v_p1_surfx5_eval_quad_node_32_r(fl); 
    fUpwindQuad_r[32] = hyb_2x3v_p1_surfx5_eval_quad_node_32_r(fc); 
  } else { 
    fUpwindQuad_l[32] = hyb_2x3v_p1_surfx5_eval_quad_node_32_l(fc); 
    fUpwindQuad_r[32] = hyb_2x3v_p1_surfx5_eval_quad_node_32_l(fr); 
  } 
  if (-(0.33541019662496846*alpha[8])+0.33541019662496846*alpha[7]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[33] = hyb_2x3v_p1_surfx5_eval_quad_node_33_r(fl); 
    fUpwindQuad_r[33] = hyb_2x3v_p1_surfx5_eval_quad_node_33_r(fc); 
  } else { 
    fUpwindQuad_l[33] = hyb_2x3v_p1_surfx5_eval_quad_node_33_l(fc); 
    fUpwindQuad_r[33] = hyb_2x3v_p1_surfx5_eval_quad_node_33_l(fr); 
  } 
  if (0.33541019662496846*(alpha[7]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[34] = hyb_2x3v_p1_surfx5_eval_quad_node_34_r(fl); 
    fUpwindQuad_r[34] = hyb_2x3v_p1_surfx5_eval_quad_node_34_r(fc); 
  } else { 
    fUpwindQuad_l[34] = hyb_2x3v_p1_surfx5_eval_quad_node_34_l(fc); 
    fUpwindQuad_r[34] = hyb_2x3v_p1_surfx5_eval_quad_node_34_l(fr); 
  } 
  if (0.33541019662496846*(alpha[8]+alpha[7]+alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[35] = hyb_2x3v_p1_surfx5_eval_quad_node_35_r(fl); 
    fUpwindQuad_r[35] = hyb_2x3v_p1_surfx5_eval_quad_node_35_r(fc); 
  } else { 
    fUpwindQuad_l[35] = hyb_2x3v_p1_surfx5_eval_quad_node_35_l(fc); 
    fUpwindQuad_r[35] = hyb_2x3v_p1_surfx5_eval_quad_node_35_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.25*alpha[8]*fUpwind_l[8]+0.25*alpha[7]*fUpwind_l[7]+0.25*alpha[4]*fUpwind_l[4]+0.25*alpha[3]*fUpwind_l[3]+0.25*alpha[2]*fUpwind_l[2]+0.25*alpha[1]*fUpwind_l[1]+0.25*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.25*alpha[7]*fUpwind_l[11]+0.25*alpha[4]*fUpwind_l[8]+0.25*fUpwind_l[4]*alpha[8]+0.25*alpha[3]*fUpwind_l[6]+0.25*alpha[2]*fUpwind_l[5]+0.25*alpha[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.25*alpha[8]*fUpwind_l[12]+0.25*alpha[4]*fUpwind_l[9]+0.25*alpha[3]*fUpwind_l[7]+0.25*fUpwind_l[3]*alpha[7]+0.25*alpha[1]*fUpwind_l[5]+0.25*alpha[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.22360679774997902*alpha[7]*fUpwind_l[18]+0.22360679774997896*alpha[3]*fUpwind_l[16]+0.25*alpha[8]*fUpwind_l[13]+0.25*alpha[4]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[7]+0.25*fUpwind_l[2]*alpha[7]+0.25*alpha[1]*fUpwind_l[6]+0.25*alpha[0]*fUpwind_l[3]+0.25*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.22360679774997902*alpha[8]*fUpwind_l[25]+0.22360679774997896*alpha[4]*fUpwind_l[24]+0.25*alpha[7]*fUpwind_l[14]+0.25*alpha[3]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[9]+0.25*alpha[1]*fUpwind_l[8]+0.25*fUpwind_l[1]*alpha[8]+0.25*alpha[0]*fUpwind_l[4]+0.25*fUpwind_l[0]*alpha[4]; 
  Ghat_l[5] = 0.25*alpha[4]*fUpwind_l[12]+0.25*alpha[3]*fUpwind_l[11]+0.25*alpha[8]*fUpwind_l[9]+0.25*fUpwind_l[6]*alpha[7]+0.25*alpha[0]*fUpwind_l[5]+0.25*alpha[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alpha[2]; 
  Ghat_l[6] = 0.22360679774997896*alpha[7]*fUpwind_l[20]+0.22360679774997902*alpha[3]*fUpwind_l[17]+0.25*alpha[4]*fUpwind_l[13]+0.25*alpha[2]*fUpwind_l[11]+0.25*alpha[8]*fUpwind_l[10]+0.25*fUpwind_l[5]*alpha[7]+0.25*alpha[0]*fUpwind_l[6]+0.25*alpha[1]*fUpwind_l[3]+0.25*fUpwind_l[1]*alpha[3]; 
  Ghat_l[7] = 0.22360679774997902*alpha[3]*fUpwind_l[18]+0.22360679774997896*alpha[7]*fUpwind_l[16]+0.25*alpha[8]*fUpwind_l[15]+0.25*alpha[4]*fUpwind_l[14]+0.25*alpha[1]*fUpwind_l[11]+0.25*alpha[0]*fUpwind_l[7]+0.25*fUpwind_l[0]*alpha[7]+0.25*alpha[2]*fUpwind_l[3]+0.25*fUpwind_l[2]*alpha[3]; 
  Ghat_l[8] = 0.22360679774997902*alpha[4]*fUpwind_l[25]+0.22360679774997896*alpha[8]*fUpwind_l[24]+0.25*alpha[7]*fUpwind_l[15]+0.25*alpha[3]*fUpwind_l[13]+0.25*alpha[2]*fUpwind_l[12]+0.25*alpha[0]*fUpwind_l[8]+0.25*fUpwind_l[0]*alpha[8]+0.25*alpha[1]*fUpwind_l[4]+0.25*fUpwind_l[1]*alpha[4]; 
  Ghat_l[9] = 0.22360679774997896*alpha[8]*fUpwind_l[28]+0.22360679774997902*alpha[4]*fUpwind_l[26]+0.25*alpha[3]*fUpwind_l[14]+0.25*alpha[1]*fUpwind_l[12]+0.25*alpha[7]*fUpwind_l[10]+0.25*alpha[0]*fUpwind_l[9]+0.25*fUpwind_l[5]*alpha[8]+0.25*alpha[2]*fUpwind_l[4]+0.25*fUpwind_l[2]*alpha[4]; 
  Ghat_l[10] = 0.22360679774997896*alpha[8]*fUpwind_l[29]+0.22360679774997902*alpha[4]*fUpwind_l[27]+0.22360679774997896*alpha[7]*fUpwind_l[22]+0.22360679774997902*alpha[3]*fUpwind_l[19]+0.25*alpha[2]*fUpwind_l[14]+0.25*alpha[1]*fUpwind_l[13]+0.25*alpha[0]*fUpwind_l[10]+0.25*alpha[7]*fUpwind_l[9]+0.25*fUpwind_l[6]*alpha[8]+0.25*alpha[3]*fUpwind_l[4]+0.25*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.22360679774997896*alpha[3]*fUpwind_l[20]+0.22360679774997902*alpha[7]*fUpwind_l[17]+0.25*alpha[4]*fUpwind_l[15]+0.25*alpha[8]*fUpwind_l[14]+0.25*alpha[0]*fUpwind_l[11]+0.25*alpha[1]*fUpwind_l[7]+0.25*fUpwind_l[1]*alpha[7]+0.25*alpha[2]*fUpwind_l[6]+0.25*alpha[3]*fUpwind_l[5]; 
  Ghat_l[12] = 0.22360679774997896*alpha[4]*fUpwind_l[28]+0.22360679774997902*alpha[8]*fUpwind_l[26]+0.25*alpha[3]*fUpwind_l[15]+0.25*alpha[7]*fUpwind_l[13]+0.25*alpha[0]*fUpwind_l[12]+0.25*alpha[1]*fUpwind_l[9]+0.25*alpha[2]*fUpwind_l[8]+0.25*fUpwind_l[2]*alpha[8]+0.25*alpha[4]*fUpwind_l[5]; 
  Ghat_l[13] = 0.22360679774997896*alpha[4]*fUpwind_l[29]+0.22360679774997902*alpha[8]*fUpwind_l[27]+0.22360679774997902*alpha[7]*fUpwind_l[23]+0.22360679774997896*alpha[3]*fUpwind_l[21]+0.25*alpha[2]*fUpwind_l[15]+0.25*alpha[0]*fUpwind_l[13]+0.25*alpha[7]*fUpwind_l[12]+0.25*alpha[1]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[8]+0.25*fUpwind_l[3]*alpha[8]+0.25*alpha[4]*fUpwind_l[6]; 
  Ghat_l[14] = 0.22360679774997902*alpha[8]*fUpwind_l[31]+0.22360679774997896*alpha[4]*fUpwind_l[30]+0.22360679774997896*alpha[3]*fUpwind_l[22]+0.22360679774997902*alpha[7]*fUpwind_l[19]+0.25*alpha[1]*fUpwind_l[15]+0.25*alpha[0]*fUpwind_l[14]+0.25*alpha[8]*fUpwind_l[11]+0.25*alpha[2]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[9]+0.25*alpha[4]*fUpwind_l[7]+0.25*fUpwind_l[4]*alpha[7]; 
  Ghat_l[15] = 0.22360679774997902*alpha[4]*fUpwind_l[31]+0.22360679774997896*alpha[8]*fUpwind_l[30]+0.22360679774997902*alpha[3]*fUpwind_l[23]+0.22360679774997896*alpha[7]*fUpwind_l[21]+0.25*alpha[0]*fUpwind_l[15]+0.25*alpha[1]*fUpwind_l[14]+0.25*alpha[2]*fUpwind_l[13]+0.25*alpha[3]*fUpwind_l[12]+0.25*alpha[4]*fUpwind_l[11]+0.25*alpha[7]*fUpwind_l[8]+0.25*fUpwind_l[7]*alpha[8]; 
  Ghat_l[16] = 0.25*alpha[8]*fUpwind_l[21]+0.25000000000000006*alpha[4]*fUpwind_l[19]+0.25000000000000006*alpha[2]*fUpwind_l[18]+0.25000000000000006*alpha[1]*fUpwind_l[17]+0.25*alpha[0]*fUpwind_l[16]+0.22360679774997896*alpha[7]*fUpwind_l[7]+0.22360679774997896*alpha[3]*fUpwind_l[3]; 
  Ghat_l[17] = 0.25000000000000006*alpha[4]*fUpwind_l[21]+0.25000000000000006*alpha[2]*fUpwind_l[20]+0.25*alpha[8]*fUpwind_l[19]+0.25*alpha[0]*fUpwind_l[17]+0.25000000000000006*alpha[1]*fUpwind_l[16]+0.22360679774997902*alpha[7]*fUpwind_l[11]+0.22360679774997902*alpha[3]*fUpwind_l[6]; 
  Ghat_l[18] = 0.25*alpha[8]*fUpwind_l[23]+0.25000000000000006*alpha[4]*fUpwind_l[22]+0.25000000000000006*alpha[1]*fUpwind_l[20]+0.25*alpha[0]*fUpwind_l[18]+0.25000000000000006*alpha[2]*fUpwind_l[16]+0.22360679774997902*alpha[3]*fUpwind_l[7]+0.22360679774997902*fUpwind_l[3]*alpha[7]; 
  Ghat_l[19] = 0.25000000000000006*alpha[2]*fUpwind_l[22]+0.25000000000000006*alpha[1]*fUpwind_l[21]+0.25*alpha[0]*fUpwind_l[19]+0.25*alpha[8]*fUpwind_l[17]+0.25000000000000006*alpha[4]*fUpwind_l[16]+0.22360679774997902*alpha[7]*fUpwind_l[14]+0.22360679774997902*alpha[3]*fUpwind_l[10]; 
  Ghat_l[20] = 0.25000000000000006*alpha[4]*fUpwind_l[23]+0.25*alpha[8]*fUpwind_l[22]+0.25*alpha[0]*fUpwind_l[20]+0.25000000000000006*alpha[1]*fUpwind_l[18]+0.25000000000000006*alpha[2]*fUpwind_l[17]+0.22360679774997896*alpha[3]*fUpwind_l[11]+0.22360679774997896*fUpwind_l[6]*alpha[7]; 
  Ghat_l[21] = 0.25000000000000006*alpha[2]*fUpwind_l[23]+0.25*alpha[0]*fUpwind_l[21]+0.25000000000000006*alpha[1]*fUpwind_l[19]+0.25000000000000006*alpha[4]*fUpwind_l[17]+0.25*alpha[8]*fUpwind_l[16]+0.22360679774997896*alpha[7]*fUpwind_l[15]+0.22360679774997896*alpha[3]*fUpwind_l[13]; 
  Ghat_l[22] = 0.25000000000000006*alpha[1]*fUpwind_l[23]+0.25*alpha[0]*fUpwind_l[22]+0.25*alpha[8]*fUpwind_l[20]+0.25000000000000006*alpha[2]*fUpwind_l[19]+0.25000000000000006*alpha[4]*fUpwind_l[18]+0.22360679774997896*alpha[3]*fUpwind_l[14]+0.22360679774997896*alpha[7]*fUpwind_l[10]; 
  Ghat_l[23] = 0.25*alpha[0]*fUpwind_l[23]+0.25000000000000006*alpha[1]*fUpwind_l[22]+0.25000000000000006*alpha[2]*fUpwind_l[21]+0.25000000000000006*alpha[4]*fUpwind_l[20]+0.25*alpha[8]*fUpwind_l[18]+0.22360679774997902*alpha[3]*fUpwind_l[15]+0.22360679774997902*alpha[7]*fUpwind_l[13]; 
  Ghat_l[24] = 0.25*alpha[7]*fUpwind_l[30]+0.25000000000000006*alpha[3]*fUpwind_l[27]+0.25000000000000006*alpha[2]*fUpwind_l[26]+0.25000000000000006*alpha[1]*fUpwind_l[25]+0.25*alpha[0]*fUpwind_l[24]+0.22360679774997896*alpha[8]*fUpwind_l[8]+0.22360679774997896*alpha[4]*fUpwind_l[4]; 
  Ghat_l[25] = 0.25*alpha[7]*fUpwind_l[31]+0.25000000000000006*alpha[3]*fUpwind_l[29]+0.25000000000000006*alpha[2]*fUpwind_l[28]+0.25*alpha[0]*fUpwind_l[25]+0.25000000000000006*alpha[1]*fUpwind_l[24]+0.22360679774997902*alpha[4]*fUpwind_l[8]+0.22360679774997902*fUpwind_l[4]*alpha[8]; 
  Ghat_l[26] = 0.25000000000000006*alpha[3]*fUpwind_l[30]+0.25000000000000006*alpha[1]*fUpwind_l[28]+0.25*alpha[7]*fUpwind_l[27]+0.25*alpha[0]*fUpwind_l[26]+0.25000000000000006*alpha[2]*fUpwind_l[24]+0.22360679774997902*alpha[8]*fUpwind_l[12]+0.22360679774997902*alpha[4]*fUpwind_l[9]; 
  Ghat_l[27] = 0.25000000000000006*alpha[2]*fUpwind_l[30]+0.25000000000000006*alpha[1]*fUpwind_l[29]+0.25*alpha[0]*fUpwind_l[27]+0.25*alpha[7]*fUpwind_l[26]+0.25000000000000006*alpha[3]*fUpwind_l[24]+0.22360679774997902*alpha[8]*fUpwind_l[13]+0.22360679774997902*alpha[4]*fUpwind_l[10]; 
  Ghat_l[28] = 0.25000000000000006*alpha[3]*fUpwind_l[31]+0.25*alpha[7]*fUpwind_l[29]+0.25*alpha[0]*fUpwind_l[28]+0.25000000000000006*alpha[1]*fUpwind_l[26]+0.25000000000000006*alpha[2]*fUpwind_l[25]+0.22360679774997896*alpha[4]*fUpwind_l[12]+0.22360679774997896*alpha[8]*fUpwind_l[9]; 
  Ghat_l[29] = 0.25000000000000006*alpha[2]*fUpwind_l[31]+0.25*alpha[0]*fUpwind_l[29]+0.25*alpha[7]*fUpwind_l[28]+0.25000000000000006*alpha[1]*fUpwind_l[27]+0.25000000000000006*alpha[3]*fUpwind_l[25]+0.22360679774997896*alpha[4]*fUpwind_l[13]+0.22360679774997896*alpha[8]*fUpwind_l[10]; 
  Ghat_l[30] = 0.25000000000000006*alpha[1]*fUpwind_l[31]+0.25*alpha[0]*fUpwind_l[30]+0.25000000000000006*alpha[2]*fUpwind_l[27]+0.25000000000000006*alpha[3]*fUpwind_l[26]+0.25*alpha[7]*fUpwind_l[24]+0.22360679774997896*alpha[8]*fUpwind_l[15]+0.22360679774997896*alpha[4]*fUpwind_l[14]; 
  Ghat_l[31] = 0.25*alpha[0]*fUpwind_l[31]+0.25000000000000006*alpha[1]*fUpwind_l[30]+0.25000000000000006*alpha[2]*fUpwind_l[29]+0.25000000000000006*alpha[3]*fUpwind_l[28]+0.25*alpha[7]*fUpwind_l[25]+0.22360679774997902*alpha[4]*fUpwind_l[15]+0.22360679774997902*alpha[8]*fUpwind_l[14]; 

  Ghat_r[0] = 0.25*alpha[8]*fUpwind_r[8]+0.25*alpha[7]*fUpwind_r[7]+0.25*alpha[4]*fUpwind_r[4]+0.25*alpha[3]*fUpwind_r[3]+0.25*alpha[2]*fUpwind_r[2]+0.25*alpha[1]*fUpwind_r[1]+0.25*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.25*alpha[7]*fUpwind_r[11]+0.25*alpha[4]*fUpwind_r[8]+0.25*fUpwind_r[4]*alpha[8]+0.25*alpha[3]*fUpwind_r[6]+0.25*alpha[2]*fUpwind_r[5]+0.25*alpha[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.25*alpha[8]*fUpwind_r[12]+0.25*alpha[4]*fUpwind_r[9]+0.25*alpha[3]*fUpwind_r[7]+0.25*fUpwind_r[3]*alpha[7]+0.25*alpha[1]*fUpwind_r[5]+0.25*alpha[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.22360679774997902*alpha[7]*fUpwind_r[18]+0.22360679774997896*alpha[3]*fUpwind_r[16]+0.25*alpha[8]*fUpwind_r[13]+0.25*alpha[4]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[7]+0.25*fUpwind_r[2]*alpha[7]+0.25*alpha[1]*fUpwind_r[6]+0.25*alpha[0]*fUpwind_r[3]+0.25*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.22360679774997902*alpha[8]*fUpwind_r[25]+0.22360679774997896*alpha[4]*fUpwind_r[24]+0.25*alpha[7]*fUpwind_r[14]+0.25*alpha[3]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[9]+0.25*alpha[1]*fUpwind_r[8]+0.25*fUpwind_r[1]*alpha[8]+0.25*alpha[0]*fUpwind_r[4]+0.25*fUpwind_r[0]*alpha[4]; 
  Ghat_r[5] = 0.25*alpha[4]*fUpwind_r[12]+0.25*alpha[3]*fUpwind_r[11]+0.25*alpha[8]*fUpwind_r[9]+0.25*fUpwind_r[6]*alpha[7]+0.25*alpha[0]*fUpwind_r[5]+0.25*alpha[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alpha[2]; 
  Ghat_r[6] = 0.22360679774997896*alpha[7]*fUpwind_r[20]+0.22360679774997902*alpha[3]*fUpwind_r[17]+0.25*alpha[4]*fUpwind_r[13]+0.25*alpha[2]*fUpwind_r[11]+0.25*alpha[8]*fUpwind_r[10]+0.25*fUpwind_r[5]*alpha[7]+0.25*alpha[0]*fUpwind_r[6]+0.25*alpha[1]*fUpwind_r[3]+0.25*fUpwind_r[1]*alpha[3]; 
  Ghat_r[7] = 0.22360679774997902*alpha[3]*fUpwind_r[18]+0.22360679774997896*alpha[7]*fUpwind_r[16]+0.25*alpha[8]*fUpwind_r[15]+0.25*alpha[4]*fUpwind_r[14]+0.25*alpha[1]*fUpwind_r[11]+0.25*alpha[0]*fUpwind_r[7]+0.25*fUpwind_r[0]*alpha[7]+0.25*alpha[2]*fUpwind_r[3]+0.25*fUpwind_r[2]*alpha[3]; 
  Ghat_r[8] = 0.22360679774997902*alpha[4]*fUpwind_r[25]+0.22360679774997896*alpha[8]*fUpwind_r[24]+0.25*alpha[7]*fUpwind_r[15]+0.25*alpha[3]*fUpwind_r[13]+0.25*alpha[2]*fUpwind_r[12]+0.25*alpha[0]*fUpwind_r[8]+0.25*fUpwind_r[0]*alpha[8]+0.25*alpha[1]*fUpwind_r[4]+0.25*fUpwind_r[1]*alpha[4]; 
  Ghat_r[9] = 0.22360679774997896*alpha[8]*fUpwind_r[28]+0.22360679774997902*alpha[4]*fUpwind_r[26]+0.25*alpha[3]*fUpwind_r[14]+0.25*alpha[1]*fUpwind_r[12]+0.25*alpha[7]*fUpwind_r[10]+0.25*alpha[0]*fUpwind_r[9]+0.25*fUpwind_r[5]*alpha[8]+0.25*alpha[2]*fUpwind_r[4]+0.25*fUpwind_r[2]*alpha[4]; 
  Ghat_r[10] = 0.22360679774997896*alpha[8]*fUpwind_r[29]+0.22360679774997902*alpha[4]*fUpwind_r[27]+0.22360679774997896*alpha[7]*fUpwind_r[22]+0.22360679774997902*alpha[3]*fUpwind_r[19]+0.25*alpha[2]*fUpwind_r[14]+0.25*alpha[1]*fUpwind_r[13]+0.25*alpha[0]*fUpwind_r[10]+0.25*alpha[7]*fUpwind_r[9]+0.25*fUpwind_r[6]*alpha[8]+0.25*alpha[3]*fUpwind_r[4]+0.25*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.22360679774997896*alpha[3]*fUpwind_r[20]+0.22360679774997902*alpha[7]*fUpwind_r[17]+0.25*alpha[4]*fUpwind_r[15]+0.25*alpha[8]*fUpwind_r[14]+0.25*alpha[0]*fUpwind_r[11]+0.25*alpha[1]*fUpwind_r[7]+0.25*fUpwind_r[1]*alpha[7]+0.25*alpha[2]*fUpwind_r[6]+0.25*alpha[3]*fUpwind_r[5]; 
  Ghat_r[12] = 0.22360679774997896*alpha[4]*fUpwind_r[28]+0.22360679774997902*alpha[8]*fUpwind_r[26]+0.25*alpha[3]*fUpwind_r[15]+0.25*alpha[7]*fUpwind_r[13]+0.25*alpha[0]*fUpwind_r[12]+0.25*alpha[1]*fUpwind_r[9]+0.25*alpha[2]*fUpwind_r[8]+0.25*fUpwind_r[2]*alpha[8]+0.25*alpha[4]*fUpwind_r[5]; 
  Ghat_r[13] = 0.22360679774997896*alpha[4]*fUpwind_r[29]+0.22360679774997902*alpha[8]*fUpwind_r[27]+0.22360679774997902*alpha[7]*fUpwind_r[23]+0.22360679774997896*alpha[3]*fUpwind_r[21]+0.25*alpha[2]*fUpwind_r[15]+0.25*alpha[0]*fUpwind_r[13]+0.25*alpha[7]*fUpwind_r[12]+0.25*alpha[1]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[8]+0.25*fUpwind_r[3]*alpha[8]+0.25*alpha[4]*fUpwind_r[6]; 
  Ghat_r[14] = 0.22360679774997902*alpha[8]*fUpwind_r[31]+0.22360679774997896*alpha[4]*fUpwind_r[30]+0.22360679774997896*alpha[3]*fUpwind_r[22]+0.22360679774997902*alpha[7]*fUpwind_r[19]+0.25*alpha[1]*fUpwind_r[15]+0.25*alpha[0]*fUpwind_r[14]+0.25*alpha[8]*fUpwind_r[11]+0.25*alpha[2]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[9]+0.25*alpha[4]*fUpwind_r[7]+0.25*fUpwind_r[4]*alpha[7]; 
  Ghat_r[15] = 0.22360679774997902*alpha[4]*fUpwind_r[31]+0.22360679774997896*alpha[8]*fUpwind_r[30]+0.22360679774997902*alpha[3]*fUpwind_r[23]+0.22360679774997896*alpha[7]*fUpwind_r[21]+0.25*alpha[0]*fUpwind_r[15]+0.25*alpha[1]*fUpwind_r[14]+0.25*alpha[2]*fUpwind_r[13]+0.25*alpha[3]*fUpwind_r[12]+0.25*alpha[4]*fUpwind_r[11]+0.25*alpha[7]*fUpwind_r[8]+0.25*fUpwind_r[7]*alpha[8]; 
  Ghat_r[16] = 0.25*alpha[8]*fUpwind_r[21]+0.25000000000000006*alpha[4]*fUpwind_r[19]+0.25000000000000006*alpha[2]*fUpwind_r[18]+0.25000000000000006*alpha[1]*fUpwind_r[17]+0.25*alpha[0]*fUpwind_r[16]+0.22360679774997896*alpha[7]*fUpwind_r[7]+0.22360679774997896*alpha[3]*fUpwind_r[3]; 
  Ghat_r[17] = 0.25000000000000006*alpha[4]*fUpwind_r[21]+0.25000000000000006*alpha[2]*fUpwind_r[20]+0.25*alpha[8]*fUpwind_r[19]+0.25*alpha[0]*fUpwind_r[17]+0.25000000000000006*alpha[1]*fUpwind_r[16]+0.22360679774997902*alpha[7]*fUpwind_r[11]+0.22360679774997902*alpha[3]*fUpwind_r[6]; 
  Ghat_r[18] = 0.25*alpha[8]*fUpwind_r[23]+0.25000000000000006*alpha[4]*fUpwind_r[22]+0.25000000000000006*alpha[1]*fUpwind_r[20]+0.25*alpha[0]*fUpwind_r[18]+0.25000000000000006*alpha[2]*fUpwind_r[16]+0.22360679774997902*alpha[3]*fUpwind_r[7]+0.22360679774997902*fUpwind_r[3]*alpha[7]; 
  Ghat_r[19] = 0.25000000000000006*alpha[2]*fUpwind_r[22]+0.25000000000000006*alpha[1]*fUpwind_r[21]+0.25*alpha[0]*fUpwind_r[19]+0.25*alpha[8]*fUpwind_r[17]+0.25000000000000006*alpha[4]*fUpwind_r[16]+0.22360679774997902*alpha[7]*fUpwind_r[14]+0.22360679774997902*alpha[3]*fUpwind_r[10]; 
  Ghat_r[20] = 0.25000000000000006*alpha[4]*fUpwind_r[23]+0.25*alpha[8]*fUpwind_r[22]+0.25*alpha[0]*fUpwind_r[20]+0.25000000000000006*alpha[1]*fUpwind_r[18]+0.25000000000000006*alpha[2]*fUpwind_r[17]+0.22360679774997896*alpha[3]*fUpwind_r[11]+0.22360679774997896*fUpwind_r[6]*alpha[7]; 
  Ghat_r[21] = 0.25000000000000006*alpha[2]*fUpwind_r[23]+0.25*alpha[0]*fUpwind_r[21]+0.25000000000000006*alpha[1]*fUpwind_r[19]+0.25000000000000006*alpha[4]*fUpwind_r[17]+0.25*alpha[8]*fUpwind_r[16]+0.22360679774997896*alpha[7]*fUpwind_r[15]+0.22360679774997896*alpha[3]*fUpwind_r[13]; 
  Ghat_r[22] = 0.25000000000000006*alpha[1]*fUpwind_r[23]+0.25*alpha[0]*fUpwind_r[22]+0.25*alpha[8]*fUpwind_r[20]+0.25000000000000006*alpha[2]*fUpwind_r[19]+0.25000000000000006*alpha[4]*fUpwind_r[18]+0.22360679774997896*alpha[3]*fUpwind_r[14]+0.22360679774997896*alpha[7]*fUpwind_r[10]; 
  Ghat_r[23] = 0.25*alpha[0]*fUpwind_r[23]+0.25000000000000006*alpha[1]*fUpwind_r[22]+0.25000000000000006*alpha[2]*fUpwind_r[21]+0.25000000000000006*alpha[4]*fUpwind_r[20]+0.25*alpha[8]*fUpwind_r[18]+0.22360679774997902*alpha[3]*fUpwind_r[15]+0.22360679774997902*alpha[7]*fUpwind_r[13]; 
  Ghat_r[24] = 0.25*alpha[7]*fUpwind_r[30]+0.25000000000000006*alpha[3]*fUpwind_r[27]+0.25000000000000006*alpha[2]*fUpwind_r[26]+0.25000000000000006*alpha[1]*fUpwind_r[25]+0.25*alpha[0]*fUpwind_r[24]+0.22360679774997896*alpha[8]*fUpwind_r[8]+0.22360679774997896*alpha[4]*fUpwind_r[4]; 
  Ghat_r[25] = 0.25*alpha[7]*fUpwind_r[31]+0.25000000000000006*alpha[3]*fUpwind_r[29]+0.25000000000000006*alpha[2]*fUpwind_r[28]+0.25*alpha[0]*fUpwind_r[25]+0.25000000000000006*alpha[1]*fUpwind_r[24]+0.22360679774997902*alpha[4]*fUpwind_r[8]+0.22360679774997902*fUpwind_r[4]*alpha[8]; 
  Ghat_r[26] = 0.25000000000000006*alpha[3]*fUpwind_r[30]+0.25000000000000006*alpha[1]*fUpwind_r[28]+0.25*alpha[7]*fUpwind_r[27]+0.25*alpha[0]*fUpwind_r[26]+0.25000000000000006*alpha[2]*fUpwind_r[24]+0.22360679774997902*alpha[8]*fUpwind_r[12]+0.22360679774997902*alpha[4]*fUpwind_r[9]; 
  Ghat_r[27] = 0.25000000000000006*alpha[2]*fUpwind_r[30]+0.25000000000000006*alpha[1]*fUpwind_r[29]+0.25*alpha[0]*fUpwind_r[27]+0.25*alpha[7]*fUpwind_r[26]+0.25000000000000006*alpha[3]*fUpwind_r[24]+0.22360679774997902*alpha[8]*fUpwind_r[13]+0.22360679774997902*alpha[4]*fUpwind_r[10]; 
  Ghat_r[28] = 0.25000000000000006*alpha[3]*fUpwind_r[31]+0.25*alpha[7]*fUpwind_r[29]+0.25*alpha[0]*fUpwind_r[28]+0.25000000000000006*alpha[1]*fUpwind_r[26]+0.25000000000000006*alpha[2]*fUpwind_r[25]+0.22360679774997896*alpha[4]*fUpwind_r[12]+0.22360679774997896*alpha[8]*fUpwind_r[9]; 
  Ghat_r[29] = 0.25000000000000006*alpha[2]*fUpwind_r[31]+0.25*alpha[0]*fUpwind_r[29]+0.25*alpha[7]*fUpwind_r[28]+0.25000000000000006*alpha[1]*fUpwind_r[27]+0.25000000000000006*alpha[3]*fUpwind_r[25]+0.22360679774997896*alpha[4]*fUpwind_r[13]+0.22360679774997896*alpha[8]*fUpwind_r[10]; 
  Ghat_r[30] = 0.25000000000000006*alpha[1]*fUpwind_r[31]+0.25*alpha[0]*fUpwind_r[30]+0.25000000000000006*alpha[2]*fUpwind_r[27]+0.25000000000000006*alpha[3]*fUpwind_r[26]+0.25*alpha[7]*fUpwind_r[24]+0.22360679774997896*alpha[8]*fUpwind_r[15]+0.22360679774997896*alpha[4]*fUpwind_r[14]; 
  Ghat_r[31] = 0.25*alpha[0]*fUpwind_r[31]+0.25000000000000006*alpha[1]*fUpwind_r[30]+0.25000000000000006*alpha[2]*fUpwind_r[29]+0.25000000000000006*alpha[3]*fUpwind_r[28]+0.25*alpha[7]*fUpwind_r[25]+0.22360679774997902*alpha[4]*fUpwind_r[15]+0.22360679774997902*alpha[8]*fUpwind_r[14]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv12; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv12; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv12; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv12; 
  out[4] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv12; 
  out[5] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv12); 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv12; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv12; 
  out[8] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv12; 
  out[9] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv12; 
  out[10] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv12; 
  out[11] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv12; 
  out[12] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv12); 
  out[13] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv12); 
  out[14] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv12); 
  out[15] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv12); 
  out[16] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv12; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv12; 
  out[18] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv12; 
  out[19] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv12; 
  out[20] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv12); 
  out[21] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv12); 
  out[22] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv12); 
  out[23] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv12); 
  out[24] += -(1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv12); 
  out[25] += -(1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv12); 
  out[26] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv12; 
  out[27] += -(1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv12); 
  out[28] += -(1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv12); 
  out[29] += -(1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv12); 
  out[30] += -(1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv12); 
  out[31] += -(1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv12); 
  out[32] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv12; 
  out[33] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv12; 
  out[34] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv12; 
  out[35] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv12; 
  out[36] += -(1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv12); 
  out[37] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dv12; 
  out[38] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dv12; 
  out[39] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dv12; 
  out[40] += -(1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv12); 
  out[41] += -(1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv12); 
  out[42] += -(1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv12); 
  out[43] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dv12; 
  out[44] += -(1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dv12); 
  out[45] += -(1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dv12); 
  out[46] += -(1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dv12); 
  out[47] += -(1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dv12); 
  out[48] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dv12; 
  out[49] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dv12; 
  out[50] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dv12; 
  out[51] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dv12; 
  out[52] += -(1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dv12); 
  out[53] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dv12; 
  out[54] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dv12; 
  out[55] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dv12; 
  out[56] += -(1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dv12); 
  out[57] += -(1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dv12); 
  out[58] += -(1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dv12); 
  out[59] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dv12; 
  out[60] += -(1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dv12); 
  out[61] += -(1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dv12); 
  out[62] += -(1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dv12); 
  out[63] += -(1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dv12); 
  out[64] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv12; 
  out[65] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv12; 
  out[66] += (1.5811388300841898*Ghat_l[2]-1.5811388300841898*Ghat_r[2])*dv12; 
  out[67] += (1.5811388300841898*Ghat_l[3]-1.5811388300841898*Ghat_r[3])*dv12; 
  out[68] += (1.5811388300841898*Ghat_l[4]-1.5811388300841898*Ghat_r[4])*dv12; 
  out[69] += (1.5811388300841895*Ghat_l[5]-1.5811388300841895*Ghat_r[5])*dv12; 
  out[70] += (1.5811388300841895*Ghat_l[6]-1.5811388300841895*Ghat_r[6])*dv12; 
  out[71] += (1.5811388300841895*Ghat_l[7]-1.5811388300841895*Ghat_r[7])*dv12; 
  out[72] += (1.5811388300841895*Ghat_l[8]-1.5811388300841895*Ghat_r[8])*dv12; 
  out[73] += (1.5811388300841895*Ghat_l[9]-1.5811388300841895*Ghat_r[9])*dv12; 
  out[74] += (1.5811388300841895*Ghat_l[10]-1.5811388300841895*Ghat_r[10])*dv12; 
  out[75] += (1.5811388300841898*Ghat_l[11]-1.5811388300841898*Ghat_r[11])*dv12; 
  out[76] += (1.5811388300841898*Ghat_l[12]-1.5811388300841898*Ghat_r[12])*dv12; 
  out[77] += (1.5811388300841898*Ghat_l[13]-1.5811388300841898*Ghat_r[13])*dv12; 
  out[78] += (1.5811388300841898*Ghat_l[14]-1.5811388300841898*Ghat_r[14])*dv12; 
  out[79] += (1.5811388300841895*Ghat_l[15]-1.5811388300841895*Ghat_r[15])*dv12; 

  return 0.;

} 
