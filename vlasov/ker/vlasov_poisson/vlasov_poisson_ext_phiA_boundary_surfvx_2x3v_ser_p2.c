#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // pots:        potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:       external E and B fields (scaled by q/m).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 

  const double dv10 = 2/dxv[2]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  double alpha[48] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ax = &pots[8]; 
  const double *Ay = &pots[16]; 
  const double *Az = &pots[24]; 

  alpha[0] = 3.4641016151377544*Az[1]*dx10*wv3-3.4641016151377544*Ax[2]*dx11*wv2+3.4641016151377544*Ay[1]*dx10*wv2-3.4641016151377544*phi[1]*dx10; 
  alpha[1] = 7.745966692414834*Az[4]*dx10*wv3-3.4641016151377544*Ax[3]*dx11*wv2+7.745966692414834*Ay[4]*dx10*wv2-7.745966692414834*phi[4]*dx10; 
  alpha[2] = 3.4641016151377544*Az[3]*dx10*wv3-7.745966692414834*Ax[5]*dx11*wv2+3.4641016151377544*Ay[3]*dx10*wv2-3.4641016151377544*phi[3]*dx10; 
  alpha[3] = Ay[1]*dv2*dx10-1.0*Ax[2]*dv2*dx11; 
  alpha[4] = Az[1]*dv3*dx10; 
  alpha[5] = 7.745966692414834*Az[6]*dx10*wv3-7.745966692414834*Ax[7]*dx11*wv2+7.745966692414834*Ay[6]*dx10*wv2-7.745966692414834*phi[6]*dx10; 
  alpha[6] = 2.23606797749979*Ay[4]*dv2*dx10-1.0*Ax[3]*dv2*dx11; 
  alpha[7] = Ay[3]*dv2*dx10-2.23606797749979*Ax[5]*dv2*dx11; 
  alpha[8] = 2.23606797749979*Az[4]*dv3*dx10; 
  alpha[9] = Az[3]*dv3*dx10; 
  alpha[11] = -(3.464101615137755*Ax[6]*dx11*wv2); 
  alpha[12] = 3.464101615137755*Az[7]*dx10*wv3+3.464101615137755*Ay[7]*dx10*wv2-3.464101615137755*phi[7]*dx10; 
  alpha[15] = 2.2360679774997902*Ay[6]*dv2*dx10-2.2360679774997902*Ax[7]*dv2*dx11; 
  alpha[16] = 2.2360679774997902*Az[6]*dv3*dx10; 
  alpha[21] = -(1.0*Ax[6]*dv2*dx11); 
  alpha[22] = Ay[7]*dv2*dx10; 
  alpha[26] = Az[7]*dv3*dx10; 

  double fUpwindQuad[81] = {0.0};
  double fUpwind[48] = {0.0};
  double Ghat[48] = {0.0}; 

  if (edge == -1) { 

  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[5]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.33541019662496846*alpha[3]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if (0.375*(alpha[26]+alpha[22])-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*(alpha[8]+alpha[6])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_l(fEdge); 
  } 
  if (0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[6]-0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_l(fEdge); 
  } 
  if (-(0.375*alpha[26])+0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]+0.45*alpha[6]+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_l(fEdge); 
  } 
  if (0.375*alpha[26]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]-0.33541019662496846*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_l(fEdge); 
  } 
  if (-(0.2795084971874737*alpha[12])+0.22360679774997896*alpha[11]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_l(fEdge); 
  } 
  if (-(0.375*alpha[26])-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_l(fEdge); 
  } 
  if (0.375*alpha[26]-0.375*alpha[22]+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]-0.45*alpha[6]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_l(fEdge); 
  } 
  if (-(0.375*alpha[22])+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[6]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_l(fEdge); 
  } 
  if (-(0.375*(alpha[26]+alpha[22]))+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*(alpha[8]+alpha[6])+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_l(fEdge); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[5]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.33541019662496846*(alpha[3]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]))+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_r(fSkin); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[22])+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]-0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_r(fSkin); 
  } else { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*alpha[22]+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_r(fSkin); 
  } else { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.33541019662496846*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_r(fSkin); 
  } else { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_l(fEdge); 
  } 
  if (0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_r(fSkin); 
  } else { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_r(fSkin); 
  } else { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_r(fSkin); 
  } else { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_r(fSkin); 
  } else { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22])-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_r(fSkin); 
  } else { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_l(fEdge); 
  } 
  if (0.375*(alpha[26]+alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_r(fSkin); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_l(fEdge); 
  } 
  if (0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_r(fSkin); 
  } else { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_l(fEdge); 
  } 
  if (-(0.375*alpha[26])+0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_r(fSkin); 
  } else { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_l(fEdge); 
  } 
  if (0.375*alpha[26]-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_r(fSkin); 
  } else { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_l(fEdge); 
  } 
  if (0.25*alpha[0]-0.2795084971874737*(alpha[12]+alpha[11]) > 0) { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_r(fSkin); 
  } else { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_l(fEdge); 
  } 
  if (-(0.375*alpha[26])-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_r(fSkin); 
  } else { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_l(fEdge); 
  } 
  if (0.375*alpha[26]-0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_r(fSkin); 
  } else { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_l(fEdge); 
  } 
  if (-(0.375*(alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_r(fSkin); 
  } else { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_l(fEdge); 
  } 
  if (-(0.375*(alpha[26]+alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_r(fSkin); 
  } else { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]))+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_r(fSkin); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[22])+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_r(fSkin); 
  } else { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*alpha[22]+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_r(fSkin); 
  } else { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_r(fSkin); 
  } else { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_l(fEdge); 
  } 
  if (0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_r(fSkin); 
  } else { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]+0.33541019662496846*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_r(fSkin); 
  } else { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_r(fSkin); 
  } else { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]+0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_r(fSkin); 
  } else { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22])-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_r(fSkin); 
  } else { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_r(fSkin); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.33541019662496846*(alpha[3]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_r(fSkin); 
  } else { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_r(fSkin); 
  } else { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_r(fSkin); 
  } else { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_l(fEdge); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[5]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_r(fSkin); 
  } else { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_r(fSkin); 
  } else { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_r(fSkin); 
  } else { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_r(fSkin); 
  } else { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_r(fSkin); 
  } else { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_l(fEdge); 
  } 
  if (0.375*(alpha[26]+alpha[22])-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*(alpha[8]+alpha[6])-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_r(fSkin); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_l(fEdge); 
  } 
  if (0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[6]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_r(fSkin); 
  } else { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_l(fEdge); 
  } 
  if (-(0.375*alpha[26])+0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]-0.45*alpha[6]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_r(fSkin); 
  } else { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_l(fEdge); 
  } 
  if (0.375*alpha[26]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_r(fSkin); 
  } else { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_l(fEdge); 
  } 
  if (-(0.2795084971874737*alpha[12])+0.22360679774997896*alpha[11]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_r(fSkin); 
  } else { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_l(fEdge); 
  } 
  if (-(0.375*alpha[26])-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]+0.33541019662496846*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_r(fSkin); 
  } else { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_l(fEdge); 
  } 
  if (0.375*alpha[26]-0.375*alpha[22]+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]+0.45*alpha[6]-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_r(fSkin); 
  } else { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_l(fEdge); 
  } 
  if (-(0.375*alpha[22])+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[6]+0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_r(fSkin); 
  } else { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_l(fEdge); 
  } 
  if (-(0.375*(alpha[26]+alpha[22]))+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*(alpha[8]+alpha[6])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_r(fSkin); 
  } else { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_r(fSkin); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_l(fEdge); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.33541019662496846*alpha[3]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_r(fSkin); 
  } else { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_r(fSkin); 
  } else { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_r(fSkin); 
  } else { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_l(fEdge); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[5]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_r(fSkin); 
  } else { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_l(fEdge); 
  } 
  if (0.29999999999999993*alpha[26]+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_r(fSkin); 
  } else { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_l(fEdge); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_r(fSkin); 
  } else { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_r(fSkin); 
  } else { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_l(fEdge); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_r(fSkin); 
  } else { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[26]*fUpwind[26]+0.25*alpha[22]*fUpwind[22]+0.25*alpha[21]*fUpwind[21]+0.25*alpha[16]*fUpwind[16]+0.25*alpha[15]*fUpwind[15]+0.25*alpha[12]*fUpwind[12]+0.25*alpha[11]*fUpwind[11]+0.25*alpha[9]*fUpwind[9]+0.25*alpha[8]*fUpwind[8]+0.25*alpha[7]*fUpwind[7]+0.25*alpha[6]*fUpwind[6]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[4]*fUpwind[4]+0.25*alpha[3]*fUpwind[3]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.25000000000000006*alpha[26]*fUpwind[36]+0.22360679774997896*alpha[16]*fUpwind[35]+0.25000000000000006*alpha[22]*fUpwind[33]+0.22360679774997896*alpha[15]*fUpwind[32]+0.22360679774997902*alpha[8]*fUpwind[25]+0.22360679774997902*alpha[6]*fUpwind[21]+0.22360679774997902*fUpwind[6]*alpha[21]+0.25000000000000006*alpha[12]*fUpwind[20]+0.22360679774997902*alpha[5]*fUpwind[19]+0.25*alpha[9]*fUpwind[16]+0.25*fUpwind[9]*alpha[16]+0.25*alpha[7]*fUpwind[15]+0.25*fUpwind[7]*alpha[15]+0.22360679774997896*alpha[1]*fUpwind[11]+0.22360679774997896*fUpwind[1]*alpha[11]+0.25*alpha[4]*fUpwind[8]+0.25*fUpwind[4]*alpha[8]+0.25*alpha[3]*fUpwind[6]+0.25*fUpwind[3]*alpha[6]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.22360679774997896*alpha[16]*fUpwind[36]+0.22360679774997896*alpha[15]*fUpwind[33]+0.25000000000000006*alpha[21]*fUpwind[32]+0.22360679774997902*alpha[9]*fUpwind[26]+0.22360679774997902*fUpwind[9]*alpha[26]+0.22360679774997902*alpha[7]*fUpwind[22]+0.22360679774997902*fUpwind[7]*alpha[22]+0.22360679774997902*alpha[5]*fUpwind[20]+0.25000000000000006*alpha[11]*fUpwind[19]+0.25*alpha[8]*fUpwind[16]+0.25*fUpwind[8]*alpha[16]+0.25*alpha[6]*fUpwind[15]+0.25*fUpwind[6]*alpha[15]+0.22360679774997896*alpha[2]*fUpwind[12]+0.22360679774997896*fUpwind[2]*alpha[12]+0.25*alpha[4]*fUpwind[9]+0.25*fUpwind[4]*alpha[9]+0.25*alpha[3]*fUpwind[7]+0.25*fUpwind[3]*alpha[7]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.25000000000000006*alpha[26]*fUpwind[38]+0.22360679774997896*alpha[15]*fUpwind[34]+0.25*alpha[16]*fUpwind[31]+0.22360679774997902*alpha[7]*fUpwind[24]+0.22360679774997902*alpha[6]*fUpwind[23]+0.25000000000000006*alpha[12]*fUpwind[22]+0.25000000000000006*fUpwind[12]*alpha[22]+0.25000000000000006*alpha[11]*fUpwind[21]+0.25000000000000006*fUpwind[11]*alpha[21]+0.25*alpha[9]*fUpwind[18]+0.25*alpha[8]*fUpwind[17]+0.25*alpha[5]*fUpwind[15]+0.25*fUpwind[5]*alpha[15]+0.22360679774997896*alpha[3]*fUpwind[13]+0.25*alpha[4]*fUpwind[10]+0.25*alpha[2]*fUpwind[7]+0.25*fUpwind[2]*alpha[7]+0.25*alpha[1]*fUpwind[6]+0.25*fUpwind[1]*alpha[6]+0.25*alpha[0]*fUpwind[3]+0.25*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.22360679774997896*alpha[16]*fUpwind[41]+0.25000000000000006*alpha[22]*fUpwind[38]+0.25000000000000006*alpha[21]*fUpwind[37]+0.25*alpha[15]*fUpwind[31]+0.22360679774997902*alpha[9]*fUpwind[29]+0.22360679774997902*alpha[8]*fUpwind[28]+0.25000000000000006*alpha[12]*fUpwind[26]+0.25000000000000006*fUpwind[12]*alpha[26]+0.25000000000000006*alpha[11]*fUpwind[25]+0.25*alpha[7]*fUpwind[18]+0.25*alpha[6]*fUpwind[17]+0.25*alpha[5]*fUpwind[16]+0.25*fUpwind[5]*alpha[16]+0.22360679774997896*alpha[4]*fUpwind[14]+0.25*alpha[3]*fUpwind[10]+0.25*alpha[2]*fUpwind[9]+0.25*fUpwind[2]*alpha[9]+0.25*alpha[1]*fUpwind[8]+0.25*fUpwind[1]*alpha[8]+0.25*alpha[0]*fUpwind[4]+0.25*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.22360679774997896*alpha[9]*fUpwind[36]+0.22360679774997896*alpha[8]*fUpwind[35]+0.22360679774997896*alpha[7]*fUpwind[33]+0.22360679774997896*alpha[6]*fUpwind[32]+0.22360679774997902*alpha[16]*fUpwind[26]+0.22360679774997902*fUpwind[16]*alpha[26]+0.22360679774997902*alpha[16]*fUpwind[25]+0.22360679774997902*alpha[15]*fUpwind[22]+0.22360679774997902*fUpwind[15]*alpha[22]+0.22360679774997902*alpha[15]*fUpwind[21]+0.22360679774997902*fUpwind[15]*alpha[21]+0.22360679774997902*alpha[2]*fUpwind[20]+0.22360679774997902*alpha[1]*fUpwind[19]+0.25*alpha[4]*fUpwind[16]+0.25*fUpwind[4]*alpha[16]+0.25*alpha[3]*fUpwind[15]+0.25*fUpwind[3]*alpha[15]+0.22360679774997896*alpha[5]*fUpwind[12]+0.22360679774997896*fUpwind[5]*alpha[12]+0.22360679774997896*alpha[5]*fUpwind[11]+0.22360679774997896*fUpwind[5]*alpha[11]+0.25*alpha[8]*fUpwind[9]+0.25*fUpwind[8]*alpha[9]+0.25*alpha[6]*fUpwind[7]+0.25*fUpwind[6]*alpha[7]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.25*alpha[26]*fUpwind[45]+0.22360679774997902*alpha[16]*fUpwind[44]+0.22360679774997896*alpha[8]*fUpwind[37]+0.22360679774997896*alpha[7]*fUpwind[34]+0.25*alpha[12]*fUpwind[33]+0.22360679774997896*alpha[5]*fUpwind[32]+0.25*alpha[9]*fUpwind[31]+0.22360679774997902*alpha[15]*fUpwind[24]+0.2*alpha[21]*fUpwind[23]+0.22360679774997902*alpha[3]*fUpwind[23]+0.25*fUpwind[20]*alpha[22]+0.22360679774997902*alpha[1]*fUpwind[21]+0.22360679774997902*fUpwind[1]*alpha[21]+0.22360679774997902*alpha[15]*fUpwind[19]+0.25*alpha[16]*fUpwind[18]+0.25*alpha[4]*fUpwind[17]+0.25*alpha[2]*fUpwind[15]+0.25*fUpwind[2]*alpha[15]+0.22360679774997896*alpha[6]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[11]+0.22360679774997896*fUpwind[6]*alpha[11]+0.25*alpha[8]*fUpwind[10]+0.25*alpha[5]*fUpwind[7]+0.25*fUpwind[5]*alpha[7]+0.25*alpha[0]*fUpwind[6]+0.25*fUpwind[0]*alpha[6]+0.25*alpha[1]*fUpwind[3]+0.25*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.22360679774997902*alpha[16]*fUpwind[45]+0.22360679774997896*alpha[9]*fUpwind[38]+0.22360679774997896*alpha[6]*fUpwind[34]+0.22360679774997896*alpha[5]*fUpwind[33]+0.25*alpha[11]*fUpwind[32]+0.25*alpha[8]*fUpwind[31]+0.22360679774997902*fUpwind[18]*alpha[26]+0.2*alpha[22]*fUpwind[24]+0.22360679774997902*alpha[3]*fUpwind[24]+0.22360679774997902*alpha[15]*fUpwind[23]+0.22360679774997902*alpha[2]*fUpwind[22]+0.22360679774997902*fUpwind[2]*alpha[22]+0.25*fUpwind[19]*alpha[21]+0.22360679774997902*alpha[15]*fUpwind[20]+0.25*alpha[4]*fUpwind[18]+0.25*alpha[16]*fUpwind[17]+0.25*alpha[1]*fUpwind[15]+0.25*fUpwind[1]*alpha[15]+0.22360679774997896*alpha[7]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[12]+0.22360679774997896*fUpwind[7]*alpha[12]+0.25*alpha[9]*fUpwind[10]+0.25*alpha[0]*fUpwind[7]+0.25*fUpwind[0]*alpha[7]+0.25*alpha[5]*fUpwind[6]+0.25*fUpwind[5]*alpha[6]+0.25*alpha[2]*fUpwind[3]+0.25*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.25*alpha[22]*fUpwind[45]+0.22360679774997902*alpha[15]*fUpwind[44]+0.22360679774997896*alpha[9]*fUpwind[41]+0.22360679774997896*alpha[6]*fUpwind[37]+0.25*alpha[12]*fUpwind[36]+0.22360679774997896*alpha[5]*fUpwind[35]+0.25*alpha[7]*fUpwind[31]+0.22360679774997902*alpha[16]*fUpwind[29]+0.22360679774997902*alpha[4]*fUpwind[28]+0.25*fUpwind[20]*alpha[26]+0.22360679774997902*alpha[1]*fUpwind[25]+0.22360679774997902*fUpwind[17]*alpha[21]+0.22360679774997902*alpha[16]*fUpwind[19]+0.25*alpha[15]*fUpwind[18]+0.25*alpha[3]*fUpwind[17]+0.25*alpha[2]*fUpwind[16]+0.25*fUpwind[2]*alpha[16]+0.22360679774997896*alpha[8]*fUpwind[14]+0.22360679774997896*alpha[8]*fUpwind[11]+0.22360679774997896*fUpwind[8]*alpha[11]+0.25*alpha[6]*fUpwind[10]+0.25*alpha[5]*fUpwind[9]+0.25*fUpwind[5]*alpha[9]+0.25*alpha[0]*fUpwind[8]+0.25*fUpwind[0]*alpha[8]+0.25*alpha[1]*fUpwind[4]+0.25*fUpwind[1]*alpha[4]; 
  Ghat[9] = 0.22360679774997902*alpha[15]*fUpwind[45]+0.25*alpha[21]*fUpwind[44]+0.22360679774997896*alpha[8]*fUpwind[41]+0.22360679774997896*alpha[7]*fUpwind[38]+0.22360679774997896*alpha[5]*fUpwind[36]+0.25*alpha[11]*fUpwind[35]+0.25*alpha[6]*fUpwind[31]+0.2*alpha[26]*fUpwind[29]+0.22360679774997902*alpha[4]*fUpwind[29]+0.22360679774997902*alpha[16]*fUpwind[28]+0.22360679774997902*alpha[2]*fUpwind[26]+0.22360679774997902*fUpwind[2]*alpha[26]+0.22360679774997902*fUpwind[18]*alpha[22]+0.22360679774997902*alpha[16]*fUpwind[20]+0.25*alpha[3]*fUpwind[18]+0.25*alpha[15]*fUpwind[17]+0.25*alpha[1]*fUpwind[16]+0.25*fUpwind[1]*alpha[16]+0.22360679774997896*alpha[9]*fUpwind[14]+0.22360679774997896*alpha[9]*fUpwind[12]+0.22360679774997896*fUpwind[9]*alpha[12]+0.25*alpha[7]*fUpwind[10]+0.25*alpha[0]*fUpwind[9]+0.25*fUpwind[0]*alpha[9]+0.25*alpha[5]*fUpwind[8]+0.25*fUpwind[5]*alpha[8]+0.25*alpha[2]*fUpwind[4]+0.25*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.22360679774997902*alpha[16]*fUpwind[47]+0.22360679774997902*alpha[15]*fUpwind[46]+0.22360679774997896*alpha[9]*fUpwind[43]+0.22360679774997896*alpha[8]*fUpwind[42]+0.22360679774997896*alpha[7]*fUpwind[40]+0.22360679774997896*alpha[6]*fUpwind[39]+0.25*alpha[12]*fUpwind[38]+0.25*alpha[11]*fUpwind[37]+0.25*alpha[5]*fUpwind[31]+0.22360679774997902*alpha[4]*fUpwind[30]+0.22360679774997902*alpha[3]*fUpwind[27]+0.25*alpha[22]*fUpwind[26]+0.25*fUpwind[22]*alpha[26]+0.25*alpha[21]*fUpwind[25]+0.25*alpha[2]*fUpwind[18]+0.25*alpha[1]*fUpwind[17]+0.25*alpha[15]*fUpwind[16]+0.25*fUpwind[15]*alpha[16]+0.25*alpha[0]*fUpwind[10]+0.25*alpha[7]*fUpwind[9]+0.25*fUpwind[7]*alpha[9]+0.25*alpha[6]*fUpwind[8]+0.25*fUpwind[6]*alpha[8]+0.25*alpha[3]*fUpwind[4]+0.25*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.25*alpha[9]*fUpwind[35]+0.25*alpha[7]*fUpwind[32]+0.25000000000000006*alpha[4]*fUpwind[25]+0.15971914124998499*alpha[21]*fUpwind[21]+0.25000000000000006*alpha[3]*fUpwind[21]+0.25000000000000006*fUpwind[3]*alpha[21]+0.25000000000000006*alpha[2]*fUpwind[19]+0.22360679774997896*alpha[16]*fUpwind[16]+0.22360679774997896*alpha[15]*fUpwind[15]+0.15971914124998499*alpha[11]*fUpwind[11]+0.25*alpha[0]*fUpwind[11]+0.25*fUpwind[0]*alpha[11]+0.22360679774997896*alpha[8]*fUpwind[8]+0.22360679774997896*alpha[6]*fUpwind[6]+0.22360679774997896*alpha[5]*fUpwind[5]+0.22360679774997896*alpha[1]*fUpwind[1]; 
  Ghat[12] = 0.25*alpha[8]*fUpwind[36]+0.25*alpha[6]*fUpwind[33]+0.15971914124998499*alpha[26]*fUpwind[26]+0.25000000000000006*alpha[4]*fUpwind[26]+0.25000000000000006*fUpwind[4]*alpha[26]+0.15971914124998499*alpha[22]*fUpwind[22]+0.25000000000000006*alpha[3]*fUpwind[22]+0.25000000000000006*fUpwind[3]*alpha[22]+0.25000000000000006*alpha[1]*fUpwind[20]+0.22360679774997896*alpha[16]*fUpwind[16]+0.22360679774997896*alpha[15]*fUpwind[15]+0.15971914124998499*alpha[12]*fUpwind[12]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.22360679774997896*alpha[9]*fUpwind[9]+0.22360679774997896*alpha[7]*fUpwind[7]+0.22360679774997896*alpha[5]*fUpwind[5]+0.22360679774997896*alpha[2]*fUpwind[2]; 
  Ghat[13] = 0.25000000000000006*alpha[16]*fUpwind[46]+0.25*alpha[9]*fUpwind[40]+0.25*alpha[8]*fUpwind[39]+0.25*alpha[5]*fUpwind[34]+0.25000000000000006*alpha[4]*fUpwind[27]+0.25000000000000006*alpha[2]*fUpwind[24]+0.25000000000000006*alpha[1]*fUpwind[23]+0.22360679774997896*alpha[22]*fUpwind[22]+0.22360679774997896*alpha[21]*fUpwind[21]+0.22360679774997896*alpha[15]*fUpwind[15]+0.25*alpha[0]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[7]+0.22360679774997896*alpha[6]*fUpwind[6]+0.22360679774997896*alpha[3]*fUpwind[3]; 
  Ghat[14] = 0.25000000000000006*alpha[15]*fUpwind[47]+0.25*alpha[7]*fUpwind[43]+0.25*alpha[6]*fUpwind[42]+0.25*alpha[5]*fUpwind[41]+0.25000000000000006*alpha[3]*fUpwind[30]+0.25000000000000006*alpha[2]*fUpwind[29]+0.25000000000000006*alpha[1]*fUpwind[28]+0.22360679774997896*alpha[26]*fUpwind[26]+0.22360679774997896*alpha[16]*fUpwind[16]+0.25*alpha[0]*fUpwind[14]+0.22360679774997896*alpha[9]*fUpwind[9]+0.22360679774997896*alpha[8]*fUpwind[8]+0.22360679774997896*alpha[4]*fUpwind[4]; 
  Ghat[15] = 0.22360679774997902*alpha[9]*fUpwind[45]+0.22360679774997902*alpha[8]*fUpwind[44]+0.22360679774997896*alpha[16]*fUpwind[38]+0.22360679774997896*alpha[16]*fUpwind[37]+0.2*alpha[22]*fUpwind[34]+0.2*alpha[21]*fUpwind[34]+0.22360679774997896*alpha[3]*fUpwind[34]+0.22360679774997896*alpha[2]*fUpwind[33]+0.22360679774997896*alpha[1]*fUpwind[32]+0.22360679774997902*alpha[26]*fUpwind[31]+0.25*alpha[4]*fUpwind[31]+0.22360679774997902*alpha[6]*fUpwind[24]+0.22360679774997902*alpha[7]*fUpwind[23]+0.22360679774997902*alpha[5]*fUpwind[22]+0.22360679774997902*fUpwind[5]*alpha[22]+0.22360679774997902*alpha[5]*fUpwind[21]+0.22360679774997902*fUpwind[5]*alpha[21]+0.22360679774997902*alpha[7]*fUpwind[20]+0.22360679774997902*alpha[6]*fUpwind[19]+0.25*alpha[8]*fUpwind[18]+0.25*alpha[9]*fUpwind[17]+0.25*fUpwind[10]*alpha[16]+0.22360679774997896*alpha[12]*fUpwind[15]+0.22360679774997896*alpha[11]*fUpwind[15]+0.25*alpha[0]*fUpwind[15]+0.22360679774997896*fUpwind[13]*alpha[15]+0.22360679774997896*fUpwind[12]*alpha[15]+0.22360679774997896*fUpwind[11]*alpha[15]+0.25*fUpwind[0]*alpha[15]+0.25*alpha[1]*fUpwind[7]+0.25*fUpwind[1]*alpha[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[2]*alpha[6]+0.25*alpha[3]*fUpwind[5]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[16] = 0.22360679774997902*alpha[7]*fUpwind[45]+0.22360679774997902*alpha[6]*fUpwind[44]+0.2*alpha[26]*fUpwind[41]+0.22360679774997896*alpha[4]*fUpwind[41]+0.22360679774997896*alpha[15]*fUpwind[38]+0.22360679774997896*alpha[15]*fUpwind[37]+0.22360679774997896*alpha[2]*fUpwind[36]+0.22360679774997896*alpha[1]*fUpwind[35]+0.22360679774997902*alpha[22]*fUpwind[31]+0.22360679774997902*alpha[21]*fUpwind[31]+0.25*alpha[3]*fUpwind[31]+0.22360679774997902*alpha[8]*fUpwind[29]+0.22360679774997902*alpha[9]*fUpwind[28]+0.22360679774997902*alpha[5]*fUpwind[26]+0.22360679774997902*fUpwind[5]*alpha[26]+0.22360679774997902*alpha[5]*fUpwind[25]+0.22360679774997902*alpha[9]*fUpwind[20]+0.22360679774997902*alpha[8]*fUpwind[19]+0.25*alpha[6]*fUpwind[18]+0.25*alpha[7]*fUpwind[17]+0.22360679774997896*alpha[12]*fUpwind[16]+0.22360679774997896*alpha[11]*fUpwind[16]+0.25*alpha[0]*fUpwind[16]+0.22360679774997896*fUpwind[14]*alpha[16]+0.22360679774997896*fUpwind[12]*alpha[16]+0.22360679774997896*fUpwind[11]*alpha[16]+0.25*fUpwind[0]*alpha[16]+0.25*fUpwind[10]*alpha[15]+0.25*alpha[1]*fUpwind[9]+0.25*fUpwind[1]*alpha[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[2]*alpha[8]+0.25*alpha[4]*fUpwind[5]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[17] = 0.22360679774997902*alpha[9]*fUpwind[47]+0.22360679774997902*alpha[7]*fUpwind[46]+0.25000000000000006*alpha[12]*fUpwind[45]+0.22360679774997902*alpha[5]*fUpwind[44]+0.22360679774997896*alpha[16]*fUpwind[43]+0.22360679774997896*alpha[4]*fUpwind[42]+0.22360679774997896*alpha[15]*fUpwind[40]+0.2*alpha[21]*fUpwind[39]+0.22360679774997896*alpha[3]*fUpwind[39]+0.22360679774997896*alpha[1]*fUpwind[37]+0.25000000000000006*alpha[22]*fUpwind[36]+0.22360679774997896*alpha[15]*fUpwind[35]+0.25000000000000006*alpha[26]*fUpwind[33]+0.22360679774997896*alpha[16]*fUpwind[32]+0.25*alpha[2]*fUpwind[31]+0.22360679774997902*alpha[8]*fUpwind[30]+0.22360679774997902*alpha[6]*fUpwind[27]+0.22360679774997902*alpha[6]*fUpwind[25]+0.22360679774997902*alpha[8]*fUpwind[21]+0.22360679774997902*fUpwind[8]*alpha[21]+0.25*alpha[5]*fUpwind[18]+0.22360679774997896*alpha[11]*fUpwind[17]+0.25*alpha[0]*fUpwind[17]+0.25*alpha[7]*fUpwind[16]+0.25*fUpwind[7]*alpha[16]+0.25*alpha[9]*fUpwind[15]+0.25*fUpwind[9]*alpha[15]+0.25*alpha[1]*fUpwind[10]+0.25*alpha[3]*fUpwind[8]+0.25*fUpwind[3]*alpha[8]+0.25*alpha[4]*fUpwind[6]+0.25*fUpwind[4]*alpha[6]; 
  Ghat[18] = 0.22360679774997902*alpha[8]*fUpwind[47]+0.22360679774997902*alpha[6]*fUpwind[46]+0.22360679774997902*alpha[5]*fUpwind[45]+0.25000000000000006*alpha[11]*fUpwind[44]+0.2*alpha[26]*fUpwind[43]+0.22360679774997896*alpha[4]*fUpwind[43]+0.22360679774997896*alpha[16]*fUpwind[42]+0.2*alpha[22]*fUpwind[40]+0.22360679774997896*alpha[3]*fUpwind[40]+0.22360679774997896*alpha[15]*fUpwind[39]+0.22360679774997896*alpha[2]*fUpwind[38]+0.22360679774997896*alpha[15]*fUpwind[36]+0.25000000000000006*alpha[21]*fUpwind[35]+0.22360679774997896*alpha[16]*fUpwind[33]+0.25*alpha[1]*fUpwind[31]+0.22360679774997902*alpha[9]*fUpwind[30]+0.22360679774997902*alpha[7]*fUpwind[27]+0.22360679774997902*alpha[7]*fUpwind[26]+0.22360679774997902*fUpwind[7]*alpha[26]+0.22360679774997902*alpha[9]*fUpwind[22]+0.22360679774997902*fUpwind[9]*alpha[22]+0.22360679774997896*alpha[12]*fUpwind[18]+0.25*alpha[0]*fUpwind[18]+0.25*alpha[5]*fUpwind[17]+0.25*alpha[6]*fUpwind[16]+0.25*fUpwind[6]*alpha[16]+0.25*alpha[8]*fUpwind[15]+0.25*fUpwind[8]*alpha[15]+0.25*alpha[2]*fUpwind[10]+0.25*alpha[3]*fUpwind[9]+0.25*fUpwind[3]*alpha[9]+0.25*alpha[4]*fUpwind[7]+0.25*fUpwind[4]*alpha[7]; 
  Ghat[19] = 0.2*alpha[16]*fUpwind[36]+0.22360679774997896*alpha[26]*fUpwind[35]+0.25000000000000006*alpha[4]*fUpwind[35]+0.2*alpha[15]*fUpwind[33]+0.22360679774997896*alpha[22]*fUpwind[32]+0.15971914124998499*alpha[21]*fUpwind[32]+0.25000000000000006*alpha[3]*fUpwind[32]+0.25*alpha[9]*fUpwind[25]+0.25*alpha[7]*fUpwind[21]+0.25*fUpwind[7]*alpha[21]+0.2*alpha[5]*fUpwind[20]+0.22360679774997896*alpha[12]*fUpwind[19]+0.15971914124998499*alpha[11]*fUpwind[19]+0.25*alpha[0]*fUpwind[19]+0.22360679774997902*alpha[8]*fUpwind[16]+0.22360679774997902*fUpwind[8]*alpha[16]+0.22360679774997902*alpha[6]*fUpwind[15]+0.22360679774997902*fUpwind[6]*alpha[15]+0.25000000000000006*alpha[2]*fUpwind[11]+0.25000000000000006*fUpwind[2]*alpha[11]+0.22360679774997902*alpha[1]*fUpwind[5]+0.22360679774997902*fUpwind[1]*alpha[5]; 
  Ghat[20] = 0.15971914124998499*alpha[26]*fUpwind[36]+0.25000000000000006*alpha[4]*fUpwind[36]+0.2*alpha[16]*fUpwind[35]+0.15971914124998499*alpha[22]*fUpwind[33]+0.22360679774997896*alpha[21]*fUpwind[33]+0.25000000000000006*alpha[3]*fUpwind[33]+0.2*alpha[15]*fUpwind[32]+0.25*alpha[8]*fUpwind[26]+0.25*fUpwind[8]*alpha[26]+0.25*alpha[6]*fUpwind[22]+0.25*fUpwind[6]*alpha[22]+0.15971914124998499*alpha[12]*fUpwind[20]+0.22360679774997896*alpha[11]*fUpwind[20]+0.25*alpha[0]*fUpwind[20]+0.2*alpha[5]*fUpwind[19]+0.22360679774997902*alpha[9]*fUpwind[16]+0.22360679774997902*fUpwind[9]*alpha[16]+0.22360679774997902*alpha[7]*fUpwind[15]+0.22360679774997902*fUpwind[7]*alpha[15]+0.25000000000000006*alpha[1]*fUpwind[12]+0.25000000000000006*fUpwind[1]*alpha[12]+0.22360679774997902*alpha[2]*fUpwind[5]+0.22360679774997902*fUpwind[2]*alpha[5]; 
  Ghat[21] = 0.25*alpha[9]*fUpwind[44]+0.25000000000000006*alpha[4]*fUpwind[37]+0.2*alpha[15]*fUpwind[34]+0.25000000000000006*alpha[2]*fUpwind[32]+0.22360679774997902*alpha[16]*fUpwind[31]+0.2*alpha[6]*fUpwind[23]+0.15971914124998499*alpha[11]*fUpwind[21]+0.25*alpha[0]*fUpwind[21]+0.22360679774997896*fUpwind[13]*alpha[21]+0.15971914124998499*fUpwind[11]*alpha[21]+0.25*fUpwind[0]*alpha[21]+0.25*alpha[7]*fUpwind[19]+0.22360679774997902*alpha[8]*fUpwind[17]+0.22360679774997902*alpha[5]*fUpwind[15]+0.22360679774997902*fUpwind[5]*alpha[15]+0.25000000000000006*alpha[3]*fUpwind[11]+0.25000000000000006*fUpwind[3]*alpha[11]+0.22360679774997902*alpha[1]*fUpwind[6]+0.22360679774997902*fUpwind[1]*alpha[6]; 
  Ghat[22] = 0.25*alpha[8]*fUpwind[45]+0.15971914124998499*alpha[26]*fUpwind[38]+0.25000000000000006*alpha[4]*fUpwind[38]+0.2*alpha[15]*fUpwind[34]+0.25000000000000006*alpha[1]*fUpwind[33]+0.22360679774997902*alpha[16]*fUpwind[31]+0.25*fUpwind[10]*alpha[26]+0.2*alpha[7]*fUpwind[24]+0.15971914124998499*alpha[12]*fUpwind[22]+0.25*alpha[0]*fUpwind[22]+0.22360679774997896*fUpwind[13]*alpha[22]+0.15971914124998499*fUpwind[12]*alpha[22]+0.25*fUpwind[0]*alpha[22]+0.25*alpha[6]*fUpwind[20]+0.22360679774997902*alpha[9]*fUpwind[18]+0.22360679774997902*alpha[5]*fUpwind[15]+0.22360679774997902*fUpwind[5]*alpha[15]+0.25000000000000006*alpha[3]*fUpwind[12]+0.25000000000000006*fUpwind[3]*alpha[12]+0.22360679774997902*alpha[2]*fUpwind[7]+0.22360679774997902*fUpwind[2]*alpha[7]; 
  Ghat[23] = 0.25*alpha[9]*fUpwind[46]+0.25000000000000006*alpha[16]*fUpwind[40]+0.25000000000000006*alpha[4]*fUpwind[39]+0.25000000000000006*alpha[2]*fUpwind[34]+0.22360679774997896*alpha[22]*fUpwind[33]+0.2*alpha[15]*fUpwind[32]+0.25*alpha[8]*fUpwind[27]+0.25*alpha[5]*fUpwind[24]+0.22360679774997896*alpha[11]*fUpwind[23]+0.25*alpha[0]*fUpwind[23]+0.2*alpha[6]*fUpwind[21]+0.2*fUpwind[6]*alpha[21]+0.22360679774997902*alpha[7]*fUpwind[15]+0.22360679774997902*fUpwind[7]*alpha[15]+0.25000000000000006*alpha[1]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[6]+0.22360679774997902*fUpwind[3]*alpha[6]; 
  Ghat[24] = 0.25*alpha[8]*fUpwind[46]+0.22360679774997896*alpha[26]*fUpwind[40]+0.25000000000000006*alpha[4]*fUpwind[40]+0.25000000000000006*alpha[16]*fUpwind[39]+0.25000000000000006*alpha[1]*fUpwind[34]+0.2*alpha[15]*fUpwind[33]+0.22360679774997896*alpha[21]*fUpwind[32]+0.25*alpha[9]*fUpwind[27]+0.22360679774997896*alpha[12]*fUpwind[24]+0.25*alpha[0]*fUpwind[24]+0.25*alpha[5]*fUpwind[23]+0.2*alpha[7]*fUpwind[22]+0.2*fUpwind[7]*alpha[22]+0.22360679774997902*alpha[6]*fUpwind[15]+0.22360679774997902*fUpwind[6]*alpha[15]+0.25000000000000006*alpha[2]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[7]+0.22360679774997902*fUpwind[3]*alpha[7]; 
  Ghat[25] = 0.25*alpha[7]*fUpwind[44]+0.2*alpha[16]*fUpwind[41]+0.15971914124998499*alpha[21]*fUpwind[37]+0.25000000000000006*alpha[3]*fUpwind[37]+0.25000000000000006*alpha[2]*fUpwind[35]+0.22360679774997902*alpha[15]*fUpwind[31]+0.2*alpha[8]*fUpwind[28]+0.15971914124998499*alpha[11]*fUpwind[25]+0.25*alpha[0]*fUpwind[25]+0.25*fUpwind[10]*alpha[21]+0.25*alpha[9]*fUpwind[19]+0.22360679774997902*alpha[6]*fUpwind[17]+0.22360679774997902*alpha[5]*fUpwind[16]+0.22360679774997902*fUpwind[5]*alpha[16]+0.25000000000000006*alpha[4]*fUpwind[11]+0.25000000000000006*fUpwind[4]*alpha[11]+0.22360679774997902*alpha[1]*fUpwind[8]+0.22360679774997902*fUpwind[1]*alpha[8]; 
  Ghat[26] = 0.25*alpha[6]*fUpwind[45]+0.2*alpha[16]*fUpwind[41]+0.15971914124998499*alpha[22]*fUpwind[38]+0.25000000000000006*alpha[3]*fUpwind[38]+0.25000000000000006*alpha[1]*fUpwind[36]+0.22360679774997902*alpha[15]*fUpwind[31]+0.2*alpha[9]*fUpwind[29]+0.15971914124998499*alpha[12]*fUpwind[26]+0.25*alpha[0]*fUpwind[26]+0.22360679774997896*fUpwind[14]*alpha[26]+0.15971914124998499*fUpwind[12]*alpha[26]+0.25*fUpwind[0]*alpha[26]+0.25*fUpwind[10]*alpha[22]+0.25*alpha[8]*fUpwind[20]+0.22360679774997902*alpha[7]*fUpwind[18]+0.22360679774997902*alpha[5]*fUpwind[16]+0.22360679774997902*fUpwind[5]*alpha[16]+0.25000000000000006*alpha[4]*fUpwind[12]+0.25000000000000006*fUpwind[4]*alpha[12]+0.22360679774997902*alpha[2]*fUpwind[9]+0.22360679774997902*fUpwind[2]*alpha[9]; 
  Ghat[27] = 0.25*alpha[5]*fUpwind[46]+0.25000000000000006*alpha[2]*fUpwind[40]+0.25000000000000006*alpha[1]*fUpwind[39]+0.22360679774997896*alpha[22]*fUpwind[38]+0.22360679774997896*alpha[21]*fUpwind[37]+0.25000000000000006*alpha[16]*fUpwind[34]+0.22360679774997902*alpha[15]*fUpwind[31]+0.25*alpha[0]*fUpwind[27]+0.25*alpha[9]*fUpwind[24]+0.25*alpha[8]*fUpwind[23]+0.22360679774997902*alpha[7]*fUpwind[18]+0.22360679774997902*alpha[6]*fUpwind[17]+0.25000000000000006*alpha[4]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[10]; 
  Ghat[28] = 0.25*alpha[7]*fUpwind[47]+0.25000000000000006*alpha[15]*fUpwind[43]+0.22360679774997896*alpha[21]*fUpwind[42]+0.25000000000000006*alpha[3]*fUpwind[42]+0.25000000000000006*alpha[2]*fUpwind[41]+0.22360679774997896*alpha[26]*fUpwind[36]+0.2*alpha[16]*fUpwind[35]+0.25*alpha[6]*fUpwind[30]+0.25*alpha[5]*fUpwind[29]+0.22360679774997896*alpha[11]*fUpwind[28]+0.25*alpha[0]*fUpwind[28]+0.2*alpha[8]*fUpwind[25]+0.22360679774997902*alpha[9]*fUpwind[16]+0.22360679774997902*fUpwind[9]*alpha[16]+0.25000000000000006*alpha[1]*fUpwind[14]+0.22360679774997902*alpha[4]*fUpwind[8]+0.22360679774997902*fUpwind[4]*alpha[8]; 
  Ghat[29] = 0.25*alpha[6]*fUpwind[47]+0.22360679774997896*alpha[22]*fUpwind[43]+0.25000000000000006*alpha[3]*fUpwind[43]+0.25000000000000006*alpha[15]*fUpwind[42]+0.25000000000000006*alpha[1]*fUpwind[41]+0.2*alpha[16]*fUpwind[36]+0.25*alpha[7]*fUpwind[30]+0.22360679774997896*alpha[12]*fUpwind[29]+0.25*alpha[0]*fUpwind[29]+0.25*alpha[5]*fUpwind[28]+0.2*alpha[9]*fUpwind[26]+0.2*fUpwind[9]*alpha[26]+0.22360679774997902*alpha[8]*fUpwind[16]+0.22360679774997902*fUpwind[8]*alpha[16]+0.25000000000000006*alpha[2]*fUpwind[14]+0.22360679774997902*alpha[4]*fUpwind[9]+0.22360679774997902*fUpwind[4]*alpha[9]; 
  Ghat[30] = 0.25*alpha[5]*fUpwind[47]+0.25000000000000006*alpha[2]*fUpwind[43]+0.25000000000000006*alpha[1]*fUpwind[42]+0.25000000000000006*alpha[15]*fUpwind[41]+0.22360679774997896*alpha[26]*fUpwind[38]+0.22360679774997902*alpha[16]*fUpwind[31]+0.25*alpha[0]*fUpwind[30]+0.25*alpha[7]*fUpwind[29]+0.25*alpha[6]*fUpwind[28]+0.22360679774997902*alpha[9]*fUpwind[18]+0.22360679774997902*alpha[8]*fUpwind[17]+0.25000000000000006*alpha[3]*fUpwind[14]+0.22360679774997902*alpha[4]*fUpwind[10]; 
  Ghat[31] = 0.2*alpha[26]*fUpwind[47]+0.22360679774997902*alpha[4]*fUpwind[47]+0.2*alpha[22]*fUpwind[46]+0.2*alpha[21]*fUpwind[46]+0.22360679774997902*alpha[3]*fUpwind[46]+0.22360679774997902*alpha[2]*fUpwind[45]+0.22360679774997902*alpha[1]*fUpwind[44]+0.22360679774997896*alpha[8]*fUpwind[43]+0.22360679774997896*alpha[9]*fUpwind[42]+0.22360679774997896*alpha[6]*fUpwind[40]+0.22360679774997896*alpha[7]*fUpwind[39]+0.22360679774997896*alpha[5]*fUpwind[38]+0.22360679774997896*alpha[5]*fUpwind[37]+0.22360679774997896*alpha[7]*fUpwind[36]+0.22360679774997896*alpha[6]*fUpwind[35]+0.22360679774997896*alpha[9]*fUpwind[33]+0.22360679774997896*alpha[8]*fUpwind[32]+0.22360679774997896*alpha[12]*fUpwind[31]+0.22360679774997896*alpha[11]*fUpwind[31]+0.25*alpha[0]*fUpwind[31]+0.22360679774997902*alpha[16]*fUpwind[30]+0.22360679774997902*alpha[15]*fUpwind[27]+0.22360679774997902*alpha[15]*fUpwind[26]+0.22360679774997902*fUpwind[15]*alpha[26]+0.22360679774997902*alpha[15]*fUpwind[25]+0.22360679774997902*alpha[16]*fUpwind[22]+0.22360679774997902*fUpwind[16]*alpha[22]+0.22360679774997902*alpha[16]*fUpwind[21]+0.22360679774997902*fUpwind[16]*alpha[21]+0.25*alpha[1]*fUpwind[18]+0.25*alpha[2]*fUpwind[17]+0.25*alpha[3]*fUpwind[16]+0.25*fUpwind[3]*alpha[16]+0.25*alpha[4]*fUpwind[15]+0.25*fUpwind[4]*alpha[15]+0.25*alpha[5]*fUpwind[10]+0.25*alpha[6]*fUpwind[9]+0.25*fUpwind[6]*alpha[9]+0.25*alpha[7]*fUpwind[8]+0.25*fUpwind[7]*alpha[8]; 
  Ghat[32] = 0.2*alpha[16]*fUpwind[45]+0.22360679774997896*alpha[26]*fUpwind[44]+0.25000000000000006*alpha[4]*fUpwind[44]+0.25*alpha[9]*fUpwind[37]+0.2*alpha[6]*fUpwind[34]+0.2*alpha[5]*fUpwind[33]+0.22360679774997896*alpha[12]*fUpwind[32]+0.15971914124998499*alpha[11]*fUpwind[32]+0.25*alpha[0]*fUpwind[32]+0.22360679774997896*alpha[8]*fUpwind[31]+0.22360679774997896*alpha[21]*fUpwind[24]+0.2*alpha[15]*fUpwind[23]+0.22360679774997896*fUpwind[19]*alpha[22]+0.25000000000000006*alpha[2]*fUpwind[21]+0.15971914124998499*fUpwind[19]*alpha[21]+0.25000000000000006*fUpwind[2]*alpha[21]+0.2*alpha[15]*fUpwind[20]+0.25000000000000006*alpha[3]*fUpwind[19]+0.22360679774997896*alpha[16]*fUpwind[17]+0.22360679774997896*alpha[1]*fUpwind[15]+0.22360679774997896*fUpwind[1]*alpha[15]+0.25*alpha[7]*fUpwind[11]+0.25*fUpwind[7]*alpha[11]+0.22360679774997896*alpha[5]*fUpwind[6]+0.22360679774997896*fUpwind[5]*alpha[6]; 
  Ghat[33] = 0.15971914124998499*alpha[26]*fUpwind[45]+0.25000000000000006*alpha[4]*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+0.25*alpha[8]*fUpwind[38]+0.2*alpha[7]*fUpwind[34]+0.15971914124998499*alpha[12]*fUpwind[33]+0.22360679774997896*alpha[11]*fUpwind[33]+0.25*alpha[0]*fUpwind[33]+0.2*alpha[5]*fUpwind[32]+0.22360679774997896*alpha[9]*fUpwind[31]+0.25000000000000006*fUpwind[17]*alpha[26]+0.2*alpha[15]*fUpwind[24]+0.22360679774997896*alpha[22]*fUpwind[23]+0.25000000000000006*alpha[1]*fUpwind[22]+0.15971914124998499*fUpwind[20]*alpha[22]+0.25000000000000006*fUpwind[1]*alpha[22]+0.22360679774997896*fUpwind[20]*alpha[21]+0.25000000000000006*alpha[3]*fUpwind[20]+0.2*alpha[15]*fUpwind[19]+0.22360679774997896*alpha[16]*fUpwind[18]+0.22360679774997896*alpha[2]*fUpwind[15]+0.22360679774997896*fUpwind[2]*alpha[15]+0.25*alpha[6]*fUpwind[12]+0.25*fUpwind[6]*alpha[12]+0.22360679774997896*alpha[5]*fUpwind[7]+0.22360679774997896*fUpwind[5]*alpha[7]; 
  Ghat[34] = 0.22360679774997896*alpha[26]*fUpwind[46]+0.25000000000000006*alpha[4]*fUpwind[46]+0.25*alpha[8]*fUpwind[40]+0.25*alpha[9]*fUpwind[39]+0.22360679774997896*alpha[12]*fUpwind[34]+0.22360679774997896*alpha[11]*fUpwind[34]+0.25*alpha[0]*fUpwind[34]+0.2*alpha[7]*fUpwind[33]+0.2*alpha[6]*fUpwind[32]+0.25000000000000006*alpha[16]*fUpwind[27]+0.25000000000000006*alpha[1]*fUpwind[24]+0.25000000000000006*alpha[2]*fUpwind[23]+0.2*alpha[15]*fUpwind[22]+0.2*fUpwind[15]*alpha[22]+0.2*alpha[15]*fUpwind[21]+0.2*fUpwind[15]*alpha[21]+0.22360679774997896*alpha[3]*fUpwind[15]+0.22360679774997896*fUpwind[3]*alpha[15]+0.25*alpha[5]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[7]+0.22360679774997896*fUpwind[6]*alpha[7]; 
  Ghat[35] = 0.2*alpha[15]*fUpwind[45]+0.22360679774997896*alpha[22]*fUpwind[44]+0.15971914124998499*alpha[21]*fUpwind[44]+0.25000000000000006*alpha[3]*fUpwind[44]+0.2*alpha[8]*fUpwind[41]+0.25*alpha[7]*fUpwind[37]+0.2*alpha[5]*fUpwind[36]+0.22360679774997896*alpha[12]*fUpwind[35]+0.15971914124998499*alpha[11]*fUpwind[35]+0.25*alpha[0]*fUpwind[35]+0.22360679774997896*alpha[6]*fUpwind[31]+0.2*alpha[16]*fUpwind[28]+0.22360679774997896*fUpwind[19]*alpha[26]+0.25000000000000006*alpha[2]*fUpwind[25]+0.25000000000000006*fUpwind[18]*alpha[21]+0.2*alpha[16]*fUpwind[20]+0.25000000000000006*alpha[4]*fUpwind[19]+0.22360679774997896*alpha[15]*fUpwind[17]+0.22360679774997896*alpha[1]*fUpwind[16]+0.22360679774997896*fUpwind[1]*alpha[16]+0.25*alpha[9]*fUpwind[11]+0.25*fUpwind[9]*alpha[11]+0.22360679774997896*alpha[5]*fUpwind[8]+0.22360679774997896*fUpwind[5]*alpha[8]; 
  Ghat[36] = 0.15971914124998499*alpha[22]*fUpwind[45]+0.22360679774997896*alpha[21]*fUpwind[45]+0.25000000000000006*alpha[3]*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+0.2*alpha[9]*fUpwind[41]+0.25*alpha[6]*fUpwind[38]+0.15971914124998499*alpha[12]*fUpwind[36]+0.22360679774997896*alpha[11]*fUpwind[36]+0.25*alpha[0]*fUpwind[36]+0.2*alpha[5]*fUpwind[35]+0.22360679774997896*alpha[7]*fUpwind[31]+0.2*alpha[16]*fUpwind[29]+0.22360679774997896*alpha[26]*fUpwind[28]+0.25000000000000006*alpha[1]*fUpwind[26]+0.15971914124998499*fUpwind[20]*alpha[26]+0.25000000000000006*fUpwind[1]*alpha[26]+0.25000000000000006*fUpwind[17]*alpha[22]+0.25000000000000006*alpha[4]*fUpwind[20]+0.2*alpha[16]*fUpwind[19]+0.22360679774997896*alpha[15]*fUpwind[18]+0.22360679774997896*alpha[2]*fUpwind[16]+0.22360679774997896*fUpwind[2]*alpha[16]+0.25*alpha[8]*fUpwind[12]+0.25*fUpwind[8]*alpha[12]+0.22360679774997896*alpha[5]*fUpwind[9]+0.22360679774997896*fUpwind[5]*alpha[9]; 
  Ghat[37] = 0.2*alpha[16]*fUpwind[47]+0.2*alpha[15]*fUpwind[46]+0.25000000000000006*alpha[2]*fUpwind[44]+0.2*alpha[8]*fUpwind[42]+0.2*alpha[6]*fUpwind[39]+0.15971914124998499*alpha[11]*fUpwind[37]+0.25*alpha[0]*fUpwind[37]+0.25*alpha[7]*fUpwind[35]+0.25*alpha[9]*fUpwind[32]+0.22360679774997896*alpha[5]*fUpwind[31]+0.22360679774997896*alpha[21]*fUpwind[27]+0.15971914124998499*alpha[21]*fUpwind[25]+0.25000000000000006*alpha[3]*fUpwind[25]+0.25000000000000006*alpha[4]*fUpwind[21]+0.25000000000000006*fUpwind[4]*alpha[21]+0.22360679774997896*alpha[1]*fUpwind[17]+0.22360679774997896*alpha[15]*fUpwind[16]+0.22360679774997896*fUpwind[15]*alpha[16]+0.25*fUpwind[10]*alpha[11]+0.22360679774997896*alpha[6]*fUpwind[8]+0.22360679774997896*fUpwind[6]*alpha[8]; 
  Ghat[38] = 0.2*alpha[16]*fUpwind[47]+0.2*alpha[15]*fUpwind[46]+0.25000000000000006*alpha[1]*fUpwind[45]+0.2*alpha[9]*fUpwind[43]+0.2*alpha[7]*fUpwind[40]+0.15971914124998499*alpha[12]*fUpwind[38]+0.25*alpha[0]*fUpwind[38]+0.25*alpha[6]*fUpwind[36]+0.25*alpha[8]*fUpwind[33]+0.22360679774997896*alpha[5]*fUpwind[31]+0.22360679774997896*alpha[26]*fUpwind[30]+0.22360679774997896*alpha[22]*fUpwind[27]+0.15971914124998499*alpha[22]*fUpwind[26]+0.25000000000000006*alpha[3]*fUpwind[26]+0.15971914124998499*fUpwind[22]*alpha[26]+0.25000000000000006*fUpwind[3]*alpha[26]+0.25000000000000006*alpha[4]*fUpwind[22]+0.25000000000000006*fUpwind[4]*alpha[22]+0.22360679774997896*alpha[2]*fUpwind[18]+0.22360679774997896*alpha[15]*fUpwind[16]+0.22360679774997896*fUpwind[15]*alpha[16]+0.25*fUpwind[10]*alpha[12]+0.22360679774997896*alpha[7]*fUpwind[9]+0.22360679774997896*fUpwind[7]*alpha[9]; 
  Ghat[39] = 0.25000000000000006*alpha[2]*fUpwind[46]+0.22360679774997896*alpha[22]*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+0.25*alpha[5]*fUpwind[40]+0.22360679774997896*alpha[11]*fUpwind[39]+0.25*alpha[0]*fUpwind[39]+0.2*alpha[6]*fUpwind[37]+0.25*alpha[9]*fUpwind[34]+0.22360679774997896*alpha[7]*fUpwind[31]+0.25000000000000006*alpha[1]*fUpwind[27]+0.25000000000000006*alpha[16]*fUpwind[24]+0.25000000000000006*alpha[4]*fUpwind[23]+0.2*fUpwind[17]*alpha[21]+0.22360679774997896*alpha[15]*fUpwind[18]+0.22360679774997896*alpha[3]*fUpwind[17]+0.25*alpha[8]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[10]; 
  Ghat[40] = 0.25000000000000006*alpha[1]*fUpwind[46]+0.2*alpha[15]*fUpwind[45]+0.22360679774997896*alpha[21]*fUpwind[44]+0.22360679774997896*alpha[12]*fUpwind[40]+0.25*alpha[0]*fUpwind[40]+0.25*alpha[5]*fUpwind[39]+0.2*alpha[7]*fUpwind[38]+0.25*alpha[8]*fUpwind[34]+0.22360679774997896*alpha[6]*fUpwind[31]+0.25000000000000006*alpha[2]*fUpwind[27]+0.22360679774997896*fUpwind[24]*alpha[26]+0.25000000000000006*alpha[4]*fUpwind[24]+0.25000000000000006*alpha[16]*fUpwind[23]+0.2*fUpwind[18]*alpha[22]+0.22360679774997896*alpha[3]*fUpwind[18]+0.22360679774997896*alpha[15]*fUpwind[17]+0.25*alpha[9]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[10]; 
  Ghat[41] = 0.22360679774997896*alpha[22]*fUpwind[47]+0.22360679774997896*alpha[21]*fUpwind[47]+0.25000000000000006*alpha[3]*fUpwind[47]+0.25*alpha[6]*fUpwind[43]+0.25*alpha[7]*fUpwind[42]+0.22360679774997896*alpha[12]*fUpwind[41]+0.22360679774997896*alpha[11]*fUpwind[41]+0.25*alpha[0]*fUpwind[41]+0.2*alpha[9]*fUpwind[36]+0.2*alpha[8]*fUpwind[35]+0.25000000000000006*alpha[15]*fUpwind[30]+0.25000000000000006*alpha[1]*fUpwind[29]+0.25000000000000006*alpha[2]*fUpwind[28]+0.2*alpha[16]*fUpwind[26]+0.2*fUpwind[16]*alpha[26]+0.2*alpha[16]*fUpwind[25]+0.22360679774997896*alpha[4]*fUpwind[16]+0.22360679774997896*fUpwind[4]*alpha[16]+0.25*alpha[5]*fUpwind[14]+0.22360679774997896*alpha[8]*fUpwind[9]+0.22360679774997896*fUpwind[8]*alpha[9]; 
  Ghat[42] = 0.25000000000000006*alpha[2]*fUpwind[47]+0.22360679774997896*alpha[26]*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+0.25*alpha[5]*fUpwind[43]+0.22360679774997896*alpha[11]*fUpwind[42]+0.25*alpha[0]*fUpwind[42]+0.25*alpha[7]*fUpwind[41]+0.2*alpha[8]*fUpwind[37]+0.22360679774997896*alpha[9]*fUpwind[31]+0.25000000000000006*alpha[1]*fUpwind[30]+0.25000000000000006*alpha[15]*fUpwind[29]+0.22360679774997896*alpha[21]*fUpwind[28]+0.25000000000000006*alpha[3]*fUpwind[28]+0.22360679774997896*alpha[16]*fUpwind[18]+0.22360679774997896*alpha[4]*fUpwind[17]+0.25*alpha[6]*fUpwind[14]+0.22360679774997896*alpha[8]*fUpwind[10]; 
  Ghat[43] = 0.25000000000000006*alpha[1]*fUpwind[47]+0.2*alpha[16]*fUpwind[45]+0.22360679774997896*alpha[12]*fUpwind[43]+0.25*alpha[0]*fUpwind[43]+0.25*alpha[5]*fUpwind[42]+0.25*alpha[6]*fUpwind[41]+0.2*alpha[9]*fUpwind[38]+0.22360679774997896*alpha[8]*fUpwind[31]+0.25000000000000006*alpha[2]*fUpwind[30]+0.22360679774997896*alpha[22]*fUpwind[29]+0.25000000000000006*alpha[3]*fUpwind[29]+0.25000000000000006*alpha[15]*fUpwind[28]+0.2*fUpwind[18]*alpha[26]+0.22360679774997896*alpha[4]*fUpwind[18]+0.22360679774997896*alpha[16]*fUpwind[17]+0.25*alpha[7]*fUpwind[14]+0.22360679774997896*alpha[9]*fUpwind[10]; 
  Ghat[44] = 0.2*alpha[8]*fUpwind[47]+0.2*alpha[6]*fUpwind[46]+0.2*alpha[5]*fUpwind[45]+0.22360679774997896*alpha[12]*fUpwind[44]+0.15971914124998499*alpha[11]*fUpwind[44]+0.25*alpha[0]*fUpwind[44]+0.2*alpha[16]*fUpwind[42]+0.22360679774997896*alpha[21]*fUpwind[40]+0.2*alpha[15]*fUpwind[39]+0.25000000000000006*alpha[2]*fUpwind[37]+0.2*alpha[15]*fUpwind[36]+0.22360679774997896*alpha[22]*fUpwind[35]+0.15971914124998499*alpha[21]*fUpwind[35]+0.25000000000000006*alpha[3]*fUpwind[35]+0.2*alpha[16]*fUpwind[33]+0.22360679774997896*alpha[26]*fUpwind[32]+0.25000000000000006*alpha[4]*fUpwind[32]+0.22360679774997902*alpha[1]*fUpwind[31]+0.25*alpha[7]*fUpwind[25]+0.25*alpha[9]*fUpwind[21]+0.25*fUpwind[9]*alpha[21]+0.25000000000000006*alpha[11]*fUpwind[18]+0.22360679774997902*alpha[5]*fUpwind[17]+0.22360679774997902*alpha[6]*fUpwind[16]+0.22360679774997902*fUpwind[6]*alpha[16]+0.22360679774997902*alpha[8]*fUpwind[15]+0.22360679774997902*fUpwind[8]*alpha[15]; 
  Ghat[45] = 0.2*alpha[9]*fUpwind[47]+0.2*alpha[7]*fUpwind[46]+0.15971914124998499*alpha[12]*fUpwind[45]+0.22360679774997896*alpha[11]*fUpwind[45]+0.25*alpha[0]*fUpwind[45]+0.2*alpha[5]*fUpwind[44]+0.2*alpha[16]*fUpwind[43]+0.22360679774997896*alpha[26]*fUpwind[42]+0.2*alpha[15]*fUpwind[40]+0.22360679774997896*alpha[22]*fUpwind[39]+0.25000000000000006*alpha[1]*fUpwind[38]+0.15971914124998499*alpha[22]*fUpwind[36]+0.22360679774997896*alpha[21]*fUpwind[36]+0.25000000000000006*alpha[3]*fUpwind[36]+0.2*alpha[15]*fUpwind[35]+0.15971914124998499*alpha[26]*fUpwind[33]+0.25000000000000006*alpha[4]*fUpwind[33]+0.2*alpha[16]*fUpwind[32]+0.22360679774997902*alpha[2]*fUpwind[31]+0.25*alpha[6]*fUpwind[26]+0.25*fUpwind[6]*alpha[26]+0.25*alpha[8]*fUpwind[22]+0.25*fUpwind[8]*alpha[22]+0.22360679774997902*alpha[5]*fUpwind[18]+0.25000000000000006*alpha[12]*fUpwind[17]+0.22360679774997902*alpha[7]*fUpwind[16]+0.22360679774997902*fUpwind[7]*alpha[16]+0.22360679774997902*alpha[9]*fUpwind[15]+0.22360679774997902*fUpwind[9]*alpha[15]; 
  Ghat[46] = 0.22360679774997896*alpha[12]*fUpwind[46]+0.22360679774997896*alpha[11]*fUpwind[46]+0.25*alpha[0]*fUpwind[46]+0.2*alpha[7]*fUpwind[45]+0.2*alpha[6]*fUpwind[44]+0.25000000000000006*alpha[1]*fUpwind[40]+0.25000000000000006*alpha[2]*fUpwind[39]+0.2*alpha[15]*fUpwind[38]+0.2*alpha[15]*fUpwind[37]+0.22360679774997896*alpha[26]*fUpwind[34]+0.25000000000000006*alpha[4]*fUpwind[34]+0.2*alpha[22]*fUpwind[31]+0.2*alpha[21]*fUpwind[31]+0.22360679774997902*alpha[3]*fUpwind[31]+0.25*alpha[5]*fUpwind[27]+0.25*alpha[8]*fUpwind[24]+0.25*alpha[9]*fUpwind[23]+0.22360679774997902*alpha[6]*fUpwind[18]+0.22360679774997902*alpha[7]*fUpwind[17]+0.25000000000000006*fUpwind[13]*alpha[16]+0.22360679774997902*fUpwind[10]*alpha[15]; 
  Ghat[47] = 0.22360679774997896*alpha[12]*fUpwind[47]+0.22360679774997896*alpha[11]*fUpwind[47]+0.25*alpha[0]*fUpwind[47]+0.2*alpha[9]*fUpwind[45]+0.2*alpha[8]*fUpwind[44]+0.25000000000000006*alpha[1]*fUpwind[43]+0.25000000000000006*alpha[2]*fUpwind[42]+0.22360679774997896*alpha[22]*fUpwind[41]+0.22360679774997896*alpha[21]*fUpwind[41]+0.25000000000000006*alpha[3]*fUpwind[41]+0.2*alpha[16]*fUpwind[38]+0.2*alpha[16]*fUpwind[37]+0.2*alpha[26]*fUpwind[31]+0.22360679774997902*alpha[4]*fUpwind[31]+0.25*alpha[5]*fUpwind[30]+0.25*alpha[6]*fUpwind[29]+0.25*alpha[7]*fUpwind[28]+0.22360679774997902*alpha[8]*fUpwind[18]+0.22360679774997902*alpha[9]*fUpwind[17]+0.22360679774997902*fUpwind[10]*alpha[16]+0.25000000000000006*fUpwind[14]*alpha[15]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv10); 
  out[3] += -(1.224744871391589*Ghat[0]*dv10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv10); 
  out[6] += -(0.7071067811865475*Ghat[5]*dv10); 
  out[7] += -(1.224744871391589*Ghat[1]*dv10); 
  out[8] += -(1.224744871391589*Ghat[2]*dv10); 
  out[9] += -(0.7071067811865475*Ghat[6]*dv10); 
  out[10] += -(0.7071067811865475*Ghat[7]*dv10); 
  out[11] += -(1.224744871391589*Ghat[3]*dv10); 
  out[12] += -(0.7071067811865475*Ghat[8]*dv10); 
  out[13] += -(0.7071067811865475*Ghat[9]*dv10); 
  out[14] += -(1.224744871391589*Ghat[4]*dv10); 
  out[15] += -(0.7071067811865475*Ghat[10]*dv10); 
  out[16] += -(0.7071067811865475*Ghat[11]*dv10); 
  out[17] += -(0.7071067811865475*Ghat[12]*dv10); 
  out[18] += -(1.5811388300841895*Ghat[0]*dv10); 
  out[19] += -(0.7071067811865475*Ghat[13]*dv10); 
  out[20] += -(0.7071067811865475*Ghat[14]*dv10); 
  out[21] += -(1.224744871391589*Ghat[5]*dv10); 
  out[22] += -(0.7071067811865475*Ghat[15]*dv10); 
  out[23] += -(1.224744871391589*Ghat[6]*dv10); 
  out[24] += -(1.224744871391589*Ghat[7]*dv10); 
  out[25] += -(0.7071067811865475*Ghat[16]*dv10); 
  out[26] += -(1.224744871391589*Ghat[8]*dv10); 
  out[27] += -(1.224744871391589*Ghat[9]*dv10); 
  out[28] += -(0.7071067811865475*Ghat[17]*dv10); 
  out[29] += -(0.7071067811865475*Ghat[18]*dv10); 
  out[30] += -(1.224744871391589*Ghat[10]*dv10); 
  out[31] += -(0.7071067811865475*Ghat[19]*dv10); 
  out[32] += -(0.7071067811865475*Ghat[20]*dv10); 
  out[33] += -(1.224744871391589*Ghat[11]*dv10); 
  out[34] += -(1.224744871391589*Ghat[12]*dv10); 
  out[35] += -(1.5811388300841898*Ghat[1]*dv10); 
  out[36] += -(1.5811388300841898*Ghat[2]*dv10); 
  out[37] += -(0.7071067811865475*Ghat[21]*dv10); 
  out[38] += -(0.7071067811865475*Ghat[22]*dv10); 
  out[39] += -(1.5811388300841898*Ghat[3]*dv10); 
  out[40] += -(0.7071067811865475*Ghat[23]*dv10); 
  out[41] += -(0.7071067811865475*Ghat[24]*dv10); 
  out[42] += -(1.224744871391589*Ghat[13]*dv10); 
  out[43] += -(0.7071067811865475*Ghat[25]*dv10); 
  out[44] += -(0.7071067811865475*Ghat[26]*dv10); 
  out[45] += -(1.5811388300841898*Ghat[4]*dv10); 
  out[46] += -(0.7071067811865475*Ghat[27]*dv10); 
  out[47] += -(0.7071067811865475*Ghat[28]*dv10); 
  out[48] += -(0.7071067811865475*Ghat[29]*dv10); 
  out[49] += -(1.224744871391589*Ghat[14]*dv10); 
  out[50] += -(0.7071067811865475*Ghat[30]*dv10); 
  out[51] += -(1.224744871391589*Ghat[15]*dv10); 
  out[52] += -(1.224744871391589*Ghat[16]*dv10); 
  out[53] += -(0.7071067811865475*Ghat[31]*dv10); 
  out[54] += -(1.224744871391589*Ghat[17]*dv10); 
  out[55] += -(1.224744871391589*Ghat[18]*dv10); 
  out[56] += -(1.224744871391589*Ghat[19]*dv10); 
  out[57] += -(1.224744871391589*Ghat[20]*dv10); 
  out[58] += -(1.5811388300841895*Ghat[5]*dv10); 
  out[59] += -(0.7071067811865475*Ghat[32]*dv10); 
  out[60] += -(0.7071067811865475*Ghat[33]*dv10); 
  out[61] += -(1.224744871391589*Ghat[21]*dv10); 
  out[62] += -(1.224744871391589*Ghat[22]*dv10); 
  out[63] += -(1.5811388300841895*Ghat[6]*dv10); 
  out[64] += -(1.5811388300841895*Ghat[7]*dv10); 
  out[65] += -(0.7071067811865475*Ghat[34]*dv10); 
  out[66] += -(1.224744871391589*Ghat[23]*dv10); 
  out[67] += -(1.224744871391589*Ghat[24]*dv10); 
  out[68] += -(0.7071067811865475*Ghat[35]*dv10); 
  out[69] += -(0.7071067811865475*Ghat[36]*dv10); 
  out[70] += -(1.224744871391589*Ghat[25]*dv10); 
  out[71] += -(1.224744871391589*Ghat[26]*dv10); 
  out[72] += -(1.5811388300841895*Ghat[8]*dv10); 
  out[73] += -(1.5811388300841895*Ghat[9]*dv10); 
  out[74] += -(0.7071067811865475*Ghat[37]*dv10); 
  out[75] += -(0.7071067811865475*Ghat[38]*dv10); 
  out[76] += -(1.5811388300841895*Ghat[10]*dv10); 
  out[77] += -(0.7071067811865475*Ghat[39]*dv10); 
  out[78] += -(0.7071067811865475*Ghat[40]*dv10); 
  out[79] += -(1.224744871391589*Ghat[27]*dv10); 
  out[80] += -(0.7071067811865475*Ghat[41]*dv10); 
  out[81] += -(1.224744871391589*Ghat[28]*dv10); 
  out[82] += -(1.224744871391589*Ghat[29]*dv10); 
  out[83] += -(0.7071067811865475*Ghat[42]*dv10); 
  out[84] += -(0.7071067811865475*Ghat[43]*dv10); 
  out[85] += -(1.224744871391589*Ghat[30]*dv10); 
  out[86] += -(1.224744871391589*Ghat[31]*dv10); 
  out[87] += -(1.224744871391589*Ghat[32]*dv10); 
  out[88] += -(1.224744871391589*Ghat[33]*dv10); 
  out[89] += -(1.5811388300841898*Ghat[15]*dv10); 
  out[90] += -(1.224744871391589*Ghat[34]*dv10); 
  out[91] += -(1.224744871391589*Ghat[35]*dv10); 
  out[92] += -(1.224744871391589*Ghat[36]*dv10); 
  out[93] += -(1.5811388300841898*Ghat[16]*dv10); 
  out[94] += -(0.7071067811865475*Ghat[44]*dv10); 
  out[95] += -(0.7071067811865475*Ghat[45]*dv10); 
  out[96] += -(1.224744871391589*Ghat[37]*dv10); 
  out[97] += -(1.224744871391589*Ghat[38]*dv10); 
  out[98] += -(1.5811388300841898*Ghat[17]*dv10); 
  out[99] += -(1.5811388300841898*Ghat[18]*dv10); 
  out[100] += -(0.7071067811865475*Ghat[46]*dv10); 
  out[101] += -(1.224744871391589*Ghat[39]*dv10); 
  out[102] += -(1.224744871391589*Ghat[40]*dv10); 
  out[103] += -(1.224744871391589*Ghat[41]*dv10); 
  out[104] += -(0.7071067811865475*Ghat[47]*dv10); 
  out[105] += -(1.224744871391589*Ghat[42]*dv10); 
  out[106] += -(1.224744871391589*Ghat[43]*dv10); 
  out[107] += -(1.224744871391589*Ghat[44]*dv10); 
  out[108] += -(1.224744871391589*Ghat[45]*dv10); 
  out[109] += -(1.5811388300841895*Ghat[31]*dv10); 
  out[110] += -(1.224744871391589*Ghat[46]*dv10); 
  out[111] += -(1.224744871391589*Ghat[47]*dv10); 

  } else { 

  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p2_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p2_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[5]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_5x_p2_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_5x_p2_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.33541019662496846*alpha[3]-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_5x_p2_surfx3_eval_quad_node_7_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_5x_p2_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if (0.375*(alpha[26]+alpha[22])-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*(alpha[8]+alpha[6])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx3_eval_quad_node_9_l(fSkin); 
  } 
  if (0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[6]-0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_5x_p2_surfx3_eval_quad_node_10_l(fSkin); 
  } 
  if (-(0.375*alpha[26])+0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]+0.45*alpha[6]+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_5x_p2_surfx3_eval_quad_node_11_l(fSkin); 
  } 
  if (0.375*alpha[26]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]-0.33541019662496846*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_5x_p2_surfx3_eval_quad_node_12_l(fSkin); 
  } 
  if (-(0.2795084971874737*alpha[12])+0.22360679774997896*alpha[11]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_5x_p2_surfx3_eval_quad_node_13_l(fSkin); 
  } 
  if (-(0.375*alpha[26])-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_5x_p2_surfx3_eval_quad_node_14_l(fSkin); 
  } 
  if (0.375*alpha[26]-0.375*alpha[22]+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]-0.45*alpha[6]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_5x_p2_surfx3_eval_quad_node_15_l(fSkin); 
  } 
  if (-(0.375*alpha[22])+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[6]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_5x_p2_surfx3_eval_quad_node_16_l(fSkin); 
  } 
  if (-(0.375*(alpha[26]+alpha[22]))+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*(alpha[8]+alpha[6])+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = ser_5x_p2_surfx3_eval_quad_node_17_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx3_eval_quad_node_18_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = ser_5x_p2_surfx3_eval_quad_node_19_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_5x_p2_surfx3_eval_quad_node_20_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_5x_p2_surfx3_eval_quad_node_21_l(fSkin); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[5]+0.33541019662496846*alpha[2]-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = ser_5x_p2_surfx3_eval_quad_node_22_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = ser_5x_p2_surfx3_eval_quad_node_23_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_5x_p2_surfx3_eval_quad_node_24_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.33541019662496846*(alpha[3]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = ser_5x_p2_surfx3_eval_quad_node_25_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])-0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = ser_5x_p2_surfx3_eval_quad_node_26_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]))+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_r(fEdge); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx3_eval_quad_node_27_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[22])+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]-0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_r(fEdge); 
  } else { 
    fUpwindQuad[28] = ser_5x_p2_surfx3_eval_quad_node_28_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*alpha[22]+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_r(fEdge); 
  } else { 
    fUpwindQuad[29] = ser_5x_p2_surfx3_eval_quad_node_29_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.33541019662496846*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_r(fEdge); 
  } else { 
    fUpwindQuad[30] = ser_5x_p2_surfx3_eval_quad_node_30_l(fSkin); 
  } 
  if (0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_r(fEdge); 
  } else { 
    fUpwindQuad[31] = ser_5x_p2_surfx3_eval_quad_node_31_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_r(fEdge); 
  } else { 
    fUpwindQuad[32] = ser_5x_p2_surfx3_eval_quad_node_32_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_r(fEdge); 
  } else { 
    fUpwindQuad[33] = ser_5x_p2_surfx3_eval_quad_node_33_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_r(fEdge); 
  } else { 
    fUpwindQuad[34] = ser_5x_p2_surfx3_eval_quad_node_34_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22])-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_r(fEdge); 
  } else { 
    fUpwindQuad[35] = ser_5x_p2_surfx3_eval_quad_node_35_l(fSkin); 
  } 
  if (0.375*(alpha[26]+alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_r(fEdge); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx3_eval_quad_node_36_l(fSkin); 
  } 
  if (0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_r(fEdge); 
  } else { 
    fUpwindQuad[37] = ser_5x_p2_surfx3_eval_quad_node_37_l(fSkin); 
  } 
  if (-(0.375*alpha[26])+0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_r(fEdge); 
  } else { 
    fUpwindQuad[38] = ser_5x_p2_surfx3_eval_quad_node_38_l(fSkin); 
  } 
  if (0.375*alpha[26]-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_r(fEdge); 
  } else { 
    fUpwindQuad[39] = ser_5x_p2_surfx3_eval_quad_node_39_l(fSkin); 
  } 
  if (0.25*alpha[0]-0.2795084971874737*(alpha[12]+alpha[11]) > 0) { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_r(fEdge); 
  } else { 
    fUpwindQuad[40] = ser_5x_p2_surfx3_eval_quad_node_40_l(fSkin); 
  } 
  if (-(0.375*alpha[26])-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*alpha[4]+0.25*alpha[0] > 0) { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_r(fEdge); 
  } else { 
    fUpwindQuad[41] = ser_5x_p2_surfx3_eval_quad_node_41_l(fSkin); 
  } 
  if (0.375*alpha[26]-0.375*(alpha[22]+alpha[21])-0.2795084971874737*(alpha[12]+alpha[11])-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_r(fEdge); 
  } else { 
    fUpwindQuad[42] = ser_5x_p2_surfx3_eval_quad_node_42_l(fSkin); 
  } 
  if (-(0.375*(alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*alpha[3]+0.25*alpha[0] > 0) { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_r(fEdge); 
  } else { 
    fUpwindQuad[43] = ser_5x_p2_surfx3_eval_quad_node_43_l(fSkin); 
  } 
  if (-(0.375*(alpha[26]+alpha[22]+alpha[21]))-0.2795084971874737*(alpha[12]+alpha[11])+0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[0] > 0) { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_r(fEdge); 
  } else { 
    fUpwindQuad[44] = ser_5x_p2_surfx3_eval_quad_node_44_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]))+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*(alpha[9]+alpha[7])-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_r(fEdge); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx3_eval_quad_node_45_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[22])+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[7]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_r(fEdge); 
  } else { 
    fUpwindQuad[46] = ser_5x_p2_surfx3_eval_quad_node_46_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*alpha[22]+0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]-0.45*alpha[7]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_r(fEdge); 
  } else { 
    fUpwindQuad[47] = ser_5x_p2_surfx3_eval_quad_node_47_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_r(fEdge); 
  } else { 
    fUpwindQuad[48] = ser_5x_p2_surfx3_eval_quad_node_48_l(fSkin); 
  } 
  if (0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.33541019662496846*alpha[2]+0.25*alpha[0] > 0) { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_r(fEdge); 
  } else { 
    fUpwindQuad[49] = ser_5x_p2_surfx3_eval_quad_node_49_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[9]+0.33541019662496846*(alpha[4]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_r(fEdge); 
  } else { 
    fUpwindQuad[50] = ser_5x_p2_surfx3_eval_quad_node_50_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]-0.45*alpha[9]+0.45*alpha[7]-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_r(fEdge); 
  } else { 
    fUpwindQuad[51] = ser_5x_p2_surfx3_eval_quad_node_51_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[22]-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*alpha[7]+0.33541019662496846*(alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_r(fEdge); 
  } else { 
    fUpwindQuad[52] = ser_5x_p2_surfx3_eval_quad_node_52_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22])-0.375*alpha[21]+0.22360679774997896*alpha[12]-0.2795084971874737*alpha[11]+0.45*(alpha[9]+alpha[7])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])+0.25*alpha[0] > 0) { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_r(fEdge); 
  } else { 
    fUpwindQuad[53] = ser_5x_p2_surfx3_eval_quad_node_53_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*alpha[8]+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[3]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_r(fEdge); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx3_eval_quad_node_54_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[7]-0.45*(alpha[6]+alpha[5])-0.33541019662496846*(alpha[3]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_r(fEdge); 
  } else { 
    fUpwindQuad[55] = ser_5x_p2_surfx3_eval_quad_node_55_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*(alpha[8]+alpha[7])-0.45*(alpha[6]+alpha[5])+0.33541019662496846*alpha[4]-0.33541019662496846*(alpha[3]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_r(fEdge); 
  } else { 
    fUpwindQuad[56] = ser_5x_p2_surfx3_eval_quad_node_56_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[5])-0.33541019662496846*(alpha[4]+alpha[2])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_r(fEdge); 
  } else { 
    fUpwindQuad[57] = ser_5x_p2_surfx3_eval_quad_node_57_l(fSkin); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[5]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_r(fEdge); 
  } else { 
    fUpwindQuad[58] = ser_5x_p2_surfx3_eval_quad_node_58_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_r(fEdge); 
  } else { 
    fUpwindQuad[59] = ser_5x_p2_surfx3_eval_quad_node_59_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[9]-0.45*(alpha[8]+alpha[7])+0.45*alpha[6]-0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_r(fEdge); 
  } else { 
    fUpwindQuad[60] = ser_5x_p2_surfx3_eval_quad_node_60_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.33541019662496846*alpha[3]-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_r(fEdge); 
  } else { 
    fUpwindQuad[61] = ser_5x_p2_surfx3_eval_quad_node_61_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*alpha[9]+0.45*alpha[8]-0.45*alpha[7]+0.45*alpha[6]-0.45*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.33541019662496846*alpha[2]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_r(fEdge); 
  } else { 
    fUpwindQuad[62] = ser_5x_p2_surfx3_eval_quad_node_62_l(fSkin); 
  } 
  if (0.375*(alpha[26]+alpha[22])-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*(alpha[8]+alpha[6])-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_r(fEdge); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx3_eval_quad_node_63_l(fSkin); 
  } 
  if (0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[6]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_r(fEdge); 
  } else { 
    fUpwindQuad[64] = ser_5x_p2_surfx3_eval_quad_node_64_l(fSkin); 
  } 
  if (-(0.375*alpha[26])+0.375*alpha[22]-0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]-0.45*alpha[6]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_r(fEdge); 
  } else { 
    fUpwindQuad[65] = ser_5x_p2_surfx3_eval_quad_node_65_l(fSkin); 
  } 
  if (0.375*alpha[26]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_r(fEdge); 
  } else { 
    fUpwindQuad[66] = ser_5x_p2_surfx3_eval_quad_node_66_l(fSkin); 
  } 
  if (-(0.2795084971874737*alpha[12])+0.22360679774997896*alpha[11]+0.33541019662496846*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_r(fEdge); 
  } else { 
    fUpwindQuad[67] = ser_5x_p2_surfx3_eval_quad_node_67_l(fSkin); 
  } 
  if (-(0.375*alpha[26])-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[8]+0.33541019662496846*(alpha[4]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_r(fEdge); 
  } else { 
    fUpwindQuad[68] = ser_5x_p2_surfx3_eval_quad_node_68_l(fSkin); 
  } 
  if (0.375*alpha[26]-0.375*alpha[22]+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]-0.45*alpha[8]+0.45*alpha[6]-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_r(fEdge); 
  } else { 
    fUpwindQuad[69] = ser_5x_p2_surfx3_eval_quad_node_69_l(fSkin); 
  } 
  if (-(0.375*alpha[22])+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*alpha[6]+0.33541019662496846*(alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_r(fEdge); 
  } else { 
    fUpwindQuad[70] = ser_5x_p2_surfx3_eval_quad_node_70_l(fSkin); 
  } 
  if (-(0.375*(alpha[26]+alpha[22]))+0.29999999999999993*alpha[21]-0.2795084971874737*alpha[12]+0.22360679774997896*alpha[11]+0.45*(alpha[8]+alpha[6])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_r(fEdge); 
  } else { 
    fUpwindQuad[71] = ser_5x_p2_surfx3_eval_quad_node_71_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[26]+alpha[22]+alpha[21]))-0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.45*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_r(fEdge); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx3_eval_quad_node_72_l(fSkin); 
  } 
  if (-(0.29999999999999993*(alpha[22]+alpha[21]))-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]-0.33541019662496846*alpha[3]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_r(fEdge); 
  } else { 
    fUpwindQuad[73] = ser_5x_p2_surfx3_eval_quad_node_73_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]-0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[16]-0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8])-0.45*(alpha[7]+alpha[6])+0.45*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_r(fEdge); 
  } else { 
    fUpwindQuad[74] = ser_5x_p2_surfx3_eval_quad_node_74_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])-0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_r(fEdge); 
  } else { 
    fUpwindQuad[75] = ser_5x_p2_surfx3_eval_quad_node_75_l(fSkin); 
  } 
  if (0.22360679774997896*(alpha[12]+alpha[11])+0.45*alpha[5]+0.33541019662496846*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_r(fEdge); 
  } else { 
    fUpwindQuad[76] = ser_5x_p2_surfx3_eval_quad_node_76_l(fSkin); 
  } 
  if (0.29999999999999993*alpha[26]+0.6037383539249431*alpha[16]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_r(fEdge); 
  } else { 
    fUpwindQuad[77] = ser_5x_p2_surfx3_eval_quad_node_77_l(fSkin); 
  } 
  if (-(0.29999999999999993*alpha[26])+0.29999999999999993*(alpha[22]+alpha[21])-0.6037383539249431*alpha[16]+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])-0.45*(alpha[9]+alpha[8])+0.45*(alpha[7]+alpha[6]+alpha[5])-0.33541019662496846*alpha[4]+0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_r(fEdge); 
  } else { 
    fUpwindQuad[78] = ser_5x_p2_surfx3_eval_quad_node_78_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[22]+alpha[21])+0.6037383539249431*alpha[15]+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[7]+alpha[6]+alpha[5])+0.33541019662496846*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_r(fEdge); 
  } else { 
    fUpwindQuad[79] = ser_5x_p2_surfx3_eval_quad_node_79_l(fSkin); 
  } 
  if (0.29999999999999993*(alpha[26]+alpha[22]+alpha[21])+0.6037383539249431*(alpha[16]+alpha[15])+0.22360679774997896*(alpha[12]+alpha[11])+0.45*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])+0.33541019662496846*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_r(fEdge); 
  } else { 
    fUpwindQuad[80] = ser_5x_p2_surfx3_eval_quad_node_80_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[26]*fUpwind[26]+0.25*alpha[22]*fUpwind[22]+0.25*alpha[21]*fUpwind[21]+0.25*alpha[16]*fUpwind[16]+0.25*alpha[15]*fUpwind[15]+0.25*alpha[12]*fUpwind[12]+0.25*alpha[11]*fUpwind[11]+0.25*alpha[9]*fUpwind[9]+0.25*alpha[8]*fUpwind[8]+0.25*alpha[7]*fUpwind[7]+0.25*alpha[6]*fUpwind[6]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[4]*fUpwind[4]+0.25*alpha[3]*fUpwind[3]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.25000000000000006*alpha[26]*fUpwind[36]+0.22360679774997896*alpha[16]*fUpwind[35]+0.25000000000000006*alpha[22]*fUpwind[33]+0.22360679774997896*alpha[15]*fUpwind[32]+0.22360679774997902*alpha[8]*fUpwind[25]+0.22360679774997902*alpha[6]*fUpwind[21]+0.22360679774997902*fUpwind[6]*alpha[21]+0.25000000000000006*alpha[12]*fUpwind[20]+0.22360679774997902*alpha[5]*fUpwind[19]+0.25*alpha[9]*fUpwind[16]+0.25*fUpwind[9]*alpha[16]+0.25*alpha[7]*fUpwind[15]+0.25*fUpwind[7]*alpha[15]+0.22360679774997896*alpha[1]*fUpwind[11]+0.22360679774997896*fUpwind[1]*alpha[11]+0.25*alpha[4]*fUpwind[8]+0.25*fUpwind[4]*alpha[8]+0.25*alpha[3]*fUpwind[6]+0.25*fUpwind[3]*alpha[6]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.22360679774997896*alpha[16]*fUpwind[36]+0.22360679774997896*alpha[15]*fUpwind[33]+0.25000000000000006*alpha[21]*fUpwind[32]+0.22360679774997902*alpha[9]*fUpwind[26]+0.22360679774997902*fUpwind[9]*alpha[26]+0.22360679774997902*alpha[7]*fUpwind[22]+0.22360679774997902*fUpwind[7]*alpha[22]+0.22360679774997902*alpha[5]*fUpwind[20]+0.25000000000000006*alpha[11]*fUpwind[19]+0.25*alpha[8]*fUpwind[16]+0.25*fUpwind[8]*alpha[16]+0.25*alpha[6]*fUpwind[15]+0.25*fUpwind[6]*alpha[15]+0.22360679774997896*alpha[2]*fUpwind[12]+0.22360679774997896*fUpwind[2]*alpha[12]+0.25*alpha[4]*fUpwind[9]+0.25*fUpwind[4]*alpha[9]+0.25*alpha[3]*fUpwind[7]+0.25*fUpwind[3]*alpha[7]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.25000000000000006*alpha[26]*fUpwind[38]+0.22360679774997896*alpha[15]*fUpwind[34]+0.25*alpha[16]*fUpwind[31]+0.22360679774997902*alpha[7]*fUpwind[24]+0.22360679774997902*alpha[6]*fUpwind[23]+0.25000000000000006*alpha[12]*fUpwind[22]+0.25000000000000006*fUpwind[12]*alpha[22]+0.25000000000000006*alpha[11]*fUpwind[21]+0.25000000000000006*fUpwind[11]*alpha[21]+0.25*alpha[9]*fUpwind[18]+0.25*alpha[8]*fUpwind[17]+0.25*alpha[5]*fUpwind[15]+0.25*fUpwind[5]*alpha[15]+0.22360679774997896*alpha[3]*fUpwind[13]+0.25*alpha[4]*fUpwind[10]+0.25*alpha[2]*fUpwind[7]+0.25*fUpwind[2]*alpha[7]+0.25*alpha[1]*fUpwind[6]+0.25*fUpwind[1]*alpha[6]+0.25*alpha[0]*fUpwind[3]+0.25*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.22360679774997896*alpha[16]*fUpwind[41]+0.25000000000000006*alpha[22]*fUpwind[38]+0.25000000000000006*alpha[21]*fUpwind[37]+0.25*alpha[15]*fUpwind[31]+0.22360679774997902*alpha[9]*fUpwind[29]+0.22360679774997902*alpha[8]*fUpwind[28]+0.25000000000000006*alpha[12]*fUpwind[26]+0.25000000000000006*fUpwind[12]*alpha[26]+0.25000000000000006*alpha[11]*fUpwind[25]+0.25*alpha[7]*fUpwind[18]+0.25*alpha[6]*fUpwind[17]+0.25*alpha[5]*fUpwind[16]+0.25*fUpwind[5]*alpha[16]+0.22360679774997896*alpha[4]*fUpwind[14]+0.25*alpha[3]*fUpwind[10]+0.25*alpha[2]*fUpwind[9]+0.25*fUpwind[2]*alpha[9]+0.25*alpha[1]*fUpwind[8]+0.25*fUpwind[1]*alpha[8]+0.25*alpha[0]*fUpwind[4]+0.25*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.22360679774997896*alpha[9]*fUpwind[36]+0.22360679774997896*alpha[8]*fUpwind[35]+0.22360679774997896*alpha[7]*fUpwind[33]+0.22360679774997896*alpha[6]*fUpwind[32]+0.22360679774997902*alpha[16]*fUpwind[26]+0.22360679774997902*fUpwind[16]*alpha[26]+0.22360679774997902*alpha[16]*fUpwind[25]+0.22360679774997902*alpha[15]*fUpwind[22]+0.22360679774997902*fUpwind[15]*alpha[22]+0.22360679774997902*alpha[15]*fUpwind[21]+0.22360679774997902*fUpwind[15]*alpha[21]+0.22360679774997902*alpha[2]*fUpwind[20]+0.22360679774997902*alpha[1]*fUpwind[19]+0.25*alpha[4]*fUpwind[16]+0.25*fUpwind[4]*alpha[16]+0.25*alpha[3]*fUpwind[15]+0.25*fUpwind[3]*alpha[15]+0.22360679774997896*alpha[5]*fUpwind[12]+0.22360679774997896*fUpwind[5]*alpha[12]+0.22360679774997896*alpha[5]*fUpwind[11]+0.22360679774997896*fUpwind[5]*alpha[11]+0.25*alpha[8]*fUpwind[9]+0.25*fUpwind[8]*alpha[9]+0.25*alpha[6]*fUpwind[7]+0.25*fUpwind[6]*alpha[7]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.25*alpha[26]*fUpwind[45]+0.22360679774997902*alpha[16]*fUpwind[44]+0.22360679774997896*alpha[8]*fUpwind[37]+0.22360679774997896*alpha[7]*fUpwind[34]+0.25*alpha[12]*fUpwind[33]+0.22360679774997896*alpha[5]*fUpwind[32]+0.25*alpha[9]*fUpwind[31]+0.22360679774997902*alpha[15]*fUpwind[24]+0.2*alpha[21]*fUpwind[23]+0.22360679774997902*alpha[3]*fUpwind[23]+0.25*fUpwind[20]*alpha[22]+0.22360679774997902*alpha[1]*fUpwind[21]+0.22360679774997902*fUpwind[1]*alpha[21]+0.22360679774997902*alpha[15]*fUpwind[19]+0.25*alpha[16]*fUpwind[18]+0.25*alpha[4]*fUpwind[17]+0.25*alpha[2]*fUpwind[15]+0.25*fUpwind[2]*alpha[15]+0.22360679774997896*alpha[6]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[11]+0.22360679774997896*fUpwind[6]*alpha[11]+0.25*alpha[8]*fUpwind[10]+0.25*alpha[5]*fUpwind[7]+0.25*fUpwind[5]*alpha[7]+0.25*alpha[0]*fUpwind[6]+0.25*fUpwind[0]*alpha[6]+0.25*alpha[1]*fUpwind[3]+0.25*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.22360679774997902*alpha[16]*fUpwind[45]+0.22360679774997896*alpha[9]*fUpwind[38]+0.22360679774997896*alpha[6]*fUpwind[34]+0.22360679774997896*alpha[5]*fUpwind[33]+0.25*alpha[11]*fUpwind[32]+0.25*alpha[8]*fUpwind[31]+0.22360679774997902*fUpwind[18]*alpha[26]+0.2*alpha[22]*fUpwind[24]+0.22360679774997902*alpha[3]*fUpwind[24]+0.22360679774997902*alpha[15]*fUpwind[23]+0.22360679774997902*alpha[2]*fUpwind[22]+0.22360679774997902*fUpwind[2]*alpha[22]+0.25*fUpwind[19]*alpha[21]+0.22360679774997902*alpha[15]*fUpwind[20]+0.25*alpha[4]*fUpwind[18]+0.25*alpha[16]*fUpwind[17]+0.25*alpha[1]*fUpwind[15]+0.25*fUpwind[1]*alpha[15]+0.22360679774997896*alpha[7]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[12]+0.22360679774997896*fUpwind[7]*alpha[12]+0.25*alpha[9]*fUpwind[10]+0.25*alpha[0]*fUpwind[7]+0.25*fUpwind[0]*alpha[7]+0.25*alpha[5]*fUpwind[6]+0.25*fUpwind[5]*alpha[6]+0.25*alpha[2]*fUpwind[3]+0.25*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.25*alpha[22]*fUpwind[45]+0.22360679774997902*alpha[15]*fUpwind[44]+0.22360679774997896*alpha[9]*fUpwind[41]+0.22360679774997896*alpha[6]*fUpwind[37]+0.25*alpha[12]*fUpwind[36]+0.22360679774997896*alpha[5]*fUpwind[35]+0.25*alpha[7]*fUpwind[31]+0.22360679774997902*alpha[16]*fUpwind[29]+0.22360679774997902*alpha[4]*fUpwind[28]+0.25*fUpwind[20]*alpha[26]+0.22360679774997902*alpha[1]*fUpwind[25]+0.22360679774997902*fUpwind[17]*alpha[21]+0.22360679774997902*alpha[16]*fUpwind[19]+0.25*alpha[15]*fUpwind[18]+0.25*alpha[3]*fUpwind[17]+0.25*alpha[2]*fUpwind[16]+0.25*fUpwind[2]*alpha[16]+0.22360679774997896*alpha[8]*fUpwind[14]+0.22360679774997896*alpha[8]*fUpwind[11]+0.22360679774997896*fUpwind[8]*alpha[11]+0.25*alpha[6]*fUpwind[10]+0.25*alpha[5]*fUpwind[9]+0.25*fUpwind[5]*alpha[9]+0.25*alpha[0]*fUpwind[8]+0.25*fUpwind[0]*alpha[8]+0.25*alpha[1]*fUpwind[4]+0.25*fUpwind[1]*alpha[4]; 
  Ghat[9] = 0.22360679774997902*alpha[15]*fUpwind[45]+0.25*alpha[21]*fUpwind[44]+0.22360679774997896*alpha[8]*fUpwind[41]+0.22360679774997896*alpha[7]*fUpwind[38]+0.22360679774997896*alpha[5]*fUpwind[36]+0.25*alpha[11]*fUpwind[35]+0.25*alpha[6]*fUpwind[31]+0.2*alpha[26]*fUpwind[29]+0.22360679774997902*alpha[4]*fUpwind[29]+0.22360679774997902*alpha[16]*fUpwind[28]+0.22360679774997902*alpha[2]*fUpwind[26]+0.22360679774997902*fUpwind[2]*alpha[26]+0.22360679774997902*fUpwind[18]*alpha[22]+0.22360679774997902*alpha[16]*fUpwind[20]+0.25*alpha[3]*fUpwind[18]+0.25*alpha[15]*fUpwind[17]+0.25*alpha[1]*fUpwind[16]+0.25*fUpwind[1]*alpha[16]+0.22360679774997896*alpha[9]*fUpwind[14]+0.22360679774997896*alpha[9]*fUpwind[12]+0.22360679774997896*fUpwind[9]*alpha[12]+0.25*alpha[7]*fUpwind[10]+0.25*alpha[0]*fUpwind[9]+0.25*fUpwind[0]*alpha[9]+0.25*alpha[5]*fUpwind[8]+0.25*fUpwind[5]*alpha[8]+0.25*alpha[2]*fUpwind[4]+0.25*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.22360679774997902*alpha[16]*fUpwind[47]+0.22360679774997902*alpha[15]*fUpwind[46]+0.22360679774997896*alpha[9]*fUpwind[43]+0.22360679774997896*alpha[8]*fUpwind[42]+0.22360679774997896*alpha[7]*fUpwind[40]+0.22360679774997896*alpha[6]*fUpwind[39]+0.25*alpha[12]*fUpwind[38]+0.25*alpha[11]*fUpwind[37]+0.25*alpha[5]*fUpwind[31]+0.22360679774997902*alpha[4]*fUpwind[30]+0.22360679774997902*alpha[3]*fUpwind[27]+0.25*alpha[22]*fUpwind[26]+0.25*fUpwind[22]*alpha[26]+0.25*alpha[21]*fUpwind[25]+0.25*alpha[2]*fUpwind[18]+0.25*alpha[1]*fUpwind[17]+0.25*alpha[15]*fUpwind[16]+0.25*fUpwind[15]*alpha[16]+0.25*alpha[0]*fUpwind[10]+0.25*alpha[7]*fUpwind[9]+0.25*fUpwind[7]*alpha[9]+0.25*alpha[6]*fUpwind[8]+0.25*fUpwind[6]*alpha[8]+0.25*alpha[3]*fUpwind[4]+0.25*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.25*alpha[9]*fUpwind[35]+0.25*alpha[7]*fUpwind[32]+0.25000000000000006*alpha[4]*fUpwind[25]+0.15971914124998499*alpha[21]*fUpwind[21]+0.25000000000000006*alpha[3]*fUpwind[21]+0.25000000000000006*fUpwind[3]*alpha[21]+0.25000000000000006*alpha[2]*fUpwind[19]+0.22360679774997896*alpha[16]*fUpwind[16]+0.22360679774997896*alpha[15]*fUpwind[15]+0.15971914124998499*alpha[11]*fUpwind[11]+0.25*alpha[0]*fUpwind[11]+0.25*fUpwind[0]*alpha[11]+0.22360679774997896*alpha[8]*fUpwind[8]+0.22360679774997896*alpha[6]*fUpwind[6]+0.22360679774997896*alpha[5]*fUpwind[5]+0.22360679774997896*alpha[1]*fUpwind[1]; 
  Ghat[12] = 0.25*alpha[8]*fUpwind[36]+0.25*alpha[6]*fUpwind[33]+0.15971914124998499*alpha[26]*fUpwind[26]+0.25000000000000006*alpha[4]*fUpwind[26]+0.25000000000000006*fUpwind[4]*alpha[26]+0.15971914124998499*alpha[22]*fUpwind[22]+0.25000000000000006*alpha[3]*fUpwind[22]+0.25000000000000006*fUpwind[3]*alpha[22]+0.25000000000000006*alpha[1]*fUpwind[20]+0.22360679774997896*alpha[16]*fUpwind[16]+0.22360679774997896*alpha[15]*fUpwind[15]+0.15971914124998499*alpha[12]*fUpwind[12]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.22360679774997896*alpha[9]*fUpwind[9]+0.22360679774997896*alpha[7]*fUpwind[7]+0.22360679774997896*alpha[5]*fUpwind[5]+0.22360679774997896*alpha[2]*fUpwind[2]; 
  Ghat[13] = 0.25000000000000006*alpha[16]*fUpwind[46]+0.25*alpha[9]*fUpwind[40]+0.25*alpha[8]*fUpwind[39]+0.25*alpha[5]*fUpwind[34]+0.25000000000000006*alpha[4]*fUpwind[27]+0.25000000000000006*alpha[2]*fUpwind[24]+0.25000000000000006*alpha[1]*fUpwind[23]+0.22360679774997896*alpha[22]*fUpwind[22]+0.22360679774997896*alpha[21]*fUpwind[21]+0.22360679774997896*alpha[15]*fUpwind[15]+0.25*alpha[0]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[7]+0.22360679774997896*alpha[6]*fUpwind[6]+0.22360679774997896*alpha[3]*fUpwind[3]; 
  Ghat[14] = 0.25000000000000006*alpha[15]*fUpwind[47]+0.25*alpha[7]*fUpwind[43]+0.25*alpha[6]*fUpwind[42]+0.25*alpha[5]*fUpwind[41]+0.25000000000000006*alpha[3]*fUpwind[30]+0.25000000000000006*alpha[2]*fUpwind[29]+0.25000000000000006*alpha[1]*fUpwind[28]+0.22360679774997896*alpha[26]*fUpwind[26]+0.22360679774997896*alpha[16]*fUpwind[16]+0.25*alpha[0]*fUpwind[14]+0.22360679774997896*alpha[9]*fUpwind[9]+0.22360679774997896*alpha[8]*fUpwind[8]+0.22360679774997896*alpha[4]*fUpwind[4]; 
  Ghat[15] = 0.22360679774997902*alpha[9]*fUpwind[45]+0.22360679774997902*alpha[8]*fUpwind[44]+0.22360679774997896*alpha[16]*fUpwind[38]+0.22360679774997896*alpha[16]*fUpwind[37]+0.2*alpha[22]*fUpwind[34]+0.2*alpha[21]*fUpwind[34]+0.22360679774997896*alpha[3]*fUpwind[34]+0.22360679774997896*alpha[2]*fUpwind[33]+0.22360679774997896*alpha[1]*fUpwind[32]+0.22360679774997902*alpha[26]*fUpwind[31]+0.25*alpha[4]*fUpwind[31]+0.22360679774997902*alpha[6]*fUpwind[24]+0.22360679774997902*alpha[7]*fUpwind[23]+0.22360679774997902*alpha[5]*fUpwind[22]+0.22360679774997902*fUpwind[5]*alpha[22]+0.22360679774997902*alpha[5]*fUpwind[21]+0.22360679774997902*fUpwind[5]*alpha[21]+0.22360679774997902*alpha[7]*fUpwind[20]+0.22360679774997902*alpha[6]*fUpwind[19]+0.25*alpha[8]*fUpwind[18]+0.25*alpha[9]*fUpwind[17]+0.25*fUpwind[10]*alpha[16]+0.22360679774997896*alpha[12]*fUpwind[15]+0.22360679774997896*alpha[11]*fUpwind[15]+0.25*alpha[0]*fUpwind[15]+0.22360679774997896*fUpwind[13]*alpha[15]+0.22360679774997896*fUpwind[12]*alpha[15]+0.22360679774997896*fUpwind[11]*alpha[15]+0.25*fUpwind[0]*alpha[15]+0.25*alpha[1]*fUpwind[7]+0.25*fUpwind[1]*alpha[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[2]*alpha[6]+0.25*alpha[3]*fUpwind[5]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[16] = 0.22360679774997902*alpha[7]*fUpwind[45]+0.22360679774997902*alpha[6]*fUpwind[44]+0.2*alpha[26]*fUpwind[41]+0.22360679774997896*alpha[4]*fUpwind[41]+0.22360679774997896*alpha[15]*fUpwind[38]+0.22360679774997896*alpha[15]*fUpwind[37]+0.22360679774997896*alpha[2]*fUpwind[36]+0.22360679774997896*alpha[1]*fUpwind[35]+0.22360679774997902*alpha[22]*fUpwind[31]+0.22360679774997902*alpha[21]*fUpwind[31]+0.25*alpha[3]*fUpwind[31]+0.22360679774997902*alpha[8]*fUpwind[29]+0.22360679774997902*alpha[9]*fUpwind[28]+0.22360679774997902*alpha[5]*fUpwind[26]+0.22360679774997902*fUpwind[5]*alpha[26]+0.22360679774997902*alpha[5]*fUpwind[25]+0.22360679774997902*alpha[9]*fUpwind[20]+0.22360679774997902*alpha[8]*fUpwind[19]+0.25*alpha[6]*fUpwind[18]+0.25*alpha[7]*fUpwind[17]+0.22360679774997896*alpha[12]*fUpwind[16]+0.22360679774997896*alpha[11]*fUpwind[16]+0.25*alpha[0]*fUpwind[16]+0.22360679774997896*fUpwind[14]*alpha[16]+0.22360679774997896*fUpwind[12]*alpha[16]+0.22360679774997896*fUpwind[11]*alpha[16]+0.25*fUpwind[0]*alpha[16]+0.25*fUpwind[10]*alpha[15]+0.25*alpha[1]*fUpwind[9]+0.25*fUpwind[1]*alpha[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[2]*alpha[8]+0.25*alpha[4]*fUpwind[5]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[17] = 0.22360679774997902*alpha[9]*fUpwind[47]+0.22360679774997902*alpha[7]*fUpwind[46]+0.25000000000000006*alpha[12]*fUpwind[45]+0.22360679774997902*alpha[5]*fUpwind[44]+0.22360679774997896*alpha[16]*fUpwind[43]+0.22360679774997896*alpha[4]*fUpwind[42]+0.22360679774997896*alpha[15]*fUpwind[40]+0.2*alpha[21]*fUpwind[39]+0.22360679774997896*alpha[3]*fUpwind[39]+0.22360679774997896*alpha[1]*fUpwind[37]+0.25000000000000006*alpha[22]*fUpwind[36]+0.22360679774997896*alpha[15]*fUpwind[35]+0.25000000000000006*alpha[26]*fUpwind[33]+0.22360679774997896*alpha[16]*fUpwind[32]+0.25*alpha[2]*fUpwind[31]+0.22360679774997902*alpha[8]*fUpwind[30]+0.22360679774997902*alpha[6]*fUpwind[27]+0.22360679774997902*alpha[6]*fUpwind[25]+0.22360679774997902*alpha[8]*fUpwind[21]+0.22360679774997902*fUpwind[8]*alpha[21]+0.25*alpha[5]*fUpwind[18]+0.22360679774997896*alpha[11]*fUpwind[17]+0.25*alpha[0]*fUpwind[17]+0.25*alpha[7]*fUpwind[16]+0.25*fUpwind[7]*alpha[16]+0.25*alpha[9]*fUpwind[15]+0.25*fUpwind[9]*alpha[15]+0.25*alpha[1]*fUpwind[10]+0.25*alpha[3]*fUpwind[8]+0.25*fUpwind[3]*alpha[8]+0.25*alpha[4]*fUpwind[6]+0.25*fUpwind[4]*alpha[6]; 
  Ghat[18] = 0.22360679774997902*alpha[8]*fUpwind[47]+0.22360679774997902*alpha[6]*fUpwind[46]+0.22360679774997902*alpha[5]*fUpwind[45]+0.25000000000000006*alpha[11]*fUpwind[44]+0.2*alpha[26]*fUpwind[43]+0.22360679774997896*alpha[4]*fUpwind[43]+0.22360679774997896*alpha[16]*fUpwind[42]+0.2*alpha[22]*fUpwind[40]+0.22360679774997896*alpha[3]*fUpwind[40]+0.22360679774997896*alpha[15]*fUpwind[39]+0.22360679774997896*alpha[2]*fUpwind[38]+0.22360679774997896*alpha[15]*fUpwind[36]+0.25000000000000006*alpha[21]*fUpwind[35]+0.22360679774997896*alpha[16]*fUpwind[33]+0.25*alpha[1]*fUpwind[31]+0.22360679774997902*alpha[9]*fUpwind[30]+0.22360679774997902*alpha[7]*fUpwind[27]+0.22360679774997902*alpha[7]*fUpwind[26]+0.22360679774997902*fUpwind[7]*alpha[26]+0.22360679774997902*alpha[9]*fUpwind[22]+0.22360679774997902*fUpwind[9]*alpha[22]+0.22360679774997896*alpha[12]*fUpwind[18]+0.25*alpha[0]*fUpwind[18]+0.25*alpha[5]*fUpwind[17]+0.25*alpha[6]*fUpwind[16]+0.25*fUpwind[6]*alpha[16]+0.25*alpha[8]*fUpwind[15]+0.25*fUpwind[8]*alpha[15]+0.25*alpha[2]*fUpwind[10]+0.25*alpha[3]*fUpwind[9]+0.25*fUpwind[3]*alpha[9]+0.25*alpha[4]*fUpwind[7]+0.25*fUpwind[4]*alpha[7]; 
  Ghat[19] = 0.2*alpha[16]*fUpwind[36]+0.22360679774997896*alpha[26]*fUpwind[35]+0.25000000000000006*alpha[4]*fUpwind[35]+0.2*alpha[15]*fUpwind[33]+0.22360679774997896*alpha[22]*fUpwind[32]+0.15971914124998499*alpha[21]*fUpwind[32]+0.25000000000000006*alpha[3]*fUpwind[32]+0.25*alpha[9]*fUpwind[25]+0.25*alpha[7]*fUpwind[21]+0.25*fUpwind[7]*alpha[21]+0.2*alpha[5]*fUpwind[20]+0.22360679774997896*alpha[12]*fUpwind[19]+0.15971914124998499*alpha[11]*fUpwind[19]+0.25*alpha[0]*fUpwind[19]+0.22360679774997902*alpha[8]*fUpwind[16]+0.22360679774997902*fUpwind[8]*alpha[16]+0.22360679774997902*alpha[6]*fUpwind[15]+0.22360679774997902*fUpwind[6]*alpha[15]+0.25000000000000006*alpha[2]*fUpwind[11]+0.25000000000000006*fUpwind[2]*alpha[11]+0.22360679774997902*alpha[1]*fUpwind[5]+0.22360679774997902*fUpwind[1]*alpha[5]; 
  Ghat[20] = 0.15971914124998499*alpha[26]*fUpwind[36]+0.25000000000000006*alpha[4]*fUpwind[36]+0.2*alpha[16]*fUpwind[35]+0.15971914124998499*alpha[22]*fUpwind[33]+0.22360679774997896*alpha[21]*fUpwind[33]+0.25000000000000006*alpha[3]*fUpwind[33]+0.2*alpha[15]*fUpwind[32]+0.25*alpha[8]*fUpwind[26]+0.25*fUpwind[8]*alpha[26]+0.25*alpha[6]*fUpwind[22]+0.25*fUpwind[6]*alpha[22]+0.15971914124998499*alpha[12]*fUpwind[20]+0.22360679774997896*alpha[11]*fUpwind[20]+0.25*alpha[0]*fUpwind[20]+0.2*alpha[5]*fUpwind[19]+0.22360679774997902*alpha[9]*fUpwind[16]+0.22360679774997902*fUpwind[9]*alpha[16]+0.22360679774997902*alpha[7]*fUpwind[15]+0.22360679774997902*fUpwind[7]*alpha[15]+0.25000000000000006*alpha[1]*fUpwind[12]+0.25000000000000006*fUpwind[1]*alpha[12]+0.22360679774997902*alpha[2]*fUpwind[5]+0.22360679774997902*fUpwind[2]*alpha[5]; 
  Ghat[21] = 0.25*alpha[9]*fUpwind[44]+0.25000000000000006*alpha[4]*fUpwind[37]+0.2*alpha[15]*fUpwind[34]+0.25000000000000006*alpha[2]*fUpwind[32]+0.22360679774997902*alpha[16]*fUpwind[31]+0.2*alpha[6]*fUpwind[23]+0.15971914124998499*alpha[11]*fUpwind[21]+0.25*alpha[0]*fUpwind[21]+0.22360679774997896*fUpwind[13]*alpha[21]+0.15971914124998499*fUpwind[11]*alpha[21]+0.25*fUpwind[0]*alpha[21]+0.25*alpha[7]*fUpwind[19]+0.22360679774997902*alpha[8]*fUpwind[17]+0.22360679774997902*alpha[5]*fUpwind[15]+0.22360679774997902*fUpwind[5]*alpha[15]+0.25000000000000006*alpha[3]*fUpwind[11]+0.25000000000000006*fUpwind[3]*alpha[11]+0.22360679774997902*alpha[1]*fUpwind[6]+0.22360679774997902*fUpwind[1]*alpha[6]; 
  Ghat[22] = 0.25*alpha[8]*fUpwind[45]+0.15971914124998499*alpha[26]*fUpwind[38]+0.25000000000000006*alpha[4]*fUpwind[38]+0.2*alpha[15]*fUpwind[34]+0.25000000000000006*alpha[1]*fUpwind[33]+0.22360679774997902*alpha[16]*fUpwind[31]+0.25*fUpwind[10]*alpha[26]+0.2*alpha[7]*fUpwind[24]+0.15971914124998499*alpha[12]*fUpwind[22]+0.25*alpha[0]*fUpwind[22]+0.22360679774997896*fUpwind[13]*alpha[22]+0.15971914124998499*fUpwind[12]*alpha[22]+0.25*fUpwind[0]*alpha[22]+0.25*alpha[6]*fUpwind[20]+0.22360679774997902*alpha[9]*fUpwind[18]+0.22360679774997902*alpha[5]*fUpwind[15]+0.22360679774997902*fUpwind[5]*alpha[15]+0.25000000000000006*alpha[3]*fUpwind[12]+0.25000000000000006*fUpwind[3]*alpha[12]+0.22360679774997902*alpha[2]*fUpwind[7]+0.22360679774997902*fUpwind[2]*alpha[7]; 
  Ghat[23] = 0.25*alpha[9]*fUpwind[46]+0.25000000000000006*alpha[16]*fUpwind[40]+0.25000000000000006*alpha[4]*fUpwind[39]+0.25000000000000006*alpha[2]*fUpwind[34]+0.22360679774997896*alpha[22]*fUpwind[33]+0.2*alpha[15]*fUpwind[32]+0.25*alpha[8]*fUpwind[27]+0.25*alpha[5]*fUpwind[24]+0.22360679774997896*alpha[11]*fUpwind[23]+0.25*alpha[0]*fUpwind[23]+0.2*alpha[6]*fUpwind[21]+0.2*fUpwind[6]*alpha[21]+0.22360679774997902*alpha[7]*fUpwind[15]+0.22360679774997902*fUpwind[7]*alpha[15]+0.25000000000000006*alpha[1]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[6]+0.22360679774997902*fUpwind[3]*alpha[6]; 
  Ghat[24] = 0.25*alpha[8]*fUpwind[46]+0.22360679774997896*alpha[26]*fUpwind[40]+0.25000000000000006*alpha[4]*fUpwind[40]+0.25000000000000006*alpha[16]*fUpwind[39]+0.25000000000000006*alpha[1]*fUpwind[34]+0.2*alpha[15]*fUpwind[33]+0.22360679774997896*alpha[21]*fUpwind[32]+0.25*alpha[9]*fUpwind[27]+0.22360679774997896*alpha[12]*fUpwind[24]+0.25*alpha[0]*fUpwind[24]+0.25*alpha[5]*fUpwind[23]+0.2*alpha[7]*fUpwind[22]+0.2*fUpwind[7]*alpha[22]+0.22360679774997902*alpha[6]*fUpwind[15]+0.22360679774997902*fUpwind[6]*alpha[15]+0.25000000000000006*alpha[2]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[7]+0.22360679774997902*fUpwind[3]*alpha[7]; 
  Ghat[25] = 0.25*alpha[7]*fUpwind[44]+0.2*alpha[16]*fUpwind[41]+0.15971914124998499*alpha[21]*fUpwind[37]+0.25000000000000006*alpha[3]*fUpwind[37]+0.25000000000000006*alpha[2]*fUpwind[35]+0.22360679774997902*alpha[15]*fUpwind[31]+0.2*alpha[8]*fUpwind[28]+0.15971914124998499*alpha[11]*fUpwind[25]+0.25*alpha[0]*fUpwind[25]+0.25*fUpwind[10]*alpha[21]+0.25*alpha[9]*fUpwind[19]+0.22360679774997902*alpha[6]*fUpwind[17]+0.22360679774997902*alpha[5]*fUpwind[16]+0.22360679774997902*fUpwind[5]*alpha[16]+0.25000000000000006*alpha[4]*fUpwind[11]+0.25000000000000006*fUpwind[4]*alpha[11]+0.22360679774997902*alpha[1]*fUpwind[8]+0.22360679774997902*fUpwind[1]*alpha[8]; 
  Ghat[26] = 0.25*alpha[6]*fUpwind[45]+0.2*alpha[16]*fUpwind[41]+0.15971914124998499*alpha[22]*fUpwind[38]+0.25000000000000006*alpha[3]*fUpwind[38]+0.25000000000000006*alpha[1]*fUpwind[36]+0.22360679774997902*alpha[15]*fUpwind[31]+0.2*alpha[9]*fUpwind[29]+0.15971914124998499*alpha[12]*fUpwind[26]+0.25*alpha[0]*fUpwind[26]+0.22360679774997896*fUpwind[14]*alpha[26]+0.15971914124998499*fUpwind[12]*alpha[26]+0.25*fUpwind[0]*alpha[26]+0.25*fUpwind[10]*alpha[22]+0.25*alpha[8]*fUpwind[20]+0.22360679774997902*alpha[7]*fUpwind[18]+0.22360679774997902*alpha[5]*fUpwind[16]+0.22360679774997902*fUpwind[5]*alpha[16]+0.25000000000000006*alpha[4]*fUpwind[12]+0.25000000000000006*fUpwind[4]*alpha[12]+0.22360679774997902*alpha[2]*fUpwind[9]+0.22360679774997902*fUpwind[2]*alpha[9]; 
  Ghat[27] = 0.25*alpha[5]*fUpwind[46]+0.25000000000000006*alpha[2]*fUpwind[40]+0.25000000000000006*alpha[1]*fUpwind[39]+0.22360679774997896*alpha[22]*fUpwind[38]+0.22360679774997896*alpha[21]*fUpwind[37]+0.25000000000000006*alpha[16]*fUpwind[34]+0.22360679774997902*alpha[15]*fUpwind[31]+0.25*alpha[0]*fUpwind[27]+0.25*alpha[9]*fUpwind[24]+0.25*alpha[8]*fUpwind[23]+0.22360679774997902*alpha[7]*fUpwind[18]+0.22360679774997902*alpha[6]*fUpwind[17]+0.25000000000000006*alpha[4]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[10]; 
  Ghat[28] = 0.25*alpha[7]*fUpwind[47]+0.25000000000000006*alpha[15]*fUpwind[43]+0.22360679774997896*alpha[21]*fUpwind[42]+0.25000000000000006*alpha[3]*fUpwind[42]+0.25000000000000006*alpha[2]*fUpwind[41]+0.22360679774997896*alpha[26]*fUpwind[36]+0.2*alpha[16]*fUpwind[35]+0.25*alpha[6]*fUpwind[30]+0.25*alpha[5]*fUpwind[29]+0.22360679774997896*alpha[11]*fUpwind[28]+0.25*alpha[0]*fUpwind[28]+0.2*alpha[8]*fUpwind[25]+0.22360679774997902*alpha[9]*fUpwind[16]+0.22360679774997902*fUpwind[9]*alpha[16]+0.25000000000000006*alpha[1]*fUpwind[14]+0.22360679774997902*alpha[4]*fUpwind[8]+0.22360679774997902*fUpwind[4]*alpha[8]; 
  Ghat[29] = 0.25*alpha[6]*fUpwind[47]+0.22360679774997896*alpha[22]*fUpwind[43]+0.25000000000000006*alpha[3]*fUpwind[43]+0.25000000000000006*alpha[15]*fUpwind[42]+0.25000000000000006*alpha[1]*fUpwind[41]+0.2*alpha[16]*fUpwind[36]+0.25*alpha[7]*fUpwind[30]+0.22360679774997896*alpha[12]*fUpwind[29]+0.25*alpha[0]*fUpwind[29]+0.25*alpha[5]*fUpwind[28]+0.2*alpha[9]*fUpwind[26]+0.2*fUpwind[9]*alpha[26]+0.22360679774997902*alpha[8]*fUpwind[16]+0.22360679774997902*fUpwind[8]*alpha[16]+0.25000000000000006*alpha[2]*fUpwind[14]+0.22360679774997902*alpha[4]*fUpwind[9]+0.22360679774997902*fUpwind[4]*alpha[9]; 
  Ghat[30] = 0.25*alpha[5]*fUpwind[47]+0.25000000000000006*alpha[2]*fUpwind[43]+0.25000000000000006*alpha[1]*fUpwind[42]+0.25000000000000006*alpha[15]*fUpwind[41]+0.22360679774997896*alpha[26]*fUpwind[38]+0.22360679774997902*alpha[16]*fUpwind[31]+0.25*alpha[0]*fUpwind[30]+0.25*alpha[7]*fUpwind[29]+0.25*alpha[6]*fUpwind[28]+0.22360679774997902*alpha[9]*fUpwind[18]+0.22360679774997902*alpha[8]*fUpwind[17]+0.25000000000000006*alpha[3]*fUpwind[14]+0.22360679774997902*alpha[4]*fUpwind[10]; 
  Ghat[31] = 0.2*alpha[26]*fUpwind[47]+0.22360679774997902*alpha[4]*fUpwind[47]+0.2*alpha[22]*fUpwind[46]+0.2*alpha[21]*fUpwind[46]+0.22360679774997902*alpha[3]*fUpwind[46]+0.22360679774997902*alpha[2]*fUpwind[45]+0.22360679774997902*alpha[1]*fUpwind[44]+0.22360679774997896*alpha[8]*fUpwind[43]+0.22360679774997896*alpha[9]*fUpwind[42]+0.22360679774997896*alpha[6]*fUpwind[40]+0.22360679774997896*alpha[7]*fUpwind[39]+0.22360679774997896*alpha[5]*fUpwind[38]+0.22360679774997896*alpha[5]*fUpwind[37]+0.22360679774997896*alpha[7]*fUpwind[36]+0.22360679774997896*alpha[6]*fUpwind[35]+0.22360679774997896*alpha[9]*fUpwind[33]+0.22360679774997896*alpha[8]*fUpwind[32]+0.22360679774997896*alpha[12]*fUpwind[31]+0.22360679774997896*alpha[11]*fUpwind[31]+0.25*alpha[0]*fUpwind[31]+0.22360679774997902*alpha[16]*fUpwind[30]+0.22360679774997902*alpha[15]*fUpwind[27]+0.22360679774997902*alpha[15]*fUpwind[26]+0.22360679774997902*fUpwind[15]*alpha[26]+0.22360679774997902*alpha[15]*fUpwind[25]+0.22360679774997902*alpha[16]*fUpwind[22]+0.22360679774997902*fUpwind[16]*alpha[22]+0.22360679774997902*alpha[16]*fUpwind[21]+0.22360679774997902*fUpwind[16]*alpha[21]+0.25*alpha[1]*fUpwind[18]+0.25*alpha[2]*fUpwind[17]+0.25*alpha[3]*fUpwind[16]+0.25*fUpwind[3]*alpha[16]+0.25*alpha[4]*fUpwind[15]+0.25*fUpwind[4]*alpha[15]+0.25*alpha[5]*fUpwind[10]+0.25*alpha[6]*fUpwind[9]+0.25*fUpwind[6]*alpha[9]+0.25*alpha[7]*fUpwind[8]+0.25*fUpwind[7]*alpha[8]; 
  Ghat[32] = 0.2*alpha[16]*fUpwind[45]+0.22360679774997896*alpha[26]*fUpwind[44]+0.25000000000000006*alpha[4]*fUpwind[44]+0.25*alpha[9]*fUpwind[37]+0.2*alpha[6]*fUpwind[34]+0.2*alpha[5]*fUpwind[33]+0.22360679774997896*alpha[12]*fUpwind[32]+0.15971914124998499*alpha[11]*fUpwind[32]+0.25*alpha[0]*fUpwind[32]+0.22360679774997896*alpha[8]*fUpwind[31]+0.22360679774997896*alpha[21]*fUpwind[24]+0.2*alpha[15]*fUpwind[23]+0.22360679774997896*fUpwind[19]*alpha[22]+0.25000000000000006*alpha[2]*fUpwind[21]+0.15971914124998499*fUpwind[19]*alpha[21]+0.25000000000000006*fUpwind[2]*alpha[21]+0.2*alpha[15]*fUpwind[20]+0.25000000000000006*alpha[3]*fUpwind[19]+0.22360679774997896*alpha[16]*fUpwind[17]+0.22360679774997896*alpha[1]*fUpwind[15]+0.22360679774997896*fUpwind[1]*alpha[15]+0.25*alpha[7]*fUpwind[11]+0.25*fUpwind[7]*alpha[11]+0.22360679774997896*alpha[5]*fUpwind[6]+0.22360679774997896*fUpwind[5]*alpha[6]; 
  Ghat[33] = 0.15971914124998499*alpha[26]*fUpwind[45]+0.25000000000000006*alpha[4]*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+0.25*alpha[8]*fUpwind[38]+0.2*alpha[7]*fUpwind[34]+0.15971914124998499*alpha[12]*fUpwind[33]+0.22360679774997896*alpha[11]*fUpwind[33]+0.25*alpha[0]*fUpwind[33]+0.2*alpha[5]*fUpwind[32]+0.22360679774997896*alpha[9]*fUpwind[31]+0.25000000000000006*fUpwind[17]*alpha[26]+0.2*alpha[15]*fUpwind[24]+0.22360679774997896*alpha[22]*fUpwind[23]+0.25000000000000006*alpha[1]*fUpwind[22]+0.15971914124998499*fUpwind[20]*alpha[22]+0.25000000000000006*fUpwind[1]*alpha[22]+0.22360679774997896*fUpwind[20]*alpha[21]+0.25000000000000006*alpha[3]*fUpwind[20]+0.2*alpha[15]*fUpwind[19]+0.22360679774997896*alpha[16]*fUpwind[18]+0.22360679774997896*alpha[2]*fUpwind[15]+0.22360679774997896*fUpwind[2]*alpha[15]+0.25*alpha[6]*fUpwind[12]+0.25*fUpwind[6]*alpha[12]+0.22360679774997896*alpha[5]*fUpwind[7]+0.22360679774997896*fUpwind[5]*alpha[7]; 
  Ghat[34] = 0.22360679774997896*alpha[26]*fUpwind[46]+0.25000000000000006*alpha[4]*fUpwind[46]+0.25*alpha[8]*fUpwind[40]+0.25*alpha[9]*fUpwind[39]+0.22360679774997896*alpha[12]*fUpwind[34]+0.22360679774997896*alpha[11]*fUpwind[34]+0.25*alpha[0]*fUpwind[34]+0.2*alpha[7]*fUpwind[33]+0.2*alpha[6]*fUpwind[32]+0.25000000000000006*alpha[16]*fUpwind[27]+0.25000000000000006*alpha[1]*fUpwind[24]+0.25000000000000006*alpha[2]*fUpwind[23]+0.2*alpha[15]*fUpwind[22]+0.2*fUpwind[15]*alpha[22]+0.2*alpha[15]*fUpwind[21]+0.2*fUpwind[15]*alpha[21]+0.22360679774997896*alpha[3]*fUpwind[15]+0.22360679774997896*fUpwind[3]*alpha[15]+0.25*alpha[5]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[7]+0.22360679774997896*fUpwind[6]*alpha[7]; 
  Ghat[35] = 0.2*alpha[15]*fUpwind[45]+0.22360679774997896*alpha[22]*fUpwind[44]+0.15971914124998499*alpha[21]*fUpwind[44]+0.25000000000000006*alpha[3]*fUpwind[44]+0.2*alpha[8]*fUpwind[41]+0.25*alpha[7]*fUpwind[37]+0.2*alpha[5]*fUpwind[36]+0.22360679774997896*alpha[12]*fUpwind[35]+0.15971914124998499*alpha[11]*fUpwind[35]+0.25*alpha[0]*fUpwind[35]+0.22360679774997896*alpha[6]*fUpwind[31]+0.2*alpha[16]*fUpwind[28]+0.22360679774997896*fUpwind[19]*alpha[26]+0.25000000000000006*alpha[2]*fUpwind[25]+0.25000000000000006*fUpwind[18]*alpha[21]+0.2*alpha[16]*fUpwind[20]+0.25000000000000006*alpha[4]*fUpwind[19]+0.22360679774997896*alpha[15]*fUpwind[17]+0.22360679774997896*alpha[1]*fUpwind[16]+0.22360679774997896*fUpwind[1]*alpha[16]+0.25*alpha[9]*fUpwind[11]+0.25*fUpwind[9]*alpha[11]+0.22360679774997896*alpha[5]*fUpwind[8]+0.22360679774997896*fUpwind[5]*alpha[8]; 
  Ghat[36] = 0.15971914124998499*alpha[22]*fUpwind[45]+0.22360679774997896*alpha[21]*fUpwind[45]+0.25000000000000006*alpha[3]*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+0.2*alpha[9]*fUpwind[41]+0.25*alpha[6]*fUpwind[38]+0.15971914124998499*alpha[12]*fUpwind[36]+0.22360679774997896*alpha[11]*fUpwind[36]+0.25*alpha[0]*fUpwind[36]+0.2*alpha[5]*fUpwind[35]+0.22360679774997896*alpha[7]*fUpwind[31]+0.2*alpha[16]*fUpwind[29]+0.22360679774997896*alpha[26]*fUpwind[28]+0.25000000000000006*alpha[1]*fUpwind[26]+0.15971914124998499*fUpwind[20]*alpha[26]+0.25000000000000006*fUpwind[1]*alpha[26]+0.25000000000000006*fUpwind[17]*alpha[22]+0.25000000000000006*alpha[4]*fUpwind[20]+0.2*alpha[16]*fUpwind[19]+0.22360679774997896*alpha[15]*fUpwind[18]+0.22360679774997896*alpha[2]*fUpwind[16]+0.22360679774997896*fUpwind[2]*alpha[16]+0.25*alpha[8]*fUpwind[12]+0.25*fUpwind[8]*alpha[12]+0.22360679774997896*alpha[5]*fUpwind[9]+0.22360679774997896*fUpwind[5]*alpha[9]; 
  Ghat[37] = 0.2*alpha[16]*fUpwind[47]+0.2*alpha[15]*fUpwind[46]+0.25000000000000006*alpha[2]*fUpwind[44]+0.2*alpha[8]*fUpwind[42]+0.2*alpha[6]*fUpwind[39]+0.15971914124998499*alpha[11]*fUpwind[37]+0.25*alpha[0]*fUpwind[37]+0.25*alpha[7]*fUpwind[35]+0.25*alpha[9]*fUpwind[32]+0.22360679774997896*alpha[5]*fUpwind[31]+0.22360679774997896*alpha[21]*fUpwind[27]+0.15971914124998499*alpha[21]*fUpwind[25]+0.25000000000000006*alpha[3]*fUpwind[25]+0.25000000000000006*alpha[4]*fUpwind[21]+0.25000000000000006*fUpwind[4]*alpha[21]+0.22360679774997896*alpha[1]*fUpwind[17]+0.22360679774997896*alpha[15]*fUpwind[16]+0.22360679774997896*fUpwind[15]*alpha[16]+0.25*fUpwind[10]*alpha[11]+0.22360679774997896*alpha[6]*fUpwind[8]+0.22360679774997896*fUpwind[6]*alpha[8]; 
  Ghat[38] = 0.2*alpha[16]*fUpwind[47]+0.2*alpha[15]*fUpwind[46]+0.25000000000000006*alpha[1]*fUpwind[45]+0.2*alpha[9]*fUpwind[43]+0.2*alpha[7]*fUpwind[40]+0.15971914124998499*alpha[12]*fUpwind[38]+0.25*alpha[0]*fUpwind[38]+0.25*alpha[6]*fUpwind[36]+0.25*alpha[8]*fUpwind[33]+0.22360679774997896*alpha[5]*fUpwind[31]+0.22360679774997896*alpha[26]*fUpwind[30]+0.22360679774997896*alpha[22]*fUpwind[27]+0.15971914124998499*alpha[22]*fUpwind[26]+0.25000000000000006*alpha[3]*fUpwind[26]+0.15971914124998499*fUpwind[22]*alpha[26]+0.25000000000000006*fUpwind[3]*alpha[26]+0.25000000000000006*alpha[4]*fUpwind[22]+0.25000000000000006*fUpwind[4]*alpha[22]+0.22360679774997896*alpha[2]*fUpwind[18]+0.22360679774997896*alpha[15]*fUpwind[16]+0.22360679774997896*fUpwind[15]*alpha[16]+0.25*fUpwind[10]*alpha[12]+0.22360679774997896*alpha[7]*fUpwind[9]+0.22360679774997896*fUpwind[7]*alpha[9]; 
  Ghat[39] = 0.25000000000000006*alpha[2]*fUpwind[46]+0.22360679774997896*alpha[22]*fUpwind[45]+0.2*alpha[15]*fUpwind[44]+0.25*alpha[5]*fUpwind[40]+0.22360679774997896*alpha[11]*fUpwind[39]+0.25*alpha[0]*fUpwind[39]+0.2*alpha[6]*fUpwind[37]+0.25*alpha[9]*fUpwind[34]+0.22360679774997896*alpha[7]*fUpwind[31]+0.25000000000000006*alpha[1]*fUpwind[27]+0.25000000000000006*alpha[16]*fUpwind[24]+0.25000000000000006*alpha[4]*fUpwind[23]+0.2*fUpwind[17]*alpha[21]+0.22360679774997896*alpha[15]*fUpwind[18]+0.22360679774997896*alpha[3]*fUpwind[17]+0.25*alpha[8]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[10]; 
  Ghat[40] = 0.25000000000000006*alpha[1]*fUpwind[46]+0.2*alpha[15]*fUpwind[45]+0.22360679774997896*alpha[21]*fUpwind[44]+0.22360679774997896*alpha[12]*fUpwind[40]+0.25*alpha[0]*fUpwind[40]+0.25*alpha[5]*fUpwind[39]+0.2*alpha[7]*fUpwind[38]+0.25*alpha[8]*fUpwind[34]+0.22360679774997896*alpha[6]*fUpwind[31]+0.25000000000000006*alpha[2]*fUpwind[27]+0.22360679774997896*fUpwind[24]*alpha[26]+0.25000000000000006*alpha[4]*fUpwind[24]+0.25000000000000006*alpha[16]*fUpwind[23]+0.2*fUpwind[18]*alpha[22]+0.22360679774997896*alpha[3]*fUpwind[18]+0.22360679774997896*alpha[15]*fUpwind[17]+0.25*alpha[9]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[10]; 
  Ghat[41] = 0.22360679774997896*alpha[22]*fUpwind[47]+0.22360679774997896*alpha[21]*fUpwind[47]+0.25000000000000006*alpha[3]*fUpwind[47]+0.25*alpha[6]*fUpwind[43]+0.25*alpha[7]*fUpwind[42]+0.22360679774997896*alpha[12]*fUpwind[41]+0.22360679774997896*alpha[11]*fUpwind[41]+0.25*alpha[0]*fUpwind[41]+0.2*alpha[9]*fUpwind[36]+0.2*alpha[8]*fUpwind[35]+0.25000000000000006*alpha[15]*fUpwind[30]+0.25000000000000006*alpha[1]*fUpwind[29]+0.25000000000000006*alpha[2]*fUpwind[28]+0.2*alpha[16]*fUpwind[26]+0.2*fUpwind[16]*alpha[26]+0.2*alpha[16]*fUpwind[25]+0.22360679774997896*alpha[4]*fUpwind[16]+0.22360679774997896*fUpwind[4]*alpha[16]+0.25*alpha[5]*fUpwind[14]+0.22360679774997896*alpha[8]*fUpwind[9]+0.22360679774997896*fUpwind[8]*alpha[9]; 
  Ghat[42] = 0.25000000000000006*alpha[2]*fUpwind[47]+0.22360679774997896*alpha[26]*fUpwind[45]+0.2*alpha[16]*fUpwind[44]+0.25*alpha[5]*fUpwind[43]+0.22360679774997896*alpha[11]*fUpwind[42]+0.25*alpha[0]*fUpwind[42]+0.25*alpha[7]*fUpwind[41]+0.2*alpha[8]*fUpwind[37]+0.22360679774997896*alpha[9]*fUpwind[31]+0.25000000000000006*alpha[1]*fUpwind[30]+0.25000000000000006*alpha[15]*fUpwind[29]+0.22360679774997896*alpha[21]*fUpwind[28]+0.25000000000000006*alpha[3]*fUpwind[28]+0.22360679774997896*alpha[16]*fUpwind[18]+0.22360679774997896*alpha[4]*fUpwind[17]+0.25*alpha[6]*fUpwind[14]+0.22360679774997896*alpha[8]*fUpwind[10]; 
  Ghat[43] = 0.25000000000000006*alpha[1]*fUpwind[47]+0.2*alpha[16]*fUpwind[45]+0.22360679774997896*alpha[12]*fUpwind[43]+0.25*alpha[0]*fUpwind[43]+0.25*alpha[5]*fUpwind[42]+0.25*alpha[6]*fUpwind[41]+0.2*alpha[9]*fUpwind[38]+0.22360679774997896*alpha[8]*fUpwind[31]+0.25000000000000006*alpha[2]*fUpwind[30]+0.22360679774997896*alpha[22]*fUpwind[29]+0.25000000000000006*alpha[3]*fUpwind[29]+0.25000000000000006*alpha[15]*fUpwind[28]+0.2*fUpwind[18]*alpha[26]+0.22360679774997896*alpha[4]*fUpwind[18]+0.22360679774997896*alpha[16]*fUpwind[17]+0.25*alpha[7]*fUpwind[14]+0.22360679774997896*alpha[9]*fUpwind[10]; 
  Ghat[44] = 0.2*alpha[8]*fUpwind[47]+0.2*alpha[6]*fUpwind[46]+0.2*alpha[5]*fUpwind[45]+0.22360679774997896*alpha[12]*fUpwind[44]+0.15971914124998499*alpha[11]*fUpwind[44]+0.25*alpha[0]*fUpwind[44]+0.2*alpha[16]*fUpwind[42]+0.22360679774997896*alpha[21]*fUpwind[40]+0.2*alpha[15]*fUpwind[39]+0.25000000000000006*alpha[2]*fUpwind[37]+0.2*alpha[15]*fUpwind[36]+0.22360679774997896*alpha[22]*fUpwind[35]+0.15971914124998499*alpha[21]*fUpwind[35]+0.25000000000000006*alpha[3]*fUpwind[35]+0.2*alpha[16]*fUpwind[33]+0.22360679774997896*alpha[26]*fUpwind[32]+0.25000000000000006*alpha[4]*fUpwind[32]+0.22360679774997902*alpha[1]*fUpwind[31]+0.25*alpha[7]*fUpwind[25]+0.25*alpha[9]*fUpwind[21]+0.25*fUpwind[9]*alpha[21]+0.25000000000000006*alpha[11]*fUpwind[18]+0.22360679774997902*alpha[5]*fUpwind[17]+0.22360679774997902*alpha[6]*fUpwind[16]+0.22360679774997902*fUpwind[6]*alpha[16]+0.22360679774997902*alpha[8]*fUpwind[15]+0.22360679774997902*fUpwind[8]*alpha[15]; 
  Ghat[45] = 0.2*alpha[9]*fUpwind[47]+0.2*alpha[7]*fUpwind[46]+0.15971914124998499*alpha[12]*fUpwind[45]+0.22360679774997896*alpha[11]*fUpwind[45]+0.25*alpha[0]*fUpwind[45]+0.2*alpha[5]*fUpwind[44]+0.2*alpha[16]*fUpwind[43]+0.22360679774997896*alpha[26]*fUpwind[42]+0.2*alpha[15]*fUpwind[40]+0.22360679774997896*alpha[22]*fUpwind[39]+0.25000000000000006*alpha[1]*fUpwind[38]+0.15971914124998499*alpha[22]*fUpwind[36]+0.22360679774997896*alpha[21]*fUpwind[36]+0.25000000000000006*alpha[3]*fUpwind[36]+0.2*alpha[15]*fUpwind[35]+0.15971914124998499*alpha[26]*fUpwind[33]+0.25000000000000006*alpha[4]*fUpwind[33]+0.2*alpha[16]*fUpwind[32]+0.22360679774997902*alpha[2]*fUpwind[31]+0.25*alpha[6]*fUpwind[26]+0.25*fUpwind[6]*alpha[26]+0.25*alpha[8]*fUpwind[22]+0.25*fUpwind[8]*alpha[22]+0.22360679774997902*alpha[5]*fUpwind[18]+0.25000000000000006*alpha[12]*fUpwind[17]+0.22360679774997902*alpha[7]*fUpwind[16]+0.22360679774997902*fUpwind[7]*alpha[16]+0.22360679774997902*alpha[9]*fUpwind[15]+0.22360679774997902*fUpwind[9]*alpha[15]; 
  Ghat[46] = 0.22360679774997896*alpha[12]*fUpwind[46]+0.22360679774997896*alpha[11]*fUpwind[46]+0.25*alpha[0]*fUpwind[46]+0.2*alpha[7]*fUpwind[45]+0.2*alpha[6]*fUpwind[44]+0.25000000000000006*alpha[1]*fUpwind[40]+0.25000000000000006*alpha[2]*fUpwind[39]+0.2*alpha[15]*fUpwind[38]+0.2*alpha[15]*fUpwind[37]+0.22360679774997896*alpha[26]*fUpwind[34]+0.25000000000000006*alpha[4]*fUpwind[34]+0.2*alpha[22]*fUpwind[31]+0.2*alpha[21]*fUpwind[31]+0.22360679774997902*alpha[3]*fUpwind[31]+0.25*alpha[5]*fUpwind[27]+0.25*alpha[8]*fUpwind[24]+0.25*alpha[9]*fUpwind[23]+0.22360679774997902*alpha[6]*fUpwind[18]+0.22360679774997902*alpha[7]*fUpwind[17]+0.25000000000000006*fUpwind[13]*alpha[16]+0.22360679774997902*fUpwind[10]*alpha[15]; 
  Ghat[47] = 0.22360679774997896*alpha[12]*fUpwind[47]+0.22360679774997896*alpha[11]*fUpwind[47]+0.25*alpha[0]*fUpwind[47]+0.2*alpha[9]*fUpwind[45]+0.2*alpha[8]*fUpwind[44]+0.25000000000000006*alpha[1]*fUpwind[43]+0.25000000000000006*alpha[2]*fUpwind[42]+0.22360679774997896*alpha[22]*fUpwind[41]+0.22360679774997896*alpha[21]*fUpwind[41]+0.25000000000000006*alpha[3]*fUpwind[41]+0.2*alpha[16]*fUpwind[38]+0.2*alpha[16]*fUpwind[37]+0.2*alpha[26]*fUpwind[31]+0.22360679774997902*alpha[4]*fUpwind[31]+0.25*alpha[5]*fUpwind[30]+0.25*alpha[6]*fUpwind[29]+0.25*alpha[7]*fUpwind[28]+0.22360679774997902*alpha[8]*fUpwind[18]+0.22360679774997902*alpha[9]*fUpwind[17]+0.22360679774997902*fUpwind[10]*alpha[16]+0.25000000000000006*fUpwind[14]*alpha[15]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -(1.224744871391589*Ghat[0]*dv10); 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -(1.224744871391589*Ghat[1]*dv10); 
  out[8] += -(1.224744871391589*Ghat[2]*dv10); 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += 0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -(1.224744871391589*Ghat[3]*dv10); 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -(1.224744871391589*Ghat[4]*dv10); 
  out[15] += 0.7071067811865475*Ghat[10]*dv10; 
  out[16] += 0.7071067811865475*Ghat[11]*dv10; 
  out[17] += 0.7071067811865475*Ghat[12]*dv10; 
  out[18] += 1.5811388300841895*Ghat[0]*dv10; 
  out[19] += 0.7071067811865475*Ghat[13]*dv10; 
  out[20] += 0.7071067811865475*Ghat[14]*dv10; 
  out[21] += -(1.224744871391589*Ghat[5]*dv10); 
  out[22] += 0.7071067811865475*Ghat[15]*dv10; 
  out[23] += -(1.224744871391589*Ghat[6]*dv10); 
  out[24] += -(1.224744871391589*Ghat[7]*dv10); 
  out[25] += 0.7071067811865475*Ghat[16]*dv10; 
  out[26] += -(1.224744871391589*Ghat[8]*dv10); 
  out[27] += -(1.224744871391589*Ghat[9]*dv10); 
  out[28] += 0.7071067811865475*Ghat[17]*dv10; 
  out[29] += 0.7071067811865475*Ghat[18]*dv10; 
  out[30] += -(1.224744871391589*Ghat[10]*dv10); 
  out[31] += 0.7071067811865475*Ghat[19]*dv10; 
  out[32] += 0.7071067811865475*Ghat[20]*dv10; 
  out[33] += -(1.224744871391589*Ghat[11]*dv10); 
  out[34] += -(1.224744871391589*Ghat[12]*dv10); 
  out[35] += 1.5811388300841898*Ghat[1]*dv10; 
  out[36] += 1.5811388300841898*Ghat[2]*dv10; 
  out[37] += 0.7071067811865475*Ghat[21]*dv10; 
  out[38] += 0.7071067811865475*Ghat[22]*dv10; 
  out[39] += 1.5811388300841898*Ghat[3]*dv10; 
  out[40] += 0.7071067811865475*Ghat[23]*dv10; 
  out[41] += 0.7071067811865475*Ghat[24]*dv10; 
  out[42] += -(1.224744871391589*Ghat[13]*dv10); 
  out[43] += 0.7071067811865475*Ghat[25]*dv10; 
  out[44] += 0.7071067811865475*Ghat[26]*dv10; 
  out[45] += 1.5811388300841898*Ghat[4]*dv10; 
  out[46] += 0.7071067811865475*Ghat[27]*dv10; 
  out[47] += 0.7071067811865475*Ghat[28]*dv10; 
  out[48] += 0.7071067811865475*Ghat[29]*dv10; 
  out[49] += -(1.224744871391589*Ghat[14]*dv10); 
  out[50] += 0.7071067811865475*Ghat[30]*dv10; 
  out[51] += -(1.224744871391589*Ghat[15]*dv10); 
  out[52] += -(1.224744871391589*Ghat[16]*dv10); 
  out[53] += 0.7071067811865475*Ghat[31]*dv10; 
  out[54] += -(1.224744871391589*Ghat[17]*dv10); 
  out[55] += -(1.224744871391589*Ghat[18]*dv10); 
  out[56] += -(1.224744871391589*Ghat[19]*dv10); 
  out[57] += -(1.224744871391589*Ghat[20]*dv10); 
  out[58] += 1.5811388300841895*Ghat[5]*dv10; 
  out[59] += 0.7071067811865475*Ghat[32]*dv10; 
  out[60] += 0.7071067811865475*Ghat[33]*dv10; 
  out[61] += -(1.224744871391589*Ghat[21]*dv10); 
  out[62] += -(1.224744871391589*Ghat[22]*dv10); 
  out[63] += 1.5811388300841895*Ghat[6]*dv10; 
  out[64] += 1.5811388300841895*Ghat[7]*dv10; 
  out[65] += 0.7071067811865475*Ghat[34]*dv10; 
  out[66] += -(1.224744871391589*Ghat[23]*dv10); 
  out[67] += -(1.224744871391589*Ghat[24]*dv10); 
  out[68] += 0.7071067811865475*Ghat[35]*dv10; 
  out[69] += 0.7071067811865475*Ghat[36]*dv10; 
  out[70] += -(1.224744871391589*Ghat[25]*dv10); 
  out[71] += -(1.224744871391589*Ghat[26]*dv10); 
  out[72] += 1.5811388300841895*Ghat[8]*dv10; 
  out[73] += 1.5811388300841895*Ghat[9]*dv10; 
  out[74] += 0.7071067811865475*Ghat[37]*dv10; 
  out[75] += 0.7071067811865475*Ghat[38]*dv10; 
  out[76] += 1.5811388300841895*Ghat[10]*dv10; 
  out[77] += 0.7071067811865475*Ghat[39]*dv10; 
  out[78] += 0.7071067811865475*Ghat[40]*dv10; 
  out[79] += -(1.224744871391589*Ghat[27]*dv10); 
  out[80] += 0.7071067811865475*Ghat[41]*dv10; 
  out[81] += -(1.224744871391589*Ghat[28]*dv10); 
  out[82] += -(1.224744871391589*Ghat[29]*dv10); 
  out[83] += 0.7071067811865475*Ghat[42]*dv10; 
  out[84] += 0.7071067811865475*Ghat[43]*dv10; 
  out[85] += -(1.224744871391589*Ghat[30]*dv10); 
  out[86] += -(1.224744871391589*Ghat[31]*dv10); 
  out[87] += -(1.224744871391589*Ghat[32]*dv10); 
  out[88] += -(1.224744871391589*Ghat[33]*dv10); 
  out[89] += 1.5811388300841898*Ghat[15]*dv10; 
  out[90] += -(1.224744871391589*Ghat[34]*dv10); 
  out[91] += -(1.224744871391589*Ghat[35]*dv10); 
  out[92] += -(1.224744871391589*Ghat[36]*dv10); 
  out[93] += 1.5811388300841898*Ghat[16]*dv10; 
  out[94] += 0.7071067811865475*Ghat[44]*dv10; 
  out[95] += 0.7071067811865475*Ghat[45]*dv10; 
  out[96] += -(1.224744871391589*Ghat[37]*dv10); 
  out[97] += -(1.224744871391589*Ghat[38]*dv10); 
  out[98] += 1.5811388300841898*Ghat[17]*dv10; 
  out[99] += 1.5811388300841898*Ghat[18]*dv10; 
  out[100] += 0.7071067811865475*Ghat[46]*dv10; 
  out[101] += -(1.224744871391589*Ghat[39]*dv10); 
  out[102] += -(1.224744871391589*Ghat[40]*dv10); 
  out[103] += -(1.224744871391589*Ghat[41]*dv10); 
  out[104] += 0.7071067811865475*Ghat[47]*dv10; 
  out[105] += -(1.224744871391589*Ghat[42]*dv10); 
  out[106] += -(1.224744871391589*Ghat[43]*dv10); 
  out[107] += -(1.224744871391589*Ghat[44]*dv10); 
  out[108] += -(1.224744871391589*Ghat[45]*dv10); 
  out[109] += 1.5811388300841895*Ghat[31]*dv10; 
  out[110] += -(1.224744871391589*Ghat[46]*dv10); 
  out[111] += -(1.224744871391589*Ghat[47]*dv10); 

  } 
  return 0.;

} 
