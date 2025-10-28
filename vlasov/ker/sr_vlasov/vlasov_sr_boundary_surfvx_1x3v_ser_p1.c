#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double *E0 = &qmem[0]; 
  double p1_over_gamma_l[8] = {0.0}; 
  double p1_over_gamma_r[8] = {0.0}; 
  p1_over_gamma_l[0] = 2.738612787525831*gamma[11]*dv11-2.121320343559642*gamma[4]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_l[1] = 2.738612787525831*gamma[8]*dv11-4.743416490252569*gamma[12]*dv11; 
  p1_over_gamma_l[2] = 2.738612787525831*gamma[17]*dv11-2.121320343559642*gamma[10]*dv11+1.224744871391589*gamma[6]*dv11; 
  p1_over_gamma_l[3] = 2.738612787525831*gamma[14]*dv11-4.743416490252569*gamma[18]*dv11; 
  p1_over_gamma_l[5] = 1.224744871391589*gamma[16]*dv11-2.121320343559642*gamma[19]*dv11; 
  p1_over_gamma_r[0] = 2.738612787525831*gamma[11]*dv11+2.121320343559642*gamma[4]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_r[1] = 4.743416490252569*gamma[12]*dv11+2.738612787525831*gamma[8]*dv11; 
  p1_over_gamma_r[2] = 2.738612787525831*gamma[17]*dv11+2.121320343559642*gamma[10]*dv11+1.224744871391589*gamma[6]*dv11; 
  p1_over_gamma_r[3] = 4.743416490252569*gamma[18]*dv11+2.738612787525831*gamma[14]*dv11; 
  p1_over_gamma_r[5] = 2.121320343559642*gamma[19]*dv11+1.224744871391589*gamma[16]*dv11; 

  double p2_over_gamma_l[8] = {0.0}; 
  double p2_over_gamma_r[8] = {0.0}; 
  p2_over_gamma_l[0] = 2.738612787525831*gamma[13]*dv12-2.121320343559642*gamma[5]*dv12+1.224744871391589*gamma[3]*dv12; 
  p2_over_gamma_l[1] = 2.738612787525831*gamma[17]*dv12-2.121320343559642*gamma[10]*dv12+1.224744871391589*gamma[6]*dv12; 
  p2_over_gamma_l[2] = 2.738612787525831*gamma[9]*dv12-4.743416490252569*gamma[15]*dv12; 
  p2_over_gamma_l[3] = 2.738612787525831*gamma[16]*dv12-4.743416490252569*gamma[19]*dv12; 
  p2_over_gamma_l[4] = 1.224744871391589*gamma[14]*dv12-2.121320343559642*gamma[18]*dv12; 
  p2_over_gamma_r[0] = 2.738612787525831*gamma[13]*dv12+2.121320343559642*gamma[5]*dv12+1.224744871391589*gamma[3]*dv12; 
  p2_over_gamma_r[1] = 2.738612787525831*gamma[17]*dv12+2.121320343559642*gamma[10]*dv12+1.224744871391589*gamma[6]*dv12; 
  p2_over_gamma_r[2] = 4.743416490252569*gamma[15]*dv12+2.738612787525831*gamma[9]*dv12; 
  p2_over_gamma_r[3] = 4.743416490252569*gamma[19]*dv12+2.738612787525831*gamma[16]*dv12; 
  p2_over_gamma_r[4] = 2.121320343559642*gamma[18]*dv12+1.224744871391589*gamma[14]*dv12; 

  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double alpha[16] = {0.0}; 

  double fUpwindQuad[18] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  alpha[0] = (-1.0*B1[0]*p2_over_gamma_r[0])+B2[0]*p1_over_gamma_r[0]+2.0*E0[0]; 
  alpha[1] = 2.0*E0[1]+p1_over_gamma_r[0]*B2[1]-1.0*p2_over_gamma_r[0]*B1[1]; 
  alpha[2] = B2[0]*p1_over_gamma_r[1]-1.0*B1[0]*p2_over_gamma_r[1]; 
  alpha[3] = B2[0]*p1_over_gamma_r[2]-1.0*B1[0]*p2_over_gamma_r[2]; 
  alpha[4] = B2[1]*p1_over_gamma_r[1]-1.0*B1[1]*p2_over_gamma_r[1]; 
  alpha[5] = B2[1]*p1_over_gamma_r[2]-1.0*B1[1]*p2_over_gamma_r[2]; 
  alpha[6] = B2[0]*p1_over_gamma_r[3]-1.0*B1[0]*p2_over_gamma_r[3]; 
  alpha[7] = B2[1]*p1_over_gamma_r[3]-1.0*B1[1]*p2_over_gamma_r[3]; 
  alpha[8] = -1.0*B1[0]*p2_over_gamma_r[4]; 
  alpha[9] = -1.0*B1[1]*p2_over_gamma_r[4]; 
  alpha[12] = B2[0]*p1_over_gamma_r[5]; 
  alpha[13] = 1.0*B2[1]*p1_over_gamma_r[5]; 

  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.3952847075210473*alpha[13]-0.3952847075210473*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.4743416490252568*alpha[4]-0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[5]+0.4743416490252568*(alpha[4]+alpha[3])-0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fEdge); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]+0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]+0.4743416490252568*alpha[5]-0.4743416490252568*alpha[3]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fEdge); 
  } 
  if (0.3952847075210473*alpha[13]-0.3952847075210473*alpha[12]+0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fEdge); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]+0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]-0.4743416490252568*alpha[5]+0.4743416490252568*alpha[3]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[5]-0.4743416490252568*(alpha[4]+alpha[3])+0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fEdge); 
  } 
  if (0.3952847075210473*alpha[13]-0.3952847075210473*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.4743416490252568*alpha[4]+0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fEdge); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*(alpha[7]+alpha[6])-0.4743416490252568*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fEdge); 
  } 
  if ((-0.3952847075210473*alpha[13])-0.3952847075210473*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.4743416490252568*(alpha[4]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*(alpha[7]+alpha[6])+0.4743416490252568*alpha[5]-0.4743416490252568*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]-0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]-0.4743416490252568*(alpha[5]+alpha[3])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fEdge); 
  } 
  if ((-0.3952847075210473*alpha[13])-0.3952847075210473*alpha[12]-0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]-0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]+0.4743416490252568*(alpha[5]+alpha[3])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*(alpha[7]+alpha[6])-0.4743416490252568*alpha[5]+0.4743416490252568*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fEdge); 
  } 
  if ((-0.3952847075210473*alpha[13])-0.3952847075210473*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.4743416490252568*(alpha[4]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fEdge); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*(alpha[7]+alpha[6])+0.4743416490252568*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[13]*fUpwind[13]+0.3535533905932737*alpha[12]*fUpwind[12]+0.3535533905932737*alpha[9]*fUpwind[9]+0.3535533905932737*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[7]*fUpwind[7]+0.3535533905932737*alpha[6]*fUpwind[6]+0.3535533905932737*alpha[5]*fUpwind[5]+0.3535533905932737*alpha[4]*fUpwind[4]+0.3535533905932737*alpha[3]*fUpwind[3]+0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[12]*fUpwind[13]+0.3535533905932737*fUpwind[12]*alpha[13]+0.3535533905932737*alpha[8]*fUpwind[9]+0.3535533905932737*fUpwind[8]*alpha[9]+0.3535533905932737*alpha[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.3535533905932737*alpha[13]*fUpwind[15]+0.3535533905932737*alpha[12]*fUpwind[14]+0.3162277660168379*alpha[7]*fUpwind[11]+0.3162277660168379*alpha[6]*fUpwind[10]+0.3162277660168379*alpha[4]*fUpwind[9]+0.3162277660168379*fUpwind[4]*alpha[9]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alpha[8]+0.3535533905932737*alpha[5]*fUpwind[7]+0.3535533905932737*fUpwind[5]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3162277660168379*alpha[7]*fUpwind[15]+0.3162277660168379*alpha[6]*fUpwind[14]+0.3162277660168379*alpha[5]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind[12]+0.3162277660168379*fUpwind[3]*alpha[12]+0.3535533905932737*alpha[9]*fUpwind[11]+0.3535533905932737*alpha[8]*fUpwind[10]+0.3535533905932737*alpha[4]*fUpwind[7]+0.3535533905932737*fUpwind[4]*alpha[7]+0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.3535533905932737*alpha[12]*fUpwind[15]+0.3535533905932737*alpha[13]*fUpwind[14]+0.3162277660168379*alpha[6]*fUpwind[11]+0.3162277660168379*alpha[7]*fUpwind[10]+0.3162277660168379*alpha[2]*fUpwind[9]+0.3162277660168379*fUpwind[2]*alpha[9]+0.3162277660168379*alpha[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alpha[8]+0.3535533905932737*alpha[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3162277660168379*alpha[6]*fUpwind[15]+0.3162277660168379*alpha[7]*fUpwind[14]+0.3162277660168379*alpha[3]*fUpwind[13]+0.3162277660168379*fUpwind[3]*alpha[13]+0.3162277660168379*alpha[5]*fUpwind[12]+0.3162277660168379*fUpwind[5]*alpha[12]+0.3535533905932737*alpha[8]*fUpwind[11]+0.3535533905932737*alpha[9]*fUpwind[10]+0.3535533905932737*alpha[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alpha[3]; 
  Ghat[6] = 0.3162277660168379*alpha[5]*fUpwind[15]+0.3162277660168379*alpha[3]*fUpwind[14]+0.3162277660168379*alpha[7]*fUpwind[13]+0.3162277660168379*fUpwind[7]*alpha[13]+0.3162277660168379*alpha[6]*fUpwind[12]+0.3162277660168379*fUpwind[6]*alpha[12]+0.3162277660168379*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[2]*fUpwind[10]+0.3162277660168379*alpha[7]*fUpwind[9]+0.3162277660168379*fUpwind[7]*alpha[9]+0.3162277660168379*alpha[6]*fUpwind[8]+0.3162277660168379*fUpwind[6]*alpha[8]+0.3535533905932737*alpha[1]*fUpwind[7]+0.3535533905932737*fUpwind[1]*alpha[7]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alpha[3]; 
  Ghat[7] = 0.3162277660168379*alpha[3]*fUpwind[15]+0.3162277660168379*alpha[5]*fUpwind[14]+0.3162277660168379*alpha[6]*fUpwind[13]+0.3162277660168379*fUpwind[6]*alpha[13]+0.3162277660168379*alpha[7]*fUpwind[12]+0.3162277660168379*fUpwind[7]*alpha[12]+0.3162277660168379*alpha[2]*fUpwind[11]+0.3162277660168379*alpha[4]*fUpwind[10]+0.3162277660168379*alpha[6]*fUpwind[9]+0.3162277660168379*fUpwind[6]*alpha[9]+0.3162277660168379*alpha[7]*fUpwind[8]+0.3162277660168379*fUpwind[7]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alpha[7]+0.3535533905932737*alpha[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alpha[4]; 
  Ghat[8] = 0.3535533905932737*alpha[5]*fUpwind[11]+0.3535533905932737*alpha[3]*fUpwind[10]+0.2258769757263128*alpha[9]*fUpwind[9]+0.3535533905932737*alpha[1]*fUpwind[9]+0.3535533905932737*fUpwind[1]*alpha[9]+0.2258769757263128*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alpha[8]+0.3162277660168379*alpha[7]*fUpwind[7]+0.3162277660168379*alpha[6]*fUpwind[6]+0.3162277660168379*alpha[4]*fUpwind[4]+0.3162277660168379*alpha[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alpha[3]*fUpwind[11]+0.3535533905932737*alpha[5]*fUpwind[10]+0.2258769757263128*alpha[8]*fUpwind[9]+0.3535533905932737*alpha[0]*fUpwind[9]+0.2258769757263128*fUpwind[8]*alpha[9]+0.3535533905932737*fUpwind[0]*alpha[9]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alpha[8]+0.3162277660168379*alpha[6]*fUpwind[7]+0.3162277660168379*fUpwind[6]*alpha[7]+0.3162277660168379*alpha[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.282842712474619*alpha[7]*fUpwind[15]+0.2828427124746191*alpha[6]*fUpwind[14]+0.3162277660168379*fUpwind[11]*alpha[13]+0.3162277660168379*fUpwind[10]*alpha[12]+0.2258769757263128*alpha[9]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[11]+0.2258769757263128*alpha[8]*fUpwind[10]+0.3535533905932737*alpha[0]*fUpwind[10]+0.3535533905932737*alpha[5]*fUpwind[9]+0.3535533905932737*fUpwind[5]*alpha[9]+0.3535533905932737*alpha[3]*fUpwind[8]+0.3535533905932737*fUpwind[3]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alpha[7]+0.3162277660168379*alpha[2]*fUpwind[6]+0.3162277660168379*fUpwind[2]*alpha[6]; 
  Ghat[11] = 0.2828427124746191*alpha[6]*fUpwind[15]+0.282842712474619*alpha[7]*fUpwind[14]+0.3162277660168379*fUpwind[10]*alpha[13]+0.3162277660168379*fUpwind[11]*alpha[12]+0.2258769757263128*alpha[8]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[11]+0.2258769757263128*alpha[9]*fUpwind[10]+0.3535533905932737*alpha[1]*fUpwind[10]+0.3535533905932737*alpha[3]*fUpwind[9]+0.3535533905932737*fUpwind[3]*alpha[9]+0.3535533905932737*alpha[5]*fUpwind[8]+0.3535533905932737*fUpwind[5]*alpha[8]+0.3162277660168379*alpha[2]*fUpwind[7]+0.3162277660168379*fUpwind[2]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind[6]+0.3162277660168379*fUpwind[4]*alpha[6]; 
  Ghat[12] = 0.3535533905932737*alpha[4]*fUpwind[15]+0.3535533905932737*alpha[2]*fUpwind[14]+0.2258769757263128*alpha[13]*fUpwind[13]+0.3535533905932737*alpha[1]*fUpwind[13]+0.3535533905932737*fUpwind[1]*alpha[13]+0.2258769757263128*alpha[12]*fUpwind[12]+0.3535533905932737*alpha[0]*fUpwind[12]+0.3535533905932737*fUpwind[0]*alpha[12]+0.3162277660168379*alpha[7]*fUpwind[7]+0.3162277660168379*alpha[6]*fUpwind[6]+0.3162277660168379*alpha[5]*fUpwind[5]+0.3162277660168379*alpha[3]*fUpwind[3]; 
  Ghat[13] = 0.3535533905932737*alpha[2]*fUpwind[15]+0.3535533905932737*alpha[4]*fUpwind[14]+0.2258769757263128*alpha[12]*fUpwind[13]+0.3535533905932737*alpha[0]*fUpwind[13]+0.2258769757263128*fUpwind[12]*alpha[13]+0.3535533905932737*fUpwind[0]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alpha[12]+0.3162277660168379*alpha[6]*fUpwind[7]+0.3162277660168379*fUpwind[6]*alpha[7]+0.3162277660168379*alpha[3]*fUpwind[5]+0.3162277660168379*fUpwind[3]*alpha[5]; 
  Ghat[14] = 0.2258769757263128*alpha[13]*fUpwind[15]+0.3162277660168379*alpha[9]*fUpwind[15]+0.3535533905932737*alpha[1]*fUpwind[15]+0.2258769757263128*alpha[12]*fUpwind[14]+0.3162277660168379*alpha[8]*fUpwind[14]+0.3535533905932737*alpha[0]*fUpwind[14]+0.3535533905932737*alpha[4]*fUpwind[13]+0.3535533905932737*fUpwind[4]*alpha[13]+0.3535533905932737*alpha[2]*fUpwind[12]+0.3535533905932737*fUpwind[2]*alpha[12]+0.282842712474619*alpha[7]*fUpwind[11]+0.2828427124746191*alpha[6]*fUpwind[10]+0.3162277660168379*alpha[5]*fUpwind[7]+0.3162277660168379*fUpwind[5]*alpha[7]+0.3162277660168379*alpha[3]*fUpwind[6]+0.3162277660168379*fUpwind[3]*alpha[6]; 
  Ghat[15] = 0.2258769757263128*alpha[12]*fUpwind[15]+0.3162277660168379*alpha[8]*fUpwind[15]+0.3535533905932737*alpha[0]*fUpwind[15]+0.2258769757263128*alpha[13]*fUpwind[14]+0.3162277660168379*alpha[9]*fUpwind[14]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3535533905932737*alpha[2]*fUpwind[13]+0.3535533905932737*fUpwind[2]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind[12]+0.3535533905932737*fUpwind[4]*alpha[12]+0.2828427124746191*alpha[6]*fUpwind[11]+0.282842712474619*alpha[7]*fUpwind[10]+0.3162277660168379*alpha[3]*fUpwind[7]+0.3162277660168379*fUpwind[3]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind[6]+0.3162277660168379*fUpwind[5]*alpha[6]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -1.224744871391589*Ghat[1]*dv10; 
  out[6] += -0.7071067811865475*Ghat[4]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += -0.7071067811865475*Ghat[5]*dv10; 
  out[9] += -1.224744871391589*Ghat[3]*dv10; 
  out[10] += -0.7071067811865475*Ghat[6]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -1.224744871391589*Ghat[5]*dv10; 
  out[13] += -0.7071067811865475*Ghat[7]*dv10; 
  out[14] += -1.224744871391589*Ghat[6]*dv10; 
  out[15] += -1.224744871391589*Ghat[7]*dv10; 
  out[16] += -1.58113883008419*Ghat[0]*dv10; 
  out[17] += -1.58113883008419*Ghat[1]*dv10; 
  out[18] += -1.58113883008419*Ghat[2]*dv10; 
  out[19] += -1.58113883008419*Ghat[3]*dv10; 
  out[20] += -1.58113883008419*Ghat[4]*dv10; 
  out[21] += -1.58113883008419*Ghat[5]*dv10; 
  out[22] += -1.58113883008419*Ghat[6]*dv10; 
  out[23] += -1.58113883008419*Ghat[7]*dv10; 
  out[24] += -0.7071067811865475*Ghat[8]*dv10; 
  out[25] += -0.7071067811865475*Ghat[9]*dv10; 
  out[26] += -1.224744871391589*Ghat[8]*dv10; 
  out[27] += -0.7071067811865475*Ghat[10]*dv10; 
  out[28] += -1.224744871391589*Ghat[9]*dv10; 
  out[29] += -0.7071067811865475*Ghat[11]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 
  out[32] += -0.7071067811865475*Ghat[12]*dv10; 
  out[33] += -0.7071067811865475*Ghat[13]*dv10; 
  out[34] += -1.224744871391589*Ghat[12]*dv10; 
  out[35] += -0.7071067811865475*Ghat[14]*dv10; 
  out[36] += -1.224744871391589*Ghat[13]*dv10; 
  out[37] += -0.7071067811865475*Ghat[15]*dv10; 
  out[38] += -1.224744871391589*Ghat[14]*dv10; 
  out[39] += -1.224744871391589*Ghat[15]*dv10; 

  } else { 

  alpha[0] = (-1.0*B1[0]*p2_over_gamma_l[0])+B2[0]*p1_over_gamma_l[0]+2.0*E0[0]; 
  alpha[1] = 2.0*E0[1]+p1_over_gamma_l[0]*B2[1]-1.0*p2_over_gamma_l[0]*B1[1]; 
  alpha[2] = B2[0]*p1_over_gamma_l[1]-1.0*B1[0]*p2_over_gamma_l[1]; 
  alpha[3] = B2[0]*p1_over_gamma_l[2]-1.0*B1[0]*p2_over_gamma_l[2]; 
  alpha[4] = B2[1]*p1_over_gamma_l[1]-1.0*B1[1]*p2_over_gamma_l[1]; 
  alpha[5] = B2[1]*p1_over_gamma_l[2]-1.0*B1[1]*p2_over_gamma_l[2]; 
  alpha[6] = B2[0]*p1_over_gamma_l[3]-1.0*B1[0]*p2_over_gamma_l[3]; 
  alpha[7] = B2[1]*p1_over_gamma_l[3]-1.0*B1[1]*p2_over_gamma_l[3]; 
  alpha[8] = -1.0*B1[0]*p2_over_gamma_l[4]; 
  alpha[9] = -1.0*B1[1]*p2_over_gamma_l[4]; 
  alpha[12] = B2[0]*p1_over_gamma_l[5]; 
  alpha[13] = 1.0*B2[1]*p1_over_gamma_l[5]; 

  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.3952847075210473*alpha[13]-0.3952847075210473*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.4743416490252568*alpha[4]-0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[5]+0.4743416490252568*(alpha[4]+alpha[3])-0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fSkin); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]+0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]+0.4743416490252568*alpha[5]-0.4743416490252568*alpha[3]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fSkin); 
  } 
  if (0.3952847075210473*alpha[13]-0.3952847075210473*alpha[12]+0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fSkin); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]+0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]-0.4743416490252568*alpha[5]+0.4743416490252568*alpha[3]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[5]-0.4743416490252568*(alpha[4]+alpha[3])+0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fSkin); 
  } 
  if (0.3952847075210473*alpha[13]-0.3952847075210473*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.4743416490252568*alpha[4]+0.4743416490252568*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fSkin); 
  } 
  if ((-0.3162277660168379*alpha[13])+0.3162277660168379*alpha[12]-0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*(alpha[7]+alpha[6])-0.4743416490252568*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fSkin); 
  } 
  if ((-0.3952847075210473*alpha[13])-0.3952847075210473*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.4743416490252568*(alpha[4]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*(alpha[7]+alpha[6])+0.4743416490252568*alpha[5]-0.4743416490252568*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]-0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]-0.4743416490252568*(alpha[5]+alpha[3])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fSkin); 
  } 
  if ((-0.3952847075210473*alpha[13])-0.3952847075210473*alpha[12]-0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]-0.3952847075210473*alpha[9]-0.3952847075210473*alpha[8]+0.4743416490252568*(alpha[5]+alpha[3])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]-0.6363961030678926*(alpha[7]+alpha[6])-0.4743416490252568*alpha[5]+0.4743416490252568*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fSkin); 
  } 
  if ((-0.3952847075210473*alpha[13])-0.3952847075210473*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.4743416490252568*(alpha[4]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fSkin); 
  } 
  if (0.3162277660168379*alpha[13]+0.3162277660168379*alpha[12]+0.3162277660168379*alpha[9]+0.3162277660168379*alpha[8]+0.6363961030678926*(alpha[7]+alpha[6])+0.4743416490252568*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[13]*fUpwind[13]+0.3535533905932737*alpha[12]*fUpwind[12]+0.3535533905932737*alpha[9]*fUpwind[9]+0.3535533905932737*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[7]*fUpwind[7]+0.3535533905932737*alpha[6]*fUpwind[6]+0.3535533905932737*alpha[5]*fUpwind[5]+0.3535533905932737*alpha[4]*fUpwind[4]+0.3535533905932737*alpha[3]*fUpwind[3]+0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[1]*fUpwind[1]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[12]*fUpwind[13]+0.3535533905932737*fUpwind[12]*alpha[13]+0.3535533905932737*alpha[8]*fUpwind[9]+0.3535533905932737*fUpwind[8]*alpha[9]+0.3535533905932737*alpha[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.3535533905932737*alpha[13]*fUpwind[15]+0.3535533905932737*alpha[12]*fUpwind[14]+0.3162277660168379*alpha[7]*fUpwind[11]+0.3162277660168379*alpha[6]*fUpwind[10]+0.3162277660168379*alpha[4]*fUpwind[9]+0.3162277660168379*fUpwind[4]*alpha[9]+0.3162277660168379*alpha[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alpha[8]+0.3535533905932737*alpha[5]*fUpwind[7]+0.3535533905932737*fUpwind[5]*alpha[7]+0.3535533905932737*alpha[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3162277660168379*alpha[7]*fUpwind[15]+0.3162277660168379*alpha[6]*fUpwind[14]+0.3162277660168379*alpha[5]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alpha[13]+0.3162277660168379*alpha[3]*fUpwind[12]+0.3162277660168379*fUpwind[3]*alpha[12]+0.3535533905932737*alpha[9]*fUpwind[11]+0.3535533905932737*alpha[8]*fUpwind[10]+0.3535533905932737*alpha[4]*fUpwind[7]+0.3535533905932737*fUpwind[4]*alpha[7]+0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.3535533905932737*alpha[12]*fUpwind[15]+0.3535533905932737*alpha[13]*fUpwind[14]+0.3162277660168379*alpha[6]*fUpwind[11]+0.3162277660168379*alpha[7]*fUpwind[10]+0.3162277660168379*alpha[2]*fUpwind[9]+0.3162277660168379*fUpwind[2]*alpha[9]+0.3162277660168379*alpha[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alpha[8]+0.3535533905932737*alpha[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alpha[7]+0.3535533905932737*alpha[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3162277660168379*alpha[6]*fUpwind[15]+0.3162277660168379*alpha[7]*fUpwind[14]+0.3162277660168379*alpha[3]*fUpwind[13]+0.3162277660168379*fUpwind[3]*alpha[13]+0.3162277660168379*alpha[5]*fUpwind[12]+0.3162277660168379*fUpwind[5]*alpha[12]+0.3535533905932737*alpha[8]*fUpwind[11]+0.3535533905932737*alpha[9]*fUpwind[10]+0.3535533905932737*alpha[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alpha[3]; 
  Ghat[6] = 0.3162277660168379*alpha[5]*fUpwind[15]+0.3162277660168379*alpha[3]*fUpwind[14]+0.3162277660168379*alpha[7]*fUpwind[13]+0.3162277660168379*fUpwind[7]*alpha[13]+0.3162277660168379*alpha[6]*fUpwind[12]+0.3162277660168379*fUpwind[6]*alpha[12]+0.3162277660168379*alpha[4]*fUpwind[11]+0.3162277660168379*alpha[2]*fUpwind[10]+0.3162277660168379*alpha[7]*fUpwind[9]+0.3162277660168379*fUpwind[7]*alpha[9]+0.3162277660168379*alpha[6]*fUpwind[8]+0.3162277660168379*fUpwind[6]*alpha[8]+0.3535533905932737*alpha[1]*fUpwind[7]+0.3535533905932737*fUpwind[1]*alpha[7]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alpha[6]+0.3535533905932737*alpha[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alpha[3]; 
  Ghat[7] = 0.3162277660168379*alpha[3]*fUpwind[15]+0.3162277660168379*alpha[5]*fUpwind[14]+0.3162277660168379*alpha[6]*fUpwind[13]+0.3162277660168379*fUpwind[6]*alpha[13]+0.3162277660168379*alpha[7]*fUpwind[12]+0.3162277660168379*fUpwind[7]*alpha[12]+0.3162277660168379*alpha[2]*fUpwind[11]+0.3162277660168379*alpha[4]*fUpwind[10]+0.3162277660168379*alpha[6]*fUpwind[9]+0.3162277660168379*fUpwind[6]*alpha[9]+0.3162277660168379*alpha[7]*fUpwind[8]+0.3162277660168379*fUpwind[7]*alpha[8]+0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alpha[7]+0.3535533905932737*alpha[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alpha[4]; 
  Ghat[8] = 0.3535533905932737*alpha[5]*fUpwind[11]+0.3535533905932737*alpha[3]*fUpwind[10]+0.2258769757263128*alpha[9]*fUpwind[9]+0.3535533905932737*alpha[1]*fUpwind[9]+0.3535533905932737*fUpwind[1]*alpha[9]+0.2258769757263128*alpha[8]*fUpwind[8]+0.3535533905932737*alpha[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alpha[8]+0.3162277660168379*alpha[7]*fUpwind[7]+0.3162277660168379*alpha[6]*fUpwind[6]+0.3162277660168379*alpha[4]*fUpwind[4]+0.3162277660168379*alpha[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alpha[3]*fUpwind[11]+0.3535533905932737*alpha[5]*fUpwind[10]+0.2258769757263128*alpha[8]*fUpwind[9]+0.3535533905932737*alpha[0]*fUpwind[9]+0.2258769757263128*fUpwind[8]*alpha[9]+0.3535533905932737*fUpwind[0]*alpha[9]+0.3535533905932737*alpha[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alpha[8]+0.3162277660168379*alpha[6]*fUpwind[7]+0.3162277660168379*fUpwind[6]*alpha[7]+0.3162277660168379*alpha[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.282842712474619*alpha[7]*fUpwind[15]+0.2828427124746191*alpha[6]*fUpwind[14]+0.3162277660168379*fUpwind[11]*alpha[13]+0.3162277660168379*fUpwind[10]*alpha[12]+0.2258769757263128*alpha[9]*fUpwind[11]+0.3535533905932737*alpha[1]*fUpwind[11]+0.2258769757263128*alpha[8]*fUpwind[10]+0.3535533905932737*alpha[0]*fUpwind[10]+0.3535533905932737*alpha[5]*fUpwind[9]+0.3535533905932737*fUpwind[5]*alpha[9]+0.3535533905932737*alpha[3]*fUpwind[8]+0.3535533905932737*fUpwind[3]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alpha[7]+0.3162277660168379*alpha[2]*fUpwind[6]+0.3162277660168379*fUpwind[2]*alpha[6]; 
  Ghat[11] = 0.2828427124746191*alpha[6]*fUpwind[15]+0.282842712474619*alpha[7]*fUpwind[14]+0.3162277660168379*fUpwind[10]*alpha[13]+0.3162277660168379*fUpwind[11]*alpha[12]+0.2258769757263128*alpha[8]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[11]+0.2258769757263128*alpha[9]*fUpwind[10]+0.3535533905932737*alpha[1]*fUpwind[10]+0.3535533905932737*alpha[3]*fUpwind[9]+0.3535533905932737*fUpwind[3]*alpha[9]+0.3535533905932737*alpha[5]*fUpwind[8]+0.3535533905932737*fUpwind[5]*alpha[8]+0.3162277660168379*alpha[2]*fUpwind[7]+0.3162277660168379*fUpwind[2]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind[6]+0.3162277660168379*fUpwind[4]*alpha[6]; 
  Ghat[12] = 0.3535533905932737*alpha[4]*fUpwind[15]+0.3535533905932737*alpha[2]*fUpwind[14]+0.2258769757263128*alpha[13]*fUpwind[13]+0.3535533905932737*alpha[1]*fUpwind[13]+0.3535533905932737*fUpwind[1]*alpha[13]+0.2258769757263128*alpha[12]*fUpwind[12]+0.3535533905932737*alpha[0]*fUpwind[12]+0.3535533905932737*fUpwind[0]*alpha[12]+0.3162277660168379*alpha[7]*fUpwind[7]+0.3162277660168379*alpha[6]*fUpwind[6]+0.3162277660168379*alpha[5]*fUpwind[5]+0.3162277660168379*alpha[3]*fUpwind[3]; 
  Ghat[13] = 0.3535533905932737*alpha[2]*fUpwind[15]+0.3535533905932737*alpha[4]*fUpwind[14]+0.2258769757263128*alpha[12]*fUpwind[13]+0.3535533905932737*alpha[0]*fUpwind[13]+0.2258769757263128*fUpwind[12]*alpha[13]+0.3535533905932737*fUpwind[0]*alpha[13]+0.3535533905932737*alpha[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alpha[12]+0.3162277660168379*alpha[6]*fUpwind[7]+0.3162277660168379*fUpwind[6]*alpha[7]+0.3162277660168379*alpha[3]*fUpwind[5]+0.3162277660168379*fUpwind[3]*alpha[5]; 
  Ghat[14] = 0.2258769757263128*alpha[13]*fUpwind[15]+0.3162277660168379*alpha[9]*fUpwind[15]+0.3535533905932737*alpha[1]*fUpwind[15]+0.2258769757263128*alpha[12]*fUpwind[14]+0.3162277660168379*alpha[8]*fUpwind[14]+0.3535533905932737*alpha[0]*fUpwind[14]+0.3535533905932737*alpha[4]*fUpwind[13]+0.3535533905932737*fUpwind[4]*alpha[13]+0.3535533905932737*alpha[2]*fUpwind[12]+0.3535533905932737*fUpwind[2]*alpha[12]+0.282842712474619*alpha[7]*fUpwind[11]+0.2828427124746191*alpha[6]*fUpwind[10]+0.3162277660168379*alpha[5]*fUpwind[7]+0.3162277660168379*fUpwind[5]*alpha[7]+0.3162277660168379*alpha[3]*fUpwind[6]+0.3162277660168379*fUpwind[3]*alpha[6]; 
  Ghat[15] = 0.2258769757263128*alpha[12]*fUpwind[15]+0.3162277660168379*alpha[8]*fUpwind[15]+0.3535533905932737*alpha[0]*fUpwind[15]+0.2258769757263128*alpha[13]*fUpwind[14]+0.3162277660168379*alpha[9]*fUpwind[14]+0.3535533905932737*alpha[1]*fUpwind[14]+0.3535533905932737*alpha[2]*fUpwind[13]+0.3535533905932737*fUpwind[2]*alpha[13]+0.3535533905932737*alpha[4]*fUpwind[12]+0.3535533905932737*fUpwind[4]*alpha[12]+0.2828427124746191*alpha[6]*fUpwind[11]+0.282842712474619*alpha[7]*fUpwind[10]+0.3162277660168379*alpha[3]*fUpwind[7]+0.3162277660168379*fUpwind[3]*alpha[7]+0.3162277660168379*alpha[5]*fUpwind[6]+0.3162277660168379*fUpwind[5]*alpha[6]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -1.224744871391589*Ghat[1]*dv10; 
  out[6] += 0.7071067811865475*Ghat[4]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += 0.7071067811865475*Ghat[5]*dv10; 
  out[9] += -1.224744871391589*Ghat[3]*dv10; 
  out[10] += 0.7071067811865475*Ghat[6]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -1.224744871391589*Ghat[5]*dv10; 
  out[13] += 0.7071067811865475*Ghat[7]*dv10; 
  out[14] += -1.224744871391589*Ghat[6]*dv10; 
  out[15] += -1.224744871391589*Ghat[7]*dv10; 
  out[16] += 1.58113883008419*Ghat[0]*dv10; 
  out[17] += 1.58113883008419*Ghat[1]*dv10; 
  out[18] += 1.58113883008419*Ghat[2]*dv10; 
  out[19] += 1.58113883008419*Ghat[3]*dv10; 
  out[20] += 1.58113883008419*Ghat[4]*dv10; 
  out[21] += 1.58113883008419*Ghat[5]*dv10; 
  out[22] += 1.58113883008419*Ghat[6]*dv10; 
  out[23] += 1.58113883008419*Ghat[7]*dv10; 
  out[24] += 0.7071067811865475*Ghat[8]*dv10; 
  out[25] += 0.7071067811865475*Ghat[9]*dv10; 
  out[26] += -1.224744871391589*Ghat[8]*dv10; 
  out[27] += 0.7071067811865475*Ghat[10]*dv10; 
  out[28] += -1.224744871391589*Ghat[9]*dv10; 
  out[29] += 0.7071067811865475*Ghat[11]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 
  out[32] += 0.7071067811865475*Ghat[12]*dv10; 
  out[33] += 0.7071067811865475*Ghat[13]*dv10; 
  out[34] += -1.224744871391589*Ghat[12]*dv10; 
  out[35] += 0.7071067811865475*Ghat[14]*dv10; 
  out[36] += -1.224744871391589*Ghat[13]*dv10; 
  out[37] += 0.7071067811865475*Ghat[15]*dv10; 
  out[38] += -1.224744871391589*Ghat[14]*dv10; 
  out[39] += -1.224744871391589*Ghat[15]*dv10; 

  } 
  return 0.;

} 
