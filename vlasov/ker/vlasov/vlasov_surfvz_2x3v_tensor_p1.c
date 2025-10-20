#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_tensor_5x_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_tensor_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_surfvz_2x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv12 = 2/dxv[4]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E2 = &field[8]; 
  const double *B0 = &field[12]; 
  const double *B1 = &field[16]; 
  const double *B2 = &field[20]; 

  double alpha[16] = {0.0}; 

  alpha[0] = (-2.0*B0[0]*wv2)+2.0*B1[0]*wv1+2.0*E2[0]; 
  alpha[1] = (-2.0*B0[1]*wv2)+2.0*B1[1]*wv1+2.0*E2[1]; 
  alpha[2] = (-2.0*B0[2]*wv2)+2.0*B1[2]*wv1+2.0*E2[2]; 
  alpha[3] = 0.5773502691896258*B1[0]*dv1; 
  alpha[4] = -0.5773502691896258*B0[0]*dv2; 
  alpha[5] = (-2.0*B0[3]*wv2)+2.0*B1[3]*wv1+2.0*E2[3]; 
  alpha[6] = 0.5773502691896258*B1[1]*dv1; 
  alpha[7] = 0.5773502691896258*B1[2]*dv1; 
  alpha[8] = -0.5773502691896258*B0[1]*dv2; 
  alpha[9] = -0.5773502691896258*B0[2]*dv2; 
  alpha[11] = 0.5773502691896258*B1[3]*dv1; 
  alpha[12] = -0.5773502691896258*B0[3]*dv2; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if ((-0.25*(alpha[12]+alpha[11]))+0.25*(alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.25*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[0] = tensor_5x_p1_surfx5_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = tensor_5x_p1_surfx5_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = tensor_5x_p1_surfx5_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = tensor_5x_p1_surfx5_eval_quad_node_0_l(fr); 
  } 
  if (0.25*alpha[12]-0.25*(alpha[11]+alpha[9]+alpha[8])+0.25*(alpha[7]+alpha[6]+alpha[5]+alpha[4])-0.25*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[1] = tensor_5x_p1_surfx5_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = tensor_5x_p1_surfx5_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = tensor_5x_p1_surfx5_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = tensor_5x_p1_surfx5_eval_quad_node_1_l(fr); 
  } 
  if ((-0.25*alpha[12])+0.25*(alpha[11]+alpha[9]+alpha[8])-0.25*(alpha[7]+alpha[6])+0.25*alpha[5]-0.25*alpha[4]+0.25*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[2] = tensor_5x_p1_surfx5_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = tensor_5x_p1_surfx5_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = tensor_5x_p1_surfx5_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = tensor_5x_p1_surfx5_eval_quad_node_2_l(fr); 
  } 
  if (0.25*(alpha[12]+alpha[11])-0.25*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[3] = tensor_5x_p1_surfx5_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = tensor_5x_p1_surfx5_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = tensor_5x_p1_surfx5_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = tensor_5x_p1_surfx5_eval_quad_node_3_l(fr); 
  } 
  if (0.25*(alpha[12]+alpha[11])-0.25*alpha[9]+0.25*alpha[8]-0.25*alpha[7]+0.25*alpha[6]-0.25*(alpha[5]+alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[4] = tensor_5x_p1_surfx5_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = tensor_5x_p1_surfx5_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = tensor_5x_p1_surfx5_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = tensor_5x_p1_surfx5_eval_quad_node_4_l(fr); 
  } 
  if ((-0.25*alpha[12])+0.25*(alpha[11]+alpha[9])-0.25*(alpha[8]+alpha[7])+0.25*alpha[6]-0.25*alpha[5]+0.25*alpha[4]-0.25*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[5] = tensor_5x_p1_surfx5_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = tensor_5x_p1_surfx5_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = tensor_5x_p1_surfx5_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = tensor_5x_p1_surfx5_eval_quad_node_5_l(fr); 
  } 
  if (0.25*alpha[12]-0.25*(alpha[11]+alpha[9])+0.25*(alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5]+alpha[4])+0.25*(alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[6] = tensor_5x_p1_surfx5_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = tensor_5x_p1_surfx5_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = tensor_5x_p1_surfx5_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = tensor_5x_p1_surfx5_eval_quad_node_6_l(fr); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]))+0.25*alpha[9]-0.25*alpha[8]+0.25*alpha[7]-0.25*(alpha[6]+alpha[5])+0.25*(alpha[4]+alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[7] = tensor_5x_p1_surfx5_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = tensor_5x_p1_surfx5_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = tensor_5x_p1_surfx5_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = tensor_5x_p1_surfx5_eval_quad_node_7_l(fr); 
  } 
  if (0.25*(alpha[12]+alpha[11]+alpha[9])-0.25*alpha[8]+0.25*alpha[7]-0.25*(alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[8] = tensor_5x_p1_surfx5_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = tensor_5x_p1_surfx5_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = tensor_5x_p1_surfx5_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = tensor_5x_p1_surfx5_eval_quad_node_8_l(fr); 
  } 
  if ((-0.25*alpha[12])+0.25*alpha[11]-0.25*alpha[9]+0.25*(alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5])+0.25*alpha[4]-0.25*(alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[9] = tensor_5x_p1_surfx5_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = tensor_5x_p1_surfx5_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = tensor_5x_p1_surfx5_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = tensor_5x_p1_surfx5_eval_quad_node_9_l(fr); 
  } 
  if (0.25*alpha[12]-0.25*alpha[11]+0.25*alpha[9]-0.25*(alpha[8]+alpha[7])+0.25*alpha[6]-0.25*(alpha[5]+alpha[4])+0.25*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[10] = tensor_5x_p1_surfx5_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = tensor_5x_p1_surfx5_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = tensor_5x_p1_surfx5_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = tensor_5x_p1_surfx5_eval_quad_node_10_l(fr); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]+alpha[9]))+0.25*alpha[8]-0.25*alpha[7]+0.25*alpha[6]-0.25*alpha[5]+0.25*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[11] = tensor_5x_p1_surfx5_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = tensor_5x_p1_surfx5_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = tensor_5x_p1_surfx5_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = tensor_5x_p1_surfx5_eval_quad_node_11_l(fr); 
  } 
  if ((-0.25*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.25*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[12] = tensor_5x_p1_surfx5_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = tensor_5x_p1_surfx5_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = tensor_5x_p1_surfx5_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = tensor_5x_p1_surfx5_eval_quad_node_12_l(fr); 
  } 
  if (0.25*alpha[12]-0.25*alpha[11]+0.25*(alpha[9]+alpha[8])-0.25*(alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4])-0.25*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[13] = tensor_5x_p1_surfx5_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = tensor_5x_p1_surfx5_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = tensor_5x_p1_surfx5_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = tensor_5x_p1_surfx5_eval_quad_node_13_l(fr); 
  } 
  if ((-0.25*alpha[12])+0.25*alpha[11]-0.25*(alpha[9]+alpha[8])+0.25*(alpha[7]+alpha[6]+alpha[5])-0.25*alpha[4]+0.25*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[14] = tensor_5x_p1_surfx5_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = tensor_5x_p1_surfx5_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = tensor_5x_p1_surfx5_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = tensor_5x_p1_surfx5_eval_quad_node_14_l(fr); 
  } 
  if (0.25*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[15] = tensor_5x_p1_surfx5_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = tensor_5x_p1_surfx5_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = tensor_5x_p1_surfx5_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = tensor_5x_p1_surfx5_eval_quad_node_15_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_5x_p1_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  tensor_5x_p1_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.25*alpha[12]*fUpwind_l[12]+0.25*alpha[11]*fUpwind_l[11]+0.25*alpha[9]*fUpwind_l[9]+0.25*alpha[8]*fUpwind_l[8]+0.25*alpha[7]*fUpwind_l[7]+0.25*alpha[6]*fUpwind_l[6]+0.25*alpha[5]*fUpwind_l[5]+0.25*alpha[4]*fUpwind_l[4]+0.25*alpha[3]*fUpwind_l[3]+0.25*alpha[2]*fUpwind_l[2]+0.25*alpha[1]*fUpwind_l[1]+0.25*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.25*alpha[9]*fUpwind_l[12]+0.25*fUpwind_l[9]*alpha[12]+0.25*alpha[7]*fUpwind_l[11]+0.25*fUpwind_l[7]*alpha[11]+0.25*alpha[4]*fUpwind_l[8]+0.25*fUpwind_l[4]*alpha[8]+0.25*alpha[3]*fUpwind_l[6]+0.25*fUpwind_l[3]*alpha[6]+0.25*alpha[2]*fUpwind_l[5]+0.25*fUpwind_l[2]*alpha[5]+0.25*alpha[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.25*alpha[8]*fUpwind_l[12]+0.25*fUpwind_l[8]*alpha[12]+0.25*alpha[6]*fUpwind_l[11]+0.25*fUpwind_l[6]*alpha[11]+0.25*alpha[4]*fUpwind_l[9]+0.25*fUpwind_l[4]*alpha[9]+0.25*alpha[3]*fUpwind_l[7]+0.25*fUpwind_l[3]*alpha[7]+0.25*alpha[1]*fUpwind_l[5]+0.25*fUpwind_l[1]*alpha[5]+0.25*alpha[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.25*alpha[12]*fUpwind_l[15]+0.25*alpha[9]*fUpwind_l[14]+0.25*alpha[8]*fUpwind_l[13]+0.25*alpha[5]*fUpwind_l[11]+0.25*fUpwind_l[5]*alpha[11]+0.25*alpha[4]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[7]+0.25*fUpwind_l[2]*alpha[7]+0.25*alpha[1]*fUpwind_l[6]+0.25*fUpwind_l[1]*alpha[6]+0.25*alpha[0]*fUpwind_l[3]+0.25*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.25*alpha[11]*fUpwind_l[15]+0.25*alpha[7]*fUpwind_l[14]+0.25*alpha[6]*fUpwind_l[13]+0.25*alpha[5]*fUpwind_l[12]+0.25*fUpwind_l[5]*alpha[12]+0.25*alpha[3]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[9]+0.25*fUpwind_l[2]*alpha[9]+0.25*alpha[1]*fUpwind_l[8]+0.25*fUpwind_l[1]*alpha[8]+0.25*alpha[0]*fUpwind_l[4]+0.25*fUpwind_l[0]*alpha[4]; 
  Ghat_l[5] = 0.25*alpha[4]*fUpwind_l[12]+0.25*fUpwind_l[4]*alpha[12]+0.25*alpha[3]*fUpwind_l[11]+0.25*fUpwind_l[3]*alpha[11]+0.25*alpha[8]*fUpwind_l[9]+0.25*fUpwind_l[8]*alpha[9]+0.25*alpha[6]*fUpwind_l[7]+0.25*fUpwind_l[6]*alpha[7]+0.25*alpha[0]*fUpwind_l[5]+0.25*fUpwind_l[0]*alpha[5]+0.25*alpha[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alpha[2]; 
  Ghat_l[6] = 0.25*alpha[9]*fUpwind_l[15]+0.25*alpha[12]*fUpwind_l[14]+0.25*alpha[4]*fUpwind_l[13]+0.25*alpha[2]*fUpwind_l[11]+0.25*fUpwind_l[2]*alpha[11]+0.25*alpha[8]*fUpwind_l[10]+0.25*alpha[5]*fUpwind_l[7]+0.25*fUpwind_l[5]*alpha[7]+0.25*alpha[0]*fUpwind_l[6]+0.25*fUpwind_l[0]*alpha[6]+0.25*alpha[1]*fUpwind_l[3]+0.25*fUpwind_l[1]*alpha[3]; 
  Ghat_l[7] = 0.25*alpha[8]*fUpwind_l[15]+0.25*alpha[4]*fUpwind_l[14]+0.25*alpha[12]*fUpwind_l[13]+0.25*alpha[1]*fUpwind_l[11]+0.25*fUpwind_l[1]*alpha[11]+0.25*alpha[9]*fUpwind_l[10]+0.25*alpha[0]*fUpwind_l[7]+0.25*fUpwind_l[0]*alpha[7]+0.25*alpha[5]*fUpwind_l[6]+0.25*fUpwind_l[5]*alpha[6]+0.25*alpha[2]*fUpwind_l[3]+0.25*fUpwind_l[2]*alpha[3]; 
  Ghat_l[8] = 0.25*alpha[7]*fUpwind_l[15]+0.25*alpha[11]*fUpwind_l[14]+0.25*alpha[3]*fUpwind_l[13]+0.25*alpha[2]*fUpwind_l[12]+0.25*fUpwind_l[2]*alpha[12]+0.25*alpha[6]*fUpwind_l[10]+0.25*alpha[5]*fUpwind_l[9]+0.25*fUpwind_l[5]*alpha[9]+0.25*alpha[0]*fUpwind_l[8]+0.25*fUpwind_l[0]*alpha[8]+0.25*alpha[1]*fUpwind_l[4]+0.25*fUpwind_l[1]*alpha[4]; 
  Ghat_l[9] = 0.25*alpha[6]*fUpwind_l[15]+0.25*alpha[3]*fUpwind_l[14]+0.25*alpha[11]*fUpwind_l[13]+0.25*alpha[1]*fUpwind_l[12]+0.25*fUpwind_l[1]*alpha[12]+0.25*alpha[7]*fUpwind_l[10]+0.25*alpha[0]*fUpwind_l[9]+0.25*fUpwind_l[0]*alpha[9]+0.25*alpha[5]*fUpwind_l[8]+0.25*fUpwind_l[5]*alpha[8]+0.25*alpha[2]*fUpwind_l[4]+0.25*fUpwind_l[2]*alpha[4]; 
  Ghat_l[10] = 0.25*alpha[5]*fUpwind_l[15]+0.25*alpha[2]*fUpwind_l[14]+0.25*alpha[1]*fUpwind_l[13]+0.25*alpha[11]*fUpwind_l[12]+0.25*fUpwind_l[11]*alpha[12]+0.25*alpha[0]*fUpwind_l[10]+0.25*alpha[7]*fUpwind_l[9]+0.25*fUpwind_l[7]*alpha[9]+0.25*alpha[6]*fUpwind_l[8]+0.25*fUpwind_l[6]*alpha[8]+0.25*alpha[3]*fUpwind_l[4]+0.25*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.25*alpha[4]*fUpwind_l[15]+0.25*alpha[8]*fUpwind_l[14]+0.25*alpha[9]*fUpwind_l[13]+0.25*fUpwind_l[10]*alpha[12]+0.25*alpha[0]*fUpwind_l[11]+0.25*fUpwind_l[0]*alpha[11]+0.25*alpha[1]*fUpwind_l[7]+0.25*fUpwind_l[1]*alpha[7]+0.25*alpha[2]*fUpwind_l[6]+0.25*fUpwind_l[2]*alpha[6]+0.25*alpha[3]*fUpwind_l[5]+0.25*fUpwind_l[3]*alpha[5]; 
  Ghat_l[12] = 0.25*alpha[3]*fUpwind_l[15]+0.25*alpha[6]*fUpwind_l[14]+0.25*alpha[7]*fUpwind_l[13]+0.25*alpha[0]*fUpwind_l[12]+0.25*fUpwind_l[0]*alpha[12]+0.25*fUpwind_l[10]*alpha[11]+0.25*alpha[1]*fUpwind_l[9]+0.25*fUpwind_l[1]*alpha[9]+0.25*alpha[2]*fUpwind_l[8]+0.25*fUpwind_l[2]*alpha[8]+0.25*alpha[4]*fUpwind_l[5]+0.25*fUpwind_l[4]*alpha[5]; 
  Ghat_l[13] = 0.25*alpha[2]*fUpwind_l[15]+0.25*alpha[5]*fUpwind_l[14]+0.25*alpha[0]*fUpwind_l[13]+0.25*alpha[7]*fUpwind_l[12]+0.25*fUpwind_l[7]*alpha[12]+0.25*alpha[9]*fUpwind_l[11]+0.25*fUpwind_l[9]*alpha[11]+0.25*alpha[1]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[8]+0.25*fUpwind_l[3]*alpha[8]+0.25*alpha[4]*fUpwind_l[6]+0.25*fUpwind_l[4]*alpha[6]; 
  Ghat_l[14] = 0.25*alpha[1]*fUpwind_l[15]+0.25*alpha[0]*fUpwind_l[14]+0.25*alpha[5]*fUpwind_l[13]+0.25*alpha[6]*fUpwind_l[12]+0.25*fUpwind_l[6]*alpha[12]+0.25*alpha[8]*fUpwind_l[11]+0.25*fUpwind_l[8]*alpha[11]+0.25*alpha[2]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[9]+0.25*fUpwind_l[3]*alpha[9]+0.25*alpha[4]*fUpwind_l[7]+0.25*fUpwind_l[4]*alpha[7]; 
  Ghat_l[15] = 0.25*alpha[0]*fUpwind_l[15]+0.25*alpha[1]*fUpwind_l[14]+0.25*alpha[2]*fUpwind_l[13]+0.25*alpha[3]*fUpwind_l[12]+0.25*fUpwind_l[3]*alpha[12]+0.25*alpha[4]*fUpwind_l[11]+0.25*fUpwind_l[4]*alpha[11]+0.25*alpha[5]*fUpwind_l[10]+0.25*alpha[6]*fUpwind_l[9]+0.25*fUpwind_l[6]*alpha[9]+0.25*alpha[7]*fUpwind_l[8]+0.25*fUpwind_l[7]*alpha[8]; 

  Ghat_r[0] = 0.25*alpha[12]*fUpwind_r[12]+0.25*alpha[11]*fUpwind_r[11]+0.25*alpha[9]*fUpwind_r[9]+0.25*alpha[8]*fUpwind_r[8]+0.25*alpha[7]*fUpwind_r[7]+0.25*alpha[6]*fUpwind_r[6]+0.25*alpha[5]*fUpwind_r[5]+0.25*alpha[4]*fUpwind_r[4]+0.25*alpha[3]*fUpwind_r[3]+0.25*alpha[2]*fUpwind_r[2]+0.25*alpha[1]*fUpwind_r[1]+0.25*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.25*alpha[9]*fUpwind_r[12]+0.25*fUpwind_r[9]*alpha[12]+0.25*alpha[7]*fUpwind_r[11]+0.25*fUpwind_r[7]*alpha[11]+0.25*alpha[4]*fUpwind_r[8]+0.25*fUpwind_r[4]*alpha[8]+0.25*alpha[3]*fUpwind_r[6]+0.25*fUpwind_r[3]*alpha[6]+0.25*alpha[2]*fUpwind_r[5]+0.25*fUpwind_r[2]*alpha[5]+0.25*alpha[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.25*alpha[8]*fUpwind_r[12]+0.25*fUpwind_r[8]*alpha[12]+0.25*alpha[6]*fUpwind_r[11]+0.25*fUpwind_r[6]*alpha[11]+0.25*alpha[4]*fUpwind_r[9]+0.25*fUpwind_r[4]*alpha[9]+0.25*alpha[3]*fUpwind_r[7]+0.25*fUpwind_r[3]*alpha[7]+0.25*alpha[1]*fUpwind_r[5]+0.25*fUpwind_r[1]*alpha[5]+0.25*alpha[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.25*alpha[12]*fUpwind_r[15]+0.25*alpha[9]*fUpwind_r[14]+0.25*alpha[8]*fUpwind_r[13]+0.25*alpha[5]*fUpwind_r[11]+0.25*fUpwind_r[5]*alpha[11]+0.25*alpha[4]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[7]+0.25*fUpwind_r[2]*alpha[7]+0.25*alpha[1]*fUpwind_r[6]+0.25*fUpwind_r[1]*alpha[6]+0.25*alpha[0]*fUpwind_r[3]+0.25*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.25*alpha[11]*fUpwind_r[15]+0.25*alpha[7]*fUpwind_r[14]+0.25*alpha[6]*fUpwind_r[13]+0.25*alpha[5]*fUpwind_r[12]+0.25*fUpwind_r[5]*alpha[12]+0.25*alpha[3]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[9]+0.25*fUpwind_r[2]*alpha[9]+0.25*alpha[1]*fUpwind_r[8]+0.25*fUpwind_r[1]*alpha[8]+0.25*alpha[0]*fUpwind_r[4]+0.25*fUpwind_r[0]*alpha[4]; 
  Ghat_r[5] = 0.25*alpha[4]*fUpwind_r[12]+0.25*fUpwind_r[4]*alpha[12]+0.25*alpha[3]*fUpwind_r[11]+0.25*fUpwind_r[3]*alpha[11]+0.25*alpha[8]*fUpwind_r[9]+0.25*fUpwind_r[8]*alpha[9]+0.25*alpha[6]*fUpwind_r[7]+0.25*fUpwind_r[6]*alpha[7]+0.25*alpha[0]*fUpwind_r[5]+0.25*fUpwind_r[0]*alpha[5]+0.25*alpha[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alpha[2]; 
  Ghat_r[6] = 0.25*alpha[9]*fUpwind_r[15]+0.25*alpha[12]*fUpwind_r[14]+0.25*alpha[4]*fUpwind_r[13]+0.25*alpha[2]*fUpwind_r[11]+0.25*fUpwind_r[2]*alpha[11]+0.25*alpha[8]*fUpwind_r[10]+0.25*alpha[5]*fUpwind_r[7]+0.25*fUpwind_r[5]*alpha[7]+0.25*alpha[0]*fUpwind_r[6]+0.25*fUpwind_r[0]*alpha[6]+0.25*alpha[1]*fUpwind_r[3]+0.25*fUpwind_r[1]*alpha[3]; 
  Ghat_r[7] = 0.25*alpha[8]*fUpwind_r[15]+0.25*alpha[4]*fUpwind_r[14]+0.25*alpha[12]*fUpwind_r[13]+0.25*alpha[1]*fUpwind_r[11]+0.25*fUpwind_r[1]*alpha[11]+0.25*alpha[9]*fUpwind_r[10]+0.25*alpha[0]*fUpwind_r[7]+0.25*fUpwind_r[0]*alpha[7]+0.25*alpha[5]*fUpwind_r[6]+0.25*fUpwind_r[5]*alpha[6]+0.25*alpha[2]*fUpwind_r[3]+0.25*fUpwind_r[2]*alpha[3]; 
  Ghat_r[8] = 0.25*alpha[7]*fUpwind_r[15]+0.25*alpha[11]*fUpwind_r[14]+0.25*alpha[3]*fUpwind_r[13]+0.25*alpha[2]*fUpwind_r[12]+0.25*fUpwind_r[2]*alpha[12]+0.25*alpha[6]*fUpwind_r[10]+0.25*alpha[5]*fUpwind_r[9]+0.25*fUpwind_r[5]*alpha[9]+0.25*alpha[0]*fUpwind_r[8]+0.25*fUpwind_r[0]*alpha[8]+0.25*alpha[1]*fUpwind_r[4]+0.25*fUpwind_r[1]*alpha[4]; 
  Ghat_r[9] = 0.25*alpha[6]*fUpwind_r[15]+0.25*alpha[3]*fUpwind_r[14]+0.25*alpha[11]*fUpwind_r[13]+0.25*alpha[1]*fUpwind_r[12]+0.25*fUpwind_r[1]*alpha[12]+0.25*alpha[7]*fUpwind_r[10]+0.25*alpha[0]*fUpwind_r[9]+0.25*fUpwind_r[0]*alpha[9]+0.25*alpha[5]*fUpwind_r[8]+0.25*fUpwind_r[5]*alpha[8]+0.25*alpha[2]*fUpwind_r[4]+0.25*fUpwind_r[2]*alpha[4]; 
  Ghat_r[10] = 0.25*alpha[5]*fUpwind_r[15]+0.25*alpha[2]*fUpwind_r[14]+0.25*alpha[1]*fUpwind_r[13]+0.25*alpha[11]*fUpwind_r[12]+0.25*fUpwind_r[11]*alpha[12]+0.25*alpha[0]*fUpwind_r[10]+0.25*alpha[7]*fUpwind_r[9]+0.25*fUpwind_r[7]*alpha[9]+0.25*alpha[6]*fUpwind_r[8]+0.25*fUpwind_r[6]*alpha[8]+0.25*alpha[3]*fUpwind_r[4]+0.25*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.25*alpha[4]*fUpwind_r[15]+0.25*alpha[8]*fUpwind_r[14]+0.25*alpha[9]*fUpwind_r[13]+0.25*fUpwind_r[10]*alpha[12]+0.25*alpha[0]*fUpwind_r[11]+0.25*fUpwind_r[0]*alpha[11]+0.25*alpha[1]*fUpwind_r[7]+0.25*fUpwind_r[1]*alpha[7]+0.25*alpha[2]*fUpwind_r[6]+0.25*fUpwind_r[2]*alpha[6]+0.25*alpha[3]*fUpwind_r[5]+0.25*fUpwind_r[3]*alpha[5]; 
  Ghat_r[12] = 0.25*alpha[3]*fUpwind_r[15]+0.25*alpha[6]*fUpwind_r[14]+0.25*alpha[7]*fUpwind_r[13]+0.25*alpha[0]*fUpwind_r[12]+0.25*fUpwind_r[0]*alpha[12]+0.25*fUpwind_r[10]*alpha[11]+0.25*alpha[1]*fUpwind_r[9]+0.25*fUpwind_r[1]*alpha[9]+0.25*alpha[2]*fUpwind_r[8]+0.25*fUpwind_r[2]*alpha[8]+0.25*alpha[4]*fUpwind_r[5]+0.25*fUpwind_r[4]*alpha[5]; 
  Ghat_r[13] = 0.25*alpha[2]*fUpwind_r[15]+0.25*alpha[5]*fUpwind_r[14]+0.25*alpha[0]*fUpwind_r[13]+0.25*alpha[7]*fUpwind_r[12]+0.25*fUpwind_r[7]*alpha[12]+0.25*alpha[9]*fUpwind_r[11]+0.25*fUpwind_r[9]*alpha[11]+0.25*alpha[1]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[8]+0.25*fUpwind_r[3]*alpha[8]+0.25*alpha[4]*fUpwind_r[6]+0.25*fUpwind_r[4]*alpha[6]; 
  Ghat_r[14] = 0.25*alpha[1]*fUpwind_r[15]+0.25*alpha[0]*fUpwind_r[14]+0.25*alpha[5]*fUpwind_r[13]+0.25*alpha[6]*fUpwind_r[12]+0.25*fUpwind_r[6]*alpha[12]+0.25*alpha[8]*fUpwind_r[11]+0.25*fUpwind_r[8]*alpha[11]+0.25*alpha[2]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[9]+0.25*fUpwind_r[3]*alpha[9]+0.25*alpha[4]*fUpwind_r[7]+0.25*fUpwind_r[4]*alpha[7]; 
  Ghat_r[15] = 0.25*alpha[0]*fUpwind_r[15]+0.25*alpha[1]*fUpwind_r[14]+0.25*alpha[2]*fUpwind_r[13]+0.25*alpha[3]*fUpwind_r[12]+0.25*fUpwind_r[3]*alpha[12]+0.25*alpha[4]*fUpwind_r[11]+0.25*fUpwind_r[4]*alpha[11]+0.25*alpha[5]*fUpwind_r[10]+0.25*alpha[6]*fUpwind_r[9]+0.25*fUpwind_r[6]*alpha[9]+0.25*alpha[7]*fUpwind_r[8]+0.25*fUpwind_r[7]*alpha[8]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv12; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv12; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv12; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv12; 
  out[4] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv12; 
  out[5] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv12; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv12; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv12; 
  out[8] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv12; 
  out[9] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv12; 
  out[10] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv12; 
  out[11] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv12; 
  out[12] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv12; 
  out[13] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv12; 
  out[14] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv12; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv12; 
  out[16] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv12; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv12; 
  out[18] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv12; 
  out[19] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv12; 
  out[20] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv12; 
  out[21] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv12; 
  out[22] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv12; 
  out[23] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv12; 
  out[24] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv12; 
  out[25] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv12; 
  out[26] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv12; 
  out[27] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv12; 
  out[28] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv12; 
  out[29] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv12; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv12; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv12; 

  return 0.;

} 
