#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void
fpo_vlasov_drag_surfvz_1x3v_ser_p1(const double* w, const double* dxv, const double* hl, const double *hc, const double* hu,const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // hl/hc/hu: Rosenbluth potentials in cells 
  // fl/fc/fu: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = (-0.7654655446197428*hl[4])-0.7654655446197428*hc[4]-0.7954951288348656*hl[0]+0.7954951288348656*hc[0]; 
  alphaDrSurf_l[1] = (-0.7654655446197428*hl[8])-0.7654655446197428*hc[8]-0.7954951288348656*hl[1]+0.7954951288348656*hc[1]; 
  alphaDrSurf_l[2] = (-0.7654655446197428*hl[9])-0.7654655446197428*hc[9]-0.7954951288348656*hl[2]+0.7954951288348656*hc[2]; 
  alphaDrSurf_l[3] = (-0.7654655446197428*hl[10])-0.7654655446197428*hc[10]-0.7954951288348656*hl[3]+0.7954951288348656*hc[3]; 
  alphaDrSurf_l[4] = (-0.7654655446197428*hl[12])-0.7654655446197428*hc[12]-0.7954951288348656*hl[5]+0.7954951288348656*hc[5]; 
  alphaDrSurf_l[5] = (-0.7654655446197428*hl[13])-0.7654655446197428*hc[13]-0.7954951288348656*hl[6]+0.7954951288348656*hc[6]; 
  alphaDrSurf_l[6] = (-0.7654655446197428*hl[14])-0.7654655446197428*hc[14]-0.7954951288348656*hl[7]+0.7954951288348656*hc[7]; 
  alphaDrSurf_l[7] = (-0.7654655446197428*hl[15])-0.7654655446197428*hc[15]-0.7954951288348656*hl[11]+0.7954951288348656*hc[11]; 

  double alphaDrSurf_u[8] = {0.0}; 
  alphaDrSurf_u[0] = (-0.7654655446197428*hu[4])-0.7654655446197428*hc[4]+0.7954951288348656*hu[0]-0.7954951288348656*hc[0]; 
  alphaDrSurf_u[1] = (-0.7654655446197428*hu[8])-0.7654655446197428*hc[8]+0.7954951288348656*hu[1]-0.7954951288348656*hc[1]; 
  alphaDrSurf_u[2] = (-0.7654655446197428*hu[9])-0.7654655446197428*hc[9]+0.7954951288348656*hu[2]-0.7954951288348656*hc[2]; 
  alphaDrSurf_u[3] = (-0.7654655446197428*hu[10])-0.7654655446197428*hc[10]+0.7954951288348656*hu[3]-0.7954951288348656*hc[3]; 
  alphaDrSurf_u[4] = (-0.7654655446197428*hu[12])-0.7654655446197428*hc[12]+0.7954951288348656*hu[5]-0.7954951288348656*hc[5]; 
  alphaDrSurf_u[5] = (-0.7654655446197428*hu[13])-0.7654655446197428*hc[13]+0.7954951288348656*hu[6]-0.7954951288348656*hc[6]; 
  alphaDrSurf_u[6] = (-0.7654655446197428*hu[14])-0.7654655446197428*hc[14]+0.7954951288348656*hu[7]-0.7954951288348656*hc[7]; 
  alphaDrSurf_u[7] = (-0.7654655446197428*hu[15])-0.7654655446197428*hc[15]+0.7954951288348656*hu[11]-0.7954951288348656*hc[11]; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_u[8] = {0.0};
  double fUpwind_l[8] = {0.0};
  double fUpwind_u[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_u[8] = {0.0}; 

  if ((-0.3535533905932737*alphaDrSurf_l[7])+0.3535533905932737*(alphaDrSurf_l[6]+alphaDrSurf_l[5]+alphaDrSurf_l[4])-0.3535533905932737*(alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx4_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx4_eval_quad_node_0_l(fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_u[7])+0.3535533905932737*(alphaDrSurf_u[6]+alphaDrSurf_u[5]+alphaDrSurf_u[4])-0.3535533905932737*(alphaDrSurf_u[3]+alphaDrSurf_u[2]+alphaDrSurf_u[1])+0.3535533905932737*alphaDrSurf_u[0] > 0) { 
    fUpwindQuad_u[0] = ser_4x_p1_surfx4_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_u[0] = ser_4x_p1_surfx4_eval_quad_node_0_l(fu); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[7]-0.3535533905932737*(alphaDrSurf_l[6]+alphaDrSurf_l[5])+0.3535533905932737*(alphaDrSurf_l[4]+alphaDrSurf_l[3])-0.3535533905932737*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx4_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx4_eval_quad_node_1_l(fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_u[7]-0.3535533905932737*(alphaDrSurf_u[6]+alphaDrSurf_u[5])+0.3535533905932737*(alphaDrSurf_u[4]+alphaDrSurf_u[3])-0.3535533905932737*(alphaDrSurf_u[2]+alphaDrSurf_u[1])+0.3535533905932737*alphaDrSurf_u[0] > 0) { 
    fUpwindQuad_u[1] = ser_4x_p1_surfx4_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_u[1] = ser_4x_p1_surfx4_eval_quad_node_1_l(fu); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[7]-0.3535533905932737*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[5]-0.3535533905932737*(alphaDrSurf_l[4]+alphaDrSurf_l[3])+0.3535533905932737*alphaDrSurf_l[2]-0.3535533905932737*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx4_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx4_eval_quad_node_2_l(fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_u[7]-0.3535533905932737*alphaDrSurf_u[6]+0.3535533905932737*alphaDrSurf_u[5]-0.3535533905932737*(alphaDrSurf_u[4]+alphaDrSurf_u[3])+0.3535533905932737*alphaDrSurf_u[2]-0.3535533905932737*alphaDrSurf_u[1]+0.3535533905932737*alphaDrSurf_u[0] > 0) { 
    fUpwindQuad_u[2] = ser_4x_p1_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_u[2] = ser_4x_p1_surfx4_eval_quad_node_2_l(fu); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[7])+0.3535533905932737*alphaDrSurf_l[6]-0.3535533905932737*(alphaDrSurf_l[5]+alphaDrSurf_l[4])+0.3535533905932737*(alphaDrSurf_l[3]+alphaDrSurf_l[2])-0.3535533905932737*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx4_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx4_eval_quad_node_3_l(fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_u[7])+0.3535533905932737*alphaDrSurf_u[6]-0.3535533905932737*(alphaDrSurf_u[5]+alphaDrSurf_u[4])+0.3535533905932737*(alphaDrSurf_u[3]+alphaDrSurf_u[2])-0.3535533905932737*alphaDrSurf_u[1]+0.3535533905932737*alphaDrSurf_u[0] > 0) { 
    fUpwindQuad_u[3] = ser_4x_p1_surfx4_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_u[3] = ser_4x_p1_surfx4_eval_quad_node_3_l(fu); 
  } 
  if (0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6])-0.3535533905932737*(alphaDrSurf_l[5]+alphaDrSurf_l[4]+alphaDrSurf_l[3]+alphaDrSurf_l[2])+0.3535533905932737*(alphaDrSurf_l[1]+alphaDrSurf_l[0]) > 0) { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx4_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx4_eval_quad_node_4_l(fc); 
  } 
  if (0.3535533905932737*(alphaDrSurf_u[7]+alphaDrSurf_u[6])-0.3535533905932737*(alphaDrSurf_u[5]+alphaDrSurf_u[4]+alphaDrSurf_u[3]+alphaDrSurf_u[2])+0.3535533905932737*(alphaDrSurf_u[1]+alphaDrSurf_u[0]) > 0) { 
    fUpwindQuad_u[4] = ser_4x_p1_surfx4_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_u[4] = ser_4x_p1_surfx4_eval_quad_node_4_l(fu); 
  } 
  if ((-0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6]))+0.3535533905932737*alphaDrSurf_l[5]-0.3535533905932737*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[3]-0.3535533905932737*alphaDrSurf_l[2]+0.3535533905932737*(alphaDrSurf_l[1]+alphaDrSurf_l[0]) > 0) { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx4_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx4_eval_quad_node_5_l(fc); 
  } 
  if ((-0.3535533905932737*(alphaDrSurf_u[7]+alphaDrSurf_u[6]))+0.3535533905932737*alphaDrSurf_u[5]-0.3535533905932737*alphaDrSurf_u[4]+0.3535533905932737*alphaDrSurf_u[3]-0.3535533905932737*alphaDrSurf_u[2]+0.3535533905932737*(alphaDrSurf_u[1]+alphaDrSurf_u[0]) > 0) { 
    fUpwindQuad_u[5] = ser_4x_p1_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_u[5] = ser_4x_p1_surfx4_eval_quad_node_5_l(fu); 
  } 
  if ((-0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[5]))+0.3535533905932737*alphaDrSurf_l[4]-0.3535533905932737*alphaDrSurf_l[3]+0.3535533905932737*(alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0]) > 0) { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx4_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx4_eval_quad_node_6_l(fc); 
  } 
  if ((-0.3535533905932737*(alphaDrSurf_u[7]+alphaDrSurf_u[6]+alphaDrSurf_u[5]))+0.3535533905932737*alphaDrSurf_u[4]-0.3535533905932737*alphaDrSurf_u[3]+0.3535533905932737*(alphaDrSurf_u[2]+alphaDrSurf_u[1]+alphaDrSurf_u[0]) > 0) { 
    fUpwindQuad_u[6] = ser_4x_p1_surfx4_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_u[6] = ser_4x_p1_surfx4_eval_quad_node_6_l(fu); 
  } 
  if (0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[5]+alphaDrSurf_l[4]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0]) > 0) { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx4_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx4_eval_quad_node_7_l(fc); 
  } 
  if (0.3535533905932737*(alphaDrSurf_u[7]+alphaDrSurf_u[6]+alphaDrSurf_u[5]+alphaDrSurf_u[4]+alphaDrSurf_u[3]+alphaDrSurf_u[2]+alphaDrSurf_u[1]+alphaDrSurf_u[0]) > 0) { 
    fUpwindQuad_u[7] = ser_4x_p1_surfx4_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_u[7] = ser_4x_p1_surfx4_eval_quad_node_7_l(fu); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad_u, fUpwind_u); 

  Ghat_l[0] = 0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[6]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[5]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[3]*fUpwind_l[3]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[2]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alphaDrSurf_l[6]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[6]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[5]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.3535533905932737*alphaDrSurf_l[5]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[5]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Ghat_l[3] = 0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[4]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[3]; 
  Ghat_l[4] = 0.3535533905932737*alphaDrSurf_l[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Ghat_l[5] = 0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[4]*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[3]; 
  Ghat_l[6] = 0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alphaDrSurf_l[5]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[3]; 
  Ghat_l[7] = 0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[5]+0.3535533905932737*alphaDrSurf_l[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[4]; 

  Ghat_u[0] = 0.3535533905932737*alphaDrSurf_u[7]*fUpwind_u[7]+0.3535533905932737*alphaDrSurf_u[6]*fUpwind_u[6]+0.3535533905932737*alphaDrSurf_u[5]*fUpwind_u[5]+0.3535533905932737*alphaDrSurf_u[4]*fUpwind_u[4]+0.3535533905932737*alphaDrSurf_u[3]*fUpwind_u[3]+0.3535533905932737*alphaDrSurf_u[2]*fUpwind_u[2]+0.3535533905932737*alphaDrSurf_u[1]*fUpwind_u[1]+0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[0]; 
  Ghat_u[1] = 0.3535533905932737*alphaDrSurf_u[6]*fUpwind_u[7]+0.3535533905932737*fUpwind_u[6]*alphaDrSurf_u[7]+0.3535533905932737*alphaDrSurf_u[3]*fUpwind_u[5]+0.3535533905932737*fUpwind_u[3]*alphaDrSurf_u[5]+0.3535533905932737*alphaDrSurf_u[2]*fUpwind_u[4]+0.3535533905932737*fUpwind_u[2]*alphaDrSurf_u[4]+0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[1]+0.3535533905932737*fUpwind_u[0]*alphaDrSurf_u[1]; 
  Ghat_u[2] = 0.3535533905932737*alphaDrSurf_u[5]*fUpwind_u[7]+0.3535533905932737*fUpwind_u[5]*alphaDrSurf_u[7]+0.3535533905932737*alphaDrSurf_u[3]*fUpwind_u[6]+0.3535533905932737*fUpwind_u[3]*alphaDrSurf_u[6]+0.3535533905932737*alphaDrSurf_u[1]*fUpwind_u[4]+0.3535533905932737*fUpwind_u[1]*alphaDrSurf_u[4]+0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[2]+0.3535533905932737*fUpwind_u[0]*alphaDrSurf_u[2]; 
  Ghat_u[3] = 0.3535533905932737*alphaDrSurf_u[4]*fUpwind_u[7]+0.3535533905932737*fUpwind_u[4]*alphaDrSurf_u[7]+0.3535533905932737*alphaDrSurf_u[2]*fUpwind_u[6]+0.3535533905932737*fUpwind_u[2]*alphaDrSurf_u[6]+0.3535533905932737*alphaDrSurf_u[1]*fUpwind_u[5]+0.3535533905932737*fUpwind_u[1]*alphaDrSurf_u[5]+0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[3]+0.3535533905932737*fUpwind_u[0]*alphaDrSurf_u[3]; 
  Ghat_u[4] = 0.3535533905932737*alphaDrSurf_u[3]*fUpwind_u[7]+0.3535533905932737*fUpwind_u[3]*alphaDrSurf_u[7]+0.3535533905932737*alphaDrSurf_u[5]*fUpwind_u[6]+0.3535533905932737*fUpwind_u[5]*alphaDrSurf_u[6]+0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[4]+0.3535533905932737*fUpwind_u[0]*alphaDrSurf_u[4]+0.3535533905932737*alphaDrSurf_u[1]*fUpwind_u[2]+0.3535533905932737*fUpwind_u[1]*alphaDrSurf_u[2]; 
  Ghat_u[5] = 0.3535533905932737*alphaDrSurf_u[2]*fUpwind_u[7]+0.3535533905932737*fUpwind_u[2]*alphaDrSurf_u[7]+0.3535533905932737*alphaDrSurf_u[4]*fUpwind_u[6]+0.3535533905932737*fUpwind_u[4]*alphaDrSurf_u[6]+0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[5]+0.3535533905932737*fUpwind_u[0]*alphaDrSurf_u[5]+0.3535533905932737*alphaDrSurf_u[1]*fUpwind_u[3]+0.3535533905932737*fUpwind_u[1]*alphaDrSurf_u[3]; 
  Ghat_u[6] = 0.3535533905932737*alphaDrSurf_u[1]*fUpwind_u[7]+0.3535533905932737*fUpwind_u[1]*alphaDrSurf_u[7]+0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[6]+0.3535533905932737*fUpwind_u[0]*alphaDrSurf_u[6]+0.3535533905932737*alphaDrSurf_u[4]*fUpwind_u[5]+0.3535533905932737*fUpwind_u[4]*alphaDrSurf_u[5]+0.3535533905932737*alphaDrSurf_u[2]*fUpwind_u[3]+0.3535533905932737*fUpwind_u[2]*alphaDrSurf_u[3]; 
  Ghat_u[7] = 0.3535533905932737*alphaDrSurf_u[0]*fUpwind_u[7]+0.3535533905932737*fUpwind_u[0]*alphaDrSurf_u[7]+0.3535533905932737*alphaDrSurf_u[1]*fUpwind_u[6]+0.3535533905932737*fUpwind_u[1]*alphaDrSurf_u[6]+0.3535533905932737*alphaDrSurf_u[2]*fUpwind_u[5]+0.3535533905932737*fUpwind_u[2]*alphaDrSurf_u[5]+0.3535533905932737*alphaDrSurf_u[3]*fUpwind_u[4]+0.3535533905932737*fUpwind_u[3]*alphaDrSurf_u[4]; 

  out[0] += 0.7071067811865475*Ghat_u[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_u[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_u[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_u[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat_u[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_u[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_u[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_u[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat_u[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat_u[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat_u[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_u[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[12] += 1.224744871391589*Ghat_u[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[13] += 1.224744871391589*Ghat_u[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat_u[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat_u[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
} 
