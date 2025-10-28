#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_3x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[5]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[16]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2 = 2.0/dxv[3]; 

  double alphaDrSurf_l[16] = {0.0}; 
  alphaDrSurf_l[0] = -(1.7320508075688772*nuSum[0]*vmap[1])+nuSum[0]*vmap[0]-1.4142135623730951*nuUSum[0]; 
  alphaDrSurf_l[1] = -(1.7320508075688772*nuSum[1]*vmap[1])-1.4142135623730951*nuUSum[1]+vmap[0]*nuSum[1]; 
  alphaDrSurf_l[2] = -(1.4142135623730951*nuUSum[2])-1.7320508075688772*vmap[1]*nuSum[2]+vmap[0]*nuSum[2]; 
  alphaDrSurf_l[3] = -(1.4142135623730951*nuUSum[3])-1.7320508075688772*vmap[1]*nuSum[3]+vmap[0]*nuSum[3]; 
  alphaDrSurf_l[5] = -(1.4142135623730951*nuUSum[4])-1.7320508075688772*vmap[1]*nuSum[4]+vmap[0]*nuSum[4]; 
  alphaDrSurf_l[6] = -(1.4142135623730951*nuUSum[5])-1.7320508075688772*vmap[1]*nuSum[5]+vmap[0]*nuSum[5]; 
  alphaDrSurf_l[7] = -(1.4142135623730951*nuUSum[6])-1.7320508075688772*vmap[1]*nuSum[6]+vmap[0]*nuSum[6]; 
  alphaDrSurf_l[11] = -(1.4142135623730951*nuUSum[7])-1.7320508075688772*vmap[1]*nuSum[7]+vmap[0]*nuSum[7]; 

  double alphaDrSurf_r[16] = {0.0}; 
  alphaDrSurf_r[0] = 1.7320508075688772*nuSum[0]*vmap[1]+nuSum[0]*vmap[0]-1.4142135623730951*nuUSum[0]; 
  alphaDrSurf_r[1] = 1.7320508075688772*nuSum[1]*vmap[1]-1.4142135623730951*nuUSum[1]+vmap[0]*nuSum[1]; 
  alphaDrSurf_r[2] = -(1.4142135623730951*nuUSum[2])+1.7320508075688772*vmap[1]*nuSum[2]+vmap[0]*nuSum[2]; 
  alphaDrSurf_r[3] = -(1.4142135623730951*nuUSum[3])+1.7320508075688772*vmap[1]*nuSum[3]+vmap[0]*nuSum[3]; 
  alphaDrSurf_r[5] = -(1.4142135623730951*nuUSum[4])+1.7320508075688772*vmap[1]*nuSum[4]+vmap[0]*nuSum[4]; 
  alphaDrSurf_r[6] = -(1.4142135623730951*nuUSum[5])+1.7320508075688772*vmap[1]*nuSum[5]+vmap[0]*nuSum[5]; 
  alphaDrSurf_r[7] = -(1.4142135623730951*nuUSum[6])+1.7320508075688772*vmap[1]*nuSum[6]+vmap[0]*nuSum[6]; 
  alphaDrSurf_r[11] = -(1.4142135623730951*nuUSum[7])+1.7320508075688772*vmap[1]*nuSum[7]+vmap[0]*nuSum[7]; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if (-alphaDrSurf_l[11]+alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[5]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fc)/vmap_prime_c[0]; 
  } 
  if (-alphaDrSurf_r[11]+alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[5]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fr)/vmap_prime_r[0]; 
  } 
  if (alphaDrSurf_l[11]-alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[5]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[11]-alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[5]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fr)/vmap_prime_r[0]; 
  } 
  if (alphaDrSurf_l[11]-alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[5]-alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[11]-alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[5]-alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fr)/vmap_prime_r[0]; 
  } 
  if (-alphaDrSurf_l[11]+alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[5]+alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fc)/vmap_prime_c[0]; 
  } 
  if (-alphaDrSurf_r[11]+alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[5]+alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fr)/vmap_prime_r[0]; 
  } 
  if (alphaDrSurf_l[11]+alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[5]-alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[11]+alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[5]-alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fr)/vmap_prime_r[0]; 
  } 
  if (-alphaDrSurf_l[11]-alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[5]+alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fc)/vmap_prime_c[0]; 
  } 
  if (-alphaDrSurf_r[11]-alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[5]+alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fr)/vmap_prime_r[0]; 
  } 
  if (-alphaDrSurf_l[11]-alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[5]-alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fc)/vmap_prime_c[0]; 
  } 
  if (-alphaDrSurf_r[11]-alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[5]-alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fr)/vmap_prime_r[0]; 
  } 
  if (alphaDrSurf_l[11]+alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[5]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[11]+alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[5]+alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fr)/vmap_prime_r[0]; 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.25*(alphaDrSurf_l[11]*fUpwind_l[11]+alphaDrSurf_l[7]*fUpwind_l[7]+alphaDrSurf_l[6]*fUpwind_l[6]+alphaDrSurf_l[5]*fUpwind_l[5]+alphaDrSurf_l[3]*fUpwind_l[3]+alphaDrSurf_l[2]*fUpwind_l[2]+alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.25*(alphaDrSurf_l[7]*fUpwind_l[11]+fUpwind_l[7]*alphaDrSurf_l[11]+alphaDrSurf_l[3]*fUpwind_l[6]+fUpwind_l[3]*alphaDrSurf_l[6]+alphaDrSurf_l[2]*fUpwind_l[5]+fUpwind_l[2]*alphaDrSurf_l[5]+alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.25*(alphaDrSurf_l[6]*fUpwind_l[11]+fUpwind_l[6]*alphaDrSurf_l[11]+alphaDrSurf_l[3]*fUpwind_l[7]+fUpwind_l[3]*alphaDrSurf_l[7]+alphaDrSurf_l[1]*fUpwind_l[5]+fUpwind_l[1]*alphaDrSurf_l[5]+alphaDrSurf_l[0]*fUpwind_l[2]+fUpwind_l[0]*alphaDrSurf_l[2]); 
  Ghat_l[3] = 0.25*(alphaDrSurf_l[5]*fUpwind_l[11]+fUpwind_l[5]*alphaDrSurf_l[11]+alphaDrSurf_l[2]*fUpwind_l[7]+fUpwind_l[2]*alphaDrSurf_l[7]+alphaDrSurf_l[1]*fUpwind_l[6]+fUpwind_l[1]*alphaDrSurf_l[6]+alphaDrSurf_l[0]*fUpwind_l[3]+fUpwind_l[0]*alphaDrSurf_l[3]); 
  Ghat_l[4] = 0.25*(alphaDrSurf_l[11]*fUpwind_l[15]+alphaDrSurf_l[7]*fUpwind_l[14]+alphaDrSurf_l[6]*fUpwind_l[13]+alphaDrSurf_l[5]*fUpwind_l[12]+alphaDrSurf_l[3]*fUpwind_l[10]+alphaDrSurf_l[2]*fUpwind_l[9]+alphaDrSurf_l[1]*fUpwind_l[8]+alphaDrSurf_l[0]*fUpwind_l[4]); 
  Ghat_l[5] = 0.25*(alphaDrSurf_l[3]*fUpwind_l[11]+fUpwind_l[3]*alphaDrSurf_l[11]+alphaDrSurf_l[6]*fUpwind_l[7]+fUpwind_l[6]*alphaDrSurf_l[7]+alphaDrSurf_l[0]*fUpwind_l[5]+fUpwind_l[0]*alphaDrSurf_l[5]+alphaDrSurf_l[1]*fUpwind_l[2]+fUpwind_l[1]*alphaDrSurf_l[2]); 
  Ghat_l[6] = 0.25*(alphaDrSurf_l[2]*fUpwind_l[11]+fUpwind_l[2]*alphaDrSurf_l[11]+alphaDrSurf_l[5]*fUpwind_l[7]+fUpwind_l[5]*alphaDrSurf_l[7]+alphaDrSurf_l[0]*fUpwind_l[6]+fUpwind_l[0]*alphaDrSurf_l[6]+alphaDrSurf_l[1]*fUpwind_l[3]+fUpwind_l[1]*alphaDrSurf_l[3]); 
  Ghat_l[7] = 0.25*(alphaDrSurf_l[1]*fUpwind_l[11]+fUpwind_l[1]*alphaDrSurf_l[11]+alphaDrSurf_l[0]*fUpwind_l[7]+fUpwind_l[0]*alphaDrSurf_l[7]+alphaDrSurf_l[5]*fUpwind_l[6]+fUpwind_l[5]*alphaDrSurf_l[6]+alphaDrSurf_l[2]*fUpwind_l[3]+fUpwind_l[2]*alphaDrSurf_l[3]); 
  Ghat_l[8] = 0.25*(alphaDrSurf_l[7]*fUpwind_l[15]+alphaDrSurf_l[11]*fUpwind_l[14]+alphaDrSurf_l[3]*fUpwind_l[13]+alphaDrSurf_l[2]*fUpwind_l[12]+alphaDrSurf_l[6]*fUpwind_l[10]+alphaDrSurf_l[5]*fUpwind_l[9]+alphaDrSurf_l[0]*fUpwind_l[8]+alphaDrSurf_l[1]*fUpwind_l[4]); 
  Ghat_l[9] = 0.25*(alphaDrSurf_l[6]*fUpwind_l[15]+alphaDrSurf_l[3]*fUpwind_l[14]+alphaDrSurf_l[11]*fUpwind_l[13]+alphaDrSurf_l[1]*fUpwind_l[12]+alphaDrSurf_l[7]*fUpwind_l[10]+alphaDrSurf_l[0]*fUpwind_l[9]+alphaDrSurf_l[5]*fUpwind_l[8]+alphaDrSurf_l[2]*fUpwind_l[4]); 
  Ghat_l[10] = 0.25*(alphaDrSurf_l[5]*fUpwind_l[15]+alphaDrSurf_l[2]*fUpwind_l[14]+alphaDrSurf_l[1]*fUpwind_l[13]+alphaDrSurf_l[11]*fUpwind_l[12]+alphaDrSurf_l[0]*fUpwind_l[10]+alphaDrSurf_l[7]*fUpwind_l[9]+alphaDrSurf_l[6]*fUpwind_l[8]+alphaDrSurf_l[3]*fUpwind_l[4]); 
  Ghat_l[11] = 0.25*(alphaDrSurf_l[0]*fUpwind_l[11]+fUpwind_l[0]*alphaDrSurf_l[11]+alphaDrSurf_l[1]*fUpwind_l[7]+fUpwind_l[1]*alphaDrSurf_l[7]+alphaDrSurf_l[2]*fUpwind_l[6]+fUpwind_l[2]*alphaDrSurf_l[6]+alphaDrSurf_l[3]*fUpwind_l[5]+fUpwind_l[3]*alphaDrSurf_l[5]); 
  Ghat_l[12] = 0.25*(alphaDrSurf_l[3]*fUpwind_l[15]+alphaDrSurf_l[6]*fUpwind_l[14]+alphaDrSurf_l[7]*fUpwind_l[13]+alphaDrSurf_l[0]*fUpwind_l[12]+fUpwind_l[10]*alphaDrSurf_l[11]+alphaDrSurf_l[1]*fUpwind_l[9]+alphaDrSurf_l[2]*fUpwind_l[8]+fUpwind_l[4]*alphaDrSurf_l[5]); 
  Ghat_l[13] = 0.25*(alphaDrSurf_l[2]*fUpwind_l[15]+alphaDrSurf_l[5]*fUpwind_l[14]+alphaDrSurf_l[0]*fUpwind_l[13]+alphaDrSurf_l[7]*fUpwind_l[12]+fUpwind_l[9]*alphaDrSurf_l[11]+alphaDrSurf_l[1]*fUpwind_l[10]+alphaDrSurf_l[3]*fUpwind_l[8]+fUpwind_l[4]*alphaDrSurf_l[6]); 
  Ghat_l[14] = 0.25*(alphaDrSurf_l[1]*fUpwind_l[15]+alphaDrSurf_l[0]*fUpwind_l[14]+alphaDrSurf_l[5]*fUpwind_l[13]+alphaDrSurf_l[6]*fUpwind_l[12]+fUpwind_l[8]*alphaDrSurf_l[11]+alphaDrSurf_l[2]*fUpwind_l[10]+alphaDrSurf_l[3]*fUpwind_l[9]+fUpwind_l[4]*alphaDrSurf_l[7]); 
  Ghat_l[15] = 0.25*(alphaDrSurf_l[0]*fUpwind_l[15]+alphaDrSurf_l[1]*fUpwind_l[14]+alphaDrSurf_l[2]*fUpwind_l[13]+alphaDrSurf_l[3]*fUpwind_l[12]+fUpwind_l[4]*alphaDrSurf_l[11]+alphaDrSurf_l[5]*fUpwind_l[10]+alphaDrSurf_l[6]*fUpwind_l[9]+alphaDrSurf_l[7]*fUpwind_l[8]); 

  Ghat_r[0] = 0.25*(alphaDrSurf_r[11]*fUpwind_r[11]+alphaDrSurf_r[7]*fUpwind_r[7]+alphaDrSurf_r[6]*fUpwind_r[6]+alphaDrSurf_r[5]*fUpwind_r[5]+alphaDrSurf_r[3]*fUpwind_r[3]+alphaDrSurf_r[2]*fUpwind_r[2]+alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.25*(alphaDrSurf_r[7]*fUpwind_r[11]+fUpwind_r[7]*alphaDrSurf_r[11]+alphaDrSurf_r[3]*fUpwind_r[6]+fUpwind_r[3]*alphaDrSurf_r[6]+alphaDrSurf_r[2]*fUpwind_r[5]+fUpwind_r[2]*alphaDrSurf_r[5]+alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.25*(alphaDrSurf_r[6]*fUpwind_r[11]+fUpwind_r[6]*alphaDrSurf_r[11]+alphaDrSurf_r[3]*fUpwind_r[7]+fUpwind_r[3]*alphaDrSurf_r[7]+alphaDrSurf_r[1]*fUpwind_r[5]+fUpwind_r[1]*alphaDrSurf_r[5]+alphaDrSurf_r[0]*fUpwind_r[2]+fUpwind_r[0]*alphaDrSurf_r[2]); 
  Ghat_r[3] = 0.25*(alphaDrSurf_r[5]*fUpwind_r[11]+fUpwind_r[5]*alphaDrSurf_r[11]+alphaDrSurf_r[2]*fUpwind_r[7]+fUpwind_r[2]*alphaDrSurf_r[7]+alphaDrSurf_r[1]*fUpwind_r[6]+fUpwind_r[1]*alphaDrSurf_r[6]+alphaDrSurf_r[0]*fUpwind_r[3]+fUpwind_r[0]*alphaDrSurf_r[3]); 
  Ghat_r[4] = 0.25*(alphaDrSurf_r[11]*fUpwind_r[15]+alphaDrSurf_r[7]*fUpwind_r[14]+alphaDrSurf_r[6]*fUpwind_r[13]+alphaDrSurf_r[5]*fUpwind_r[12]+alphaDrSurf_r[3]*fUpwind_r[10]+alphaDrSurf_r[2]*fUpwind_r[9]+alphaDrSurf_r[1]*fUpwind_r[8]+alphaDrSurf_r[0]*fUpwind_r[4]); 
  Ghat_r[5] = 0.25*(alphaDrSurf_r[3]*fUpwind_r[11]+fUpwind_r[3]*alphaDrSurf_r[11]+alphaDrSurf_r[6]*fUpwind_r[7]+fUpwind_r[6]*alphaDrSurf_r[7]+alphaDrSurf_r[0]*fUpwind_r[5]+fUpwind_r[0]*alphaDrSurf_r[5]+alphaDrSurf_r[1]*fUpwind_r[2]+fUpwind_r[1]*alphaDrSurf_r[2]); 
  Ghat_r[6] = 0.25*(alphaDrSurf_r[2]*fUpwind_r[11]+fUpwind_r[2]*alphaDrSurf_r[11]+alphaDrSurf_r[5]*fUpwind_r[7]+fUpwind_r[5]*alphaDrSurf_r[7]+alphaDrSurf_r[0]*fUpwind_r[6]+fUpwind_r[0]*alphaDrSurf_r[6]+alphaDrSurf_r[1]*fUpwind_r[3]+fUpwind_r[1]*alphaDrSurf_r[3]); 
  Ghat_r[7] = 0.25*(alphaDrSurf_r[1]*fUpwind_r[11]+fUpwind_r[1]*alphaDrSurf_r[11]+alphaDrSurf_r[0]*fUpwind_r[7]+fUpwind_r[0]*alphaDrSurf_r[7]+alphaDrSurf_r[5]*fUpwind_r[6]+fUpwind_r[5]*alphaDrSurf_r[6]+alphaDrSurf_r[2]*fUpwind_r[3]+fUpwind_r[2]*alphaDrSurf_r[3]); 
  Ghat_r[8] = 0.25*(alphaDrSurf_r[7]*fUpwind_r[15]+alphaDrSurf_r[11]*fUpwind_r[14]+alphaDrSurf_r[3]*fUpwind_r[13]+alphaDrSurf_r[2]*fUpwind_r[12]+alphaDrSurf_r[6]*fUpwind_r[10]+alphaDrSurf_r[5]*fUpwind_r[9]+alphaDrSurf_r[0]*fUpwind_r[8]+alphaDrSurf_r[1]*fUpwind_r[4]); 
  Ghat_r[9] = 0.25*(alphaDrSurf_r[6]*fUpwind_r[15]+alphaDrSurf_r[3]*fUpwind_r[14]+alphaDrSurf_r[11]*fUpwind_r[13]+alphaDrSurf_r[1]*fUpwind_r[12]+alphaDrSurf_r[7]*fUpwind_r[10]+alphaDrSurf_r[0]*fUpwind_r[9]+alphaDrSurf_r[5]*fUpwind_r[8]+alphaDrSurf_r[2]*fUpwind_r[4]); 
  Ghat_r[10] = 0.25*(alphaDrSurf_r[5]*fUpwind_r[15]+alphaDrSurf_r[2]*fUpwind_r[14]+alphaDrSurf_r[1]*fUpwind_r[13]+alphaDrSurf_r[11]*fUpwind_r[12]+alphaDrSurf_r[0]*fUpwind_r[10]+alphaDrSurf_r[7]*fUpwind_r[9]+alphaDrSurf_r[6]*fUpwind_r[8]+alphaDrSurf_r[3]*fUpwind_r[4]); 
  Ghat_r[11] = 0.25*(alphaDrSurf_r[0]*fUpwind_r[11]+fUpwind_r[0]*alphaDrSurf_r[11]+alphaDrSurf_r[1]*fUpwind_r[7]+fUpwind_r[1]*alphaDrSurf_r[7]+alphaDrSurf_r[2]*fUpwind_r[6]+fUpwind_r[2]*alphaDrSurf_r[6]+alphaDrSurf_r[3]*fUpwind_r[5]+fUpwind_r[3]*alphaDrSurf_r[5]); 
  Ghat_r[12] = 0.25*(alphaDrSurf_r[3]*fUpwind_r[15]+alphaDrSurf_r[6]*fUpwind_r[14]+alphaDrSurf_r[7]*fUpwind_r[13]+alphaDrSurf_r[0]*fUpwind_r[12]+fUpwind_r[10]*alphaDrSurf_r[11]+alphaDrSurf_r[1]*fUpwind_r[9]+alphaDrSurf_r[2]*fUpwind_r[8]+fUpwind_r[4]*alphaDrSurf_r[5]); 
  Ghat_r[13] = 0.25*(alphaDrSurf_r[2]*fUpwind_r[15]+alphaDrSurf_r[5]*fUpwind_r[14]+alphaDrSurf_r[0]*fUpwind_r[13]+alphaDrSurf_r[7]*fUpwind_r[12]+fUpwind_r[9]*alphaDrSurf_r[11]+alphaDrSurf_r[1]*fUpwind_r[10]+alphaDrSurf_r[3]*fUpwind_r[8]+fUpwind_r[4]*alphaDrSurf_r[6]); 
  Ghat_r[14] = 0.25*(alphaDrSurf_r[1]*fUpwind_r[15]+alphaDrSurf_r[0]*fUpwind_r[14]+alphaDrSurf_r[5]*fUpwind_r[13]+alphaDrSurf_r[6]*fUpwind_r[12]+fUpwind_r[8]*alphaDrSurf_r[11]+alphaDrSurf_r[2]*fUpwind_r[10]+alphaDrSurf_r[3]*fUpwind_r[9]+fUpwind_r[4]*alphaDrSurf_r[7]); 
  Ghat_r[15] = 0.25*(alphaDrSurf_r[0]*fUpwind_r[15]+alphaDrSurf_r[1]*fUpwind_r[14]+alphaDrSurf_r[2]*fUpwind_r[13]+alphaDrSurf_r[3]*fUpwind_r[12]+fUpwind_r[4]*alphaDrSurf_r[11]+alphaDrSurf_r[5]*fUpwind_r[10]+alphaDrSurf_r[6]*fUpwind_r[9]+alphaDrSurf_r[7]*fUpwind_r[8]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[9] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[10] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[11] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[14] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[15] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[16] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[17] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[18] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[19] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[20] += 0.7071067811865475*Ghat_r[12]*rdv2-0.7071067811865475*Ghat_l[12]*rdv2; 
  out[21] += 0.7071067811865475*Ghat_r[13]*rdv2-0.7071067811865475*Ghat_l[13]*rdv2; 
  out[22] += 0.7071067811865475*Ghat_r[14]*rdv2-0.7071067811865475*Ghat_l[14]*rdv2; 
  out[23] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[24] += 1.224744871391589*Ghat_r[9]*rdv2+1.224744871391589*Ghat_l[9]*rdv2; 
  out[25] += 1.224744871391589*Ghat_r[10]*rdv2+1.224744871391589*Ghat_l[10]*rdv2; 
  out[26] += 1.224744871391589*Ghat_r[11]*rdv2+1.224744871391589*Ghat_l[11]*rdv2; 
  out[27] += 0.7071067811865475*Ghat_r[15]*rdv2-0.7071067811865475*Ghat_l[15]*rdv2; 
  out[28] += 1.224744871391589*Ghat_r[12]*rdv2+1.224744871391589*Ghat_l[12]*rdv2; 
  out[29] += 1.224744871391589*Ghat_r[13]*rdv2+1.224744871391589*Ghat_l[13]*rdv2; 
  out[30] += 1.224744871391589*Ghat_r[14]*rdv2+1.224744871391589*Ghat_l[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat_r[15]*rdv2+1.224744871391589*Ghat_l[15]*rdv2; 
  out[32] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[33] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 
  out[34] += 1.5811388300841898*Ghat_r[2]*rdv2-1.5811388300841898*Ghat_l[2]*rdv2; 
  out[35] += 1.5811388300841898*Ghat_r[3]*rdv2-1.5811388300841898*Ghat_l[3]*rdv2; 
  out[36] += 1.5811388300841898*Ghat_r[4]*rdv2-1.5811388300841898*Ghat_l[4]*rdv2; 
  out[37] += 1.5811388300841895*Ghat_r[5]*rdv2-1.5811388300841895*Ghat_l[5]*rdv2; 
  out[38] += 1.5811388300841895*Ghat_r[6]*rdv2-1.5811388300841895*Ghat_l[6]*rdv2; 
  out[39] += 1.5811388300841895*Ghat_r[7]*rdv2-1.5811388300841895*Ghat_l[7]*rdv2; 
  out[40] += 1.5811388300841895*Ghat_r[8]*rdv2-1.5811388300841895*Ghat_l[8]*rdv2; 
  out[41] += 1.5811388300841895*Ghat_r[9]*rdv2-1.5811388300841895*Ghat_l[9]*rdv2; 
  out[42] += 1.5811388300841895*Ghat_r[10]*rdv2-1.5811388300841895*Ghat_l[10]*rdv2; 
  out[43] += 1.5811388300841898*Ghat_r[11]*rdv2-1.5811388300841898*Ghat_l[11]*rdv2; 
  out[44] += 1.5811388300841898*Ghat_r[12]*rdv2-1.5811388300841898*Ghat_l[12]*rdv2; 
  out[45] += 1.5811388300841898*Ghat_r[13]*rdv2-1.5811388300841898*Ghat_l[13]*rdv2; 
  out[46] += 1.5811388300841898*Ghat_r[14]*rdv2-1.5811388300841898*Ghat_l[14]*rdv2; 
  out[47] += 1.5811388300841895*Ghat_r[15]*rdv2-1.5811388300841895*Ghat_l[15]*rdv2; 

  return 0.;

} 
