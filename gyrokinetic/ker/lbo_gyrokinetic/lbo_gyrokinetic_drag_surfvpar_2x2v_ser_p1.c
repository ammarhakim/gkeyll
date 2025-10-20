#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_2x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_2x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[4]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2 = 2.0/dxv[2]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = -(1.7320508075688772*nuSum[0]*vmap[1])+nuSum[0]*vmap[0]-1.4142135623730951*nuUSum[0]; 
  alphaDrSurf_l[1] = -(1.7320508075688772*nuSum[1]*vmap[1])-1.4142135623730951*nuUSum[1]+vmap[0]*nuSum[1]; 
  alphaDrSurf_l[2] = -(1.4142135623730951*nuUSum[2])-1.7320508075688772*vmap[1]*nuSum[2]+vmap[0]*nuSum[2]; 
  alphaDrSurf_l[4] = -(1.4142135623730951*nuUSum[3])-1.7320508075688772*vmap[1]*nuSum[3]+vmap[0]*nuSum[3]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = 1.7320508075688772*nuSum[0]*vmap[1]+nuSum[0]*vmap[0]-1.4142135623730951*nuUSum[0]; 
  alphaDrSurf_r[1] = 1.7320508075688772*nuSum[1]*vmap[1]-1.4142135623730951*nuUSum[1]+vmap[0]*nuSum[1]; 
  alphaDrSurf_r[2] = -(1.4142135623730951*nuUSum[2])+1.7320508075688772*vmap[1]*nuSum[2]+vmap[0]*nuSum[2]; 
  alphaDrSurf_r[4] = -(1.4142135623730951*nuUSum[3])+1.7320508075688772*vmap[1]*nuSum[3]+vmap[0]*nuSum[3]; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_r[8] = {0.0};
  double fUpwind_l[8] = {0.0};
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_l(fr)/vmap_prime_r[0]; 
  } 
  if (-alphaDrSurf_l[4]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_l(fc)/vmap_prime_c[0]; 
  } 
  if (-alphaDrSurf_r[4]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_l(fr)/vmap_prime_r[0]; 
  } 
  if (-alphaDrSurf_l[4]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_l(fc)/vmap_prime_c[0]; 
  } 
  if (-alphaDrSurf_r[4]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_l(fr)/vmap_prime_r[0]; 
  } 
  if (alphaDrSurf_l[4]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_r(fl)/vmap_prime_l[0]; 
    fUpwindQuad_l[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_l(fc)/vmap_prime_c[0]; 
    fUpwindQuad_l[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[4]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_r(fc)/vmap_prime_c[0]; 
    fUpwindQuad_r[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_l(fr)/vmap_prime_r[0]; 
    fUpwindQuad_r[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_l(fr)/vmap_prime_r[0]; 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_2x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  gkhyb_2x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*(alphaDrSurf_l[4]*fUpwind_l[4]+alphaDrSurf_l[2]*fUpwind_l[2]+alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.3535533905932737*(alphaDrSurf_l[2]*fUpwind_l[4]+fUpwind_l[2]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[4]+fUpwind_l[1]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fUpwind_l[2]+fUpwind_l[0]*alphaDrSurf_l[2]); 
  Ghat_l[3] = 0.3535533905932737*(alphaDrSurf_l[4]*fUpwind_l[7]+alphaDrSurf_l[2]*fUpwind_l[6]+alphaDrSurf_l[1]*fUpwind_l[5]+alphaDrSurf_l[0]*fUpwind_l[3]); 
  Ghat_l[4] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[4]+fUpwind_l[0]*alphaDrSurf_l[4]+alphaDrSurf_l[1]*fUpwind_l[2]+fUpwind_l[1]*alphaDrSurf_l[2]); 
  Ghat_l[5] = 0.3535533905932737*(alphaDrSurf_l[2]*fUpwind_l[7]+alphaDrSurf_l[4]*fUpwind_l[6]+alphaDrSurf_l[0]*fUpwind_l[5]+alphaDrSurf_l[1]*fUpwind_l[3]); 
  Ghat_l[6] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[7]+alphaDrSurf_l[0]*fUpwind_l[6]+alphaDrSurf_l[4]*fUpwind_l[5]+alphaDrSurf_l[2]*fUpwind_l[3]); 
  Ghat_l[7] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[7]+alphaDrSurf_l[1]*fUpwind_l[6]+alphaDrSurf_l[2]*fUpwind_l[5]+fUpwind_l[3]*alphaDrSurf_l[4]); 

  Ghat_r[0] = 0.3535533905932737*(alphaDrSurf_r[4]*fUpwind_r[4]+alphaDrSurf_r[2]*fUpwind_r[2]+alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.3535533905932737*(alphaDrSurf_r[2]*fUpwind_r[4]+fUpwind_r[2]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[4]+fUpwind_r[1]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fUpwind_r[2]+fUpwind_r[0]*alphaDrSurf_r[2]); 
  Ghat_r[3] = 0.3535533905932737*(alphaDrSurf_r[4]*fUpwind_r[7]+alphaDrSurf_r[2]*fUpwind_r[6]+alphaDrSurf_r[1]*fUpwind_r[5]+alphaDrSurf_r[0]*fUpwind_r[3]); 
  Ghat_r[4] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[4]+fUpwind_r[0]*alphaDrSurf_r[4]+alphaDrSurf_r[1]*fUpwind_r[2]+fUpwind_r[1]*alphaDrSurf_r[2]); 
  Ghat_r[5] = 0.3535533905932737*(alphaDrSurf_r[2]*fUpwind_r[7]+alphaDrSurf_r[4]*fUpwind_r[6]+alphaDrSurf_r[0]*fUpwind_r[5]+alphaDrSurf_r[1]*fUpwind_r[3]); 
  Ghat_r[6] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[7]+alphaDrSurf_r[0]*fUpwind_r[6]+alphaDrSurf_r[4]*fUpwind_r[5]+alphaDrSurf_r[2]*fUpwind_r[3]); 
  Ghat_r[7] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[7]+alphaDrSurf_r[1]*fUpwind_r[6]+alphaDrSurf_r[2]*fUpwind_r[5]+fUpwind_r[3]*alphaDrSurf_r[4]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[7] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[9] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[10] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[16] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[17] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 
  out[18] += 1.5811388300841898*Ghat_r[2]*rdv2-1.5811388300841898*Ghat_l[2]*rdv2; 
  out[19] += 1.5811388300841898*Ghat_r[3]*rdv2-1.5811388300841898*Ghat_l[3]*rdv2; 
  out[20] += 1.5811388300841895*Ghat_r[4]*rdv2-1.5811388300841895*Ghat_l[4]*rdv2; 
  out[21] += 1.5811388300841895*Ghat_r[5]*rdv2-1.5811388300841895*Ghat_l[5]*rdv2; 
  out[22] += 1.5811388300841895*Ghat_r[6]*rdv2-1.5811388300841895*Ghat_l[6]*rdv2; 
  out[23] += 1.5811388300841898*Ghat_r[7]*rdv2-1.5811388300841898*Ghat_l[7]*rdv2; 

  return 0.;

} 
