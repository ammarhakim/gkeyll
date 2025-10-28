#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_1x1v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[2]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2 = 2.0/dxv[1]; 

  double alphaDrSurf_l[2] = {0.0}; 
  alphaDrSurf_l[0] = -(1.224744871391589*nuSum[0]*vmap[1])+0.7071067811865475*nuSum[0]*vmap[0]-1.0*nuUSum[0]; 
  alphaDrSurf_l[1] = -(1.224744871391589*nuSum[1]*vmap[1])-1.0*nuUSum[1]+0.7071067811865475*vmap[0]*nuSum[1]; 

  double alphaDrSurf_r[2] = {0.0}; 
  alphaDrSurf_r[0] = 1.224744871391589*nuSum[0]*vmap[1]+0.7071067811865475*nuSum[0]*vmap[0]-1.0*nuUSum[0]; 
  alphaDrSurf_r[1] = 1.224744871391589*nuSum[1]*vmap[1]-1.0*nuUSum[1]+0.7071067811865475*vmap[0]*nuSum[1]; 

  double fUpwindQuad_l[2] = {0.0};
  double fUpwindQuad_r[2] = {0.0};
  double fUpwind_l[2] = {0.0};
  double fUpwind_r[2] = {0.0};
  double Ghat_l[2] = {0.0}; 
  double Ghat_r[2] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_l(fr)/vmap_prime_r[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_r(fl)/vmap_prime_l[0]; 
  } else { 
    fUpwindQuad_l[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_l(fc)/vmap_prime_c[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_r(fc)/vmap_prime_c[0]; 
  } else { 
    fUpwindQuad_r[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_l(fr)/vmap_prime_r[0]; 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*(alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.7071067811865475*(alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 

  Ghat_r[0] = 0.7071067811865475*(alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.7071067811865475*(alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[4] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[5] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 

  return 0.;

} 
