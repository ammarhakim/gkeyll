#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_edge, const double *vmap_prime_skin, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // dxv[3]: cell spacing. 
  // vmap: velocity space mapping in skin cell.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fEdge,fSkin: Distribution function in edge and skin cells 
  // out: Incremented distribution function in skin cell 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2 = 2.0/dxv[1]; 


  double alphaDrSurf[4] = {0.0}; 
  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};;
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(1.7320508075688772*vmap[1]+vmap[0])-1.4142135623730951*nuUSum[0]; 
  alphaDrSurf[1] = 1.7320508075688772*nuSum[1]*vmap[1]-1.4142135623730951*nuUSum[1]+vmap[0]*nuSum[1]; 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0.0) { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_r(fSkin)/vmap_prime_skin[0]; 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_r(fSkin)/vmap_prime_skin[0]; 
  } else { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l(fEdge)/vmap_prime_edge[0]; 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l(fEdge)/vmap_prime_edge[0]; 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0.0) { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_r(fSkin)/vmap_prime_skin[0]; 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_r(fSkin)/vmap_prime_skin[0]; 
  } else { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l(fEdge)/vmap_prime_edge[0]; 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l(fEdge)/vmap_prime_edge[0]; 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[4] += 1.224744871391589*Ghat[1]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[9] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[10] += 1.5811388300841898*Ghat[2]*rdv2; 
  out[11] += 1.5811388300841895*Ghat[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(vmap[0]-1.7320508075688772*vmap[1])-1.4142135623730951*nuUSum[0]; 
  alphaDrSurf[1] = -(1.7320508075688772*nuSum[1]*vmap[1])-1.4142135623730951*nuUSum[1]+vmap[0]*nuSum[1]; 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0.0) { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_r(fEdge)/vmap_prime_edge[0]; 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_r(fEdge)/vmap_prime_edge[0]; 
  } else { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l(fSkin)/vmap_prime_skin[0]; 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l(fSkin)/vmap_prime_skin[0]; 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0.0) { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_r(fEdge)/vmap_prime_edge[0]; 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_r(fEdge)/vmap_prime_edge[0]; 
  } else { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l(fSkin)/vmap_prime_skin[0]; 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l(fSkin)/vmap_prime_skin[0]; 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[4] += 1.224744871391589*Ghat[1]*rdv2; 
  out[5] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[9] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[10] += -(1.5811388300841898*Ghat[2]*rdv2); 
  out[11] += -(1.5811388300841895*Ghat[3]*rdv2); 

  } 

  return 0.;

} 
