#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_edge, const double *vmap_prime_skin, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // dxv[3]: cell spacing. 
  // vmap: velocity space mapping in skin cell.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fEdge,fSkin: Distribution function in edge and skin cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  double alphaDrSurf[6] = {0.0}; 
  double Ghat[6] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[1] = nuSum[1]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 

  Ghat[0] = -((0.25*(2.4494897427831783*(alphaDrSurf[1]*fEdge[5]+alphaDrSurf[0]*fEdge[3])-1.4142135623730951*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])))/vmap_prime_edge[1]); 
  Ghat[1] = -((0.25*(2.4494897427831783*(alphaDrSurf[0]*fEdge[5]+alphaDrSurf[1]*fEdge[3])-1.4142135623730951*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])))/vmap_prime_edge[1]); 
  Ghat[2] = -((0.25*(2.4494897427831783*(alphaDrSurf[1]*fEdge[7]+alphaDrSurf[0]*fEdge[6])-1.4142135623730951*(alphaDrSurf[1]*fEdge[4]+alphaDrSurf[0]*fEdge[2])))/vmap_prime_edge[1]); 
  Ghat[3] = -((0.25*(2.4494897427831783*(alphaDrSurf[0]*fEdge[7]+alphaDrSurf[1]*fEdge[6])-1.4142135623730951*(alphaDrSurf[0]*fEdge[4]+alphaDrSurf[1]*fEdge[2])))/vmap_prime_edge[1]); 
  Ghat[4] = -((0.016666666666666666*(36.74234614174767*alphaDrSurf[1]*fEdge[11]+36.74234614174768*alphaDrSurf[0]*fEdge[10]-21.21320343559643*alphaDrSurf[1]*fEdge[9]-21.213203435596427*alphaDrSurf[0]*fEdge[8]))/vmap_prime_edge[1]); 
  Ghat[5] = -((0.016666666666666666*(36.74234614174768*alphaDrSurf[0]*fEdge[11]+36.74234614174767*alphaDrSurf[1]*fEdge[10]-21.213203435596427*alphaDrSurf[0]*fEdge[9]-21.21320343559643*alphaDrSurf[1]*fEdge[8]))/vmap_prime_edge[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[10] += 1.224744871391589*Ghat[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat[5]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[1] = nuSum[1]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 

  Ghat[0] = -((0.25*(2.4494897427831783*(alphaDrSurf[1]*fSkin[5]+alphaDrSurf[0]*fSkin[3])-1.4142135623730951*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])))/vmap_prime_skin[1]); 
  Ghat[1] = -((0.25*(2.4494897427831783*(alphaDrSurf[0]*fSkin[5]+alphaDrSurf[1]*fSkin[3])-1.4142135623730951*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])))/vmap_prime_skin[1]); 
  Ghat[2] = -((0.25*(2.4494897427831783*(alphaDrSurf[1]*fSkin[7]+alphaDrSurf[0]*fSkin[6])-1.4142135623730951*(alphaDrSurf[1]*fSkin[4]+alphaDrSurf[0]*fSkin[2])))/vmap_prime_skin[1]); 
  Ghat[3] = -((0.25*(2.4494897427831783*(alphaDrSurf[0]*fSkin[7]+alphaDrSurf[1]*fSkin[6])-1.4142135623730951*(alphaDrSurf[0]*fSkin[4]+alphaDrSurf[1]*fSkin[2])))/vmap_prime_skin[1]); 
  Ghat[4] = -((0.016666666666666666*(36.74234614174767*alphaDrSurf[1]*fSkin[11]+36.74234614174768*alphaDrSurf[0]*fSkin[10]-21.21320343559643*alphaDrSurf[1]*fSkin[9]-21.213203435596427*alphaDrSurf[0]*fSkin[8]))/vmap_prime_skin[1]); 
  Ghat[5] = -((0.016666666666666666*(36.74234614174768*alphaDrSurf[0]*fSkin[11]+36.74234614174767*alphaDrSurf[1]*fSkin[10]-21.213203435596427*alphaDrSurf[0]*fSkin[9]-21.21320343559643*alphaDrSurf[1]*fSkin[8]))/vmap_prime_skin[1]); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[9] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[10] += 1.224744871391589*Ghat[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat[5]*rdv2; 

  } 

  return 0.;

} 
