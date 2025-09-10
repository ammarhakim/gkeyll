#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_edge, const double *vmap_prime_skin, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // dxv[4]: cell spacing. 
  // vmap: velocity space mapping in skin cell.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fEdge,fSkin: Distribution function in edge and skin cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  double alphaDrSurf[12] = {0.0}; 
  double Ghat[12] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[1] = nuSum[1]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[2] = nuSum[2]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[4] = nuSum[3]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 

  Ghat[0] = -((0.25*(1.7320508075688772*(alphaDrSurf[4]*fEdge[12]+alphaDrSurf[2]*fEdge[9]+alphaDrSurf[1]*fEdge[8])-1.0*alphaDrSurf[4]*fEdge[5]+1.7320508075688772*alphaDrSurf[0]*fEdge[4]-1.0*(alphaDrSurf[2]*fEdge[2]+alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])))/vmap_prime_edge[1]); 
  Ghat[1] = -((0.25*(1.7320508075688772*(alphaDrSurf[2]*fEdge[12]+alphaDrSurf[4]*fEdge[9]+alphaDrSurf[0]*fEdge[8])-1.0*alphaDrSurf[2]*fEdge[5]+1.7320508075688772*alphaDrSurf[1]*fEdge[4]-1.0*(fEdge[2]*alphaDrSurf[4]+alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])))/vmap_prime_edge[1]); 
  Ghat[2] = -((0.25*(1.7320508075688772*(alphaDrSurf[1]*fEdge[12]+alphaDrSurf[0]*fEdge[9]+alphaDrSurf[4]*fEdge[8])-1.0*alphaDrSurf[1]*fEdge[5]+1.7320508075688772*alphaDrSurf[2]*fEdge[4]-1.0*(fEdge[1]*alphaDrSurf[4]+alphaDrSurf[0]*fEdge[2]+fEdge[0]*alphaDrSurf[2])))/vmap_prime_edge[1]); 
  Ghat[3] = -((0.25*(1.7320508075688772*(alphaDrSurf[4]*fEdge[15]+alphaDrSurf[2]*fEdge[14]+alphaDrSurf[1]*fEdge[13])-1.0*alphaDrSurf[4]*fEdge[11]+1.7320508075688772*alphaDrSurf[0]*fEdge[10]-1.0*(alphaDrSurf[2]*fEdge[7]+alphaDrSurf[1]*fEdge[6]+alphaDrSurf[0]*fEdge[3])))/vmap_prime_edge[1]); 
  Ghat[4] = -((0.25*(1.7320508075688772*(alphaDrSurf[0]*fEdge[12]+alphaDrSurf[1]*fEdge[9]+alphaDrSurf[2]*fEdge[8])-1.0*alphaDrSurf[0]*fEdge[5]+1.7320508075688772*alphaDrSurf[4]*fEdge[4]-1.0*(fEdge[0]*alphaDrSurf[4]+alphaDrSurf[1]*fEdge[2]+fEdge[1]*alphaDrSurf[2])))/vmap_prime_edge[1]); 
  Ghat[5] = -((0.25*(1.7320508075688772*(alphaDrSurf[2]*fEdge[15]+alphaDrSurf[4]*fEdge[14]+alphaDrSurf[0]*fEdge[13])-1.0*alphaDrSurf[2]*fEdge[11]+1.7320508075688772*alphaDrSurf[1]*fEdge[10]-1.0*(alphaDrSurf[4]*fEdge[7]+alphaDrSurf[0]*fEdge[6]+alphaDrSurf[1]*fEdge[3])))/vmap_prime_edge[1]); 
  Ghat[6] = -((0.25*(1.7320508075688772*(alphaDrSurf[1]*fEdge[15]+alphaDrSurf[0]*fEdge[14]+alphaDrSurf[4]*fEdge[13])-1.0*alphaDrSurf[1]*fEdge[11]+1.7320508075688772*alphaDrSurf[2]*fEdge[10]-1.0*(alphaDrSurf[0]*fEdge[7]+alphaDrSurf[4]*fEdge[6]+alphaDrSurf[2]*fEdge[3])))/vmap_prime_edge[1]); 
  Ghat[7] = -((0.25*(1.7320508075688772*(alphaDrSurf[0]*fEdge[15]+alphaDrSurf[1]*fEdge[14]+alphaDrSurf[2]*fEdge[13])-1.0*alphaDrSurf[0]*fEdge[11]+1.7320508075688772*alphaDrSurf[4]*fEdge[10]-1.0*(alphaDrSurf[1]*fEdge[7]+alphaDrSurf[2]*fEdge[6]+fEdge[3]*alphaDrSurf[4])))/vmap_prime_edge[1]); 
  Ghat[8] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf[4]*fEdge[23]+25.980762113533157*(alphaDrSurf[2]*fEdge[22]+alphaDrSurf[1]*fEdge[21])-15.0*alphaDrSurf[4]*fEdge[20]+25.98076211353316*alphaDrSurf[0]*fEdge[19]-15.000000000000002*(alphaDrSurf[2]*fEdge[18]+alphaDrSurf[1]*fEdge[17])-15.0*alphaDrSurf[0]*fEdge[16]))/vmap_prime_edge[1]); 
  Ghat[9] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf[2]*fEdge[23]+25.98076211353316*(alphaDrSurf[4]*fEdge[22]+alphaDrSurf[0]*fEdge[21])-15.000000000000002*alphaDrSurf[2]*fEdge[20]+25.980762113533157*alphaDrSurf[1]*fEdge[19]-15.0*(alphaDrSurf[4]*fEdge[18]+alphaDrSurf[0]*fEdge[17])-15.000000000000002*alphaDrSurf[1]*fEdge[16]))/vmap_prime_edge[1]); 
  Ghat[10] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf[1]*fEdge[23]+25.98076211353316*(alphaDrSurf[0]*fEdge[22]+alphaDrSurf[4]*fEdge[21])-15.000000000000002*alphaDrSurf[1]*fEdge[20]+25.980762113533157*alphaDrSurf[2]*fEdge[19]-15.0*(alphaDrSurf[0]*fEdge[18]+alphaDrSurf[4]*fEdge[17])-15.000000000000002*alphaDrSurf[2]*fEdge[16]))/vmap_prime_edge[1]); 
  Ghat[11] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf[0]*fEdge[23]+25.980762113533157*(alphaDrSurf[1]*fEdge[22]+alphaDrSurf[2]*fEdge[21])-15.0*alphaDrSurf[0]*fEdge[20]+25.98076211353316*alphaDrSurf[4]*fEdge[19]-15.000000000000002*(alphaDrSurf[1]*fEdge[18]+alphaDrSurf[2]*fEdge[17])-15.0*alphaDrSurf[4]*fEdge[16]))/vmap_prime_edge[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 1.224744871391589*Ghat[4]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 
  out[16] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[19] += 1.224744871391589*Ghat[8]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[21] += 1.224744871391589*Ghat[9]*rdv2; 
  out[22] += 1.224744871391589*Ghat[10]*rdv2; 
  out[23] += 1.224744871391589*Ghat[11]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[1] = nuSum[1]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[2] = nuSum[2]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[4] = nuSum[3]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 

  Ghat[0] = -((0.25*(1.7320508075688772*(alphaDrSurf[4]*fSkin[12]+alphaDrSurf[2]*fSkin[9]+alphaDrSurf[1]*fSkin[8])-1.0*alphaDrSurf[4]*fSkin[5]+1.7320508075688772*alphaDrSurf[0]*fSkin[4]-1.0*(alphaDrSurf[2]*fSkin[2]+alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])))/vmap_prime_skin[1]); 
  Ghat[1] = -((0.25*(1.7320508075688772*(alphaDrSurf[2]*fSkin[12]+alphaDrSurf[4]*fSkin[9]+alphaDrSurf[0]*fSkin[8])-1.0*alphaDrSurf[2]*fSkin[5]+1.7320508075688772*alphaDrSurf[1]*fSkin[4]-1.0*(fSkin[2]*alphaDrSurf[4]+alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])))/vmap_prime_skin[1]); 
  Ghat[2] = -((0.25*(1.7320508075688772*(alphaDrSurf[1]*fSkin[12]+alphaDrSurf[0]*fSkin[9]+alphaDrSurf[4]*fSkin[8])-1.0*alphaDrSurf[1]*fSkin[5]+1.7320508075688772*alphaDrSurf[2]*fSkin[4]-1.0*(fSkin[1]*alphaDrSurf[4]+alphaDrSurf[0]*fSkin[2]+fSkin[0]*alphaDrSurf[2])))/vmap_prime_skin[1]); 
  Ghat[3] = -((0.25*(1.7320508075688772*(alphaDrSurf[4]*fSkin[15]+alphaDrSurf[2]*fSkin[14]+alphaDrSurf[1]*fSkin[13])-1.0*alphaDrSurf[4]*fSkin[11]+1.7320508075688772*alphaDrSurf[0]*fSkin[10]-1.0*(alphaDrSurf[2]*fSkin[7]+alphaDrSurf[1]*fSkin[6]+alphaDrSurf[0]*fSkin[3])))/vmap_prime_skin[1]); 
  Ghat[4] = -((0.25*(1.7320508075688772*(alphaDrSurf[0]*fSkin[12]+alphaDrSurf[1]*fSkin[9]+alphaDrSurf[2]*fSkin[8])-1.0*alphaDrSurf[0]*fSkin[5]+1.7320508075688772*alphaDrSurf[4]*fSkin[4]-1.0*(fSkin[0]*alphaDrSurf[4]+alphaDrSurf[1]*fSkin[2]+fSkin[1]*alphaDrSurf[2])))/vmap_prime_skin[1]); 
  Ghat[5] = -((0.25*(1.7320508075688772*(alphaDrSurf[2]*fSkin[15]+alphaDrSurf[4]*fSkin[14]+alphaDrSurf[0]*fSkin[13])-1.0*alphaDrSurf[2]*fSkin[11]+1.7320508075688772*alphaDrSurf[1]*fSkin[10]-1.0*(alphaDrSurf[4]*fSkin[7]+alphaDrSurf[0]*fSkin[6]+alphaDrSurf[1]*fSkin[3])))/vmap_prime_skin[1]); 
  Ghat[6] = -((0.25*(1.7320508075688772*(alphaDrSurf[1]*fSkin[15]+alphaDrSurf[0]*fSkin[14]+alphaDrSurf[4]*fSkin[13])-1.0*alphaDrSurf[1]*fSkin[11]+1.7320508075688772*alphaDrSurf[2]*fSkin[10]-1.0*(alphaDrSurf[0]*fSkin[7]+alphaDrSurf[4]*fSkin[6]+alphaDrSurf[2]*fSkin[3])))/vmap_prime_skin[1]); 
  Ghat[7] = -((0.25*(1.7320508075688772*(alphaDrSurf[0]*fSkin[15]+alphaDrSurf[1]*fSkin[14]+alphaDrSurf[2]*fSkin[13])-1.0*alphaDrSurf[0]*fSkin[11]+1.7320508075688772*alphaDrSurf[4]*fSkin[10]-1.0*(alphaDrSurf[1]*fSkin[7]+alphaDrSurf[2]*fSkin[6]+fSkin[3]*alphaDrSurf[4])))/vmap_prime_skin[1]); 
  Ghat[8] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf[4]*fSkin[23]+25.980762113533157*(alphaDrSurf[2]*fSkin[22]+alphaDrSurf[1]*fSkin[21])-15.0*alphaDrSurf[4]*fSkin[20]+25.98076211353316*alphaDrSurf[0]*fSkin[19]-15.000000000000002*(alphaDrSurf[2]*fSkin[18]+alphaDrSurf[1]*fSkin[17])-15.0*alphaDrSurf[0]*fSkin[16]))/vmap_prime_skin[1]); 
  Ghat[9] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf[2]*fSkin[23]+25.98076211353316*(alphaDrSurf[4]*fSkin[22]+alphaDrSurf[0]*fSkin[21])-15.000000000000002*alphaDrSurf[2]*fSkin[20]+25.980762113533157*alphaDrSurf[1]*fSkin[19]-15.0*(alphaDrSurf[4]*fSkin[18]+alphaDrSurf[0]*fSkin[17])-15.000000000000002*alphaDrSurf[1]*fSkin[16]))/vmap_prime_skin[1]); 
  Ghat[10] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf[1]*fSkin[23]+25.98076211353316*(alphaDrSurf[0]*fSkin[22]+alphaDrSurf[4]*fSkin[21])-15.000000000000002*alphaDrSurf[1]*fSkin[20]+25.980762113533157*alphaDrSurf[2]*fSkin[19]-15.0*(alphaDrSurf[0]*fSkin[18]+alphaDrSurf[4]*fSkin[17])-15.000000000000002*alphaDrSurf[2]*fSkin[16]))/vmap_prime_skin[1]); 
  Ghat[11] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf[0]*fSkin[23]+25.980762113533157*(alphaDrSurf[1]*fSkin[22]+alphaDrSurf[2]*fSkin[21])-15.0*alphaDrSurf[0]*fSkin[20]+25.98076211353316*alphaDrSurf[4]*fSkin[19]-15.000000000000002*(alphaDrSurf[1]*fSkin[18]+alphaDrSurf[2]*fSkin[17])-15.0*alphaDrSurf[4]*fSkin[16]))/vmap_prime_skin[1]); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[6] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[7] += -(0.7071067811865475*Ghat[6]*rdv2); 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -(0.7071067811865475*Ghat[7]*rdv2); 
  out[12] += 1.224744871391589*Ghat[4]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 
  out[16] += -(0.7071067811865475*Ghat[8]*rdv2); 
  out[17] += -(0.7071067811865475*Ghat[9]*rdv2); 
  out[18] += -(0.7071067811865475*Ghat[10]*rdv2); 
  out[19] += 1.224744871391589*Ghat[8]*rdv2; 
  out[20] += -(0.7071067811865475*Ghat[11]*rdv2); 
  out[21] += 1.224744871391589*Ghat[9]*rdv2; 
  out[22] += 1.224744871391589*Ghat[10]*rdv2; 
  out[23] += 1.224744871391589*Ghat[11]*rdv2; 

  } 

  return 0.;

} 
