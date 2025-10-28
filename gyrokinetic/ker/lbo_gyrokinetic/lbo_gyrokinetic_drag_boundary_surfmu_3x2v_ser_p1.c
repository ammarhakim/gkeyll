#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_edge, const double *vmap_prime_skin, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // dxv[5]: cell spacing. 
  // vmap: velocity space mapping in skin cell.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[16]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fEdge,fSkin: Distribution function in edge and skin cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[4]; 

  double alphaDrSurf[24] = {0.0}; 
  double Ghat[24] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[1] = nuSum[1]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[2] = nuSum[2]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[3] = nuSum[3]*(3.4641016151377544*vmap[3]+2.0*vmap[2]); 
  alphaDrSurf[5] = (3.4641016151377544*vmap[3]+2.0*vmap[2])*nuSum[4]; 
  alphaDrSurf[6] = (3.4641016151377544*vmap[3]+2.0*vmap[2])*nuSum[5]; 
  alphaDrSurf[7] = (3.4641016151377544*vmap[3]+2.0*vmap[2])*nuSum[6]; 
  alphaDrSurf[11] = (3.4641016151377544*vmap[3]+2.0*vmap[2])*nuSum[7]; 

  Ghat[0] = -((0.125*(2.4494897427831783*(alphaDrSurf[11]*fEdge[27]+alphaDrSurf[7]*fEdge[22]+alphaDrSurf[6]*fEdge[21]+alphaDrSurf[5]*fEdge[20])-1.4142135623730951*alphaDrSurf[11]*fEdge[16]+2.4494897427831783*(alphaDrSurf[3]*fEdge[14]+alphaDrSurf[2]*fEdge[13]+alphaDrSurf[1]*fEdge[12])-1.4142135623730951*(alphaDrSurf[7]*fEdge[8]+alphaDrSurf[6]*fEdge[7]+alphaDrSurf[5]*fEdge[6])+2.4494897427831783*alphaDrSurf[0]*fEdge[5]-1.4142135623730951*(alphaDrSurf[3]*fEdge[3]+alphaDrSurf[2]*fEdge[2]+alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])))/vmap_prime_edge[1]); 
  Ghat[1] = -((0.125*(2.4494897427831783*(alphaDrSurf[7]*fEdge[27]+alphaDrSurf[11]*fEdge[22]+alphaDrSurf[3]*fEdge[21]+alphaDrSurf[2]*fEdge[20])-1.4142135623730951*alphaDrSurf[7]*fEdge[16]+2.4494897427831783*(alphaDrSurf[6]*fEdge[14]+alphaDrSurf[5]*fEdge[13]+alphaDrSurf[0]*fEdge[12])-1.4142135623730951*(fEdge[8]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[7]+alphaDrSurf[2]*fEdge[6]+fEdge[3]*alphaDrSurf[6])+2.4494897427831783*alphaDrSurf[1]*fEdge[5]-1.4142135623730951*(fEdge[2]*alphaDrSurf[5]+alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])))/vmap_prime_edge[1]); 
  Ghat[2] = -((0.125*(2.4494897427831783*(alphaDrSurf[6]*fEdge[27]+alphaDrSurf[3]*fEdge[22]+alphaDrSurf[11]*fEdge[21]+alphaDrSurf[1]*fEdge[20])-1.4142135623730951*alphaDrSurf[6]*fEdge[16]+2.4494897427831783*(alphaDrSurf[7]*fEdge[14]+alphaDrSurf[0]*fEdge[13]+alphaDrSurf[5]*fEdge[12])-1.4142135623730951*(fEdge[7]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[8]+fEdge[3]*alphaDrSurf[7]+alphaDrSurf[1]*fEdge[6])+2.4494897427831783*alphaDrSurf[2]*fEdge[5]-1.4142135623730951*(fEdge[1]*alphaDrSurf[5]+alphaDrSurf[0]*fEdge[2]+fEdge[0]*alphaDrSurf[2])))/vmap_prime_edge[1]); 
  Ghat[3] = -((0.125*(2.4494897427831783*(alphaDrSurf[5]*fEdge[27]+alphaDrSurf[2]*fEdge[22]+alphaDrSurf[1]*fEdge[21]+alphaDrSurf[11]*fEdge[20])-1.4142135623730951*alphaDrSurf[5]*fEdge[16]+2.4494897427831783*(alphaDrSurf[0]*fEdge[14]+alphaDrSurf[7]*fEdge[13]+alphaDrSurf[6]*fEdge[12])-1.4142135623730951*(fEdge[6]*alphaDrSurf[11]+alphaDrSurf[2]*fEdge[8]+alphaDrSurf[1]*fEdge[7]+fEdge[2]*alphaDrSurf[7]+fEdge[1]*alphaDrSurf[6])+2.4494897427831783*alphaDrSurf[3]*fEdge[5]-1.4142135623730951*(alphaDrSurf[0]*fEdge[3]+fEdge[0]*alphaDrSurf[3])))/vmap_prime_edge[1]); 
  Ghat[4] = -((0.125*(2.4494897427831783*(alphaDrSurf[11]*fEdge[31]+alphaDrSurf[7]*fEdge[30]+alphaDrSurf[6]*fEdge[29]+alphaDrSurf[5]*fEdge[28])-1.4142135623730951*alphaDrSurf[11]*fEdge[26]+2.4494897427831783*(alphaDrSurf[3]*fEdge[25]+alphaDrSurf[2]*fEdge[24]+alphaDrSurf[1]*fEdge[23])-1.4142135623730951*(alphaDrSurf[7]*fEdge[19]+alphaDrSurf[6]*fEdge[18]+alphaDrSurf[5]*fEdge[17])+2.4494897427831783*alphaDrSurf[0]*fEdge[15]-1.4142135623730951*(alphaDrSurf[3]*fEdge[11]+alphaDrSurf[2]*fEdge[10]+alphaDrSurf[1]*fEdge[9]+alphaDrSurf[0]*fEdge[4])))/vmap_prime_edge[1]); 
  Ghat[5] = -((0.125*(2.4494897427831783*(alphaDrSurf[3]*fEdge[27]+alphaDrSurf[6]*fEdge[22]+alphaDrSurf[7]*fEdge[21]+alphaDrSurf[0]*fEdge[20])-1.4142135623730951*alphaDrSurf[3]*fEdge[16]+2.4494897427831783*(alphaDrSurf[11]*fEdge[14]+alphaDrSurf[1]*fEdge[13]+alphaDrSurf[2]*fEdge[12])-1.4142135623730951*(fEdge[3]*alphaDrSurf[11]+alphaDrSurf[6]*fEdge[8]+alphaDrSurf[7]*fEdge[7]+alphaDrSurf[0]*fEdge[6])+2.4494897427831783*alphaDrSurf[5]*fEdge[5]-1.4142135623730951*(fEdge[0]*alphaDrSurf[5]+alphaDrSurf[1]*fEdge[2]+fEdge[1]*alphaDrSurf[2])))/vmap_prime_edge[1]); 
  Ghat[6] = -((0.125*(2.4494897427831783*(alphaDrSurf[2]*fEdge[27]+alphaDrSurf[5]*fEdge[22]+alphaDrSurf[0]*fEdge[21]+alphaDrSurf[7]*fEdge[20])-1.4142135623730951*alphaDrSurf[2]*fEdge[16]+2.4494897427831783*(alphaDrSurf[1]*fEdge[14]+alphaDrSurf[11]*fEdge[13]+alphaDrSurf[3]*fEdge[12])-1.4142135623730951*(fEdge[2]*alphaDrSurf[11]+alphaDrSurf[5]*fEdge[8]+alphaDrSurf[0]*fEdge[7]+fEdge[6]*alphaDrSurf[7])+(2.4494897427831783*fEdge[5]-1.4142135623730951*fEdge[0])*alphaDrSurf[6]-1.4142135623730951*(alphaDrSurf[1]*fEdge[3]+fEdge[1]*alphaDrSurf[3])))/vmap_prime_edge[1]); 
  Ghat[7] = -((0.125*(2.4494897427831783*(alphaDrSurf[1]*fEdge[27]+alphaDrSurf[0]*fEdge[22]+alphaDrSurf[5]*fEdge[21]+alphaDrSurf[6]*fEdge[20])-1.4142135623730951*alphaDrSurf[1]*fEdge[16]+2.4494897427831783*(alphaDrSurf[2]*fEdge[14]+alphaDrSurf[3]*fEdge[13]+alphaDrSurf[11]*fEdge[12])-1.4142135623730951*(fEdge[1]*alphaDrSurf[11]+alphaDrSurf[0]*fEdge[8]+alphaDrSurf[5]*fEdge[7])+(2.4494897427831783*fEdge[5]-1.4142135623730951*fEdge[0])*alphaDrSurf[7]-1.4142135623730951*(alphaDrSurf[6]*fEdge[6]+alphaDrSurf[2]*fEdge[3]+fEdge[2]*alphaDrSurf[3])))/vmap_prime_edge[1]); 
  Ghat[8] = -((0.125*(2.4494897427831783*(alphaDrSurf[7]*fEdge[31]+alphaDrSurf[11]*fEdge[30]+alphaDrSurf[3]*fEdge[29]+alphaDrSurf[2]*fEdge[28])-1.4142135623730951*alphaDrSurf[7]*fEdge[26]+2.4494897427831783*(alphaDrSurf[6]*fEdge[25]+alphaDrSurf[5]*fEdge[24]+alphaDrSurf[0]*fEdge[23])-1.4142135623730951*(alphaDrSurf[11]*fEdge[19]+alphaDrSurf[3]*fEdge[18]+alphaDrSurf[2]*fEdge[17])+2.4494897427831783*alphaDrSurf[1]*fEdge[15]-1.4142135623730951*(alphaDrSurf[6]*fEdge[11]+alphaDrSurf[5]*fEdge[10]+alphaDrSurf[0]*fEdge[9]+alphaDrSurf[1]*fEdge[4])))/vmap_prime_edge[1]); 
  Ghat[9] = -((0.125*(2.4494897427831783*(alphaDrSurf[6]*fEdge[31]+alphaDrSurf[3]*fEdge[30]+alphaDrSurf[11]*fEdge[29]+alphaDrSurf[1]*fEdge[28])-1.4142135623730951*alphaDrSurf[6]*fEdge[26]+2.4494897427831783*(alphaDrSurf[7]*fEdge[25]+alphaDrSurf[0]*fEdge[24]+alphaDrSurf[5]*fEdge[23])-1.4142135623730951*(alphaDrSurf[3]*fEdge[19]+alphaDrSurf[11]*fEdge[18]+alphaDrSurf[1]*fEdge[17])+2.4494897427831783*alphaDrSurf[2]*fEdge[15]-1.4142135623730951*(alphaDrSurf[7]*fEdge[11]+alphaDrSurf[0]*fEdge[10]+alphaDrSurf[5]*fEdge[9]+alphaDrSurf[2]*fEdge[4])))/vmap_prime_edge[1]); 
  Ghat[10] = -((0.125*(2.4494897427831783*(alphaDrSurf[5]*fEdge[31]+alphaDrSurf[2]*fEdge[30]+alphaDrSurf[1]*fEdge[29]+alphaDrSurf[11]*fEdge[28])-1.4142135623730951*alphaDrSurf[5]*fEdge[26]+2.4494897427831783*(alphaDrSurf[0]*fEdge[25]+alphaDrSurf[7]*fEdge[24]+alphaDrSurf[6]*fEdge[23])-1.4142135623730951*(alphaDrSurf[2]*fEdge[19]+alphaDrSurf[1]*fEdge[18]+alphaDrSurf[11]*fEdge[17])+2.4494897427831783*alphaDrSurf[3]*fEdge[15]-1.4142135623730951*(alphaDrSurf[0]*fEdge[11]+alphaDrSurf[7]*fEdge[10]+alphaDrSurf[6]*fEdge[9]+alphaDrSurf[3]*fEdge[4])))/vmap_prime_edge[1]); 
  Ghat[11] = -((0.125*(2.4494897427831783*(alphaDrSurf[0]*fEdge[27]+alphaDrSurf[1]*fEdge[22]+alphaDrSurf[2]*fEdge[21]+alphaDrSurf[3]*fEdge[20])-1.4142135623730951*alphaDrSurf[0]*fEdge[16]+2.4494897427831783*(alphaDrSurf[5]*fEdge[14]+alphaDrSurf[6]*fEdge[13]+alphaDrSurf[7]*fEdge[12])+(2.4494897427831783*fEdge[5]-1.4142135623730951*fEdge[0])*alphaDrSurf[11]-1.4142135623730951*(alphaDrSurf[1]*fEdge[8]+alphaDrSurf[2]*fEdge[7]+fEdge[1]*alphaDrSurf[7]+alphaDrSurf[3]*fEdge[6]+fEdge[2]*alphaDrSurf[6]+fEdge[3]*alphaDrSurf[5])))/vmap_prime_edge[1]); 
  Ghat[12] = -((0.125*(2.4494897427831783*(alphaDrSurf[3]*fEdge[31]+alphaDrSurf[6]*fEdge[30]+alphaDrSurf[7]*fEdge[29]+alphaDrSurf[0]*fEdge[28])-1.4142135623730951*alphaDrSurf[3]*fEdge[26]+2.4494897427831783*(alphaDrSurf[11]*fEdge[25]+alphaDrSurf[1]*fEdge[24]+alphaDrSurf[2]*fEdge[23])-1.4142135623730951*(alphaDrSurf[6]*fEdge[19]+alphaDrSurf[7]*fEdge[18]+alphaDrSurf[0]*fEdge[17])+2.4494897427831783*alphaDrSurf[5]*fEdge[15]-1.4142135623730951*(alphaDrSurf[11]*fEdge[11]+alphaDrSurf[1]*fEdge[10]+alphaDrSurf[2]*fEdge[9]+fEdge[4]*alphaDrSurf[5])))/vmap_prime_edge[1]); 
  Ghat[13] = -((0.125*(2.4494897427831783*(alphaDrSurf[2]*fEdge[31]+alphaDrSurf[5]*fEdge[30]+alphaDrSurf[0]*fEdge[29]+alphaDrSurf[7]*fEdge[28])-1.4142135623730951*alphaDrSurf[2]*fEdge[26]+2.4494897427831783*(alphaDrSurf[1]*fEdge[25]+alphaDrSurf[11]*fEdge[24]+alphaDrSurf[3]*fEdge[23])-1.4142135623730951*(alphaDrSurf[5]*fEdge[19]+alphaDrSurf[0]*fEdge[18]+alphaDrSurf[7]*fEdge[17])+2.4494897427831783*alphaDrSurf[6]*fEdge[15]-1.4142135623730951*(alphaDrSurf[1]*fEdge[11]+fEdge[10]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[9]+fEdge[4]*alphaDrSurf[6])))/vmap_prime_edge[1]); 
  Ghat[14] = -((0.125*(2.4494897427831783*(alphaDrSurf[1]*fEdge[31]+alphaDrSurf[0]*fEdge[30]+alphaDrSurf[5]*fEdge[29]+alphaDrSurf[6]*fEdge[28])-1.4142135623730951*alphaDrSurf[1]*fEdge[26]+2.4494897427831783*(alphaDrSurf[2]*fEdge[25]+alphaDrSurf[3]*fEdge[24]+alphaDrSurf[11]*fEdge[23])-1.4142135623730951*(alphaDrSurf[0]*fEdge[19]+alphaDrSurf[5]*fEdge[18]+alphaDrSurf[6]*fEdge[17])+2.4494897427831783*alphaDrSurf[7]*fEdge[15]-1.4142135623730951*(alphaDrSurf[2]*fEdge[11]+fEdge[9]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[10]+fEdge[4]*alphaDrSurf[7])))/vmap_prime_edge[1]); 
  Ghat[15] = -((0.125*(2.4494897427831783*(alphaDrSurf[0]*fEdge[31]+alphaDrSurf[1]*fEdge[30]+alphaDrSurf[2]*fEdge[29]+alphaDrSurf[3]*fEdge[28])-1.4142135623730951*alphaDrSurf[0]*fEdge[26]+2.4494897427831783*(alphaDrSurf[5]*fEdge[25]+alphaDrSurf[6]*fEdge[24]+alphaDrSurf[7]*fEdge[23])-1.4142135623730951*(alphaDrSurf[1]*fEdge[19]+alphaDrSurf[2]*fEdge[18]+alphaDrSurf[3]*fEdge[17])+2.4494897427831783*alphaDrSurf[11]*fEdge[15]-1.4142135623730951*(alphaDrSurf[5]*fEdge[11]+fEdge[4]*alphaDrSurf[11]+alphaDrSurf[6]*fEdge[10]+alphaDrSurf[7]*fEdge[9])))/vmap_prime_edge[1]); 
  Ghat[16] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[11]*fEdge[47]+36.74234614174768*(alphaDrSurf[7]*fEdge[46]+alphaDrSurf[6]*fEdge[45]+alphaDrSurf[5]*fEdge[44])-21.21320343559643*alphaDrSurf[11]*fEdge[43]+36.74234614174767*(alphaDrSurf[3]*fEdge[42]+alphaDrSurf[2]*fEdge[41]+alphaDrSurf[1]*fEdge[40])-21.213203435596427*(alphaDrSurf[7]*fEdge[39]+alphaDrSurf[6]*fEdge[38]+alphaDrSurf[5]*fEdge[37])+36.74234614174768*alphaDrSurf[0]*fEdge[36]-21.21320343559643*(alphaDrSurf[3]*fEdge[35]+alphaDrSurf[2]*fEdge[34]+alphaDrSurf[1]*fEdge[33])-21.213203435596427*alphaDrSurf[0]*fEdge[32]))/vmap_prime_edge[1]); 
  Ghat[17] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[7]*fEdge[47]+36.74234614174767*(alphaDrSurf[11]*fEdge[46]+alphaDrSurf[3]*fEdge[45]+alphaDrSurf[2]*fEdge[44])-21.213203435596427*alphaDrSurf[7]*fEdge[43]+36.74234614174768*(alphaDrSurf[6]*fEdge[42]+alphaDrSurf[5]*fEdge[41]+alphaDrSurf[0]*fEdge[40])-21.21320343559643*(alphaDrSurf[11]*fEdge[39]+alphaDrSurf[3]*fEdge[38]+alphaDrSurf[2]*fEdge[37])+36.74234614174767*alphaDrSurf[1]*fEdge[36]-21.213203435596427*(alphaDrSurf[6]*fEdge[35]+alphaDrSurf[5]*fEdge[34]+alphaDrSurf[0]*fEdge[33])-21.21320343559643*alphaDrSurf[1]*fEdge[32]))/vmap_prime_edge[1]); 
  Ghat[18] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[6]*fEdge[47]+36.74234614174767*(alphaDrSurf[3]*fEdge[46]+alphaDrSurf[11]*fEdge[45]+alphaDrSurf[1]*fEdge[44])-21.213203435596427*alphaDrSurf[6]*fEdge[43]+36.74234614174768*(alphaDrSurf[7]*fEdge[42]+alphaDrSurf[0]*fEdge[41]+alphaDrSurf[5]*fEdge[40])-21.21320343559643*(alphaDrSurf[3]*fEdge[39]+alphaDrSurf[11]*fEdge[38]+alphaDrSurf[1]*fEdge[37])+36.74234614174767*alphaDrSurf[2]*fEdge[36]-21.213203435596427*(alphaDrSurf[7]*fEdge[35]+alphaDrSurf[0]*fEdge[34]+alphaDrSurf[5]*fEdge[33])-21.21320343559643*alphaDrSurf[2]*fEdge[32]))/vmap_prime_edge[1]); 
  Ghat[19] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[5]*fEdge[47]+36.74234614174767*(alphaDrSurf[2]*fEdge[46]+alphaDrSurf[1]*fEdge[45]+alphaDrSurf[11]*fEdge[44])-21.213203435596427*alphaDrSurf[5]*fEdge[43]+36.74234614174768*(alphaDrSurf[0]*fEdge[42]+alphaDrSurf[7]*fEdge[41]+alphaDrSurf[6]*fEdge[40])-21.21320343559643*(alphaDrSurf[2]*fEdge[39]+alphaDrSurf[1]*fEdge[38]+alphaDrSurf[11]*fEdge[37])+36.74234614174767*alphaDrSurf[3]*fEdge[36]-21.213203435596427*(alphaDrSurf[0]*fEdge[35]+alphaDrSurf[7]*fEdge[34]+alphaDrSurf[6]*fEdge[33])-21.21320343559643*alphaDrSurf[3]*fEdge[32]))/vmap_prime_edge[1]); 
  Ghat[20] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[3]*fEdge[47]+36.74234614174768*(alphaDrSurf[6]*fEdge[46]+alphaDrSurf[7]*fEdge[45]+alphaDrSurf[0]*fEdge[44])-21.21320343559643*alphaDrSurf[3]*fEdge[43]+36.74234614174767*(alphaDrSurf[11]*fEdge[42]+alphaDrSurf[1]*fEdge[41]+alphaDrSurf[2]*fEdge[40])-21.213203435596427*(alphaDrSurf[6]*fEdge[39]+alphaDrSurf[7]*fEdge[38]+alphaDrSurf[0]*fEdge[37])+36.74234614174768*alphaDrSurf[5]*fEdge[36]-21.21320343559643*(alphaDrSurf[11]*fEdge[35]+alphaDrSurf[1]*fEdge[34]+alphaDrSurf[2]*fEdge[33])-21.213203435596427*alphaDrSurf[5]*fEdge[32]))/vmap_prime_edge[1]); 
  Ghat[21] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[2]*fEdge[47]+36.74234614174768*(alphaDrSurf[5]*fEdge[46]+alphaDrSurf[0]*fEdge[45]+alphaDrSurf[7]*fEdge[44])-21.21320343559643*alphaDrSurf[2]*fEdge[43]+36.74234614174767*(alphaDrSurf[1]*fEdge[42]+alphaDrSurf[11]*fEdge[41]+alphaDrSurf[3]*fEdge[40])-21.213203435596427*(alphaDrSurf[5]*fEdge[39]+alphaDrSurf[0]*fEdge[38]+alphaDrSurf[7]*fEdge[37])+36.74234614174768*alphaDrSurf[6]*fEdge[36]-21.21320343559643*(alphaDrSurf[1]*fEdge[35]+alphaDrSurf[11]*fEdge[34]+alphaDrSurf[3]*fEdge[33])-21.213203435596427*alphaDrSurf[6]*fEdge[32]))/vmap_prime_edge[1]); 
  Ghat[22] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[1]*fEdge[47]+36.74234614174768*(alphaDrSurf[0]*fEdge[46]+alphaDrSurf[5]*fEdge[45]+alphaDrSurf[6]*fEdge[44])-21.21320343559643*alphaDrSurf[1]*fEdge[43]+36.74234614174767*(alphaDrSurf[2]*fEdge[42]+alphaDrSurf[3]*fEdge[41]+alphaDrSurf[11]*fEdge[40])-21.213203435596427*(alphaDrSurf[0]*fEdge[39]+alphaDrSurf[5]*fEdge[38]+alphaDrSurf[6]*fEdge[37])+36.74234614174768*alphaDrSurf[7]*fEdge[36]-21.21320343559643*(alphaDrSurf[2]*fEdge[35]+alphaDrSurf[3]*fEdge[34]+alphaDrSurf[11]*fEdge[33])-21.213203435596427*alphaDrSurf[7]*fEdge[32]))/vmap_prime_edge[1]); 
  Ghat[23] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[0]*fEdge[47]+36.74234614174767*(alphaDrSurf[1]*fEdge[46]+alphaDrSurf[2]*fEdge[45]+alphaDrSurf[3]*fEdge[44])-21.213203435596427*alphaDrSurf[0]*fEdge[43]+36.74234614174768*(alphaDrSurf[5]*fEdge[42]+alphaDrSurf[6]*fEdge[41]+alphaDrSurf[7]*fEdge[40])-21.21320343559643*(alphaDrSurf[1]*fEdge[39]+alphaDrSurf[2]*fEdge[38]+alphaDrSurf[3]*fEdge[37])+36.74234614174767*alphaDrSurf[11]*fEdge[36]-21.213203435596427*(alphaDrSurf[5]*fEdge[35]+alphaDrSurf[6]*fEdge[34]+alphaDrSurf[7]*fEdge[33])-21.21320343559643*alphaDrSurf[11]*fEdge[32]))/vmap_prime_edge[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[5] += 1.224744871391589*Ghat[0]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[10] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[12] += 1.224744871391589*Ghat[1]*rdv2; 
  out[13] += 1.224744871391589*Ghat[2]*rdv2; 
  out[14] += 1.224744871391589*Ghat[3]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[19] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[20] += 1.224744871391589*Ghat[5]*rdv2; 
  out[21] += 1.224744871391589*Ghat[6]*rdv2; 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.224744871391589*Ghat[9]*rdv2; 
  out[25] += 1.224744871391589*Ghat[10]*rdv2; 
  out[26] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[27] += 1.224744871391589*Ghat[11]*rdv2; 
  out[28] += 1.224744871391589*Ghat[12]*rdv2; 
  out[29] += 1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat[15]*rdv2; 
  out[32] += 0.7071067811865475*Ghat[16]*rdv2; 
  out[33] += 0.7071067811865475*Ghat[17]*rdv2; 
  out[34] += 0.7071067811865475*Ghat[18]*rdv2; 
  out[35] += 0.7071067811865475*Ghat[19]*rdv2; 
  out[36] += 1.224744871391589*Ghat[16]*rdv2; 
  out[37] += 0.7071067811865475*Ghat[20]*rdv2; 
  out[38] += 0.7071067811865475*Ghat[21]*rdv2; 
  out[39] += 0.7071067811865475*Ghat[22]*rdv2; 
  out[40] += 1.224744871391589*Ghat[17]*rdv2; 
  out[41] += 1.224744871391589*Ghat[18]*rdv2; 
  out[42] += 1.224744871391589*Ghat[19]*rdv2; 
  out[43] += 0.7071067811865475*Ghat[23]*rdv2; 
  out[44] += 1.224744871391589*Ghat[20]*rdv2; 
  out[45] += 1.224744871391589*Ghat[21]*rdv2; 
  out[46] += 1.224744871391589*Ghat[22]*rdv2; 
  out[47] += 1.224744871391589*Ghat[23]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[1] = nuSum[1]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[2] = nuSum[2]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[3] = nuSum[3]*(2.0*vmap[2]-3.4641016151377544*vmap[3]); 
  alphaDrSurf[5] = (2.0*vmap[2]-3.4641016151377544*vmap[3])*nuSum[4]; 
  alphaDrSurf[6] = (2.0*vmap[2]-3.4641016151377544*vmap[3])*nuSum[5]; 
  alphaDrSurf[7] = (2.0*vmap[2]-3.4641016151377544*vmap[3])*nuSum[6]; 
  alphaDrSurf[11] = (2.0*vmap[2]-3.4641016151377544*vmap[3])*nuSum[7]; 

  Ghat[0] = -((0.125*(2.4494897427831783*(alphaDrSurf[11]*fSkin[27]+alphaDrSurf[7]*fSkin[22]+alphaDrSurf[6]*fSkin[21]+alphaDrSurf[5]*fSkin[20])-1.4142135623730951*alphaDrSurf[11]*fSkin[16]+2.4494897427831783*(alphaDrSurf[3]*fSkin[14]+alphaDrSurf[2]*fSkin[13]+alphaDrSurf[1]*fSkin[12])-1.4142135623730951*(alphaDrSurf[7]*fSkin[8]+alphaDrSurf[6]*fSkin[7]+alphaDrSurf[5]*fSkin[6])+2.4494897427831783*alphaDrSurf[0]*fSkin[5]-1.4142135623730951*(alphaDrSurf[3]*fSkin[3]+alphaDrSurf[2]*fSkin[2]+alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])))/vmap_prime_skin[1]); 
  Ghat[1] = -((0.125*(2.4494897427831783*(alphaDrSurf[7]*fSkin[27]+alphaDrSurf[11]*fSkin[22]+alphaDrSurf[3]*fSkin[21]+alphaDrSurf[2]*fSkin[20])-1.4142135623730951*alphaDrSurf[7]*fSkin[16]+2.4494897427831783*(alphaDrSurf[6]*fSkin[14]+alphaDrSurf[5]*fSkin[13]+alphaDrSurf[0]*fSkin[12])-1.4142135623730951*(fSkin[8]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[7]+alphaDrSurf[2]*fSkin[6]+fSkin[3]*alphaDrSurf[6])+2.4494897427831783*alphaDrSurf[1]*fSkin[5]-1.4142135623730951*(fSkin[2]*alphaDrSurf[5]+alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])))/vmap_prime_skin[1]); 
  Ghat[2] = -((0.125*(2.4494897427831783*(alphaDrSurf[6]*fSkin[27]+alphaDrSurf[3]*fSkin[22]+alphaDrSurf[11]*fSkin[21]+alphaDrSurf[1]*fSkin[20])-1.4142135623730951*alphaDrSurf[6]*fSkin[16]+2.4494897427831783*(alphaDrSurf[7]*fSkin[14]+alphaDrSurf[0]*fSkin[13]+alphaDrSurf[5]*fSkin[12])-1.4142135623730951*(fSkin[7]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[8]+fSkin[3]*alphaDrSurf[7]+alphaDrSurf[1]*fSkin[6])+2.4494897427831783*alphaDrSurf[2]*fSkin[5]-1.4142135623730951*(fSkin[1]*alphaDrSurf[5]+alphaDrSurf[0]*fSkin[2]+fSkin[0]*alphaDrSurf[2])))/vmap_prime_skin[1]); 
  Ghat[3] = -((0.125*(2.4494897427831783*(alphaDrSurf[5]*fSkin[27]+alphaDrSurf[2]*fSkin[22]+alphaDrSurf[1]*fSkin[21]+alphaDrSurf[11]*fSkin[20])-1.4142135623730951*alphaDrSurf[5]*fSkin[16]+2.4494897427831783*(alphaDrSurf[0]*fSkin[14]+alphaDrSurf[7]*fSkin[13]+alphaDrSurf[6]*fSkin[12])-1.4142135623730951*(fSkin[6]*alphaDrSurf[11]+alphaDrSurf[2]*fSkin[8]+alphaDrSurf[1]*fSkin[7]+fSkin[2]*alphaDrSurf[7]+fSkin[1]*alphaDrSurf[6])+2.4494897427831783*alphaDrSurf[3]*fSkin[5]-1.4142135623730951*(alphaDrSurf[0]*fSkin[3]+fSkin[0]*alphaDrSurf[3])))/vmap_prime_skin[1]); 
  Ghat[4] = -((0.125*(2.4494897427831783*(alphaDrSurf[11]*fSkin[31]+alphaDrSurf[7]*fSkin[30]+alphaDrSurf[6]*fSkin[29]+alphaDrSurf[5]*fSkin[28])-1.4142135623730951*alphaDrSurf[11]*fSkin[26]+2.4494897427831783*(alphaDrSurf[3]*fSkin[25]+alphaDrSurf[2]*fSkin[24]+alphaDrSurf[1]*fSkin[23])-1.4142135623730951*(alphaDrSurf[7]*fSkin[19]+alphaDrSurf[6]*fSkin[18]+alphaDrSurf[5]*fSkin[17])+2.4494897427831783*alphaDrSurf[0]*fSkin[15]-1.4142135623730951*(alphaDrSurf[3]*fSkin[11]+alphaDrSurf[2]*fSkin[10]+alphaDrSurf[1]*fSkin[9]+alphaDrSurf[0]*fSkin[4])))/vmap_prime_skin[1]); 
  Ghat[5] = -((0.125*(2.4494897427831783*(alphaDrSurf[3]*fSkin[27]+alphaDrSurf[6]*fSkin[22]+alphaDrSurf[7]*fSkin[21]+alphaDrSurf[0]*fSkin[20])-1.4142135623730951*alphaDrSurf[3]*fSkin[16]+2.4494897427831783*(alphaDrSurf[11]*fSkin[14]+alphaDrSurf[1]*fSkin[13]+alphaDrSurf[2]*fSkin[12])-1.4142135623730951*(fSkin[3]*alphaDrSurf[11]+alphaDrSurf[6]*fSkin[8]+alphaDrSurf[7]*fSkin[7]+alphaDrSurf[0]*fSkin[6])+2.4494897427831783*alphaDrSurf[5]*fSkin[5]-1.4142135623730951*(fSkin[0]*alphaDrSurf[5]+alphaDrSurf[1]*fSkin[2]+fSkin[1]*alphaDrSurf[2])))/vmap_prime_skin[1]); 
  Ghat[6] = -((0.125*(2.4494897427831783*(alphaDrSurf[2]*fSkin[27]+alphaDrSurf[5]*fSkin[22]+alphaDrSurf[0]*fSkin[21]+alphaDrSurf[7]*fSkin[20])-1.4142135623730951*alphaDrSurf[2]*fSkin[16]+2.4494897427831783*(alphaDrSurf[1]*fSkin[14]+alphaDrSurf[11]*fSkin[13]+alphaDrSurf[3]*fSkin[12])-1.4142135623730951*(fSkin[2]*alphaDrSurf[11]+alphaDrSurf[5]*fSkin[8]+alphaDrSurf[0]*fSkin[7]+fSkin[6]*alphaDrSurf[7])+(2.4494897427831783*fSkin[5]-1.4142135623730951*fSkin[0])*alphaDrSurf[6]-1.4142135623730951*(alphaDrSurf[1]*fSkin[3]+fSkin[1]*alphaDrSurf[3])))/vmap_prime_skin[1]); 
  Ghat[7] = -((0.125*(2.4494897427831783*(alphaDrSurf[1]*fSkin[27]+alphaDrSurf[0]*fSkin[22]+alphaDrSurf[5]*fSkin[21]+alphaDrSurf[6]*fSkin[20])-1.4142135623730951*alphaDrSurf[1]*fSkin[16]+2.4494897427831783*(alphaDrSurf[2]*fSkin[14]+alphaDrSurf[3]*fSkin[13]+alphaDrSurf[11]*fSkin[12])-1.4142135623730951*(fSkin[1]*alphaDrSurf[11]+alphaDrSurf[0]*fSkin[8]+alphaDrSurf[5]*fSkin[7])+(2.4494897427831783*fSkin[5]-1.4142135623730951*fSkin[0])*alphaDrSurf[7]-1.4142135623730951*(alphaDrSurf[6]*fSkin[6]+alphaDrSurf[2]*fSkin[3]+fSkin[2]*alphaDrSurf[3])))/vmap_prime_skin[1]); 
  Ghat[8] = -((0.125*(2.4494897427831783*(alphaDrSurf[7]*fSkin[31]+alphaDrSurf[11]*fSkin[30]+alphaDrSurf[3]*fSkin[29]+alphaDrSurf[2]*fSkin[28])-1.4142135623730951*alphaDrSurf[7]*fSkin[26]+2.4494897427831783*(alphaDrSurf[6]*fSkin[25]+alphaDrSurf[5]*fSkin[24]+alphaDrSurf[0]*fSkin[23])-1.4142135623730951*(alphaDrSurf[11]*fSkin[19]+alphaDrSurf[3]*fSkin[18]+alphaDrSurf[2]*fSkin[17])+2.4494897427831783*alphaDrSurf[1]*fSkin[15]-1.4142135623730951*(alphaDrSurf[6]*fSkin[11]+alphaDrSurf[5]*fSkin[10]+alphaDrSurf[0]*fSkin[9]+alphaDrSurf[1]*fSkin[4])))/vmap_prime_skin[1]); 
  Ghat[9] = -((0.125*(2.4494897427831783*(alphaDrSurf[6]*fSkin[31]+alphaDrSurf[3]*fSkin[30]+alphaDrSurf[11]*fSkin[29]+alphaDrSurf[1]*fSkin[28])-1.4142135623730951*alphaDrSurf[6]*fSkin[26]+2.4494897427831783*(alphaDrSurf[7]*fSkin[25]+alphaDrSurf[0]*fSkin[24]+alphaDrSurf[5]*fSkin[23])-1.4142135623730951*(alphaDrSurf[3]*fSkin[19]+alphaDrSurf[11]*fSkin[18]+alphaDrSurf[1]*fSkin[17])+2.4494897427831783*alphaDrSurf[2]*fSkin[15]-1.4142135623730951*(alphaDrSurf[7]*fSkin[11]+alphaDrSurf[0]*fSkin[10]+alphaDrSurf[5]*fSkin[9]+alphaDrSurf[2]*fSkin[4])))/vmap_prime_skin[1]); 
  Ghat[10] = -((0.125*(2.4494897427831783*(alphaDrSurf[5]*fSkin[31]+alphaDrSurf[2]*fSkin[30]+alphaDrSurf[1]*fSkin[29]+alphaDrSurf[11]*fSkin[28])-1.4142135623730951*alphaDrSurf[5]*fSkin[26]+2.4494897427831783*(alphaDrSurf[0]*fSkin[25]+alphaDrSurf[7]*fSkin[24]+alphaDrSurf[6]*fSkin[23])-1.4142135623730951*(alphaDrSurf[2]*fSkin[19]+alphaDrSurf[1]*fSkin[18]+alphaDrSurf[11]*fSkin[17])+2.4494897427831783*alphaDrSurf[3]*fSkin[15]-1.4142135623730951*(alphaDrSurf[0]*fSkin[11]+alphaDrSurf[7]*fSkin[10]+alphaDrSurf[6]*fSkin[9]+alphaDrSurf[3]*fSkin[4])))/vmap_prime_skin[1]); 
  Ghat[11] = -((0.125*(2.4494897427831783*(alphaDrSurf[0]*fSkin[27]+alphaDrSurf[1]*fSkin[22]+alphaDrSurf[2]*fSkin[21]+alphaDrSurf[3]*fSkin[20])-1.4142135623730951*alphaDrSurf[0]*fSkin[16]+2.4494897427831783*(alphaDrSurf[5]*fSkin[14]+alphaDrSurf[6]*fSkin[13]+alphaDrSurf[7]*fSkin[12])+(2.4494897427831783*fSkin[5]-1.4142135623730951*fSkin[0])*alphaDrSurf[11]-1.4142135623730951*(alphaDrSurf[1]*fSkin[8]+alphaDrSurf[2]*fSkin[7]+fSkin[1]*alphaDrSurf[7]+alphaDrSurf[3]*fSkin[6]+fSkin[2]*alphaDrSurf[6]+fSkin[3]*alphaDrSurf[5])))/vmap_prime_skin[1]); 
  Ghat[12] = -((0.125*(2.4494897427831783*(alphaDrSurf[3]*fSkin[31]+alphaDrSurf[6]*fSkin[30]+alphaDrSurf[7]*fSkin[29]+alphaDrSurf[0]*fSkin[28])-1.4142135623730951*alphaDrSurf[3]*fSkin[26]+2.4494897427831783*(alphaDrSurf[11]*fSkin[25]+alphaDrSurf[1]*fSkin[24]+alphaDrSurf[2]*fSkin[23])-1.4142135623730951*(alphaDrSurf[6]*fSkin[19]+alphaDrSurf[7]*fSkin[18]+alphaDrSurf[0]*fSkin[17])+2.4494897427831783*alphaDrSurf[5]*fSkin[15]-1.4142135623730951*(alphaDrSurf[11]*fSkin[11]+alphaDrSurf[1]*fSkin[10]+alphaDrSurf[2]*fSkin[9]+fSkin[4]*alphaDrSurf[5])))/vmap_prime_skin[1]); 
  Ghat[13] = -((0.125*(2.4494897427831783*(alphaDrSurf[2]*fSkin[31]+alphaDrSurf[5]*fSkin[30]+alphaDrSurf[0]*fSkin[29]+alphaDrSurf[7]*fSkin[28])-1.4142135623730951*alphaDrSurf[2]*fSkin[26]+2.4494897427831783*(alphaDrSurf[1]*fSkin[25]+alphaDrSurf[11]*fSkin[24]+alphaDrSurf[3]*fSkin[23])-1.4142135623730951*(alphaDrSurf[5]*fSkin[19]+alphaDrSurf[0]*fSkin[18]+alphaDrSurf[7]*fSkin[17])+2.4494897427831783*alphaDrSurf[6]*fSkin[15]-1.4142135623730951*(alphaDrSurf[1]*fSkin[11]+fSkin[10]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[9]+fSkin[4]*alphaDrSurf[6])))/vmap_prime_skin[1]); 
  Ghat[14] = -((0.125*(2.4494897427831783*(alphaDrSurf[1]*fSkin[31]+alphaDrSurf[0]*fSkin[30]+alphaDrSurf[5]*fSkin[29]+alphaDrSurf[6]*fSkin[28])-1.4142135623730951*alphaDrSurf[1]*fSkin[26]+2.4494897427831783*(alphaDrSurf[2]*fSkin[25]+alphaDrSurf[3]*fSkin[24]+alphaDrSurf[11]*fSkin[23])-1.4142135623730951*(alphaDrSurf[0]*fSkin[19]+alphaDrSurf[5]*fSkin[18]+alphaDrSurf[6]*fSkin[17])+2.4494897427831783*alphaDrSurf[7]*fSkin[15]-1.4142135623730951*(alphaDrSurf[2]*fSkin[11]+fSkin[9]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[10]+fSkin[4]*alphaDrSurf[7])))/vmap_prime_skin[1]); 
  Ghat[15] = -((0.125*(2.4494897427831783*(alphaDrSurf[0]*fSkin[31]+alphaDrSurf[1]*fSkin[30]+alphaDrSurf[2]*fSkin[29]+alphaDrSurf[3]*fSkin[28])-1.4142135623730951*alphaDrSurf[0]*fSkin[26]+2.4494897427831783*(alphaDrSurf[5]*fSkin[25]+alphaDrSurf[6]*fSkin[24]+alphaDrSurf[7]*fSkin[23])-1.4142135623730951*(alphaDrSurf[1]*fSkin[19]+alphaDrSurf[2]*fSkin[18]+alphaDrSurf[3]*fSkin[17])+2.4494897427831783*alphaDrSurf[11]*fSkin[15]-1.4142135623730951*(alphaDrSurf[5]*fSkin[11]+fSkin[4]*alphaDrSurf[11]+alphaDrSurf[6]*fSkin[10]+alphaDrSurf[7]*fSkin[9])))/vmap_prime_skin[1]); 
  Ghat[16] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[11]*fSkin[47]+36.74234614174768*(alphaDrSurf[7]*fSkin[46]+alphaDrSurf[6]*fSkin[45]+alphaDrSurf[5]*fSkin[44])-21.21320343559643*alphaDrSurf[11]*fSkin[43]+36.74234614174767*(alphaDrSurf[3]*fSkin[42]+alphaDrSurf[2]*fSkin[41]+alphaDrSurf[1]*fSkin[40])-21.213203435596427*(alphaDrSurf[7]*fSkin[39]+alphaDrSurf[6]*fSkin[38]+alphaDrSurf[5]*fSkin[37])+36.74234614174768*alphaDrSurf[0]*fSkin[36]-21.21320343559643*(alphaDrSurf[3]*fSkin[35]+alphaDrSurf[2]*fSkin[34]+alphaDrSurf[1]*fSkin[33])-21.213203435596427*alphaDrSurf[0]*fSkin[32]))/vmap_prime_skin[1]); 
  Ghat[17] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[7]*fSkin[47]+36.74234614174767*(alphaDrSurf[11]*fSkin[46]+alphaDrSurf[3]*fSkin[45]+alphaDrSurf[2]*fSkin[44])-21.213203435596427*alphaDrSurf[7]*fSkin[43]+36.74234614174768*(alphaDrSurf[6]*fSkin[42]+alphaDrSurf[5]*fSkin[41]+alphaDrSurf[0]*fSkin[40])-21.21320343559643*(alphaDrSurf[11]*fSkin[39]+alphaDrSurf[3]*fSkin[38]+alphaDrSurf[2]*fSkin[37])+36.74234614174767*alphaDrSurf[1]*fSkin[36]-21.213203435596427*(alphaDrSurf[6]*fSkin[35]+alphaDrSurf[5]*fSkin[34]+alphaDrSurf[0]*fSkin[33])-21.21320343559643*alphaDrSurf[1]*fSkin[32]))/vmap_prime_skin[1]); 
  Ghat[18] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[6]*fSkin[47]+36.74234614174767*(alphaDrSurf[3]*fSkin[46]+alphaDrSurf[11]*fSkin[45]+alphaDrSurf[1]*fSkin[44])-21.213203435596427*alphaDrSurf[6]*fSkin[43]+36.74234614174768*(alphaDrSurf[7]*fSkin[42]+alphaDrSurf[0]*fSkin[41]+alphaDrSurf[5]*fSkin[40])-21.21320343559643*(alphaDrSurf[3]*fSkin[39]+alphaDrSurf[11]*fSkin[38]+alphaDrSurf[1]*fSkin[37])+36.74234614174767*alphaDrSurf[2]*fSkin[36]-21.213203435596427*(alphaDrSurf[7]*fSkin[35]+alphaDrSurf[0]*fSkin[34]+alphaDrSurf[5]*fSkin[33])-21.21320343559643*alphaDrSurf[2]*fSkin[32]))/vmap_prime_skin[1]); 
  Ghat[19] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[5]*fSkin[47]+36.74234614174767*(alphaDrSurf[2]*fSkin[46]+alphaDrSurf[1]*fSkin[45]+alphaDrSurf[11]*fSkin[44])-21.213203435596427*alphaDrSurf[5]*fSkin[43]+36.74234614174768*(alphaDrSurf[0]*fSkin[42]+alphaDrSurf[7]*fSkin[41]+alphaDrSurf[6]*fSkin[40])-21.21320343559643*(alphaDrSurf[2]*fSkin[39]+alphaDrSurf[1]*fSkin[38]+alphaDrSurf[11]*fSkin[37])+36.74234614174767*alphaDrSurf[3]*fSkin[36]-21.213203435596427*(alphaDrSurf[0]*fSkin[35]+alphaDrSurf[7]*fSkin[34]+alphaDrSurf[6]*fSkin[33])-21.21320343559643*alphaDrSurf[3]*fSkin[32]))/vmap_prime_skin[1]); 
  Ghat[20] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[3]*fSkin[47]+36.74234614174768*(alphaDrSurf[6]*fSkin[46]+alphaDrSurf[7]*fSkin[45]+alphaDrSurf[0]*fSkin[44])-21.21320343559643*alphaDrSurf[3]*fSkin[43]+36.74234614174767*(alphaDrSurf[11]*fSkin[42]+alphaDrSurf[1]*fSkin[41]+alphaDrSurf[2]*fSkin[40])-21.213203435596427*(alphaDrSurf[6]*fSkin[39]+alphaDrSurf[7]*fSkin[38]+alphaDrSurf[0]*fSkin[37])+36.74234614174768*alphaDrSurf[5]*fSkin[36]-21.21320343559643*(alphaDrSurf[11]*fSkin[35]+alphaDrSurf[1]*fSkin[34]+alphaDrSurf[2]*fSkin[33])-21.213203435596427*alphaDrSurf[5]*fSkin[32]))/vmap_prime_skin[1]); 
  Ghat[21] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[2]*fSkin[47]+36.74234614174768*(alphaDrSurf[5]*fSkin[46]+alphaDrSurf[0]*fSkin[45]+alphaDrSurf[7]*fSkin[44])-21.21320343559643*alphaDrSurf[2]*fSkin[43]+36.74234614174767*(alphaDrSurf[1]*fSkin[42]+alphaDrSurf[11]*fSkin[41]+alphaDrSurf[3]*fSkin[40])-21.213203435596427*(alphaDrSurf[5]*fSkin[39]+alphaDrSurf[0]*fSkin[38]+alphaDrSurf[7]*fSkin[37])+36.74234614174768*alphaDrSurf[6]*fSkin[36]-21.21320343559643*(alphaDrSurf[1]*fSkin[35]+alphaDrSurf[11]*fSkin[34]+alphaDrSurf[3]*fSkin[33])-21.213203435596427*alphaDrSurf[6]*fSkin[32]))/vmap_prime_skin[1]); 
  Ghat[22] = -((0.008333333333333333*(36.74234614174767*alphaDrSurf[1]*fSkin[47]+36.74234614174768*(alphaDrSurf[0]*fSkin[46]+alphaDrSurf[5]*fSkin[45]+alphaDrSurf[6]*fSkin[44])-21.21320343559643*alphaDrSurf[1]*fSkin[43]+36.74234614174767*(alphaDrSurf[2]*fSkin[42]+alphaDrSurf[3]*fSkin[41]+alphaDrSurf[11]*fSkin[40])-21.213203435596427*(alphaDrSurf[0]*fSkin[39]+alphaDrSurf[5]*fSkin[38]+alphaDrSurf[6]*fSkin[37])+36.74234614174768*alphaDrSurf[7]*fSkin[36]-21.21320343559643*(alphaDrSurf[2]*fSkin[35]+alphaDrSurf[3]*fSkin[34]+alphaDrSurf[11]*fSkin[33])-21.213203435596427*alphaDrSurf[7]*fSkin[32]))/vmap_prime_skin[1]); 
  Ghat[23] = -((0.008333333333333333*(36.74234614174768*alphaDrSurf[0]*fSkin[47]+36.74234614174767*(alphaDrSurf[1]*fSkin[46]+alphaDrSurf[2]*fSkin[45]+alphaDrSurf[3]*fSkin[44])-21.213203435596427*alphaDrSurf[0]*fSkin[43]+36.74234614174768*(alphaDrSurf[5]*fSkin[42]+alphaDrSurf[6]*fSkin[41]+alphaDrSurf[7]*fSkin[40])-21.21320343559643*(alphaDrSurf[1]*fSkin[39]+alphaDrSurf[2]*fSkin[38]+alphaDrSurf[3]*fSkin[37])+36.74234614174767*alphaDrSurf[11]*fSkin[36]-21.213203435596427*(alphaDrSurf[5]*fSkin[35]+alphaDrSurf[6]*fSkin[34]+alphaDrSurf[7]*fSkin[33])-21.21320343559643*alphaDrSurf[11]*fSkin[32]))/vmap_prime_skin[1]); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[4] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[5] += 1.224744871391589*Ghat[0]*rdv2; 
  out[6] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[7] += -(0.7071067811865475*Ghat[6]*rdv2); 
  out[8] += -(0.7071067811865475*Ghat[7]*rdv2); 
  out[9] += -(0.7071067811865475*Ghat[8]*rdv2); 
  out[10] += -(0.7071067811865475*Ghat[9]*rdv2); 
  out[11] += -(0.7071067811865475*Ghat[10]*rdv2); 
  out[12] += 1.224744871391589*Ghat[1]*rdv2; 
  out[13] += 1.224744871391589*Ghat[2]*rdv2; 
  out[14] += 1.224744871391589*Ghat[3]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += -(0.7071067811865475*Ghat[11]*rdv2); 
  out[17] += -(0.7071067811865475*Ghat[12]*rdv2); 
  out[18] += -(0.7071067811865475*Ghat[13]*rdv2); 
  out[19] += -(0.7071067811865475*Ghat[14]*rdv2); 
  out[20] += 1.224744871391589*Ghat[5]*rdv2; 
  out[21] += 1.224744871391589*Ghat[6]*rdv2; 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.224744871391589*Ghat[9]*rdv2; 
  out[25] += 1.224744871391589*Ghat[10]*rdv2; 
  out[26] += -(0.7071067811865475*Ghat[15]*rdv2); 
  out[27] += 1.224744871391589*Ghat[11]*rdv2; 
  out[28] += 1.224744871391589*Ghat[12]*rdv2; 
  out[29] += 1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat[15]*rdv2; 
  out[32] += -(0.7071067811865475*Ghat[16]*rdv2); 
  out[33] += -(0.7071067811865475*Ghat[17]*rdv2); 
  out[34] += -(0.7071067811865475*Ghat[18]*rdv2); 
  out[35] += -(0.7071067811865475*Ghat[19]*rdv2); 
  out[36] += 1.224744871391589*Ghat[16]*rdv2; 
  out[37] += -(0.7071067811865475*Ghat[20]*rdv2); 
  out[38] += -(0.7071067811865475*Ghat[21]*rdv2); 
  out[39] += -(0.7071067811865475*Ghat[22]*rdv2); 
  out[40] += 1.224744871391589*Ghat[17]*rdv2; 
  out[41] += 1.224744871391589*Ghat[18]*rdv2; 
  out[42] += 1.224744871391589*Ghat[19]*rdv2; 
  out[43] += -(0.7071067811865475*Ghat[23]*rdv2); 
  out[44] += 1.224744871391589*Ghat[20]*rdv2; 
  out[45] += 1.224744871391589*Ghat[21]*rdv2; 
  out[46] += 1.224744871391589*Ghat[22]*rdv2; 
  out[47] += 1.224744871391589*Ghat[23]*rdv2; 

  } 

  return 0.;

} 
