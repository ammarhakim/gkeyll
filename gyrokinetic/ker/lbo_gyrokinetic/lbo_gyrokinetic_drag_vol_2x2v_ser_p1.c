#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_2x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // dxv[4]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // f: input distribution function.
  // out: incremented output 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2[2]; 
  rdv2[0] = 2.0/dxv[2]; 
  rdv2[1] = 2.0/dxv[3]; 

  double alphaDrag[48] = {0.}; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = (rdv2[0]*(2.0*nuUSum[0]-1.414213562373095*nuSum[0]*vmap[0]))/vmap_prime[0]; 
  alphaDrag[1] = (rdv2[0]*(2.0*nuUSum[1]-1.414213562373095*vmap[0]*nuSum[1]))/vmap_prime[0]; 
  alphaDrag[2] = (rdv2[0]*(2.0*nuUSum[2]-1.414213562373095*vmap[0]*nuSum[2]))/vmap_prime[0]; 
  alphaDrag[3] = -(1.414213562373095*nuSum[0]*rdv2[0]*vmap[1])/vmap_prime[0]; 
  alphaDrag[5] = (rdv2[0]*(2.0*nuUSum[3]-1.414213562373095*vmap[0]*nuSum[3]))/vmap_prime[0]; 
  alphaDrag[6] = -(1.414213562373095*rdv2[0]*nuSum[1]*vmap[1])/vmap_prime[0]; 
  alphaDrag[7] = -(1.414213562373095*rdv2[0]*vmap[1]*nuSum[2])/vmap_prime[0]; 
  alphaDrag[11] = -(1.414213562373095*rdv2[0]*vmap[1]*nuSum[3])/vmap_prime[0]; 

  // Expand rdv2*nu*2*mu in phase basis.
  alphaDrag[24] = -(2.828427124746191*nuSum[0]*rdv2[1]*vmap[2])/vmap_prime[1]; 
  alphaDrag[25] = -(2.828427124746191*nuSum[1]*rdv2[1]*vmap[2])/vmap_prime[1]; 
  alphaDrag[26] = -(2.828427124746191*rdv2[1]*nuSum[2]*vmap[2])/vmap_prime[1]; 
  alphaDrag[28] = -(2.828427124746191*nuSum[0]*rdv2[1]*vmap[3])/vmap_prime[1]; 
  alphaDrag[29] = -(2.828427124746191*rdv2[1]*vmap[2]*nuSum[3])/vmap_prime[1]; 
  alphaDrag[32] = -(2.828427124746191*nuSum[1]*rdv2[1]*vmap[3])/vmap_prime[1]; 
  alphaDrag[33] = -(2.828427124746191*rdv2[1]*nuSum[2]*vmap[3])/vmap_prime[1]; 
  alphaDrag[36] = -(2.828427124746191*rdv2[1]*nuSum[3]*vmap[3])/vmap_prime[1]; 

  out[3] += 0.4330127018922193*(alphaDrag[11]*f[11]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[12]*alphaDrag[36]+f[9]*alphaDrag[33]+f[8]*alphaDrag[32]+f[5]*alphaDrag[29]+f[4]*alphaDrag[28]+f[2]*alphaDrag[26]+f[1]*alphaDrag[25]+f[0]*alphaDrag[24]); 
  out[6] += 0.4330127018922193*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[7] += 0.4330127018922193*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[8] += 0.4330127018922193*(f[9]*alphaDrag[36]+f[12]*alphaDrag[33]+f[4]*alphaDrag[32]+f[2]*alphaDrag[29]+f[8]*alphaDrag[28]+f[5]*alphaDrag[26]+f[0]*alphaDrag[25]+f[1]*alphaDrag[24]); 
  out[9] += 0.4330127018922193*(f[8]*alphaDrag[36]+f[4]*alphaDrag[33]+f[12]*alphaDrag[32]+f[1]*alphaDrag[29]+f[9]*alphaDrag[28]+f[0]*alphaDrag[26]+f[5]*alphaDrag[25]+f[2]*alphaDrag[24]); 
  out[10] += 0.4330127018922193*(f[15]*alphaDrag[36]+f[14]*alphaDrag[33]+f[13]*alphaDrag[32]+f[11]*alphaDrag[29]+f[10]*alphaDrag[28]+f[7]*alphaDrag[26]+f[6]*alphaDrag[25]+f[3]*alphaDrag[24]+alphaDrag[11]*f[15]+alphaDrag[7]*f[14]+alphaDrag[6]*f[13]+alphaDrag[5]*f[12]+alphaDrag[3]*f[10]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[11] += 0.4330127018922193*(alphaDrag[3]*f[11]+f[3]*alphaDrag[11]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[12] += 0.4330127018922193*(f[4]*alphaDrag[36]+f[8]*alphaDrag[33]+f[9]*alphaDrag[32]+f[0]*alphaDrag[29]+f[12]*alphaDrag[28]+f[1]*alphaDrag[26]+f[2]*alphaDrag[25]+f[5]*alphaDrag[24]); 
  out[13] += 0.4330127018922193*(f[14]*alphaDrag[36]+f[15]*alphaDrag[33]+f[10]*alphaDrag[32]+f[7]*alphaDrag[29]+f[13]*alphaDrag[28]+f[11]*alphaDrag[26]+f[3]*alphaDrag[25]+f[6]*alphaDrag[24]+alphaDrag[7]*f[15]+alphaDrag[11]*f[14]+alphaDrag[3]*f[13]+alphaDrag[2]*f[12]+alphaDrag[6]*f[10]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[14] += 0.4330127018922193*(f[13]*alphaDrag[36]+f[10]*alphaDrag[33]+f[15]*alphaDrag[32]+f[6]*alphaDrag[29]+f[14]*alphaDrag[28]+f[3]*alphaDrag[26]+f[11]*alphaDrag[25]+f[7]*alphaDrag[24]+alphaDrag[6]*f[15]+alphaDrag[3]*f[14]+alphaDrag[11]*f[13]+alphaDrag[1]*f[12]+alphaDrag[7]*f[10]+alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[15] += 0.4330127018922193*(f[10]*alphaDrag[36]+f[13]*alphaDrag[33]+f[14]*alphaDrag[32]+f[3]*alphaDrag[29]+f[15]*alphaDrag[28]+f[6]*alphaDrag[26]+f[7]*alphaDrag[25]+f[11]*alphaDrag[24]+alphaDrag[3]*f[15]+alphaDrag[6]*f[14]+alphaDrag[7]*f[13]+alphaDrag[0]*f[12]+f[10]*alphaDrag[11]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]+f[4]*alphaDrag[5]); 
  out[16] += 0.8660254037844386*(alphaDrag[11]*f[20]+alphaDrag[7]*f[18]+alphaDrag[6]*f[17]+alphaDrag[3]*f[16])+0.9682458365518543*(alphaDrag[5]*f[11]+f[5]*alphaDrag[11]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[17] += 0.8660254037844386*(alphaDrag[7]*f[20]+alphaDrag[11]*f[18]+alphaDrag[3]*f[17]+alphaDrag[6]*f[16])+0.9682458365518543*(alphaDrag[2]*f[11]+f[2]*alphaDrag[11]+alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[18] += 0.8660254037844386*(alphaDrag[6]*f[20]+alphaDrag[3]*f[18]+alphaDrag[11]*f[17]+alphaDrag[7]*f[16])+0.9682458365518543*(alphaDrag[1]*f[11]+f[1]*alphaDrag[11]+alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[19] += 0.4330127018922193*(f[23]*alphaDrag[36]+f[22]*alphaDrag[33]+f[21]*alphaDrag[32]+f[20]*alphaDrag[29]+f[19]*alphaDrag[28]+f[18]*alphaDrag[26]+f[17]*alphaDrag[25]+f[16]*alphaDrag[24])+0.8660254037844386*(alphaDrag[11]*f[23]+alphaDrag[7]*f[22]+alphaDrag[6]*f[21]+alphaDrag[3]*f[19])+0.9682458365518543*(alphaDrag[5]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[11]*f[12]+alphaDrag[0]*f[10]+alphaDrag[7]*f[9]+alphaDrag[6]*f[8]+alphaDrag[3]*f[4]); 
  out[20] += 0.8660254037844386*(alphaDrag[3]*f[20]+alphaDrag[6]*f[18]+alphaDrag[7]*f[17]+alphaDrag[11]*f[16])+0.9682458365518543*(alphaDrag[0]*f[11]+f[0]*alphaDrag[11]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]); 
  out[21] += 0.4330127018922193*(f[22]*alphaDrag[36]+f[23]*alphaDrag[33]+f[19]*alphaDrag[32]+f[18]*alphaDrag[29]+f[21]*alphaDrag[28]+f[20]*alphaDrag[26]+f[16]*alphaDrag[25]+f[17]*alphaDrag[24])+0.8660254037844386*(alphaDrag[7]*f[23]+alphaDrag[11]*f[22]+alphaDrag[3]*f[21]+alphaDrag[6]*f[19])+0.9682458365518543*(alphaDrag[2]*f[15]+alphaDrag[5]*f[14]+alphaDrag[0]*f[13]+alphaDrag[7]*f[12]+f[9]*alphaDrag[11]+alphaDrag[1]*f[10]+alphaDrag[3]*f[8]+f[4]*alphaDrag[6]); 
  out[22] += 0.4330127018922193*(f[21]*alphaDrag[36]+f[19]*alphaDrag[33]+f[23]*alphaDrag[32]+f[17]*alphaDrag[29]+f[22]*alphaDrag[28]+f[16]*alphaDrag[26]+f[20]*alphaDrag[25]+f[18]*alphaDrag[24])+0.8660254037844386*(alphaDrag[6]*f[23]+alphaDrag[3]*f[22]+alphaDrag[11]*f[21]+alphaDrag[7]*f[19])+0.9682458365518543*(alphaDrag[1]*f[15]+alphaDrag[0]*f[14]+alphaDrag[5]*f[13]+alphaDrag[6]*f[12]+f[8]*alphaDrag[11]+alphaDrag[2]*f[10]+alphaDrag[3]*f[9]+f[4]*alphaDrag[7]); 
  out[23] += 0.4330127018922193*(f[19]*alphaDrag[36]+f[21]*alphaDrag[33]+f[22]*alphaDrag[32]+f[16]*alphaDrag[29]+f[23]*alphaDrag[28]+f[17]*alphaDrag[26]+f[18]*alphaDrag[25]+f[20]*alphaDrag[24])+0.8660254037844386*(alphaDrag[3]*f[23]+alphaDrag[6]*f[22]+alphaDrag[7]*f[21]+alphaDrag[11]*f[19])+0.9682458365518543*(alphaDrag[0]*f[15]+alphaDrag[1]*f[14]+alphaDrag[2]*f[13]+alphaDrag[3]*f[12]+f[4]*alphaDrag[11]+alphaDrag[5]*f[10]+alphaDrag[6]*f[9]+alphaDrag[7]*f[8]); 

  return fabs(0.625*alphaDrag[0])+fabs(0.375*alphaDrag[24]); 

} 
