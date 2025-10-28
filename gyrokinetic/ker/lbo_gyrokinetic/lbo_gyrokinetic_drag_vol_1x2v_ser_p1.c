#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // dxv[3]: cell spacing. 
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
  rdv2[0] = 2.0/dxv[1]; 
  rdv2[1] = 2.0/dxv[2]; 

  double alphaDrag[24] = {0.}; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = (rdv2[0]*(2.0*nuUSum[0]-1.414213562373095*nuSum[0]*vmap[0]))/vmap_prime[0]; 
  alphaDrag[1] = (rdv2[0]*(2.0*nuUSum[1]-1.414213562373095*vmap[0]*nuSum[1]))/vmap_prime[0]; 
  alphaDrag[2] = -(1.414213562373095*nuSum[0]*rdv2[0]*vmap[1])/vmap_prime[0]; 
  alphaDrag[4] = -(1.414213562373095*rdv2[0]*nuSum[1]*vmap[1])/vmap_prime[0]; 

  // Expand rdv2*nu*2*mu in phase basis.
  alphaDrag[12] = -(2.828427124746191*nuSum[0]*rdv2[1]*vmap[2])/vmap_prime[1]; 
  alphaDrag[13] = -(2.828427124746191*nuSum[1]*rdv2[1]*vmap[2])/vmap_prime[1]; 
  alphaDrag[15] = -(2.828427124746191*nuSum[0]*rdv2[1]*vmap[3])/vmap_prime[1]; 
  alphaDrag[17] = -(2.828427124746191*nuSum[1]*rdv2[1]*vmap[3])/vmap_prime[1]; 

  out[2] += 0.6123724356957944*(alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[5]*alphaDrag[17]+f[3]*alphaDrag[15]+f[1]*alphaDrag[13]+f[0]*alphaDrag[12]); 
  out[4] += 0.6123724356957944*(alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(f[3]*alphaDrag[17]+f[5]*alphaDrag[15]+f[0]*alphaDrag[13]+f[1]*alphaDrag[12]); 
  out[6] += 0.6123724356957944*(f[7]*alphaDrag[17]+f[6]*alphaDrag[15]+f[4]*alphaDrag[13]+f[2]*alphaDrag[12]+alphaDrag[4]*f[7]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[6]*alphaDrag[17]+f[7]*alphaDrag[15]+f[2]*alphaDrag[13]+f[4]*alphaDrag[12]+alphaDrag[2]*f[7]+alphaDrag[4]*f[6]+alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[8] += 1.224744871391589*(alphaDrag[4]*f[9]+alphaDrag[2]*f[8])+1.369306393762915*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[2]*f[9]+alphaDrag[4]*f[8])+1.369306393762915*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[10] += 0.6123724356957944*(f[11]*alphaDrag[17]+f[10]*alphaDrag[15]+f[9]*alphaDrag[13]+f[8]*alphaDrag[12])+1.224744871391589*(alphaDrag[4]*f[11]+alphaDrag[2]*f[10])+1.369306393762915*(alphaDrag[1]*f[7]+alphaDrag[0]*f[6]+alphaDrag[4]*f[5]+alphaDrag[2]*f[3]); 
  out[11] += 0.6123724356957944*(f[10]*alphaDrag[17]+f[11]*alphaDrag[15]+f[8]*alphaDrag[13]+f[9]*alphaDrag[12])+1.224744871391589*(alphaDrag[2]*f[11]+alphaDrag[4]*f[10])+1.369306393762915*(alphaDrag[0]*f[7]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]+f[3]*alphaDrag[4]); 

  return fabs(0.8838834764831842*alphaDrag[0])+fabs(0.5303300858899105*alphaDrag[12]); 

} 
