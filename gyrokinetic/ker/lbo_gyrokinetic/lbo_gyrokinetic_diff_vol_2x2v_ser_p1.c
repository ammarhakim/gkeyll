#include <gkyl_lbo_gyrokinetic_kernels.h> 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_2x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // dxv[4]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fin: input distribution function.
  // out: incremented output 

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvVmapPrimeSq4[2]; 
  rdvVmapPrimeSq4[0] = 4.0/(dxv[2]*dxv[2]*vmap_prime[0]*vmap_prime[0]); 
  rdvVmapPrimeSq4[1] = 4.0/(dxv[3]*dxv[3]*vmap_prime[1]*vmap_prime[1]); 

  // Expand nuVtSqSum/vpar'^2 in conf basis.
  double facDiffVpar[4] = {0.};
  facDiffVpar[0] = nuVtSqSum[0]*rdvVmapPrimeSq4[0]; 
  facDiffVpar[1] = rdvVmapPrimeSq4[0]*nuVtSqSum[1]; 
  facDiffVpar[2] = rdvVmapPrimeSq4[0]*nuVtSqSum[2]; 
  facDiffVpar[3] = rdvVmapPrimeSq4[0]*nuVtSqSum[3]; 

  // Expand 2*m*nuVtSqSum/bmag/mu'^2 in conf basis.
  double facDiffMu[4] = {0.};
  facDiffMu[0] = (bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*rdvVmapPrimeSq4[1]*m_; 
  facDiffMu[1] = (bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*rdvVmapPrimeSq4[1]*m_; 
  facDiffMu[2] = rdvVmapPrimeSq4[1]*(bmag_inv[1]*nuVtSqSum[3]+bmag_inv[0]*nuVtSqSum[2])*m_; 
  facDiffMu[3] = rdvVmapPrimeSq4[1]*(bmag_inv[0]*nuVtSqSum[3]+bmag_inv[1]*nuVtSqSum[2])*m_; 

  out[4] += 1.060660171779821*facDiffMu[3]*vmap[3]*fin[5]+1.060660171779821*facDiffMu[2]*fin[2]*vmap[3]+1.060660171779821*facDiffMu[1]*fin[1]*vmap[3]+1.060660171779821*facDiffMu[0]*fin[0]*vmap[3]; 
  out[8] += 1.060660171779821*facDiffMu[2]*vmap[3]*fin[5]+1.060660171779821*fin[2]*facDiffMu[3]*vmap[3]+1.060660171779821*facDiffMu[0]*fin[1]*vmap[3]+1.060660171779821*fin[0]*facDiffMu[1]*vmap[3]; 
  out[9] += 1.060660171779821*facDiffMu[1]*vmap[3]*fin[5]+1.060660171779821*fin[1]*facDiffMu[3]*vmap[3]+1.060660171779821*facDiffMu[0]*fin[2]*vmap[3]+1.060660171779821*fin[0]*facDiffMu[2]*vmap[3]; 
  out[10] += 1.060660171779821*facDiffMu[3]*vmap[3]*fin[11]+1.060660171779821*facDiffMu[2]*vmap[3]*fin[7]+1.060660171779821*facDiffMu[1]*vmap[3]*fin[6]+1.060660171779821*facDiffMu[0]*fin[3]*vmap[3]; 
  out[12] += 1.060660171779821*facDiffMu[0]*vmap[3]*fin[5]+1.060660171779821*fin[0]*facDiffMu[3]*vmap[3]+1.060660171779821*facDiffMu[1]*fin[2]*vmap[3]+1.060660171779821*fin[1]*facDiffMu[2]*vmap[3]; 
  out[13] += 1.060660171779821*facDiffMu[2]*vmap[3]*fin[11]+1.060660171779821*facDiffMu[3]*vmap[3]*fin[7]+1.060660171779821*facDiffMu[0]*vmap[3]*fin[6]+1.060660171779821*facDiffMu[1]*fin[3]*vmap[3]; 
  out[14] += 1.060660171779821*facDiffMu[1]*vmap[3]*fin[11]+1.060660171779821*facDiffMu[0]*vmap[3]*fin[7]+1.060660171779821*facDiffMu[3]*vmap[3]*fin[6]+1.060660171779821*facDiffMu[2]*fin[3]*vmap[3]; 
  out[15] += 1.060660171779821*facDiffMu[0]*vmap[3]*fin[11]+1.060660171779821*facDiffMu[1]*vmap[3]*fin[7]+1.060660171779821*facDiffMu[2]*vmap[3]*fin[6]+1.060660171779821*facDiffMu[3]*fin[3]*vmap[3]; 
  out[16] += 3.354101966249685*facDiffVpar[3]*fin[5]+3.354101966249685*facDiffVpar[2]*fin[2]+3.354101966249685*facDiffVpar[1]*fin[1]+3.354101966249685*facDiffVpar[0]*fin[0]; 
  out[17] += 3.354101966249684*facDiffVpar[2]*fin[5]+3.354101966249684*fin[2]*facDiffVpar[3]+3.354101966249684*facDiffVpar[0]*fin[1]+3.354101966249684*fin[0]*facDiffVpar[1]; 
  out[18] += 3.354101966249684*facDiffVpar[1]*fin[5]+3.354101966249684*fin[1]*facDiffVpar[3]+3.354101966249684*facDiffVpar[0]*fin[2]+3.354101966249684*fin[0]*facDiffVpar[2]; 
  out[19] += 1.060660171779821*facDiffMu[3]*vmap[3]*fin[20]+1.060660171779821*facDiffMu[2]*vmap[3]*fin[18]+1.060660171779821*facDiffMu[1]*vmap[3]*fin[17]+1.060660171779821*facDiffMu[0]*vmap[3]*fin[16]+3.354101966249684*facDiffVpar[3]*fin[12]+3.354101966249684*facDiffVpar[2]*fin[9]+3.354101966249684*facDiffVpar[1]*fin[8]+3.354101966249684*facDiffVpar[0]*fin[4]; 
  out[20] += 3.354101966249685*facDiffVpar[0]*fin[5]+3.354101966249685*fin[0]*facDiffVpar[3]+3.354101966249685*facDiffVpar[1]*fin[2]+3.354101966249685*fin[1]*facDiffVpar[2]; 
  out[21] += 1.060660171779821*facDiffMu[2]*vmap[3]*fin[20]+1.060660171779821*facDiffMu[3]*vmap[3]*fin[18]+1.060660171779821*facDiffMu[0]*vmap[3]*fin[17]+1.060660171779821*facDiffMu[1]*vmap[3]*fin[16]+3.354101966249685*facDiffVpar[2]*fin[12]+3.354101966249685*facDiffVpar[3]*fin[9]+3.354101966249685*facDiffVpar[0]*fin[8]+3.354101966249685*facDiffVpar[1]*fin[4]; 
  out[22] += 1.060660171779821*facDiffMu[1]*vmap[3]*fin[20]+1.060660171779821*facDiffMu[0]*vmap[3]*fin[18]+1.060660171779821*facDiffMu[3]*vmap[3]*fin[17]+1.060660171779821*facDiffMu[2]*vmap[3]*fin[16]+3.354101966249685*facDiffVpar[1]*fin[12]+3.354101966249685*facDiffVpar[0]*fin[9]+3.354101966249685*facDiffVpar[3]*fin[8]+3.354101966249685*facDiffVpar[2]*fin[4]; 
  out[23] += 1.060660171779821*facDiffMu[0]*vmap[3]*fin[20]+1.060660171779821*facDiffMu[1]*vmap[3]*fin[18]+1.060660171779821*facDiffMu[2]*vmap[3]*fin[17]+1.060660171779821*facDiffMu[3]*vmap[3]*fin[16]+3.354101966249684*facDiffVpar[0]*fin[12]+3.354101966249684*facDiffVpar[1]*fin[9]+3.354101966249684*facDiffVpar[2]*fin[8]+3.354101966249684*facDiffVpar[3]*fin[4]; 

  return fabs(1.414213562373095*facDiffMu[0]*vmap[2]+4.5*facDiffVpar[0]); 

} 
