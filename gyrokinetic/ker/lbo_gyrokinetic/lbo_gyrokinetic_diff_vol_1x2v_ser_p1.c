#include <gkyl_lbo_gyrokinetic_kernels.h> 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // dxv[3]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fin: input distribution function.
  // out: incremented output 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvVmapPrimeSq4[2]; 
  rdvVmapPrimeSq4[0] = 4.0/(dxv[1]*dxv[1]*vmap_prime[0]*vmap_prime[0]); 
  rdvVmapPrimeSq4[1] = 4.0/(dxv[2]*dxv[2]*vmap_prime[1]*vmap_prime[1]); 

  // Expand nuVtSqSum/vpar'^2 in conf basis.
  double facDiffVpar[2] = {0.};
  facDiffVpar[0] = nuVtSqSum[0]*rdvVmapPrimeSq4[0]; 
  facDiffVpar[1] = rdvVmapPrimeSq4[0]*nuVtSqSum[1]; 

  // Expand 2*m*nuVtSqSum/bmag/mu'^2 in conf basis.
  double facDiffMu[2] = {0.};
  facDiffMu[0] = 1.414213562373095*(bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*rdvVmapPrimeSq4[1]*m_; 
  facDiffMu[1] = 1.414213562373095*(bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*rdvVmapPrimeSq4[1]*m_; 

  out[3] += 1.5*facDiffMu[1]*fin[1]*vmap[3]+1.5*facDiffMu[0]*fin[0]*vmap[3]; 
  out[5] += 1.5*facDiffMu[0]*fin[1]*vmap[3]+1.5*fin[0]*facDiffMu[1]*vmap[3]; 
  out[6] += 1.5*facDiffMu[1]*vmap[3]*fin[4]+1.5*facDiffMu[0]*fin[2]*vmap[3]; 
  out[7] += 1.5*facDiffMu[0]*vmap[3]*fin[4]+1.5*facDiffMu[1]*fin[2]*vmap[3]; 
  out[8] += 4.743416490252569*facDiffVpar[1]*fin[1]+4.743416490252569*facDiffVpar[0]*fin[0]; 
  out[9] += 4.743416490252569*facDiffVpar[0]*fin[1]+4.743416490252569*fin[0]*facDiffVpar[1]; 
  out[10] += 1.5*facDiffMu[1]*vmap[3]*fin[9]+1.5*facDiffMu[0]*vmap[3]*fin[8]+4.743416490252569*facDiffVpar[1]*fin[5]+4.743416490252569*facDiffVpar[0]*fin[3]; 
  out[11] += 1.5*facDiffMu[0]*vmap[3]*fin[9]+1.5*facDiffMu[1]*vmap[3]*fin[8]+4.743416490252569*facDiffVpar[0]*fin[5]+4.743416490252569*facDiffVpar[1]*fin[3]; 

  return fabs(2.0*facDiffMu[0]*vmap[2]+6.363961030678928*facDiffVpar[0]); 

} 
