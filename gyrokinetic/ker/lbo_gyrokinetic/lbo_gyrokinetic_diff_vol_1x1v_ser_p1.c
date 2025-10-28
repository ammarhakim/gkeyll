#include <gkyl_lbo_gyrokinetic_kernels.h> 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x1v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // dxv[2]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fin: input distribution function.
  // out: incremented output 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvVmapPrimeSq4[1]; 
  rdvVmapPrimeSq4[0] = 4.0/(dxv[1]*dxv[1]*vmap_prime[0]*vmap_prime[0]); 

  // Expand nuVtSqSum/vpar'^2 in conf basis.
  double facDiffVpar[2] = {0.};
  facDiffVpar[0] = nuVtSqSum[0]*rdvVmapPrimeSq4[0]; 
  facDiffVpar[1] = rdvVmapPrimeSq4[0]*nuVtSqSum[1]; 

  out[4] += 4.743416490252569*facDiffVpar[1]*fin[1]+4.743416490252569*facDiffVpar[0]*fin[0]; 
  out[5] += 4.743416490252569*facDiffVpar[0]*fin[1]+4.743416490252569*fin[0]*facDiffVpar[1]; 

  return fabs(6.363961030678928*facDiffVpar[0]); 

} 
