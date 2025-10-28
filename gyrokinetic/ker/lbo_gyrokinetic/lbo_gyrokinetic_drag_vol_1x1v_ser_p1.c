#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x1v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // dxv[2]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // f: input distribution function.
  // out: incremented output 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2[1]; 
  rdv2[0] = 2.0/dxv[1]; 

  double alphaDrag[6] = {0.}; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = (rdv2[0]*(1.414213562373095*nuUSum[0]-1.0*nuSum[0]*vmap[0]))/vmap_prime[0]; 
  alphaDrag[1] = (rdv2[0]*(1.414213562373095*nuUSum[1]-1.0*vmap[0]*nuSum[1]))/vmap_prime[0]; 
  alphaDrag[2] = -(1.0*nuSum[0]*rdv2[0]*vmap[1])/vmap_prime[0]; 
  alphaDrag[3] = -(1.0*rdv2[0]*nuSum[1]*vmap[1])/vmap_prime[0]; 

  out[2] += 0.8660254037844386*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[4] += 1.732050807568877*(alphaDrag[3]*f[5]+alphaDrag[2]*f[4])+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[5] += 1.732050807568877*(alphaDrag[2]*f[5]+alphaDrag[3]*f[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return fabs(1.25*alphaDrag[0]); 

} 
