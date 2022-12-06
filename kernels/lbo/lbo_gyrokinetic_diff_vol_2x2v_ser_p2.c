#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_vol_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fin: input distribution function.
  // out: incremented output 

  const double *nuVtSqSum = &nuPrimMomsSum[8];

  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[2]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[1]   = 2.0/dxv[3]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 

  double facDiffVpar[8];
  // Expand nuVtSqSum in phase basis.
  facDiffVpar[0] = nuVtSqSum[0]; 
  facDiffVpar[4] = nuVtSqSum[4]; 
  facDiffVpar[5] = nuVtSqSum[5]; 

  double facDiffMu[48];
  // Expand mu diffusion coefficient in phase basis.
  facDiffMu[0] = rdv2[1]*w[3]*(bmag_inv[4]*nuVtSqSum[4]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0]); 
  facDiffMu[4] = 0.5773502691896258*(bmag_inv[4]*nuVtSqSum[4]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0]); 
  facDiffMu[5] = rdv2[1]*w[3]*(0.8944271909999159*bmag_inv[1]*nuVtSqSum[6]+nuVtSqSum[3]*(0.8944271909999159*bmag_inv[4]+bmag_inv[0])+bmag_inv[1]*nuVtSqSum[2]); 

  return fabs(rdv2[1]*(9.0*facDiffMu[0]-10.06230589874905*(facDiffMu[5]+facDiffMu[4]))*m_+rdvSq4[0]*(4.5*facDiffVpar[0]-5.031152949374527*(facDiffVpar[5]+facDiffVpar[4]))); 

} 
