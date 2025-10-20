#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_diff_vol_3x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f:             Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)].
  // out:           Incremented output distribution functions. 
  const double rdvparSq4 = 4.0/(dxv[3]*dxv[3]); 

  const double *nuVtSqSum = &nuPrimMomsSum[8];

  return fabs(3.181980515339463*nuVtSqSum[0]*rdvparSq4); 

} 
