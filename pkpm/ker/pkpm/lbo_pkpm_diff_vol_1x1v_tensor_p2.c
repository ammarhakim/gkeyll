#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_diff_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f:             Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)].
  // out:           Incremented output distribution functions. 
  const double rdvparSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *nuVtSqSum = &nuPrimMomsSum[3];

  return fabs((6.363961030678928*nuVtSqSum[0]-7.115124735378852*nuVtSqSum[2])*rdvparSq4); 

} 
