#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_diff_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr:      Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells 
  // out:           Incremented distribution functions in center cell. 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double *nuVtSqSum = &nuPrimMomsSum[2];

  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[6]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[6]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double incr_F_0[6] = {0.0}; 
  double incr_G_1[6] = {0.0}; 

  double F_0_xx[6] = {0.0}; 
  double G_1_xx[6] = {0.0}; 
  F_0_xx[0] = 0.6708203932499369*F_0r[4]+0.6708203932499369*F_0l[4]-1.341640786499874*F_0c[4]-1.190784930203603*F_0r[2]+1.190784930203603*F_0l[2]+0.9375*F_0r[0]+0.9375*F_0l[0]-1.875*F_0c[0]; 
  F_0_xx[1] = 0.6708203932499369*F_0r[5]+0.6708203932499369*F_0l[5]-1.341640786499874*F_0c[5]-1.190784930203603*F_0r[3]+1.190784930203603*F_0l[3]+0.9375*F_0r[1]+0.9375*F_0l[1]-1.875*F_0c[1]; 
  F_0_xx[2] = 0.7382874503707888*F_0r[4]-0.7382874503707888*F_0l[4]-1.453125*F_0r[2]-1.453125*F_0l[2]-5.34375*F_0c[2]+1.190784930203603*F_0r[0]-1.190784930203603*F_0l[0]; 
  F_0_xx[3] = 0.7382874503707888*F_0r[5]-0.7382874503707888*F_0l[5]-1.453125*F_0r[3]-1.453125*F_0l[3]-5.34375*F_0c[3]+1.190784930203603*F_0r[1]-1.190784930203603*F_0l[1]; 
  F_0_xx[4] = (-0.140625*F_0r[4])-0.140625*F_0l[4]-6.28125*F_0c[4]-0.3025768239224545*F_0r[2]+0.3025768239224545*F_0l[2]+0.4192627457812106*F_0r[0]+0.4192627457812106*F_0l[0]-0.8385254915624212*F_0c[0]; 
  F_0_xx[5] = (-0.140625*F_0r[5])-0.140625*F_0l[5]-6.28125*F_0c[5]-0.3025768239224544*F_0r[3]+0.3025768239224544*F_0l[3]+0.4192627457812105*F_0r[1]+0.4192627457812105*F_0l[1]-0.8385254915624211*F_0c[1]; 
  G_1_xx[0] = 0.6708203932499369*G_1r[4]+0.6708203932499369*G_1l[4]-1.341640786499874*G_1c[4]-1.190784930203603*G_1r[2]+1.190784930203603*G_1l[2]+0.9375*G_1r[0]+0.9375*G_1l[0]-1.875*G_1c[0]; 
  G_1_xx[1] = 0.6708203932499369*G_1r[5]+0.6708203932499369*G_1l[5]-1.341640786499874*G_1c[5]-1.190784930203603*G_1r[3]+1.190784930203603*G_1l[3]+0.9375*G_1r[1]+0.9375*G_1l[1]-1.875*G_1c[1]; 
  G_1_xx[2] = 0.7382874503707888*G_1r[4]-0.7382874503707888*G_1l[4]-1.453125*G_1r[2]-1.453125*G_1l[2]-5.34375*G_1c[2]+1.190784930203603*G_1r[0]-1.190784930203603*G_1l[0]; 
  G_1_xx[3] = 0.7382874503707888*G_1r[5]-0.7382874503707888*G_1l[5]-1.453125*G_1r[3]-1.453125*G_1l[3]-5.34375*G_1c[3]+1.190784930203603*G_1r[1]-1.190784930203603*G_1l[1]; 
  G_1_xx[4] = (-0.140625*G_1r[4])-0.140625*G_1l[4]-6.28125*G_1c[4]-0.3025768239224545*G_1r[2]+0.3025768239224545*G_1l[2]+0.4192627457812106*G_1r[0]+0.4192627457812106*G_1l[0]-0.8385254915624212*G_1c[0]; 
  G_1_xx[5] = (-0.140625*G_1r[5])-0.140625*G_1l[5]-6.28125*G_1c[5]-0.3025768239224544*G_1r[3]+0.3025768239224544*G_1l[3]+0.4192627457812105*G_1r[1]+0.4192627457812105*G_1l[1]-0.8385254915624211*G_1c[1]; 

  incr_F_0[0] = 0.7071067811865475*F_0_xx[1]*nuVtSqSum[1]+0.7071067811865475*F_0_xx[0]*nuVtSqSum[0]; 
  incr_F_0[1] = 0.7071067811865475*F_0_xx[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*F_0_xx[1]; 
  incr_F_0[2] = 0.7071067811865475*nuVtSqSum[1]*F_0_xx[3]+0.7071067811865475*nuVtSqSum[0]*F_0_xx[2]; 
  incr_F_0[3] = 0.7071067811865475*nuVtSqSum[0]*F_0_xx[3]+0.7071067811865475*nuVtSqSum[1]*F_0_xx[2]; 
  incr_F_0[4] = 0.7071067811865475*nuVtSqSum[1]*F_0_xx[5]+0.7071067811865475*nuVtSqSum[0]*F_0_xx[4]; 
  incr_F_0[5] = 0.7071067811865475*nuVtSqSum[0]*F_0_xx[5]+0.7071067811865475*nuVtSqSum[1]*F_0_xx[4]; 
  incr_G_1[0] = 0.7071067811865475*G_1_xx[1]*nuVtSqSum[1]+0.7071067811865475*G_1_xx[0]*nuVtSqSum[0]; 
  incr_G_1[1] = 0.7071067811865475*G_1_xx[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*G_1_xx[1]; 
  incr_G_1[2] = 0.7071067811865475*nuVtSqSum[1]*G_1_xx[3]+0.7071067811865475*nuVtSqSum[0]*G_1_xx[2]; 
  incr_G_1[3] = 0.7071067811865475*nuVtSqSum[0]*G_1_xx[3]+0.7071067811865475*nuVtSqSum[1]*G_1_xx[2]; 
  incr_G_1[4] = 0.7071067811865475*nuVtSqSum[1]*G_1_xx[5]+0.7071067811865475*nuVtSqSum[0]*G_1_xx[4]; 
  incr_G_1[5] = 0.7071067811865475*nuVtSqSum[0]*G_1_xx[5]+0.7071067811865475*nuVtSqSum[1]*G_1_xx[4]; 

  out_F_0[0] += incr_F_0[0]*rdvSq4; 
  out_F_0[1] += incr_F_0[1]*rdvSq4; 
  out_F_0[2] += incr_F_0[2]*rdvSq4; 
  out_F_0[3] += incr_F_0[3]*rdvSq4; 
  out_F_0[4] += incr_F_0[4]*rdvSq4; 
  out_F_0[5] += incr_F_0[5]*rdvSq4; 
  out_G_1[0] += incr_G_1[0]*rdvSq4; 
  out_G_1[1] += incr_G_1[1]*rdvSq4; 
  out_G_1[2] += incr_G_1[2]*rdvSq4; 
  out_G_1[3] += incr_G_1[3]*rdvSq4; 
  out_G_1[4] += incr_G_1[4]*rdvSq4; 
  out_G_1[5] += incr_G_1[5]*rdvSq4; 

  return 0.;

} 
