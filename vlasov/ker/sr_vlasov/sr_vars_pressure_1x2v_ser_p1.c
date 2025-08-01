#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_pressure_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
{ 
  // w:           Cell-center coordinates of velocity grid.
  // dxv:         Cell spacing of velocity grid.
  // vmap:        Momentum-space nonuniform mapping (unused in uniform grid simulations).
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv:   Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
  // u_i:         Spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:      Squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:      Bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq:   Squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
  // f:           Input distribution function.
  // sr_pressure: Output relativistic pressure.
  const double volFact = dxv[1]*dxv[2]/4; 
 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double *V_1 = &u_i[2]; 
  const double *V_1_sq = &u_i_sq[2]; 
 
  double temp[8] = {0.0}; 
  double temp_sq[8] = {0.0}; 
  double p_fac[8] = {0.0}; 
  temp[0] = 2.0*V_1[0]*wx2+2.0*V_0[0]*wx1; 
  temp[1] = 2.0*V_1[1]*wx2+2.0*V_0[1]*wx1; 
  temp[2] = 0.5773502691896258*V_0[0]*dv1; 
  temp[3] = 0.5773502691896258*V_1[0]*dv2; 
  temp[4] = 0.5773502691896258*V_0[1]*dv1; 
  temp[5] = 0.5773502691896258*V_1[1]*dv2; 

  temp_sq[0] = 2.0*V_1_sq[0]*wx2_sq+2.8284271247461907*V_0[1]*V_1[1]*wx1*wx2+2.8284271247461907*V_0[0]*V_1[0]*wx1*wx2+2.0*V_0_sq[0]*wx1_sq+0.16666666666666666*V_1_sq[0]*dv2_sq+0.16666666666666666*V_0_sq[0]*dv1_sq; 
  temp_sq[1] = 2.0*V_1_sq[1]*wx2_sq+2.8284271247461907*V_0[0]*V_1[1]*wx1*wx2+2.8284271247461907*V_1[0]*V_0[1]*wx1*wx2+2.0*V_0_sq[1]*wx1_sq+0.16666666666666666*V_1_sq[1]*dv2_sq+0.16666666666666666*V_0_sq[1]*dv1_sq; 
  temp_sq[2] = 0.8164965809277261*V_0[1]*V_1[1]*dv1*wx2+0.8164965809277261*V_0[0]*V_1[0]*dv1*wx2+1.1547005383792517*V_0_sq[0]*dv1*wx1; 
  temp_sq[3] = 1.1547005383792517*V_1_sq[0]*dv2*wx2+0.8164965809277261*V_0[1]*V_1[1]*dv2*wx1+0.8164965809277261*V_0[0]*V_1[0]*dv2*wx1; 
  temp_sq[4] = 0.8164965809277261*V_0[0]*V_1[1]*dv1*wx2+0.8164965809277261*V_1[0]*V_0[1]*dv1*wx2+1.1547005383792517*V_0_sq[1]*dv1*wx1; 
  temp_sq[5] = 1.1547005383792517*V_1_sq[1]*dv2*wx2+0.8164965809277261*V_0[0]*V_1[1]*dv2*wx1+0.8164965809277261*V_1[0]*V_0[1]*dv2*wx1; 
  temp_sq[6] = 0.2357022603955158*V_0[1]*V_1[1]*dv1*dv2+0.2357022603955158*V_0[0]*V_1[0]*dv1*dv2; 
  temp_sq[7] = 0.2357022603955158*V_0[0]*V_1[1]*dv1*dv2+0.2357022603955158*V_1[0]*V_0[1]*dv1*dv2; 

  p_fac[0] = 0.5*gamma_inv[3]*temp_sq[6]+0.5*gamma_inv[2]*temp_sq[3]+0.5*gamma_inv[1]*temp_sq[2]-1.4142135623730951*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.5*gamma_inv[0]*temp_sq[0]-1.4142135623730951*GammaV[0]*temp[0]-1.4142135623730951*gamma_inv[0]; 
  p_fac[1] = 0.5*gamma_inv[3]*temp_sq[7]+0.5*gamma_inv[2]*temp_sq[5]+0.5*gamma_inv[1]*temp_sq[4]+0.5*gamma_inv[0]*temp_sq[1]-1.4142135623730951*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.4142135623730951*temp[0]*GammaV[1]; 
  p_fac[2] = 0.5*gamma_inv[2]*temp_sq[6]-1.4142135623730951*GammaV[1]*temp[4]+0.5*gamma_inv[3]*temp_sq[3]+0.5*gamma_inv[0]*temp_sq[2]-1.4142135623730951*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.5*temp_sq[0]*gamma_inv[1]-1.4142135623730951*gamma_inv[1]; 
  p_fac[3] = 0.5*gamma_inv[1]*temp_sq[6]-1.4142135623730951*GammaV[1]*temp[5]+0.5*gamma_inv[0]*temp_sq[3]-1.4142135623730951*GammaV[0]*temp[3]+0.5*temp_sq[2]*gamma_inv[3]+GammaV_sq[0]*gamma[2]+0.5*temp_sq[0]*gamma_inv[2]-1.4142135623730951*gamma_inv[2]; 
  p_fac[4] = 0.5*gamma_inv[2]*temp_sq[7]+0.5*gamma_inv[3]*temp_sq[5]+0.5*gamma_inv[0]*temp_sq[4]-1.4142135623730951*GammaV[0]*temp[4]-1.4142135623730951*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.5*gamma_inv[1]*temp_sq[1]; 
  p_fac[5] = 0.5*gamma_inv[1]*temp_sq[7]+0.5*gamma_inv[0]*temp_sq[5]-1.4142135623730951*GammaV[0]*temp[5]+0.5*gamma_inv[3]*temp_sq[4]-1.4142135623730951*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[2]+0.5*temp_sq[1]*gamma_inv[2]; 
  p_fac[6] = 0.5*gamma_inv[0]*temp_sq[6]+GammaV_sq[0]*gamma[3]+0.5*gamma_inv[1]*temp_sq[3]+0.5*temp_sq[0]*gamma_inv[3]-1.4142135623730951*gamma_inv[3]+0.5*gamma_inv[2]*temp_sq[2]; 
  p_fac[7] = 0.5*gamma_inv[0]*temp_sq[7]+0.5*gamma_inv[1]*temp_sq[5]+0.5*gamma_inv[2]*temp_sq[4]+GammaV_sq[1]*gamma[3]+0.5*temp_sq[1]*gamma_inv[3]; 

  sr_pressure[0] += (0.3535533905932737*f[7]*p_fac[7]+0.3535533905932737*f[6]*p_fac[6]+0.3535533905932737*f[5]*p_fac[5]+0.3535533905932737*f[4]*p_fac[4]+0.3535533905932737*f[3]*p_fac[3]+0.3535533905932737*f[2]*p_fac[2]+0.3535533905932737*f[1]*p_fac[1]+0.3535533905932737*f[0]*p_fac[0])*volFact; 
  sr_pressure[1] += (0.3535533905932737*f[6]*p_fac[7]+0.3535533905932737*p_fac[6]*f[7]+0.3535533905932737*f[3]*p_fac[5]+0.3535533905932737*p_fac[3]*f[5]+0.3535533905932737*f[2]*p_fac[4]+0.3535533905932737*p_fac[2]*f[4]+0.3535533905932737*f[0]*p_fac[1]+0.3535533905932737*p_fac[0]*f[1])*volFact; 
} 
