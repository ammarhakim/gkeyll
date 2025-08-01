#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_pressure_2x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
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
  const double volFact = dxv[2]/2; 
 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  double temp[8] = {0.0}; 
  double temp_sq[8] = {0.0}; 
  double p_fac[8] = {0.0}; 
  temp[0] = 1.4142135623730951*V_0[0]*wx1; 
  temp[1] = 1.4142135623730951*V_0[1]*wx1; 
  temp[2] = 1.4142135623730951*V_0[2]*wx1; 
  temp[3] = 0.408248290463863*V_0[0]*dv1; 
  temp[4] = 1.4142135623730951*V_0[3]*wx1; 
  temp[5] = 0.408248290463863*V_0[1]*dv1; 
  temp[6] = 0.408248290463863*V_0[2]*dv1; 
  temp[7] = 0.408248290463863*V_0[3]*dv1; 

  temp_sq[0] = 1.4142135623730951*V_0_sq[0]*wx1_sq+0.11785113019775789*V_0_sq[0]*dv1_sq; 
  temp_sq[1] = 1.4142135623730951*V_0_sq[1]*wx1_sq+0.11785113019775789*V_0_sq[1]*dv1_sq; 
  temp_sq[2] = 1.4142135623730951*V_0_sq[2]*wx1_sq+0.11785113019775789*V_0_sq[2]*dv1_sq; 
  temp_sq[3] = 0.8164965809277261*V_0_sq[0]*dv1*wx1; 
  temp_sq[4] = 1.4142135623730951*V_0_sq[3]*wx1_sq+0.11785113019775789*V_0_sq[3]*dv1_sq; 
  temp_sq[5] = 0.8164965809277261*V_0_sq[1]*dv1*wx1; 
  temp_sq[6] = 0.8164965809277261*V_0_sq[2]*dv1*wx1; 
  temp_sq[7] = 0.8164965809277261*V_0_sq[3]*dv1*wx1; 

  p_fac[0] = -(1.0*GammaV[3]*temp[4])+0.7071067811865475*gamma_inv[1]*temp_sq[3]-1.0*GammaV[2]*temp[2]-1.0*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.7071067811865475*gamma_inv[0]*temp_sq[0]-1.0*GammaV[0]*temp[0]-2.0*gamma_inv[0]; 
  p_fac[1] = 0.7071067811865475*gamma_inv[1]*temp_sq[5]-1.0*GammaV[2]*temp[4]-1.0*temp[2]*GammaV[3]+0.7071067811865475*gamma_inv[0]*temp_sq[1]-1.0*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.0*temp[0]*GammaV[1]; 
  p_fac[2] = 0.7071067811865475*gamma_inv[1]*temp_sq[6]-1.0*GammaV[1]*temp[4]-1.0*temp[1]*GammaV[3]+0.7071067811865475*gamma_inv[0]*temp_sq[2]-1.0*GammaV[0]*temp[2]+gamma[0]*GammaV_sq[2]-1.0*temp[0]*GammaV[2]; 
  p_fac[3] = -(1.0*GammaV[3]*temp[7])-1.0*GammaV[2]*temp[6]-1.0*GammaV[1]*temp[5]+0.7071067811865475*gamma_inv[0]*temp_sq[3]-1.0*GammaV[0]*temp[3]+GammaV_sq[0]*gamma[1]+0.7071067811865475*temp_sq[0]*gamma_inv[1]-2.0*gamma_inv[1]; 
  p_fac[4] = 0.7071067811865475*gamma_inv[1]*temp_sq[7]+0.7071067811865475*gamma_inv[0]*temp_sq[4]-1.0*GammaV[0]*temp[4]+gamma[0]*GammaV_sq[3]-1.0*temp[0]*GammaV[3]-1.0*GammaV[1]*temp[2]-1.0*temp[1]*GammaV[2]; 
  p_fac[5] = -(1.0*GammaV[2]*temp[7])-1.0*GammaV[3]*temp[6]+0.7071067811865475*gamma_inv[0]*temp_sq[5]-1.0*GammaV[0]*temp[5]-1.0*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[1]+0.7071067811865475*gamma_inv[1]*temp_sq[1]; 
  p_fac[6] = -(1.0*GammaV[1]*temp[7])+0.7071067811865475*gamma_inv[0]*temp_sq[6]-1.0*GammaV[0]*temp[6]-1.0*GammaV[3]*temp[5]-1.0*GammaV[2]*temp[3]+0.7071067811865475*gamma_inv[1]*temp_sq[2]+gamma[1]*GammaV_sq[2]; 
  p_fac[7] = 0.7071067811865475*gamma_inv[0]*temp_sq[7]-1.0*GammaV[0]*temp[7]-1.0*GammaV[1]*temp[6]-1.0*GammaV[2]*temp[5]+0.7071067811865475*gamma_inv[1]*temp_sq[4]-1.0*GammaV[3]*temp[3]+gamma[1]*GammaV_sq[3]; 

  sr_pressure[0] += (0.5*f[7]*p_fac[7]+0.5*f[6]*p_fac[6]+0.5*f[5]*p_fac[5]+0.5*f[4]*p_fac[4]+0.5*f[3]*p_fac[3]+0.5*f[2]*p_fac[2]+0.5*f[1]*p_fac[1]+0.5*f[0]*p_fac[0])*volFact; 
  sr_pressure[1] += (0.5*f[6]*p_fac[7]+0.5*p_fac[6]*f[7]+0.5*f[3]*p_fac[5]+0.5*p_fac[3]*f[5]+0.5*f[2]*p_fac[4]+0.5*p_fac[2]*f[4]+0.5*f[0]*p_fac[1]+0.5*p_fac[0]*f[1])*volFact; 
  sr_pressure[2] += (0.5*f[5]*p_fac[7]+0.5*p_fac[5]*f[7]+0.5*f[3]*p_fac[6]+0.5*p_fac[3]*f[6]+0.5*f[1]*p_fac[4]+0.5*p_fac[1]*f[4]+0.5*f[0]*p_fac[2]+0.5*p_fac[0]*f[2])*volFact; 
  sr_pressure[3] += (0.5*f[3]*p_fac[7]+0.5*p_fac[3]*f[7]+0.5*f[5]*p_fac[6]+0.5*p_fac[5]*f[6]+0.5*f[0]*p_fac[4]+0.5*p_fac[0]*f[4]+0.5*f[1]*p_fac[2]+0.5*p_fac[1]*f[2])*volFact; 
} 
