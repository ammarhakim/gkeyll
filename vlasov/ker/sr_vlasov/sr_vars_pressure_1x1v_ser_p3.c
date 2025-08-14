#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_pressure_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
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
  const double volFact = dxv[1]/2; 
 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  double temp[12] = {0.0}; 
  double temp_sq[12] = {0.0}; 
  double p_fac[12] = {0.0}; 
  temp[0] = 1.4142135623730951*V_0[0]*wx1; 
  temp[1] = 1.4142135623730951*V_0[1]*wx1; 
  temp[2] = 0.408248290463863*V_0[0]*dv1; 
  temp[3] = 0.408248290463863*V_0[1]*dv1; 
  temp[4] = 1.4142135623730951*V_0[2]*wx1; 
  temp[6] = 0.40824829046386296*V_0[2]*dv1; 
  temp[8] = 1.4142135623730951*V_0[3]*wx1; 
  temp[10] = 0.4082482904638629*V_0[3]*dv1; 

  temp_sq[0] = 1.4142135623730951*V_0_sq[0]*wx1_sq+0.11785113019775789*V_0_sq[0]*dv1_sq; 
  temp_sq[1] = 1.4142135623730951*V_0_sq[1]*wx1_sq+0.11785113019775789*V_0_sq[1]*dv1_sq; 
  temp_sq[2] = 0.8164965809277261*V_0_sq[0]*dv1*wx1; 
  temp_sq[3] = 0.8164965809277261*V_0_sq[1]*dv1*wx1; 
  temp_sq[4] = 1.4142135623730951*V_0_sq[2]*wx1_sq+0.11785113019775789*V_0_sq[2]*dv1_sq; 
  temp_sq[5] = 0.10540925533894596*V_0_sq[0]*dv1_sq; 
  temp_sq[6] = 0.816496580927726*V_0_sq[2]*dv1*wx1; 
  temp_sq[7] = 0.105409255338946*V_0_sq[1]*dv1_sq; 
  temp_sq[8] = 1.4142135623730951*V_0_sq[3]*wx1_sq+0.11785113019775789*V_0_sq[3]*dv1_sq; 
  temp_sq[10] = 0.8164965809277258*V_0_sq[3]*dv1*wx1; 

  p_fac[0] = -(1.4142135623730951*GammaV[3]*temp[8])+0.7071067811865475*gamma_inv[2]*temp_sq[5]-1.4142135623730951*GammaV[2]*temp[4]+0.7071067811865475*gamma_inv[1]*temp_sq[2]-1.4142135623730951*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.7071067811865475*gamma_inv[0]*temp_sq[0]-1.4142135623730951*GammaV[0]*temp[0]-1.4142135623730951*gamma_inv[0]; 
  p_fac[1] = -(1.2421180068162374*GammaV[2]*temp[8])+0.7071067811865475*gamma_inv[2]*temp_sq[7]-1.2421180068162374*GammaV[3]*temp[4]-1.264911064067352*GammaV[1]*temp[4]+0.7071067811865475*gamma_inv[1]*temp_sq[3]-1.264911064067352*temp[1]*GammaV[2]+0.7071067811865475*gamma_inv[0]*temp_sq[1]-1.4142135623730951*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.4142135623730951*temp[0]*GammaV[1]; 
  p_fac[2] = -(1.414213562373095*GammaV[3]*temp[10])-1.4142135623730951*GammaV[2]*temp[6]+0.6210590034081186*gamma_inv[3]*temp_sq[5]+0.6324555320336759*gamma_inv[1]*temp_sq[5]-1.4142135623730951*GammaV[1]*temp[3]+0.6324555320336759*gamma_inv[2]*temp_sq[2]+0.7071067811865475*gamma_inv[0]*temp_sq[2]-1.4142135623730951*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.7071067811865475*temp_sq[0]*gamma_inv[1]-1.4142135623730951*gamma_inv[1]; 
  p_fac[3] = -(1.2421180068162376*GammaV[2]*temp[10])+0.6210590034081187*gamma_inv[3]*temp_sq[7]+0.632455532033676*gamma_inv[1]*temp_sq[7]-1.2421180068162376*GammaV[3]*temp[6]-1.264911064067352*GammaV[1]*temp[6]+0.6324555320336759*gamma_inv[2]*temp_sq[3]+0.7071067811865475*gamma_inv[0]*temp_sq[3]-1.264911064067352*GammaV[2]*temp[3]-1.4142135623730951*GammaV[0]*temp[3]-1.4142135623730951*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.7071067811865475*gamma_inv[1]*temp_sq[1]; 
  p_fac[4] = -(0.8432740427115681*GammaV[3]*temp[8])-1.2421180068162374*GammaV[1]*temp[8]+0.7071067811865475*gamma_inv[1]*temp_sq[6]+0.7071067811865475*gamma_inv[0]*temp_sq[4]-0.9035079029052515*GammaV[2]*temp[4]-1.4142135623730951*GammaV[0]*temp[4]-1.2421180068162374*temp[1]*GammaV[3]+gamma[0]*GammaV_sq[2]-1.4142135623730951*temp[0]*GammaV[2]-1.264911064067352*GammaV[1]*temp[1]; 
  p_fac[5] = 0.45175395145262565*gamma_inv[2]*temp_sq[5]+0.7071067811865475*gamma_inv[0]*temp_sq[5]+0.6210590034081186*temp_sq[2]*gamma_inv[3]+GammaV_sq[0]*gamma[2]+0.6324555320336759*gamma_inv[1]*temp_sq[2]+0.7071067811865475*temp_sq[0]*gamma_inv[2]-1.4142135623730951*gamma_inv[2]; 
  p_fac[6] = -(0.8432740427115679*GammaV[3]*temp[10])-1.2421180068162374*GammaV[1]*temp[10]+0.6324555320336759*gamma_inv[2]*temp_sq[6]+0.7071067811865475*gamma_inv[0]*temp_sq[6]-0.9035079029052515*GammaV[2]*temp[6]-1.4142135623730951*GammaV[0]*temp[6]+0.7071067811865475*gamma_inv[1]*temp_sq[4]-1.2421180068162376*GammaV[3]*temp[3]-1.264911064067352*GammaV[1]*temp[3]-1.4142135623730951*GammaV[2]*temp[2]+1.0000000000000002*gamma[1]*GammaV_sq[2]; 
  p_fac[7] = 0.45175395145262565*gamma_inv[2]*temp_sq[7]+0.7071067811865475*gamma_inv[0]*temp_sq[7]+0.6210590034081187*gamma_inv[3]*temp_sq[3]+0.632455532033676*gamma_inv[1]*temp_sq[3]+1.0000000000000002*GammaV_sq[1]*gamma[2]+0.7071067811865475*temp_sq[1]*gamma_inv[2]; 
  p_fac[8] = 0.7071067811865474*gamma_inv[1]*temp_sq[10]+0.7071067811865475*gamma_inv[0]*temp_sq[8]-0.8432740427115681*GammaV[2]*temp[8]-1.4142135623730951*GammaV[0]*temp[8]-0.8432740427115681*GammaV[3]*temp[4]-1.2421180068162374*GammaV[1]*temp[4]+gamma[0]*GammaV_sq[3]-1.4142135623730951*temp[0]*GammaV[3]-1.2421180068162374*temp[1]*GammaV[2]; 
  p_fac[9] = 0.421637021355784*gamma_inv[3]*temp_sq[5]+0.6210590034081186*gamma_inv[1]*temp_sq[5]+GammaV_sq[0]*gamma[3]+0.7071067811865475*temp_sq[0]*gamma_inv[3]-1.4142135623730951*gamma_inv[3]+0.6210590034081186*gamma_inv[2]*temp_sq[2]; 
  p_fac[10] = 0.6324555320336759*gamma_inv[2]*temp_sq[10]+0.7071067811865475*gamma_inv[0]*temp_sq[10]-0.8432740427115681*GammaV[2]*temp[10]-1.4142135623730951*GammaV[0]*temp[10]+0.7071067811865474*gamma_inv[1]*temp_sq[8]-0.8432740427115679*GammaV[3]*temp[6]-1.2421180068162374*GammaV[1]*temp[6]-1.2421180068162376*GammaV[2]*temp[3]+1.0*gamma[1]*GammaV_sq[3]-1.414213562373095*temp[2]*GammaV[3]; 
  p_fac[11] = 0.4216370213557839*gamma_inv[3]*temp_sq[7]+0.6210590034081187*gamma_inv[1]*temp_sq[7]+1.0*GammaV_sq[1]*gamma[3]+0.6210590034081187*gamma_inv[2]*temp_sq[3]+0.7071067811865474*temp_sq[1]*gamma_inv[3]; 

  sr_pressure[0] += (0.7071067811865475*f[11]*p_fac[11]+0.7071067811865475*f[10]*p_fac[10]+0.7071067811865475*f[9]*p_fac[9]+0.7071067811865475*f[8]*p_fac[8]+0.7071067811865475*f[7]*p_fac[7]+0.7071067811865475*f[6]*p_fac[6]+0.7071067811865475*f[5]*p_fac[5]+0.7071067811865475*f[4]*p_fac[4]+0.7071067811865475*f[3]*p_fac[3]+0.7071067811865475*f[2]*p_fac[2]+0.7071067811865475*f[1]*p_fac[1]+0.7071067811865475*f[0]*p_fac[0])*volFact; 
  sr_pressure[1] += (0.7071067811865474*f[9]*p_fac[11]+0.7071067811865474*p_fac[9]*f[11]+0.6210590034081187*f[6]*p_fac[10]+0.6210590034081187*p_fac[6]*f[10]+0.6210590034081186*f[4]*p_fac[8]+0.6210590034081186*p_fac[4]*f[8]+0.7071067811865475*f[5]*p_fac[7]+0.7071067811865475*p_fac[5]*f[7]+0.632455532033676*f[3]*p_fac[6]+0.632455532033676*p_fac[3]*f[6]+0.6324555320336759*f[1]*p_fac[4]+0.6324555320336759*p_fac[1]*f[4]+0.7071067811865475*f[2]*p_fac[3]+0.7071067811865475*p_fac[2]*f[3]+0.7071067811865475*f[0]*p_fac[1]+0.7071067811865475*p_fac[0]*f[1])*volFact; 
  sr_pressure[2] += (0.6324555320336759*f[11]*p_fac[11]+0.421637021355784*f[10]*p_fac[10]+0.6210590034081187*f[3]*p_fac[10]+0.6210590034081187*p_fac[3]*f[10]+0.421637021355784*f[8]*p_fac[8]+0.6210590034081186*f[1]*p_fac[8]+0.6210590034081186*p_fac[1]*f[8]+0.6324555320336759*f[7]*p_fac[7]+0.45175395145262565*f[6]*p_fac[6]+0.7071067811865475*f[2]*p_fac[6]+0.7071067811865475*p_fac[2]*f[6]+0.45175395145262565*f[4]*p_fac[4]+0.7071067811865475*f[0]*p_fac[4]+0.7071067811865475*p_fac[0]*f[4]+0.6324555320336759*f[3]*p_fac[3]+0.6324555320336759*f[1]*p_fac[1])*volFact; 
  sr_pressure[3] += (0.4216370213557839*f[6]*p_fac[10]+0.7071067811865474*f[2]*p_fac[10]+0.4216370213557839*p_fac[6]*f[10]+0.7071067811865474*p_fac[2]*f[10]+0.421637021355784*f[4]*p_fac[8]+0.7071067811865475*f[0]*p_fac[8]+0.421637021355784*p_fac[4]*f[8]+0.7071067811865475*p_fac[0]*f[8]+0.6210590034081187*f[3]*p_fac[6]+0.6210590034081187*p_fac[3]*f[6]+0.6210590034081186*f[1]*p_fac[4]+0.6210590034081186*p_fac[1]*f[4])*volFact; 
} 
