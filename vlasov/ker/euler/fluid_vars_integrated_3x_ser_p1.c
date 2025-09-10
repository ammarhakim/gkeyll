#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_integrated_3x_ser_p1(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars) 
{ 
  // fluid:          [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // u_i:            [ux, uy, uz], Input flow velocity.
  // p_ij:           Input pressure (scalar for Euler/5 moment, tensor for 10 moment).
  // int_fluid_vars: Output integrated variables.

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[8]; 
  const double *rhouy = &fluid[16]; 
  const double *rhouz = &fluid[24]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[8]; 
  const double *uz = &u_i[16]; 

  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rhou^2, p) 
  // Note: there is a factor of 1/2^cdim to properly handle the basis function normalization when computing integrated moments 
  int_fluid_vars[0] += 0.3535533905932738*rho[0]; 
  int_fluid_vars[1] += 0.3535533905932738*rhoux[0]; 
  int_fluid_vars[2] += 0.3535533905932738*rhouy[0]; 
  int_fluid_vars[3] += 0.3535533905932738*rhouz[0]; 
  int_fluid_vars[4] += 0.125*(rhouz[7]*uz[7]+rhouy[7]*uy[7]+rhoux[7]*ux[7]+rhouz[6]*uz[6]+rhouy[6]*uy[6]+rhoux[6]*ux[6]+rhouz[5]*uz[5]+rhouy[5]*uy[5]+rhoux[5]*ux[5]+rhouz[4]*uz[4]+rhouy[4]*uy[4]+rhoux[4]*ux[4]+rhouz[3]*uz[3]+rhouy[3]*uy[3]+rhoux[3]*ux[3]+rhouz[2]*uz[2]+rhouy[2]*uy[2]+rhoux[2]*ux[2]+rhouz[1]*uz[1]+rhouy[1]*uy[1]+rhoux[1]*ux[1]+rhouz[0]*uz[0]+rhouy[0]*uy[0]+rhoux[0]*ux[0]); 
  int_fluid_vars[5] += 0.3535533905932738*p_ij[0]; 
} 
