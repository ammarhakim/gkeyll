#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_integrated_1x_ser_p1(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars) 
{ 
  // fluid:          [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // u_i:            [ux, uy, uz], Input flow velocity.
  // p_ij:           Input pressure (scalar for Euler/5 moment, tensor for 10 moment).
  // int_fluid_vars: Output integrated variables.

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[2]; 
  const double *rhouy = &fluid[4]; 
  const double *rhouz = &fluid[6]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[2]; 
  const double *uz = &u_i[4]; 

  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rhou^2, p) 
  // Note: there is a factor of 1/2^cdim to properly handle the basis function normalization when computing integrated moments 
  int_fluid_vars[0] += 0.7071067811865476*rho[0]; 
  int_fluid_vars[1] += 0.7071067811865476*rhoux[0]; 
  int_fluid_vars[2] += 0.7071067811865476*rhouy[0]; 
  int_fluid_vars[3] += 0.7071067811865476*rhouz[0]; 
  int_fluid_vars[4] += 0.5*(rhouz[1]*uz[1]+rhouy[1]*uy[1]+rhoux[1]*ux[1]+rhouz[0]*uz[0]+rhouy[0]*uy[0]+rhoux[0]*ux[0]); 
  int_fluid_vars[5] += 0.7071067811865476*p_ij[0]; 
} 
