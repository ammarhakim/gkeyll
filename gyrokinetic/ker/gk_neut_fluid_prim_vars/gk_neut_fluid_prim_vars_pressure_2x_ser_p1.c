#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_pressure_2x_ser_p1(double gas_gamma, const double *moms, const double *udrift, 
    double* GKYL_RESTRICT out) 
{ 
  // gas_gamma: Adiabatic index.
  // moms: Moments (rho, rho ux, rho uy, rho uz, totalE).
  // udrift: Input volume expansion of flow velocity: [ux, uy, uz]. 
  // out: Output volume expansion of pressure.

  const double *rho    = &moms[0]; 
  const double *rhoux  = &moms[4]; 
  const double *rhouy  = &moms[8]; 
  const double *rhouz  = &moms[12]; 
  const double *energy = &moms[16]; 

  const double *ux = &udrift[0]; 
  const double *uy = &udrift[4]; 
  const double *uz = &udrift[8]; 

  double rhoux2[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhoux, ux, rhoux2); 
 
  double rhouy2[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhouy, uy, rhouy2); 
 
  double rhouz2[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhouz, uz, rhouz2); 
 
  out[0] = (gas_gamma - 1.0)*(energy[0] - 0.5*(rhoux2[0] + rhouy2[0] + rhouz2[0])); 
  out[1] = (gas_gamma - 1.0)*(energy[1] - 0.5*(rhoux2[1] + rhouy2[1] + rhouz2[1])); 
  out[2] = (gas_gamma - 1.0)*(energy[2] - 0.5*(rhoux2[2] + rhouy2[2] + rhouz2[2])); 
  out[3] = (gas_gamma - 1.0)*(energy[3] - 0.5*(rhoux2[3] + rhouy2[3] + rhouz2[3])); 
} 
