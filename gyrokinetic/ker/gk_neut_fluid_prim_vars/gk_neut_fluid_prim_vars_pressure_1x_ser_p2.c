#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_pressure_1x_ser_p2(double gas_gamma, const double *moms, const double *udrift, 
    double* GKYL_RESTRICT out) 
{ 
  // gas_gamma: Adiabatic index.
  // moms: Moments (rho, rho ux, rho uy, rho uz, totalE).
  // udrift: Input volume expansion of flow velocity: [ux, uy, uz]. 
  // out: Output volume expansion of pressure.

  const double *rho    = &moms[0]; 
  const double *rhoux  = &moms[3]; 
  const double *rhouy  = &moms[6]; 
  const double *rhouz  = &moms[9]; 
  const double *energy = &moms[12]; 

  const double *ux = &udrift[0]; 
  const double *uy = &udrift[3]; 
  const double *uz = &udrift[6]; 

  double rhoux2[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhoux, ux, rhoux2); 
 
  double rhouy2[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouy, uy, rhouy2); 
 
  double rhouz2[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouz, uz, rhouz2); 
 
  out[0] = (gas_gamma - 1.0)*(energy[0] - 0.5*(rhoux2[0] + rhouy2[0] + rhouz2[0])); 
  out[1] = (gas_gamma - 1.0)*(energy[1] - 0.5*(rhoux2[1] + rhouy2[1] + rhouz2[1])); 
  out[2] = (gas_gamma - 1.0)*(energy[2] - 0.5*(rhoux2[2] + rhouy2[2] + rhouz2[2])); 
} 
