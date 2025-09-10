#include <gkyl_euler_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void fluid_vars_ke_3x_ser_p1(const double *fluid, const double *u_i, 
    double* GKYL_RESTRICT ke) 
{ 
  // fluid:  [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // u_i:    Input volume expansion of flow velocity: [ux, uy, uz]. 
  // ke:     Output volume expansion of kinetic energy.

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[8]; 
  const double *rhouy = &fluid[16]; 
  const double *rhouz = &fluid[24]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[8]; 
  const double *uz = &u_i[16]; 

  double rhoux2[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhoux, ux, rhoux2); 
 
  double rhouy2[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouy, uy, rhouy2); 
 
  double rhouz2[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouz, uz, rhouz2); 
 
  ke[0] = 0.5*(rhoux2[0] + rhouy2[0] + rhouz2[0]); 
  ke[1] = 0.5*(rhoux2[1] + rhouy2[1] + rhouz2[1]); 
  ke[2] = 0.5*(rhoux2[2] + rhouy2[2] + rhouz2[2]); 
  ke[3] = 0.5*(rhoux2[3] + rhouy2[3] + rhouz2[3]); 
  ke[4] = 0.5*(rhoux2[4] + rhouy2[4] + rhouz2[4]); 
  ke[5] = 0.5*(rhoux2[5] + rhouy2[5] + rhouz2[5]); 
  ke[6] = 0.5*(rhoux2[6] + rhouy2[6] + rhouz2[6]); 
  ke[7] = 0.5*(rhoux2[7] + rhouy2[7] + rhouz2[7]); 
} 
