#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_energy_1x_ser_p1(const double* ke_old, const double* ke_new, double* GKYL_RESTRICT fluid) 
{ 
  // ke_old: Kinetic energy at the old time step.
  // ke_new: Kinetic energy at the new time step.
  // fluid:  [rho, rho ux, rho uy, rho uz...], Fluid output state vector.
  //         Computes the energy update from the old and new kinetic energy. 

  double *energy = &fluid[8]; 

  energy[0] += ke_new[0]-1.0*ke_old[0]; 
  energy[1] += ke_new[1]-1.0*ke_old[1]; 
} 
