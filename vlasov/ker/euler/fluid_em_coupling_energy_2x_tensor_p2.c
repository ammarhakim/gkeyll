#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_energy_2x_tensor_p2(const double* ke_old, const double* ke_new, double* GKYL_RESTRICT fluid) 
{ 
  // ke_old: Kinetic energy at the old time step.
  // ke_new: Kinetic energy at the new time step.
  // fluid:  [rho, rho ux, rho uy, rho uz...], Fluid output state vector.
  //         Computes the energy update from the old and new kinetic energy. 

  double *energy = &fluid[36]; 

  energy[0] += ke_new[0]-1.0*ke_old[0]; 
  energy[1] += ke_new[1]-1.0*ke_old[1]; 
  energy[2] += ke_new[2]-1.0*ke_old[2]; 
  energy[3] += ke_new[3]-1.0*ke_old[3]; 
  energy[4] += ke_new[4]-1.0*ke_old[4]; 
  energy[5] += ke_new[5]-1.0*ke_old[5]; 
  energy[6] += ke_new[6]-1.0*ke_old[6]; 
  energy[7] += ke_new[7]-1.0*ke_old[7]; 
  energy[8] += ke_new[8]-1.0*ke_old[8]; 
} 
