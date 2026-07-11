#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *Bz = &qmem[10]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  double Bz_quad = 0.0;
  Bz_quad = 0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1];
  alpha_quad[0] -= 2.0/(dxv[1]*jacob_vel_surf_vx[0])*(1.060660171779821*hamil[7]-1.060660171779821*hamil[6]-0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*Bz_quad;
  alpha_quad[1] -= 2.0/(dxv[1]*jacob_vel_surf_vx[0])*(1.060660171779821*hamil[7]-1.060660171779821*hamil[6]-0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*Bz_quad;

  Bz_quad = 0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[2] -= 2.0/(dxv[1]*jacob_vel_surf_vx[0])*(-(1.060660171779821*hamil[7])-1.060660171779821*hamil[6]+0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*Bz_quad;
  alpha_quad[3] -= 2.0/(dxv[1]*jacob_vel_surf_vx[0])*(-(1.060660171779821*hamil[7])-1.060660171779821*hamil[6]+0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*Bz_quad;

} 
