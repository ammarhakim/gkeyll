#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[3]; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bx = &qmem[6]; 
  const double *By = &qmem[8]; 
  double Bx_quad = 0.0;
  double By_quad = 0.0;
  Bx_quad = 0.7071067811865475*Bx[0]-0.7071067811865475*Bx[1];
  By_quad = 0.7071067811865475*By[0]-0.7071067811865475*By[1];
  alpha_quad[0] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[1] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[2] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[3] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);

  Bx_quad = 0.7071067811865475*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = 0.7071067811865475*By[1]+0.7071067811865475*By[0];
  alpha_quad[4] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[5] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[6] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[7] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);

} 
