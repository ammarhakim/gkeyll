#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[4]; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bx = &qmem[9]; 
  const double *By = &qmem[12]; 
  double Bx_quad = 0.0;
  double By_quad = 0.0;
  Bx_quad = 0.6324555320336759*Bx[2]-0.9486832980505137*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = 0.6324555320336759*By[2]-0.9486832980505137*By[1]+0.7071067811865475*By[0];
  alpha_quad[0] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[1] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[2] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[3] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[4] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[5] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[6] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[7] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[8] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);

  Bx_quad = 0.7071067811865475*Bx[0]-0.7905694150420947*Bx[2];
  By_quad = 0.7071067811865475*By[0]-0.7905694150420947*By[2];
  alpha_quad[9] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[10] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[11] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[12] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[13] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[14] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[15] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[16] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[17] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);

  Bx_quad = 0.6324555320336759*Bx[2]+0.9486832980505137*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = 0.6324555320336759*By[2]+0.9486832980505137*By[1]+0.7071067811865475*By[0];
  alpha_quad[18] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[19] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[20] += dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[21] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[22] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[23] += dv10*(0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[24] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[25] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);
  alpha_quad[26] += dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(jacob_vel_surf_vx[0]) - dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(jacob_vel_surf_vy[0]);

} 
