#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *jacob_vel_surf_vy = &jacob_vel_surf[4]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[8]; 
  double dv11 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[3]; 
  const double *By = &qmem[12]; 
  const double *Bz = &qmem[15]; 
  double By_quad = 0.0;
  double Bz_quad = 0.0;
  By_quad = 0.6324555320336759*By[2]-0.9486832980505137*By[1]+0.7071067811865475*By[0];
  Bz_quad = 0.6324555320336759*Bz[2]-0.9486832980505137*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[0] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[1] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[2] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[3] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[4] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[5] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[6] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[7] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[8] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);

  By_quad = 0.7071067811865475*By[0]-0.7905694150420947*By[2];
  Bz_quad = 0.7071067811865475*Bz[0]-0.7905694150420947*Bz[2];
  alpha_quad[9] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[10] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[11] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[12] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[13] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[14] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[15] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[16] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[17] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);

  By_quad = 0.6324555320336759*By[2]+0.9486832980505137*By[1]+0.7071067811865475*By[0];
  Bz_quad = 0.6324555320336759*Bz[2]+0.9486832980505137*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[18] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[19] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[20] += dv11*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[21] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[22] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[23] += dv11*(0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[24] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[25] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[26] += dv11*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*By_quad/(jacob_vel_surf_vz[0]);

} 
