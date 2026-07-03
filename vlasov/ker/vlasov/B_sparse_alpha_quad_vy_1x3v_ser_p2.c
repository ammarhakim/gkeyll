#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[8]; 
  double dv10 = 2.0/dxv[1]; 
  double dv12 = 2.0/dxv[3]; 
  const double *Bx = &qmem[9]; 
  const double *Bz = &qmem[15]; 
  double Bx_quad = 0.0;
  double Bz_quad = 0.0;
  Bx_quad = 0.6324555320336759*Bx[2]-0.9486832980505137*Bx[1]+0.7071067811865475*Bx[0];
  Bz_quad = 0.6324555320336759*Bz[2]-0.9486832980505137*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[0] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[1] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[2] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[3] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[4] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[5] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[6] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[7] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[8] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);

  Bx_quad = 0.7071067811865475*Bx[0]-0.7905694150420947*Bx[2];
  Bz_quad = 0.7071067811865475*Bz[0]-0.7905694150420947*Bz[2];
  alpha_quad[9] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[10] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[11] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[12] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[13] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[14] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[15] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[16] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[17] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);

  Bx_quad = 0.6324555320336759*Bx[2]+0.9486832980505137*Bx[1]+0.7071067811865475*Bx[0];
  Bz_quad = 0.6324555320336759*Bz[2]+0.9486832980505137*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[18] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[19] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[20] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-1.837117307087383*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[21] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[22] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[23] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[24] += dv12*(0.6123724356957944*hamil[3]-1.837117307087383*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[25] += dv12*(0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[26] += dv12*(1.837117307087383*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(1.837117307087383*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);

} 
