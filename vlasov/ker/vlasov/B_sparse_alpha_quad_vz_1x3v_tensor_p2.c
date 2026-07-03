#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, const double *jacob_vel_surf, 
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
  Bx_quad = 0.9681844646844028*Bx[2]-1.0546722811938851*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = 0.9681844646844028*By[2]-1.0546722811938851*By[1]+0.7071067811865475*By[0];
  alpha_quad[0] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[1] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[2] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[3] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[4] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[5] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[6] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[7] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[8] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[9] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[10] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[11] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[12] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[13] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[14] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[15] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);

  Bx_quad = -(0.5164305132317774*Bx[2])-0.41639003950091297*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = -(0.5164305132317774*By[2])-0.41639003950091297*By[1]+0.7071067811865475*By[0];
  alpha_quad[16] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[17] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[18] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[19] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[20] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[21] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[22] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[23] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[24] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[25] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[26] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[27] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[28] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[29] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[30] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[31] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);

  Bx_quad = -(0.5164305132317774*Bx[2])+0.41639003950091297*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = -(0.5164305132317774*By[2])+0.41639003950091297*By[1]+0.7071067811865475*By[0];
  alpha_quad[32] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[33] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[34] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[35] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[36] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[37] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[38] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[39] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[40] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[41] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[42] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[43] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[44] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[45] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[46] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[47] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);

  Bx_quad = 0.9681844646844028*Bx[2]+1.0546722811938851*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = 0.9681844646844028*By[2]+1.0546722811938851*By[1]+0.7071067811865475*By[0];
  alpha_quad[48] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[49] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[50] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[51] += dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[0]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[52] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[53] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[54] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[55] += dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*By_quad/(1.0*jacob_vel_surf_vx[1]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[56] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[57] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[58] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[59] += dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[2]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);
  alpha_quad[60] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-2.0423640903852514*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[0]);
  alpha_quad[61] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.6123724356957944*hamil[2]-0.8063358442568421*hamil[8])*Bx_quad/(1.0*jacob_vel_surf_vy[1]);
  alpha_quad[62] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(0.8063358442568421*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[2]);
  alpha_quad[63] += dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*By_quad/(1.0*jacob_vel_surf_vx[3]) - dv11*(2.0423640903852514*hamil[8]+0.6123724356957944*hamil[2])*Bx_quad/(1.0*jacob_vel_surf_vy[3]);

} 
