#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_sparse_ho_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
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
  Bx_quad = 0.9681844646844028*Bx[2]-1.0546722811938851*Bx[1]+0.7071067811865475*Bx[0];
  Bz_quad = 0.9681844646844028*Bz[2]-1.0546722811938851*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[0] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[1] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[2] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[3] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[4] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[5] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[6] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[7] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[8] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[9] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[10] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[11] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[12] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[13] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[14] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[15] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);

  Bx_quad = -(0.5164305132317774*Bx[2])-0.41639003950091297*Bx[1]+0.7071067811865475*Bx[0];
  Bz_quad = -(0.5164305132317774*Bz[2])-0.41639003950091297*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[16] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[17] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[18] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[19] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[20] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[21] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[22] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[23] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[24] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[25] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[26] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[27] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[28] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[29] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[30] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[31] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);

  Bx_quad = -(0.5164305132317774*Bx[2])+0.41639003950091297*Bx[1]+0.7071067811865475*Bx[0];
  Bz_quad = -(0.5164305132317774*Bz[2])+0.41639003950091297*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[32] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[33] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[34] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[35] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[36] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[37] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[38] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[39] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[40] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[41] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[42] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[43] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[44] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[45] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[46] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[47] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);

  Bx_quad = 0.9681844646844028*Bx[2]+1.0546722811938851*Bx[1]+0.7071067811865475*Bx[0];
  Bz_quad = 0.9681844646844028*Bz[2]+1.0546722811938851*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[48] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[49] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[50] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[51] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-2.0423640903852514*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[52] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[53] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[54] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[55] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.6123724356957944*hamil[1]-0.8063358442568421*hamil[7])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[56] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[57] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[58] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[59] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(0.8063358442568421*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[60] += dv12*(0.6123724356957944*hamil[3]-2.0423640903852514*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[61] += dv12*(0.6123724356957944*hamil[3]-0.8063358442568421*hamil[9])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[62] += dv12*(0.8063358442568421*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);
  alpha_quad[63] += dv12*(2.0423640903852514*hamil[9]+0.6123724356957944*hamil[3])*Bx_quad/(jacob_vel_surf_vz[0]) - dv10*(2.0423640903852514*hamil[7]+0.6123724356957944*hamil[1])*Bz_quad/(jacob_vel_surf_vx[0]);

} 
