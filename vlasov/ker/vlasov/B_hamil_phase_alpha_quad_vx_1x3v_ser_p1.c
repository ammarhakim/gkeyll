#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *jacob_vel_surf_vy = &jacob_vel_surf[3]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[6]; 
  double dv11 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[3]; 
  const double *By = &qmem[8]; 
  const double *Bz = &qmem[10]; 
  double By_quad = 0.0;
  double Bz_quad = 0.0;
  By_quad = 0.7071067811865475*By[0]-0.7071067811865475*By[1];
  Bz_quad = 0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1];
  alpha_quad[0] += dv11*(-(0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[11]-0.4330127018922193*hamil[10]-0.75*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(-(0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[1] += dv11*(0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[11]+0.4330127018922193*hamil[10]-0.75*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(-(0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[2] += dv11*(-(0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[11]-0.4330127018922193*hamil[10]-0.75*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[3] += dv11*(0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[11]+0.4330127018922193*hamil[10]-0.75*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);

  By_quad = 0.7071067811865475*By[1]+0.7071067811865475*By[0];
  Bz_quad = 0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[4] += dv11*(0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[11]-0.4330127018922193*hamil[10]-0.75*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[5] += dv11*(-(0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[11]+0.4330127018922193*hamil[10]-0.75*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[6] += dv11*(0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[11]-0.4330127018922193*hamil[10]-0.75*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(-(0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);
  alpha_quad[7] += dv11*(-(0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[11]+0.4330127018922193*hamil[10]-0.75*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bz_quad/(jacob_vel_surf_vy[0]) - dv12*(-(0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*By_quad/(jacob_vel_surf_vz[0]);

} 
