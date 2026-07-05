#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_tensor_p2.h> 
GKYL_CU_DH void B_dense_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[16]; 
  for (int j = 0; j < 16; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvx[j] += vst_1x3v_tensor_p2_vel_dv0_v2[j*27 + b]*hamil[b]; 
  } 
  double dH_dvy[16]; 
  for (int j = 0; j < 16; ++j) { 
    dH_dvy[j] = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvy[j] += vst_1x3v_tensor_p2_vel_dv1_v2[j*27 + b]*hamil[b]; 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[4]; 
  double jac_vx_n[16]; 
  jac_vx_n[0] = 1.0*jacob_vel_surf_vx[0]; 
  jac_vx_n[1] = 1.0*jacob_vel_surf_vx[0]; 
  jac_vx_n[2] = 1.0*jacob_vel_surf_vx[0]; 
  jac_vx_n[3] = 1.0*jacob_vel_surf_vx[0]; 
  jac_vx_n[4] = 1.0*jacob_vel_surf_vx[1]; 
  jac_vx_n[5] = 1.0*jacob_vel_surf_vx[1]; 
  jac_vx_n[6] = 1.0*jacob_vel_surf_vx[1]; 
  jac_vx_n[7] = 1.0*jacob_vel_surf_vx[1]; 
  jac_vx_n[8] = 1.0*jacob_vel_surf_vx[2]; 
  jac_vx_n[9] = 1.0*jacob_vel_surf_vx[2]; 
  jac_vx_n[10] = 1.0*jacob_vel_surf_vx[2]; 
  jac_vx_n[11] = 1.0*jacob_vel_surf_vx[2]; 
  jac_vx_n[12] = 1.0*jacob_vel_surf_vx[3]; 
  jac_vx_n[13] = 1.0*jacob_vel_surf_vx[3]; 
  jac_vx_n[14] = 1.0*jacob_vel_surf_vx[3]; 
  jac_vx_n[15] = 1.0*jacob_vel_surf_vx[3]; 
  double jac_vy_n[16]; 
  jac_vy_n[0] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[1] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[2] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[3] = 1.0*jacob_vel_surf_vy[3]; 
  jac_vy_n[4] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[5] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[6] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[7] = 1.0*jacob_vel_surf_vy[3]; 
  jac_vy_n[8] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[9] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[10] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[11] = 1.0*jacob_vel_surf_vy[3]; 
  jac_vy_n[12] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[13] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[14] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[15] = 1.0*jacob_vel_surf_vy[3]; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bx = &qmem[9]; 
  const double *By = &qmem[12]; 
  for (int i = 0; i < 4; ++i) { 
    double Bx_quad = 0.0; 
    double By_quad = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      Bx_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*Bx[a]; 
      By_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*By[a]; 
    } 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += dv10*dH_dvx[j]*By_quad/jac_vx_n[j] - dv11*dH_dvy[j]*Bx_quad/jac_vy_n[j]; 
  } 
} 
