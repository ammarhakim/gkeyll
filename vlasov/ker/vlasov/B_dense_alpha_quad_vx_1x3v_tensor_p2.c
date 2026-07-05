#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_tensor_p2.h> 
GKYL_CU_DH void B_dense_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvy[16]; 
  for (int j = 0; j < 16; ++j) { 
    dH_dvy[j] = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvy[j] += vst_1x3v_tensor_p2_vel_dv1_v0[j*27 + b]*hamil[b]; 
  } 
  double dH_dvz[16]; 
  for (int j = 0; j < 16; ++j) { 
    dH_dvz[j] = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvz[j] += vst_1x3v_tensor_p2_vel_dv2_v0[j*27 + b]*hamil[b]; 
  } 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[4]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[8]; 
  double jac_vy_n[16]; 
  jac_vy_n[0] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[1] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[2] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[3] = 1.0*jacob_vel_surf_vy[0]; 
  jac_vy_n[4] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[5] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[6] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[7] = 1.0*jacob_vel_surf_vy[1]; 
  jac_vy_n[8] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[9] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[10] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[11] = 1.0*jacob_vel_surf_vy[2]; 
  jac_vy_n[12] = 1.0*jacob_vel_surf_vy[3]; 
  jac_vy_n[13] = 1.0*jacob_vel_surf_vy[3]; 
  jac_vy_n[14] = 1.0*jacob_vel_surf_vy[3]; 
  jac_vy_n[15] = 1.0*jacob_vel_surf_vy[3]; 
  double jac_vz_n[16]; 
  jac_vz_n[0] = 1.0*jacob_vel_surf_vz[0]; 
  jac_vz_n[1] = 1.0*jacob_vel_surf_vz[1]; 
  jac_vz_n[2] = 1.0*jacob_vel_surf_vz[2]; 
  jac_vz_n[3] = 1.0*jacob_vel_surf_vz[3]; 
  jac_vz_n[4] = 1.0*jacob_vel_surf_vz[0]; 
  jac_vz_n[5] = 1.0*jacob_vel_surf_vz[1]; 
  jac_vz_n[6] = 1.0*jacob_vel_surf_vz[2]; 
  jac_vz_n[7] = 1.0*jacob_vel_surf_vz[3]; 
  jac_vz_n[8] = 1.0*jacob_vel_surf_vz[0]; 
  jac_vz_n[9] = 1.0*jacob_vel_surf_vz[1]; 
  jac_vz_n[10] = 1.0*jacob_vel_surf_vz[2]; 
  jac_vz_n[11] = 1.0*jacob_vel_surf_vz[3]; 
  jac_vz_n[12] = 1.0*jacob_vel_surf_vz[0]; 
  jac_vz_n[13] = 1.0*jacob_vel_surf_vz[1]; 
  jac_vz_n[14] = 1.0*jacob_vel_surf_vz[2]; 
  jac_vz_n[15] = 1.0*jacob_vel_surf_vz[3]; 
  double dv11 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[3]; 
  const double *By = &qmem[12]; 
  const double *Bz = &qmem[15]; 
  for (int i = 0; i < 4; ++i) { 
    double By_quad = 0.0; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      By_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*By[a]; 
      Bz_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*Bz[a]; 
    } 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += dv11*dH_dvy[j]*Bz_quad/jac_vy_n[j] - dv12*dH_dvz[j]*By_quad/jac_vz_n[j]; 
  } 
} 
