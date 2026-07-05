#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_tensor_p2.h> 
GKYL_CU_DH void B_dense_alpha_quad_vy_1x3v_tensor_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[16]; 
  for (int j = 0; j < 16; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvx[j] += vst_1x3v_tensor_p2_vel_dv0_v1[j*27 + b]*hamil[b]; 
  } 
  double dH_dvz[16]; 
  for (int j = 0; j < 16; ++j) { 
    dH_dvz[j] = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvz[j] += vst_1x3v_tensor_p2_vel_dv2_v1[j*27 + b]*hamil[b]; 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[8]; 
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
  double dv10 = 2.0/dxv[1]; 
  double dv12 = 2.0/dxv[3]; 
  const double *Bx = &qmem[9]; 
  const double *Bz = &qmem[15]; 
  for (int i = 0; i < 4; ++i) { 
    double Bx_quad = 0.0; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      Bx_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*Bx[a]; 
      Bz_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*Bz[a]; 
    } 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += dv12*dH_dvz[j]*Bx_quad/jac_vz_n[j] - dv10*dH_dvx[j]*Bz_quad/jac_vx_n[j]; 
  } 
} 
