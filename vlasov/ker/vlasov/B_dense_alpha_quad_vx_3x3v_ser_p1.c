#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_3x3v_ser_p1.h> 
GKYL_CU_DH void B_dense_alpha_quad_vx_3x3v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvy[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvy[j] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvy[j] += vst_3x3v_ser_p1_vel_dv1_v0[j*8 + b]*hamil[b]; 
  } 
  double dH_dvz[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvz[j] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvz[j] += vst_3x3v_ser_p1_vel_dv2_v0[j*8 + b]*hamil[b]; 
  } 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[3]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[6]; 
  double dv11 = 2.0/dxv[4]; 
  double dv12 = 2.0/dxv[5]; 
  const double *By = &qmem[32]; 
  const double *Bz = &qmem[40]; 
  for (int i = 0; i < 8; ++i) { 
    double By_quad = 0.0; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      By_quad += vst_3x3v_ser_p1_conf_ev[i*8 + a]*By[a]; 
      Bz_quad += vst_3x3v_ser_p1_conf_ev[i*8 + a]*Bz[a]; 
    } 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += dv11*dH_dvy[j]*Bz_quad/jacob_vel_surf_vy[0] - dv12*dH_dvz[j]*By_quad/jacob_vel_surf_vz[0]; 
  } 
} 
