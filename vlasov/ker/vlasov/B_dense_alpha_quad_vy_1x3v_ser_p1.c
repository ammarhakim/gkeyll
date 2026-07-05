#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p1.h> 
GKYL_CU_DH void B_dense_alpha_quad_vy_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvx[j] += vst_1x3v_ser_p1_vel_dv0_v1[j*8 + b]*hamil[b]; 
  } 
  double dH_dvz[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvz[j] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvz[j] += vst_1x3v_ser_p1_vel_dv2_v1[j*8 + b]*hamil[b]; 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[6]; 
  double dv10 = 2.0/dxv[1]; 
  double dv12 = 2.0/dxv[3]; 
  const double *Bx = &qmem[6]; 
  const double *Bz = &qmem[10]; 
  for (int i = 0; i < 2; ++i) { 
    double Bx_quad = 0.0; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      Bx_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Bx[a]; 
      Bz_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Bz[a]; 
    } 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += dv12*dH_dvz[j]*Bx_quad/jacob_vel_surf_vz[0] - dv10*dH_dvx[j]*Bz_quad/jacob_vel_surf_vx[0]; 
  } 
} 
