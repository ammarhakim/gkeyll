#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p1.h> 
GKYL_CU_DH void B_dense_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvx[j] += vst_1x3v_ser_p1_vel_dv0_v2[j*8 + b]*hamil[b]; 
  } 
  double dH_dvy[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvy[j] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvy[j] += vst_1x3v_ser_p1_vel_dv1_v2[j*8 + b]*hamil[b]; 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[3]; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bx = &qmem[6]; 
  const double *By = &qmem[8]; 
  for (int i = 0; i < 2; ++i) { 
    double Bx_quad = 0.0; 
    double By_quad = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      Bx_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Bx[a]; 
      By_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*By[a]; 
    } 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += dv10*dH_dvx[j]*By_quad/jacob_vel_surf_vx[0] - dv11*dH_dvy[j]*Bx_quad/jacob_vel_surf_vy[0]; 
  } 
} 
