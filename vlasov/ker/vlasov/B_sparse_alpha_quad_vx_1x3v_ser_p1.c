#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p1.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vx_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvy[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvy[j] = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_1x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dvy[j] += vst_1x3v_ser_p1_vel_dv1_v0[j*8 + b]*hamil[b]; 
    } 
  } 
  double dH_dvz[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvz[j] = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_1x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dvz[j] += vst_1x3v_ser_p1_vel_dv2_v0[j*8 + b]*hamil[b]; 
    } 
  } 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[3]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[6]; 
  double dv11 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[3]; 
  const double *By = &qmem[8]; 
  const double *Bz = &qmem[10]; 
  for (int i = 0; i < 2; ++i) { 
    double By_quad = 0.0; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      By_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*By[a]; 
      Bz_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Bz[a]; 
    } 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += dv11*dH_dvy[j]*Bz_quad/jacob_vel_surf_vy[0] - dv12*dH_dvz[j]*By_quad/jacob_vel_surf_vz[0]; 
  } 
} 
