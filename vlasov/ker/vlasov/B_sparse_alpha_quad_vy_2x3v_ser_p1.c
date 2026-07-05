#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p1.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vy_2x3v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_2x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dvx[j] += vst_2x3v_ser_p1_vel_dv0_v1[j*8 + b]*hamil[b]; 
    } 
  } 
  double dH_dvz[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvz[j] = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_2x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dvz[j] += vst_2x3v_ser_p1_vel_dv2_v1[j*8 + b]*hamil[b]; 
    } 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf[6]; 
  double dv10 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[4]; 
  const double *Bx = &qmem[12]; 
  const double *Bz = &qmem[20]; 
  for (int i = 0; i < 4; ++i) { 
    double Bx_quad = 0.0; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      Bx_quad += vst_2x3v_ser_p1_conf_ev[i*4 + a]*Bx[a]; 
      Bz_quad += vst_2x3v_ser_p1_conf_ev[i*4 + a]*Bz[a]; 
    } 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += dv12*dH_dvz[j]*Bx_quad/jacob_vel_surf_vz[0] - dv10*dH_dvx[j]*Bz_quad/jacob_vel_surf_vx[0]; 
  } 
} 
