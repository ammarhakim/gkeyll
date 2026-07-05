#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void B_dense_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[3]; 
  for (int j = 0; j < 3; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvx[j] += vst_2x2v_ser_p2_vel_dv0_v1[j*8 + b]*hamil[b]; 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *Bz = &qmem[40]; 
  for (int i = 0; i < 9; ++i) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 8; ++a) Bz_quad += vst_2x2v_ser_p2_conf_ev[i*8 + a]*Bz[a]; 
    for (int j = 0; j < 3; ++j) alpha_quad[i*3 + j] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*dH_dvx[j]*Bz_quad; 
  } 
} 
