#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p1.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[2]; 
  for (int j = 0; j < 2; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int s = 0; s < 3; ++s) { 
      const int b = vst_1x2v_ser_p1_vel_sparse_idx[s]; 
      dH_dvx[j] += vst_1x2v_ser_p1_vel_dv0_v1[j*4 + b]*hamil[b]; 
    } 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *Bz = &qmem[10]; 
  for (int i = 0; i < 2; ++i) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 2; ++a) Bz_quad += vst_1x2v_ser_p1_conf_ev[i*2 + a]*Bz[a]; 
    for (int j = 0; j < 2; ++j) alpha_quad[i*2 + j] -= 2.0/(dxv[1]*jacob_vel_surf_vx[0])*dH_dvx[j]*Bz_quad; 
  } 
} 
