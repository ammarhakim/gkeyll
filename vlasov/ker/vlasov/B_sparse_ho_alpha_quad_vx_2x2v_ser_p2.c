#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void B_sparse_ho_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvy[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvy[j] = 0.0; 
    for (int s = 0; s < 5; ++s) { 
      const int b = vst_2x2v_ser_p2_ho_vel_sparse_idx[s]; 
      dH_dvy[j] += vst_2x2v_ser_p2_ho_vel_dv1_v0[j*8 + b]*hamil[b]; 
    } 
  } 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[4]; 
  const double *Bz = &qmem[40]; 
  for (int i = 0; i < 16; ++i) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 8; ++a) Bz_quad += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*Bz[a]; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += 2.0/(dxv[3]*jacob_vel_surf_vy[0])*dH_dvy[j]*Bz_quad; 
  } 
} 
