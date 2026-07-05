#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_tensor_p2.h> 
GKYL_CU_DH void B_sparse_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dH_dvx[4]; 
  for (int j = 0; j < 4; ++j) { 
    dH_dvx[j] = 0.0; 
    for (int s = 0; s < 5; ++s) { 
      const int b = vst_2x2v_tensor_p2_vel_sparse_idx[s]; 
      dH_dvx[j] += vst_2x2v_tensor_p2_vel_dv0_v1[j*9 + b]*hamil[b]; 
    } 
  } 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  const double *Bz = &qmem[45]; 
  for (int i = 0; i < 16; ++i) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 9; ++a) Bz_quad += vst_2x2v_tensor_p2_conf_ev[i*9 + a]*Bz[a]; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] -= 2.0/(dxv[2]*jacob_vel_surf_vx[j])*dH_dvx[j]*Bz_quad; 
  } 
} 
