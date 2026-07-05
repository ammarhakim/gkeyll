#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void phi_ho_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double *jacob_cy = &jacob_pos[3]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 

  for (int i = 0; i < 16; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_2x2v_ser_p2_ho_conf_dx1[i*8 + a]*phi[a]; 
    force_quad *= jacob_cy_inv; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] -= dx11*force_quad; 
  } 
} 
