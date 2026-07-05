#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p1.h> 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_ser_p1(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double *jacob_cy = &jacob_pos[2]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 

  for (int i = 0; i < 4; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 4; ++a) force_quad += vst_2x2v_ser_p1_conf_dx1[i*4 + a]*phi[a]; 
    force_quad *= jacob_cy_inv; 
    for (int j = 0; j < 2; ++j) alpha_quad[i*2 + j] -= dx11*force_quad; 
  } 
} 
