#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_tensor_p2.h> 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_tensor_p2(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 

  for (int i = 0; i < 4; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 3; ++a) force_quad += vst_1x1v_tensor_p2_conf_dx0[i*3 + a]*phi[a]; 
    force_quad *= jacob_cx_inv; 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] -= dx10*force_quad; 
  } 
} 
