#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_tensor_p2.h> 
GKYL_CU_DH void E_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  for (int i = 0; i < 16; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 9; ++a) force_quad += vst_2x2v_tensor_p2_conf_ev[i*9 + a]*Ex[a]; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += force_quad; 
  } 
} 
