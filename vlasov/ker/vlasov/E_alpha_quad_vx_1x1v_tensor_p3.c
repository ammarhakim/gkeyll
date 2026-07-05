#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_tensor_p3.h> 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  for (int i = 0; i < 5; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 4; ++a) force_quad += vst_1x1v_tensor_p3_conf_ev[i*4 + a]*Ex[a]; 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += force_quad; 
  } 
} 
