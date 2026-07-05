#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_tensor_p2.h> 
GKYL_CU_DH void E_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ez = &qmem[6]; 

  for (int i = 0; i < 4; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 3; ++a) force_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*Ez[a]; 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += force_quad; 
  } 
} 
