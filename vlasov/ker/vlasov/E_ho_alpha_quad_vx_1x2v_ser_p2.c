#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p2.h> 
GKYL_CU_DH void E_ho_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  for (int i = 0; i < 4; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 3; ++a) force_quad += vst_1x2v_ser_p2_ho_conf_ev[i*3 + a]*Ex[a]; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += force_quad; 
  } 
} 
