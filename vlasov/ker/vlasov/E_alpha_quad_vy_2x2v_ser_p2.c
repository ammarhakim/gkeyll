#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void E_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ey = &qmem[8]; 

  for (int i = 0; i < 9; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_2x2v_ser_p2_conf_ev[i*8 + a]*Ey[a]; 
    for (int j = 0; j < 3; ++j) alpha_quad[i*3 + j] += force_quad; 
  } 
} 
