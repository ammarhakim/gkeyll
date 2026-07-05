#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p1.h> 
GKYL_CU_DH void E_alpha_quad_vy_1x3v_ser_p1(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ey = &qmem[2]; 

  for (int i = 0; i < 2; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 2; ++a) force_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Ey[a]; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += force_quad; 
  } 
} 
