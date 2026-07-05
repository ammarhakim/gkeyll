#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void rad_alpha_quad_vy_2x2v_ser_p2(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vy = &rad[8]; 

  double rad_quad[3]; 
  for (int j = 0; j < 3; ++j) { 
    rad_quad[j] = 0.0; 
    for (int b = 0; b < 8; ++b) rad_quad[j] += vst_2x2v_ser_p2_vel_ev_v1[j*8 + b]*rad_vy[b]; 
  } 
  for (int i = 0; i < 9; ++i) { 
    for (int j = 0; j < 3; ++j) alpha_quad[i*3 + j] += rad_quad[j]; 
  } 
} 
