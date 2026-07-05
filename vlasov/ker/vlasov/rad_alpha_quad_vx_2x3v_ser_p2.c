#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH void rad_alpha_quad_vx_2x3v_ser_p2(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vx = &rad[0]; 

  double rad_quad[9]; 
  for (int j = 0; j < 9; ++j) { 
    rad_quad[j] = 0.0; 
    for (int b = 0; b < 20; ++b) rad_quad[j] += vst_2x3v_ser_p2_vel_ev_v0[j*20 + b]*rad_vx[b]; 
  } 
  for (int i = 0; i < 9; ++i) { 
    for (int j = 0; j < 9; ++j) alpha_quad[i*9 + j] += rad_quad[j]; 
  } 
} 
