#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p2.h> 
GKYL_CU_DH void rad_alpha_quad_vx_1x2v_ser_p2(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vx = &rad[0]; 

  double rad_quad[3]; 
  for (int j = 0; j < 3; ++j) { 
    rad_quad[j] = 0.0; 
    for (int b = 0; b < 8; ++b) rad_quad[j] += vst_1x2v_ser_p2_vel_ev_v0[j*8 + b]*rad_vx[b]; 
  } 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) alpha_quad[i*3 + j] += rad_quad[j]; 
  } 
} 
