#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_ser_p3.h> 
GKYL_CU_DH void rad_ho_alpha_quad_vx_1x1v_ser_p3(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vx = &rad[0]; 

  double rad_quad[1]; 
  for (int j = 0; j < 1; ++j) { 
    rad_quad[j] = 0.0; 
    for (int b = 0; b < 4; ++b) rad_quad[j] += vst_1x1v_ser_p3_ho_vel_ev_v0[j*4 + b]*rad_vx[b]; 
  } 
  for (int i = 0; i < 5; ++i) { 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += rad_quad[j]; 
  } 
} 
