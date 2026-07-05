#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_tensor_p2.h> 
GKYL_CU_DH void rad_alpha_quad_vz_1x3v_tensor_p2(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vz = &rad[54]; 

  double rad_quad[16]; 
  for (int j = 0; j < 16; ++j) { 
    rad_quad[j] = 0.0; 
    for (int b = 0; b < 27; ++b) rad_quad[j] += vst_1x3v_tensor_p2_vel_ev_v2[j*27 + b]*rad_vz[b]; 
  } 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += rad_quad[j]; 
  } 
} 
