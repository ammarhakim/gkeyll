#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH void rad_ho_alpha_quad_vz_2x3v_ser_p2(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vz = &rad[40]; 

  double rad_quad[16]; 
  for (int j = 0; j < 16; ++j) { 
    rad_quad[j] = 0.0; 
    for (int b = 0; b < 20; ++b) rad_quad[j] += vst_2x3v_ser_p2_ho_vel_ev_v2[j*20 + b]*rad_vz[b]; 
  } 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += rad_quad[j]; 
  } 
} 
