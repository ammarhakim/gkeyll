#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_3x3v_ser_p1.h> 
GKYL_CU_DH void phi_alpha_quad_vz_3x3v_ser_p1(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx12 = 2.0/dxv[2]; 
  const double *jacob_cz = &jacob_pos[4]; 
  const double jacob_cz_inv = 1.0/jacob_cz[0]; 

  for (int i = 0; i < 8; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_3x3v_ser_p1_conf_dx2[i*8 + a]*phi[a]; 
    force_quad *= jacob_cz_inv; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] -= dx12*force_quad; 
  } 
} 
