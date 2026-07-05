#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void phi_ho_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 

  for (int i = 0; i < 16; ++i) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_2x2v_ser_p2_ho_conf_dx0[i*8 + a]*phi[a]; 
    force_quad *= jacob_cx_inv; 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] -= dx10*force_quad; 
  } 
} 
