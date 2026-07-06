#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_ser_p2.h> 
GKYL_CU_DH double nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  return 0.0; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
