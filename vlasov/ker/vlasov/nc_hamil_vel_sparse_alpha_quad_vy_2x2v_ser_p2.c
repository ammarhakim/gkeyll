#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH double nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv0 = 0.0; 
  for (int s = 0; s < 5; ++s) { 
    const int b = vst_2x2v_ser_p2_vel_sparse_idx[s]; 
    dH_dv0 += vst_2x2v_ser_p2_vel_dv0_v1[j*8 + b]*hamil[b]; 
  } 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv10 = 2.0/dxv[2]; 
  const double vt1 = wx1 + 0.5*dv0*vst_2x2v_ser_p2_vel_nodes_v1[j*1 + 0]; 
  const double vt2 = wx2 - 0.5*dv1; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[32]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[40]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    p0_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p1_q; 
  return -om12*dH_dv0*dv10; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 9; ++i) { 
    for (int j = 0; j < 3; ++j) alpha_quad[i*3 + j] += nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
