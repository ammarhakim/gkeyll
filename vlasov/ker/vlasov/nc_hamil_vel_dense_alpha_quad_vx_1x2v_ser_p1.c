#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p1.h> 
GKYL_CU_DH double nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv1 = 0.0; 
  for (int b = 0; b < 4; ++b) dH_dv1 += vst_1x2v_ser_p1_vel_dv1_v0[j*4 + b]*hamil[b]; 
  const double wx1 = w[1]; 
  const double wx2 = w[2]; 
  const double dv0 = dxv[1]; 
  const double dv1 = dxv[2]; 
  const double dv11 = 2.0/dxv[2]; 
  const double vt1 = wx1 - 0.5*dv0; 
  const double vt2 = wx2 + 0.5*dv1*vst_1x2v_ser_p1_vel_nodes_v0[j*1 + 0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[10]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    p0_q += vst_1x2v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_1x2v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_1[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p1_q; 
  return om12*dH_dv1*dv11; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 2; ++j) alpha_quad[i*2 + j] += nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
