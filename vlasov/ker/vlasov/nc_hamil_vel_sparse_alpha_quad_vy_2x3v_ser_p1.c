#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p1.h> 
GKYL_CU_DH double nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv0 = 0.0; 
  for (int s = 0; s < 4; ++s) { 
    const int b = vst_2x3v_ser_p1_vel_sparse_idx[s]; 
    dH_dv0 += vst_2x3v_ser_p1_vel_dv0_v1[j*8 + b]*hamil[b]; 
  } 
  double dH_dv2 = 0.0; 
  for (int s = 0; s < 4; ++s) { 
    const int b = vst_2x3v_ser_p1_vel_sparse_idx[s]; 
    dH_dv2 += vst_2x3v_ser_p1_vel_dv2_v1[j*8 + b]*hamil[b]; 
  } 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double wx3 = w[4]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv2 = dxv[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[4]; 
  const double vt1 = wx1 + 0.5*dv0*vst_2x3v_ser_p1_vel_nodes_v1[j*2 + 0]; 
  const double vt2 = wx2 - 0.5*dv1; 
  const double vt3 = wx3 + 0.5*dv2*vst_2x3v_ser_p1_vel_nodes_v1[j*2 + 1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[36]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[44]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[48]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[56]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[60]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[68]; 
  double p0_q = 0.0; 
  double p2_q = 0.0; 
  double p3_q = 0.0; 
  double p5_q = 0.0; 
  double p6_q = 0.0; 
  double p8_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    p0_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
    p2_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_2[a]; 
    p3_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_3[a]; 
    p5_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_5[a]; 
    p6_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_6[a]; 
    p8_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_8[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p3_q + vt3*p6_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -om12*dH_dv0*dv10 + om23*dH_dv2*dv12; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
