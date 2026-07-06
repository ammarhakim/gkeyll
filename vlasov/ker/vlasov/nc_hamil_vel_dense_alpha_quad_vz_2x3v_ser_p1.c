#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p1.h> 
GKYL_CU_DH double nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv0 = 0.0; 
  for (int b = 0; b < 8; ++b) dH_dv0 += vst_2x3v_ser_p1_vel_dv0_v2[j*8 + b]*hamil[b]; 
  double dH_dv1 = 0.0; 
  for (int b = 0; b < 8; ++b) dH_dv1 += vst_2x3v_ser_p1_vel_dv1_v2[j*8 + b]*hamil[b]; 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double wx3 = w[4]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv2 = dxv[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double vt1 = wx1 + 0.5*dv0*vst_2x3v_ser_p1_vel_nodes_v2[j*2 + 0]; 
  const double vt2 = wx2 + 0.5*dv1*vst_2x3v_ser_p1_vel_nodes_v2[j*2 + 1]; 
  const double vt3 = wx3 - 0.5*dv2; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[40]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[44]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[52]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[56]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[64]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[68]; 
  double p1_q = 0.0; 
  double p2_q = 0.0; 
  double p4_q = 0.0; 
  double p5_q = 0.0; 
  double p7_q = 0.0; 
  double p8_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    p1_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
    p2_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_2[a]; 
    p4_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_4[a]; 
    p5_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_5[a]; 
    p7_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_7[a]; 
    p8_q += vst_2x3v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_8[a]; 
  } 
  const double om13 = vt1*p1_q + vt2*p4_q + vt3*p7_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -om13*dH_dv0*dv10 + -om23*dH_dv1*dv11; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
