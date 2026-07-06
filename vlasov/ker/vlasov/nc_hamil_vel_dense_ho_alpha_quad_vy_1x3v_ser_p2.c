#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv0 = 0.0; 
  for (int b = 0; b < 20; ++b) dH_dv0 += vst_1x3v_ser_p2_ho_vel_dv0_v1[j*20 + b]*hamil[b]; 
  double dH_dv2 = 0.0; 
  for (int b = 0; b < 20; ++b) dH_dv2 += vst_1x3v_ser_p2_ho_vel_dv2_v1[j*20 + b]*hamil[b]; 
  const double wx1 = w[1]; 
  const double wx2 = w[2]; 
  const double wx3 = w[3]; 
  const double dv0 = dxv[1]; 
  const double dv1 = dxv[2]; 
  const double dv2 = dxv[3]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv12 = 2.0/dxv[3]; 
  const double vt1 = wx1 + 0.5*dv0*vst_1x3v_ser_p2_ho_vel_nodes_v1[j*2 + 0]; 
  const double vt2 = wx2 - 0.5*dv1; 
  const double vt3 = wx3 + 0.5*dv2*vst_1x3v_ser_p2_ho_vel_nodes_v1[j*2 + 1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[27]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[33]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[36]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[42]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[45]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[51]; 
  double p0_q = 0.0; 
  double p2_q = 0.0; 
  double p3_q = 0.0; 
  double p5_q = 0.0; 
  double p6_q = 0.0; 
  double p8_q = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    p0_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_0[a]; 
    p2_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_2[a]; 
    p3_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_3[a]; 
    p5_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_5[a]; 
    p6_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_6[a]; 
    p8_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_8[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p3_q + vt3*p6_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -om12*dH_dv0*dv10 + om23*dH_dv2*dv12; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_ser_p2_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
