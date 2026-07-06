#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double hamil_vel_sparse_alpha_quad_x_2x3v_ser_p2_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[16]; 
  double dH_dv0 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_2x3v_ser_p2_vel_sparse_idx[s]; 
    dH_dv0 += vst_2x3v_ser_p2_vel_vol_dv0[m*20 + b]*hamil[b]; 
  } 
  double dH_dv1 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_2x3v_ser_p2_vel_sparse_idx[s]; 
    dH_dv1 += vst_2x3v_ser_p2_vel_vol_dv1[m*20 + b]*hamil[b]; 
  } 
  double dH_dv2 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_2x3v_ser_p2_vel_sparse_idx[s]; 
    dH_dv2 += vst_2x3v_ser_p2_vel_vol_dv2[m*20 + b]*hamil[b]; 
  } 
  if (hamil_pt_edge == -1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_2x3v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_2x3v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 8; ++a) P2 += vst_2x3v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_2[a]; 
    return P0*dH_dv0*dv10 + P1*dH_dv1*dv11 + P2*dH_dv2*dv12; 
  } 
  else if (hamil_pt_edge == 1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_2x3v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_2x3v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 8; ++a) P2 += vst_2x3v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_2[a]; 
    return P0*dH_dv0*dv10 + P1*dH_dv1*dv11 + P2*dH_dv2*dv12; 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 3; ++i) { 
    for (int m = 0; m < 27; ++m) alpha_quad[i*27 + m] += hamil_vel_sparse_alpha_quad_x_2x3v_ser_p2_node(i, m, hamil_pt_edge, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
