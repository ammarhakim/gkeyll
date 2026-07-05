#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[24]; 
  double dH_dv0[16]; 
  for (int m = 0; m < 16; ++m) { 
    dH_dv0[m] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dv0[m] += vst_2x2v_ser_p2_ho_vel_vol_dv0[m*8 + b]*hamil[b]; 
  } 
  double dH_dv1[16]; 
  for (int m = 0; m < 16; ++m) { 
    dH_dv1[m] = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dv1[m] += vst_2x2v_ser_p2_ho_vel_vol_dv1[m*8 + b]*hamil[b]; 
  } 
  if (hamil_pt_edge == -1) { 
    for (int i = 0; i < 4; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int m = 0; m < 16; ++m) alpha_quad[i*16 + m] += P0*dH_dv0[m]*dv10 + P1*dH_dv1[m]*dv11; 
    } 
  } 
  else if (hamil_pt_edge == 1) { 
    for (int i = 0; i < 4; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int m = 0; m < 16; ++m) alpha_quad[i*16 + m] += P0*dH_dv0[m]*dv10 + P1*dH_dv1[m]*dv11; 
    } 
  } 

} 
