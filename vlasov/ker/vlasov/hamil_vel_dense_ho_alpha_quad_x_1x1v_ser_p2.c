#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_ser_p2.h> 
GKYL_CU_DH void hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[1]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  double dH_dv0[4]; 
  for (int m = 0; m < 4; ++m) { 
    dH_dv0[m] = 0.0; 
    for (int b = 0; b < 3; ++b) dH_dv0[m] += vst_1x1v_ser_p2_ho_vel_vol_dv0[m*3 + b]*hamil[b]; 
  } 
  if (hamil_pt_edge == -1) { 
    for (int i = 0; i < 1; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 3; ++a) P0 += vst_1x1v_ser_p2_ho_confsurf_x0_ev_r[i*3 + a]*poisson_tensor_conf_0[a]; 
      for (int m = 0; m < 4; ++m) alpha_quad[i*4 + m] += P0*dH_dv0[m]*dv10; 
    } 
  } 
  else if (hamil_pt_edge == 1) { 
    for (int i = 0; i < 1; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 3; ++a) P0 += vst_1x1v_ser_p2_ho_confsurf_x0_ev_l[i*3 + a]*poisson_tensor_conf_0[a]; 
      for (int m = 0; m < 4; ++m) alpha_quad[i*4 + m] += P0*dH_dv0[m]*dv10; 
    } 
  } 

} 
