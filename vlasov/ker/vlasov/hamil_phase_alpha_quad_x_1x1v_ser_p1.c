#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_ser_p1.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[1]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  if (hamil_pt_edge == -1) { 
    double G0[2]; 
    for (int q = 0; q < 2; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 4; ++k) { 
      const int a = vst_1x1v_ser_p1_ph_x0_cmap[k]; 
      const double ch = vst_1x1v_ser_p1_ph_x0_dcoefr0[k]*hamil[k]; 
      for (int m = 0; m < 2; ++m) G0[m*1 + a] += vst_1x1v_ser_p1_ph_x0_Vd0[m*3 + vst_1x1v_ser_p1_ph_x0_vrd0map[k]]*ch; 
    } 
    for (int i = 0; i < 1; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x1v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
      for (int m = 0; m < 2; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 1; ++a) dH0 += vst_1x1v_ser_p1_ph_x0_Cm[i*1 + a]*G0[m*1 + a]; 
        alpha_quad[i*2 + m] += P0*dH0*dv10; 
      } 
    } 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[2]; 
    for (int q = 0; q < 2; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 4; ++k) { 
      const int a = vst_1x1v_ser_p1_ph_x0_cmap[k]; 
      const double ch = vst_1x1v_ser_p1_ph_x0_dcoefl0[k]*hamil[k]; 
      for (int m = 0; m < 2; ++m) G0[m*1 + a] += vst_1x1v_ser_p1_ph_x0_Vd0[m*3 + vst_1x1v_ser_p1_ph_x0_vld0map[k]]*ch; 
    } 
    for (int i = 0; i < 1; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x1v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
      for (int m = 0; m < 2; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 1; ++a) dH0 += vst_1x1v_ser_p1_ph_x0_Cm[i*1 + a]*G0[m*1 + a]; 
        alpha_quad[i*2 + m] += P0*dH0*dv10; 
      } 
    } 
  } 

} 
