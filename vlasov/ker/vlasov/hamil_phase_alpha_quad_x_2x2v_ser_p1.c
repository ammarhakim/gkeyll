#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p1.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[4]; 
  if (hamil_pt_edge == -1) { 
    double G0[8]; 
    for (int q = 0; q < 8; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 16; ++k) { 
      const int a = vst_2x2v_ser_p1_ph_x0_cmap[k]; 
      const double ch = vst_2x2v_ser_p1_ph_x0_dcoefr0[k]*hamil[k]; 
      for (int m = 0; m < 4; ++m) G0[m*2 + a] += vst_2x2v_ser_p1_ph_x0_Vd0[m*5 + vst_2x2v_ser_p1_ph_x0_vrd0map[k]]*ch; 
    } 
    double G1[8]; 
    for (int q = 0; q < 8; ++q) G1[q] = 0.0; 
    for (int k = 0; k < 16; ++k) { 
      const int a = vst_2x2v_ser_p1_ph_x0_cmap[k]; 
      const double ch = vst_2x2v_ser_p1_ph_x0_dcoefr1[k]*hamil[k]; 
      for (int m = 0; m < 4; ++m) G1[m*2 + a] += vst_2x2v_ser_p1_ph_x0_Vd1[m*5 + vst_2x2v_ser_p1_ph_x0_vrd1map[k]]*ch; 
    } 
    for (int i = 0; i < 2; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x2v_ser_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x2v_ser_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_1[a]; 
      for (int m = 0; m < 4; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 2; ++a) dH0 += vst_2x2v_ser_p1_ph_x0_Cm[i*2 + a]*G0[m*2 + a]; 
        double dH1 = 0.0; 
        for (int a = 0; a < 2; ++a) dH1 += vst_2x2v_ser_p1_ph_x0_Cm[i*2 + a]*G1[m*2 + a]; 
        alpha_quad[i*4 + m] += P0*dH0*dv10 + P1*dH1*dv11; 
      } 
    } 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[8]; 
    for (int q = 0; q < 8; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 16; ++k) { 
      const int a = vst_2x2v_ser_p1_ph_x0_cmap[k]; 
      const double ch = vst_2x2v_ser_p1_ph_x0_dcoefl0[k]*hamil[k]; 
      for (int m = 0; m < 4; ++m) G0[m*2 + a] += vst_2x2v_ser_p1_ph_x0_Vd0[m*5 + vst_2x2v_ser_p1_ph_x0_vld0map[k]]*ch; 
    } 
    double G1[8]; 
    for (int q = 0; q < 8; ++q) G1[q] = 0.0; 
    for (int k = 0; k < 16; ++k) { 
      const int a = vst_2x2v_ser_p1_ph_x0_cmap[k]; 
      const double ch = vst_2x2v_ser_p1_ph_x0_dcoefl1[k]*hamil[k]; 
      for (int m = 0; m < 4; ++m) G1[m*2 + a] += vst_2x2v_ser_p1_ph_x0_Vd1[m*5 + vst_2x2v_ser_p1_ph_x0_vld1map[k]]*ch; 
    } 
    for (int i = 0; i < 2; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x2v_ser_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x2v_ser_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_1[a]; 
      for (int m = 0; m < 4; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 2; ++a) dH0 += vst_2x2v_ser_p1_ph_x0_Cm[i*2 + a]*G0[m*2 + a]; 
        double dH1 = 0.0; 
        for (int a = 0; a < 2; ++a) dH1 += vst_2x2v_ser_p1_ph_x0_Cm[i*2 + a]*G1[m*2 + a]; 
        alpha_quad[i*4 + m] += P0*dH0*dv10 + P1*dH1*dv11; 
      } 
    } 
  } 

} 
