#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_x_2x3v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[16]; 
  if (hamil_pt_edge == -1) { 
    double G0[81]; 
    for (int q = 0; q < 81; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 112; ++k) { 
      const int a = vst_2x3v_ser_p2_ph_x0_cmap[k]; 
      const double ch = vst_2x3v_ser_p2_ph_x0_dcoefr0[k]*hamil[k]; 
      for (int m = 0; m < 27; ++m) G0[m*3 + a] += vst_2x3v_ser_p2_ph_x0_Vd0[m*29 + vst_2x3v_ser_p2_ph_x0_vrd0map[k]]*ch; 
    } 
    double G1[81]; 
    for (int q = 0; q < 81; ++q) G1[q] = 0.0; 
    for (int k = 0; k < 112; ++k) { 
      const int a = vst_2x3v_ser_p2_ph_x0_cmap[k]; 
      const double ch = vst_2x3v_ser_p2_ph_x0_dcoefr1[k]*hamil[k]; 
      for (int m = 0; m < 27; ++m) G1[m*3 + a] += vst_2x3v_ser_p2_ph_x0_Vd1[m*29 + vst_2x3v_ser_p2_ph_x0_vrd1map[k]]*ch; 
    } 
    double G2[81]; 
    for (int q = 0; q < 81; ++q) G2[q] = 0.0; 
    for (int k = 0; k < 112; ++k) { 
      const int a = vst_2x3v_ser_p2_ph_x0_cmap[k]; 
      const double ch = vst_2x3v_ser_p2_ph_x0_dcoefr2[k]*hamil[k]; 
      for (int m = 0; m < 27; ++m) G2[m*3 + a] += vst_2x3v_ser_p2_ph_x0_Vd2[m*29 + vst_2x3v_ser_p2_ph_x0_vrd2map[k]]*ch; 
    } 
    for (int i = 0; i < 3; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x3v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x3v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 8; ++a) P2 += vst_2x3v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_2[a]; 
      for (int m = 0; m < 27; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 3; ++a) dH0 += vst_2x3v_ser_p2_ph_x0_Cm[i*3 + a]*G0[m*3 + a]; 
        double dH1 = 0.0; 
        for (int a = 0; a < 3; ++a) dH1 += vst_2x3v_ser_p2_ph_x0_Cm[i*3 + a]*G1[m*3 + a]; 
        double dH2 = 0.0; 
        for (int a = 0; a < 3; ++a) dH2 += vst_2x3v_ser_p2_ph_x0_Cm[i*3 + a]*G2[m*3 + a]; 
        alpha_quad[i*27 + m] += P0*dH0*dv10 + P1*dH1*dv11 + P2*dH2*dv12; 
      } 
    } 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[81]; 
    for (int q = 0; q < 81; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 112; ++k) { 
      const int a = vst_2x3v_ser_p2_ph_x0_cmap[k]; 
      const double ch = vst_2x3v_ser_p2_ph_x0_dcoefl0[k]*hamil[k]; 
      for (int m = 0; m < 27; ++m) G0[m*3 + a] += vst_2x3v_ser_p2_ph_x0_Vd0[m*29 + vst_2x3v_ser_p2_ph_x0_vld0map[k]]*ch; 
    } 
    double G1[81]; 
    for (int q = 0; q < 81; ++q) G1[q] = 0.0; 
    for (int k = 0; k < 112; ++k) { 
      const int a = vst_2x3v_ser_p2_ph_x0_cmap[k]; 
      const double ch = vst_2x3v_ser_p2_ph_x0_dcoefl1[k]*hamil[k]; 
      for (int m = 0; m < 27; ++m) G1[m*3 + a] += vst_2x3v_ser_p2_ph_x0_Vd1[m*29 + vst_2x3v_ser_p2_ph_x0_vld1map[k]]*ch; 
    } 
    double G2[81]; 
    for (int q = 0; q < 81; ++q) G2[q] = 0.0; 
    for (int k = 0; k < 112; ++k) { 
      const int a = vst_2x3v_ser_p2_ph_x0_cmap[k]; 
      const double ch = vst_2x3v_ser_p2_ph_x0_dcoefl2[k]*hamil[k]; 
      for (int m = 0; m < 27; ++m) G2[m*3 + a] += vst_2x3v_ser_p2_ph_x0_Vd2[m*29 + vst_2x3v_ser_p2_ph_x0_vld2map[k]]*ch; 
    } 
    for (int i = 0; i < 3; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x3v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x3v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 8; ++a) P2 += vst_2x3v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_2[a]; 
      for (int m = 0; m < 27; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 3; ++a) dH0 += vst_2x3v_ser_p2_ph_x0_Cm[i*3 + a]*G0[m*3 + a]; 
        double dH1 = 0.0; 
        for (int a = 0; a < 3; ++a) dH1 += vst_2x3v_ser_p2_ph_x0_Cm[i*3 + a]*G1[m*3 + a]; 
        double dH2 = 0.0; 
        for (int a = 0; a < 3; ++a) dH2 += vst_2x3v_ser_p2_ph_x0_Cm[i*3 + a]*G2[m*3 + a]; 
        alpha_quad[i*27 + m] += P0*dH0*dv10 + P1*dH1*dv11 + P2*dH2*dv12; 
      } 
    } 
  } 

} 
