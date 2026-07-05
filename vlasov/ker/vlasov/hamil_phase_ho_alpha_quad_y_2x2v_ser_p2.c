#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[24]; 
  if (hamil_pt_edge == -1) { 
    double G0[48]; 
    for (int q = 0; q < 48; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      const int a = vst_2x2v_ser_p2_ho_ph_x1_cmap[k]; 
      const double ch = vst_2x2v_ser_p2_ho_ph_x1_dcoefr0[k]*hamil[k]; 
      for (int m = 0; m < 16; ++m) G0[m*3 + a] += vst_2x2v_ser_p2_ho_ph_x1_Vd0[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vrd0map[k]]*ch; 
    } 
    double G1[48]; 
    for (int q = 0; q < 48; ++q) G1[q] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      const int a = vst_2x2v_ser_p2_ho_ph_x1_cmap[k]; 
      const double ch = vst_2x2v_ser_p2_ho_ph_x1_dcoefr1[k]*hamil[k]; 
      for (int m = 0; m < 16; ++m) G1[m*3 + a] += vst_2x2v_ser_p2_ho_ph_x1_Vd1[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vrd1map[k]]*ch; 
    } 
    for (int i = 0; i < 4; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int m = 0; m < 16; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 3; ++a) dH0 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G0[m*3 + a]; 
        double dH1 = 0.0; 
        for (int a = 0; a < 3; ++a) dH1 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G1[m*3 + a]; 
        alpha_quad[i*16 + m] += P0*dH0*dv10 + P1*dH1*dv11; 
      } 
    } 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[48]; 
    for (int q = 0; q < 48; ++q) G0[q] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      const int a = vst_2x2v_ser_p2_ho_ph_x1_cmap[k]; 
      const double ch = vst_2x2v_ser_p2_ho_ph_x1_dcoefl0[k]*hamil[k]; 
      for (int m = 0; m < 16; ++m) G0[m*3 + a] += vst_2x2v_ser_p2_ho_ph_x1_Vd0[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vld0map[k]]*ch; 
    } 
    double G1[48]; 
    for (int q = 0; q < 48; ++q) G1[q] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      const int a = vst_2x2v_ser_p2_ho_ph_x1_cmap[k]; 
      const double ch = vst_2x2v_ser_p2_ho_ph_x1_dcoefl1[k]*hamil[k]; 
      for (int m = 0; m < 16; ++m) G1[m*3 + a] += vst_2x2v_ser_p2_ho_ph_x1_Vd1[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vld1map[k]]*ch; 
    } 
    for (int i = 0; i < 4; ++i) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int m = 0; m < 16; ++m) { 
        double dH0 = 0.0; 
        for (int a = 0; a < 3; ++a) dH0 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G0[m*3 + a]; 
        double dH1 = 0.0; 
        for (int a = 0; a < 3; ++a) dH1 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G1[m*3 + a]; 
        alpha_quad[i*16 + m] += P0*dH0*dv10 + P1*dH1*dv11; 
      } 
    } 
  } 

} 
