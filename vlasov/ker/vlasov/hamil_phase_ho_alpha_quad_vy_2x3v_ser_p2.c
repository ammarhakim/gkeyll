#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_2x3v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double G[128]; 
  for (int q = 0; q < 128; ++q) G[q] = 0.0; 
  for (int k = 0; k < 112; ++k) { 
    const int a = vst_2x3v_ser_p2_ho_ph_v1_cmap[k]; 
    const double ch = vst_2x3v_ser_p2_ho_ph_v1_coefr[k]*hamil[k]; 
    for (int j = 0; j < 16; ++j) G[j*8 + a] += vst_2x3v_ser_p2_ho_ph_v1_V[j*20 + vst_2x3v_ser_p2_ho_ph_v1_vrmap[k]]*ch; 
  } 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double dH_dx = 0.0; 
      for (int a = 0; a < 8; ++a) dH_dx += vst_2x3v_ser_p2_ho_ph_v1_CmD[i*8 + a]*G[j*8 + a]; 
      alpha_quad[i*16 + j] += -dx11*dH_dx; 
    } 
  } 
} 
