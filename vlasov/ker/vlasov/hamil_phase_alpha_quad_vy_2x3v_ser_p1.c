#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p1.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x3v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double G[16]; 
  for (int q = 0; q < 16; ++q) G[q] = 0.0; 
  for (int k = 0; k < 32; ++k) { 
    const int a = vst_2x3v_ser_p1_ph_v1_cmap[k]; 
    const double ch = vst_2x3v_ser_p1_ph_v1_coefr[k]*hamil[k]; 
    for (int j = 0; j < 4; ++j) G[j*4 + a] += vst_2x3v_ser_p1_ph_v1_V[j*8 + vst_2x3v_ser_p1_ph_v1_vrmap[k]]*ch; 
  } 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double dH_dx = 0.0; 
      for (int a = 0; a < 4; ++a) dH_dx += vst_2x3v_ser_p1_ph_v1_CmD[i*4 + a]*G[j*4 + a]; 
      alpha_quad[i*4 + j] += -dx11*dH_dx; 
    } 
  } 
} 
