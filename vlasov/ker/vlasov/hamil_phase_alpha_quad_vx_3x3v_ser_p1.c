#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_3x3v_ser_p1.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_3x3v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double G[32]; 
  for (int q = 0; q < 32; ++q) G[q] = 0.0; 
  for (int k = 0; k < 64; ++k) { 
    const int a = vst_3x3v_ser_p1_ph_v0_cmap[k]; 
    const double ch = vst_3x3v_ser_p1_ph_v0_coefr[k]*hamil[k]; 
    for (int j = 0; j < 4; ++j) G[j*8 + a] += vst_3x3v_ser_p1_ph_v0_V[j*8 + vst_3x3v_ser_p1_ph_v0_vrmap[k]]*ch; 
  } 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double dH_dx = 0.0; 
      for (int a = 0; a < 8; ++a) dH_dx += vst_3x3v_ser_p1_ph_v0_CmD[i*8 + a]*G[j*8 + a]; 
      alpha_quad[i*4 + j] += -dx10*dH_dx; 
    } 
  } 
} 
