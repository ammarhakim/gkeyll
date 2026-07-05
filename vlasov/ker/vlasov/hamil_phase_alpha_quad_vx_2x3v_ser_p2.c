#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x3v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double G[72]; 
  for (int q = 0; q < 72; ++q) G[q] = 0.0; 
  for (int k = 0; k < 112; ++k) { 
    const int a = vst_2x3v_ser_p2_ph_v0_cmap[k]; 
    const double ch = vst_2x3v_ser_p2_ph_v0_coefr[k]*hamil[k]; 
    for (int j = 0; j < 9; ++j) G[j*8 + a] += vst_2x3v_ser_p2_ph_v0_V[j*20 + vst_2x3v_ser_p2_ph_v0_vrmap[k]]*ch; 
  } 
  for (int i = 0; i < 9; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double dH_dx = 0.0; 
      for (int a = 0; a < 8; ++a) dH_dx += vst_2x3v_ser_p2_ph_v0_CmD[i*8 + a]*G[j*8 + a]; 
      alpha_quad[i*9 + j] += -dx10*dH_dx; 
    } 
  } 
} 
