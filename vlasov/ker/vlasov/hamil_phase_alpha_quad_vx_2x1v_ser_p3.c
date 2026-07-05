#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x1v_ser_p3.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x1v_ser_p3(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double G[12]; 
  for (int q = 0; q < 12; ++q) G[q] = 0.0; 
  for (int k = 0; k < 32; ++k) { 
    const int a = vst_2x1v_ser_p3_ph_v0_cmap[k]; 
    const double ch = vst_2x1v_ser_p3_ph_v0_coefr[k]*hamil[k]; 
    for (int j = 0; j < 1; ++j) G[j*12 + a] += vst_2x1v_ser_p3_ph_v0_V[j*4 + vst_2x1v_ser_p3_ph_v0_vrmap[k]]*ch; 
  } 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      double dH_dx = 0.0; 
      for (int a = 0; a < 12; ++a) dH_dx += vst_2x1v_ser_p3_ph_v0_CmD[i*12 + a]*G[j*12 + a]; 
      alpha_quad[i*1 + j] += -dx10*dH_dx; 
    } 
  } 
} 
