#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p2.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double G[9]; 
  for (int q = 0; q < 9; ++q) G[q] = 0.0; 
  for (int k = 0; k < 20; ++k) { 
    const int a = vst_1x2v_ser_p2_ph_v0_cmap[k]; 
    const double ch = vst_1x2v_ser_p2_ph_v0_coefr[k]*hamil[k]; 
    for (int j = 0; j < 3; ++j) G[j*3 + a] += vst_1x2v_ser_p2_ph_v0_V[j*8 + vst_1x2v_ser_p2_ph_v0_vrmap[k]]*ch; 
  } 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double dH_dx = 0.0; 
      for (int a = 0; a < 3; ++a) dH_dx += vst_1x2v_ser_p2_ph_v0_CmD[i*3 + a]*G[j*3 + a]; 
      alpha_quad[i*3 + j] += -dx10*dH_dx; 
    } 
  } 
} 
