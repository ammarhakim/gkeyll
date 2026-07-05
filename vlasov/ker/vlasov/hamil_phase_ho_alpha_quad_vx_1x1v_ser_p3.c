#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_ser_p3.h> 
GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double G[4]; 
  for (int q = 0; q < 4; ++q) G[q] = 0.0; 
  for (int k = 0; k < 12; ++k) { 
    const int a = vst_1x1v_ser_p3_ho_ph_v0_cmap[k]; 
    const double ch = vst_1x1v_ser_p3_ho_ph_v0_coefr[k]*hamil[k]; 
    for (int j = 0; j < 1; ++j) G[j*4 + a] += vst_1x1v_ser_p3_ho_ph_v0_V[j*4 + vst_1x1v_ser_p3_ho_ph_v0_vrmap[k]]*ch; 
  } 
  for (int i = 0; i < 5; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      double dH_dx = 0.0; 
      for (int a = 0; a < 4; ++a) dH_dx += vst_1x1v_ser_p3_ho_ph_v0_CmD[i*4 + a]*G[j*4 + a]; 
      alpha_quad[i*1 + j] += -dx10*dH_dx; 
    } 
  } 
} 
