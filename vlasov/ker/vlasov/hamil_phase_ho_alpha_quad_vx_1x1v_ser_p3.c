#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x1v_ser_p3.h> 
GKYL_CU_DH double hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double G[4]; 
  for (int a = 0; a < 4; ++a) G[a] = 0.0; 
  for (int k = 0; k < 12; ++k) { 
    G[vst_1x1v_ser_p3_ho_ph_v0_cmap[k]] += vst_1x1v_ser_p3_ho_ph_v0_V[j*4 + vst_1x1v_ser_p3_ho_ph_v0_vrmap[k]]*(vst_1x1v_ser_p3_ho_ph_v0_coefr[k]*hamil[k]); 
  } 
  double dH_dx = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dx += vst_1x1v_ser_p3_ho_ph_v0_CmD[i*4 + a]*G[a]; 
  return -dx10*dH_dx; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 5; ++i) { 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
