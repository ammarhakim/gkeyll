#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double hamil_phase_ho_alpha_quad_vy_2x3v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double G[8]; 
  for (int a = 0; a < 8; ++a) G[a] = 0.0; 
  for (int k = 0; k < 112; ++k) { 
    G[vst_2x3v_ser_p2_ho_ph_v1_cmap[k]] += vst_2x3v_ser_p2_ho_ph_v1_V[j*20 + vst_2x3v_ser_p2_ho_ph_v1_vrmap[k]]*(vst_2x3v_ser_p2_ho_ph_v1_coefr[k]*hamil[k]); 
  } 
  double dH_dx = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dx += vst_2x3v_ser_p2_ho_ph_v1_CmD[i*8 + a]*G[a]; 
  return -dx11*dH_dx; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_2x3v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += hamil_phase_ho_alpha_quad_vy_2x3v_ser_p2_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
