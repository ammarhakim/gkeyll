#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double G[3]; 
  for (int a = 0; a < 3; ++a) G[a] = 0.0; 
  for (int k = 0; k < 48; ++k) { 
    G[vst_1x3v_ser_p2_ho_ph_v0_cmap[k]] += vst_1x3v_ser_p2_ho_ph_v0_V[j*20 + vst_1x3v_ser_p2_ho_ph_v0_vrmap[k]]*(vst_1x3v_ser_p2_ho_ph_v0_coefr[k]*hamil[k]); 
  } 
  double dH_dx = 0.0; 
  for (int a = 0; a < 3; ++a) dH_dx += vst_1x3v_ser_p2_ho_ph_v0_CmD[i*3 + a]*G[a]; 
  return -dx10*dH_dx; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
