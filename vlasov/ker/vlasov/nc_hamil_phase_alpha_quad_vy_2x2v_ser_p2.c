#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH double nc_hamil_phase_alpha_quad_vy_2x2v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double G[8]; 
  double Gd0[8]; 
  for (int a = 0; a < 8; ++a) { G[a] = 0.0; Gd0[a] = 0.0; } 
  for (int k = 0; k < 48; ++k) { 
    const int a = vst_2x2v_ser_p2_ph_v1_cmap[k]; 
    G[a] += vst_2x2v_ser_p2_ph_v1_V[j*8 + vst_2x2v_ser_p2_ph_v1_vrmap[k]]*(vst_2x2v_ser_p2_ph_v1_coefr[k]*hamil[k]); 
    Gd0[a] += vst_2x2v_ser_p2_ph_v1_Vd0[j*6 + vst_2x2v_ser_p2_ph_v1_vrd0map[k]]*(vst_2x2v_ser_p2_ph_v1_dcoefr0[k]*hamil[k]); 
  } 
  double dH_dv0 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dv0 += vst_2x2v_ser_p2_ph_v1_Cm[i*8 + a]*Gd0[a]; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_x3 = &poisson_tensor_conf[24]; 
  double px0_q = 0.0; 
  double px1_q = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    px0_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_x1[a]; 
    px1_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_x3[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dx0 += vst_2x2v_ser_p2_ph_v1_CmDx0[i*8 + a]*G[a]; 
  double dH_dx1 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dx1 += vst_2x2v_ser_p2_ph_v1_CmDx1[i*8 + a]*G[a]; 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv10 = 2.0/dxv[2]; 
  const double vt1 = wx1 + 0.5*dv0*vst_2x2v_ser_p2_vel_nodes_v1[j*1 + 0]; 
  const double vt2 = wx2 - 0.5*dv1; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[32]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[40]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    p0_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p1_q; 
  return -(px0_q*dH_dx0*dx10 + px1_q*dH_dx1*dx11) + -om12*dH_dv0*dv10; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 9; ++i) { 
    for (int j = 0; j < 3; ++j) alpha_quad[i*3 + j] += nc_hamil_phase_alpha_quad_vy_2x2v_ser_p2_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
