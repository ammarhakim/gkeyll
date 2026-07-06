#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double nc_hamil_phase_alpha_quad_vx_1x3v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double G[3]; 
  double Gd1[3]; 
  double Gd2[3]; 
  for (int a = 0; a < 3; ++a) { G[a] = 0.0; Gd1[a] = 0.0; Gd2[a] = 0.0; } 
  for (int k = 0; k < 48; ++k) { 
    const int a = vst_1x3v_ser_p2_ph_v0_cmap[k]; 
    G[a] += vst_1x3v_ser_p2_ph_v0_V[j*20 + vst_1x3v_ser_p2_ph_v0_vrmap[k]]*(vst_1x3v_ser_p2_ph_v0_coefr[k]*hamil[k]); 
    Gd1[a] += vst_1x3v_ser_p2_ph_v0_Vd1[j*13 + vst_1x3v_ser_p2_ph_v0_vrd1map[k]]*(vst_1x3v_ser_p2_ph_v0_dcoefr1[k]*hamil[k]); 
    Gd2[a] += vst_1x3v_ser_p2_ph_v0_Vd2[j*13 + vst_1x3v_ser_p2_ph_v0_vrd2map[k]]*(vst_1x3v_ser_p2_ph_v0_dcoefr2[k]*hamil[k]); 
  } 
  double dH_dv1 = 0.0; 
  for (int a = 0; a < 3; ++a) dH_dv1 += vst_1x3v_ser_p2_ph_v0_Cm[i*3 + a]*Gd1[a]; 
  double dH_dv2 = 0.0; 
  for (int a = 0; a < 3; ++a) dH_dv2 += vst_1x3v_ser_p2_ph_v0_Cm[i*3 + a]*Gd2[a]; 
  const double dx10 = 2.0/dxv[0]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  double px0_q = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    px0_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_x0[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 3; ++a) dH_dx0 += vst_1x3v_ser_p2_ph_v0_CmDx0[i*3 + a]*G[a]; 
  const double wx1 = w[1]; 
  const double wx2 = w[2]; 
  const double wx3 = w[3]; 
  const double dv0 = dxv[1]; 
  const double dv1 = dxv[2]; 
  const double dv2 = dxv[3]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double vt1 = wx1 - 0.5*dv0; 
  const double vt2 = wx2 + 0.5*dv1*vst_1x3v_ser_p2_vel_nodes_v0[j*2 + 0]; 
  const double vt3 = wx3 + 0.5*dv2*vst_1x3v_ser_p2_vel_nodes_v0[j*2 + 1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[27]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[30]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[36]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[39]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[45]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[48]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  double p3_q = 0.0; 
  double p4_q = 0.0; 
  double p6_q = 0.0; 
  double p7_q = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    p0_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_1[a]; 
    p3_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_3[a]; 
    p4_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_4[a]; 
    p6_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_6[a]; 
    p7_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_7[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p3_q + vt3*p6_q; 
  const double om13 = vt1*p1_q + vt2*p4_q + vt3*p7_q; 
  return -(px0_q*dH_dx0*dx10) + om12*dH_dv1*dv11 + om13*dH_dv2*dv12; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) alpha_quad[i*9 + j] += nc_hamil_phase_alpha_quad_vx_1x3v_ser_p2_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
