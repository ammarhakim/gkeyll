#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1_node(int i, int j, const double *w, const double *dxv,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double G[8]; 
  double Gd0[8]; 
  double Gd1[8]; 
  for (int a = 0; a < 8; ++a) { G[a] = 0.0; Gd0[a] = 0.0; Gd1[a] = 0.0; } 
  for (int k = 0; k < 64; ++k) { 
    const int a = vst_3x3v_ser_p1_ph_v2_cmap[k]; 
    G[a] += vst_3x3v_ser_p1_ph_v2_V[j*8 + vst_3x3v_ser_p1_ph_v2_vrmap[k]]*(vst_3x3v_ser_p1_ph_v2_coefr[k]*hamil[k]); 
    Gd0[a] += vst_3x3v_ser_p1_ph_v2_Vd0[j*5 + vst_3x3v_ser_p1_ph_v2_vrd0map[k]]*(vst_3x3v_ser_p1_ph_v2_dcoefr0[k]*hamil[k]); 
    Gd1[a] += vst_3x3v_ser_p1_ph_v2_Vd1[j*5 + vst_3x3v_ser_p1_ph_v2_vrd1map[k]]*(vst_3x3v_ser_p1_ph_v2_dcoefr1[k]*hamil[k]); 
  } 
  double dH_dv0 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dv0 += vst_3x3v_ser_p1_ph_v2_Cm[i*8 + a]*Gd0[a]; 
  double dH_dv1 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dv1 += vst_3x3v_ser_p1_ph_v2_Cm[i*8 + a]*Gd1[a]; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dx12 = 2.0/dxv[2]; 
  const double *poisson_tensor_conf_x2 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_x5 = &poisson_tensor_conf[40]; 
  const double *poisson_tensor_conf_x8 = &poisson_tensor_conf[64]; 
  double px0_q = 0.0; 
  double px1_q = 0.0; 
  double px2_q = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    px0_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_x2[a]; 
    px1_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_x5[a]; 
    px2_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_x8[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dx0 += vst_3x3v_ser_p1_ph_v2_CmDx0[i*8 + a]*G[a]; 
  double dH_dx1 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dx1 += vst_3x3v_ser_p1_ph_v2_CmDx1[i*8 + a]*G[a]; 
  double dH_dx2 = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dx2 += vst_3x3v_ser_p1_ph_v2_CmDx2[i*8 + a]*G[a]; 
  const double wx1 = w[3]; 
  const double wx2 = w[4]; 
  const double wx3 = w[5]; 
  const double dv0 = dxv[3]; 
  const double dv1 = dxv[4]; 
  const double dv2 = dxv[5]; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double vt1 = wx1 + 0.5*dv0*vst_3x3v_ser_p1_vel_nodes_v2[j*2 + 0]; 
  const double vt2 = wx2 + 0.5*dv1*vst_3x3v_ser_p1_vel_nodes_v2[j*2 + 1]; 
  const double vt3 = wx3 - 0.5*dv2; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[80]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[88]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[104]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[112]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[128]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[136]; 
  double p1_q = 0.0; 
  double p2_q = 0.0; 
  double p4_q = 0.0; 
  double p5_q = 0.0; 
  double p7_q = 0.0; 
  double p8_q = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    p1_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
    p2_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_2[a]; 
    p4_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_4[a]; 
    p5_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_5[a]; 
    p7_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_7[a]; 
    p8_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_8[a]; 
  } 
  const double om13 = vt1*p1_q + vt2*p4_q + vt3*p7_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -(px0_q*dH_dx0*dx10 + px1_q*dH_dx1*dx11 + px2_q*dH_dx2*dx12) + -om13*dH_dv0*dv10 + -om23*dH_dv1*dv11; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1_node(i, j, w, dxv, poisson_tensor_conf, hamil); 
  } 
} 
