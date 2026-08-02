#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p2.h> 
GKYL_CU_DH void vlasov_surfx_1x2v_ser_p2_mode(int m, double dx10,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  const int a = vst_1x2v_ser_p2_prj_x0_pm_a[m]; 
  const int b = vst_1x2v_ser_p2_prj_x0_pm_b[m]; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int i = 0; i < 1; ++i) { 
    double t_l = 0.0; 
    double t_r = 0.0; 
    for (int j = 0; j < 9; ++j) { 
      t_l += vst_1x2v_ser_p2_prj_x0_Vw[j*8 + b]*Fhat_l_nodal[0 + i*9 + j]; 
      t_r += vst_1x2v_ser_p2_prj_x0_Vw[j*8 + b]*Fhat_r_nodal[0 + i*9 + j]; 
    } 
    g_l += vst_1x2v_ser_p2_prj_x0_Cw[i*1 + a]*t_l; 
    g_r += vst_1x2v_ser_p2_prj_x0_Cw[i*1 + a]*t_r; 
  } 
  out[m] += dx10*(vst_1x2v_ser_p2_prj_x0_pm_cl[m]*g_l + vst_1x2v_ser_p2_prj_x0_pm_cr[m]*g_r); 
} 

GKYL_CU_DH double vlasov_surfx_1x2v_ser_p2(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  for (int k = 0; k < 8; ++k) { 
  const int a = vst_1x2v_ser_p2_prj_x0_kamap[k]; 
  const int b = vst_1x2v_ser_p2_prj_x0_kbmap[k]; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int i = 0; i < 1; ++i) { 
    double t_l = 0.0; 
    double t_r = 0.0; 
    for (int j = 0; j < 9; ++j) { 
      t_l += vst_1x2v_ser_p2_prj_x0_Vw[j*8 + b]*Fhat_l_nodal[0 + i*9 + j]; 
      t_r += vst_1x2v_ser_p2_prj_x0_Vw[j*8 + b]*Fhat_r_nodal[0 + i*9 + j]; 
    } 
    g_l += vst_1x2v_ser_p2_prj_x0_Cw[i*1 + a]*t_l; 
    g_r += vst_1x2v_ser_p2_prj_x0_Cw[i*1 + a]*t_r; 
  } 
  for (int q = vst_1x2v_ser_p2_prj_x0_out_off[k]; q < vst_1x2v_ser_p2_prj_x0_out_off[k+1]; ++q) { 
    out[vst_1x2v_ser_p2_prj_x0_out_mode[q]] += dx10*(vst_1x2v_ser_p2_prj_x0_out_cl[q]*g_l + vst_1x2v_ser_p2_prj_x0_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
