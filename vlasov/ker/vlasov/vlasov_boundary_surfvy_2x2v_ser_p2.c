#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void vlasov_boundary_surfvy_2x2v_ser_p2_mode(int m, int edge, double dv11,
  const double *flux, double* GKYL_RESTRICT out) 
{ 
  const int a = vst_2x2v_ser_p2_prj_v1_pm_a[m]; 
  const int b = vst_2x2v_ser_p2_prj_v1_pm_b[m]; 
  double g = 0.0; 
  for (int i = 0; i < 9; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 3; ++j) { 
      t += vst_2x2v_ser_p2_prj_v1_Vw[j*3 + b]*flux[27 + i*3 + j]; 
    } 
    g += vst_2x2v_ser_p2_prj_v1_Cw[i*8 + a]*t; 
  } 
  if (edge == -1) { 
    out[m] += dv11*vst_2x2v_ser_p2_prj_v1_pm_cr[m]*g; 
  } else { 
    out[m] += dv11*vst_2x2v_ser_p2_prj_v1_pm_cl[m]*g; 
  } 
} 

GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[3]; 
  for (int k = 0; k < 20; ++k) { 
  const int a = vst_2x2v_ser_p2_prj_v1_kamap[k]; 
  const int b = vst_2x2v_ser_p2_prj_v1_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 9; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 3; ++j) { 
      t += vst_2x2v_ser_p2_prj_v1_Vw[j*3 + b]*flux[27 + i*3 + j]; 
    } 
    g += vst_2x2v_ser_p2_prj_v1_Cw[i*8 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_2x2v_ser_p2_prj_v1_out_off[k]; q < vst_2x2v_ser_p2_prj_v1_out_off[k+1]; ++q) { 
      out[vst_2x2v_ser_p2_prj_v1_out_mode[q]] += dv11*vst_2x2v_ser_p2_prj_v1_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x2v_ser_p2_prj_v1_out_off[k]; q < vst_2x2v_ser_p2_prj_v1_out_off[k+1]; ++q) { 
      out[vst_2x2v_ser_p2_prj_v1_out_mode[q]] += dv11*vst_2x2v_ser_p2_prj_v1_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
