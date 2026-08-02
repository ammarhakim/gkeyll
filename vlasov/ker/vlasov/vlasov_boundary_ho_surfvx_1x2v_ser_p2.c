#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p2.h> 
GKYL_CU_DH void vlasov_boundary_ho_surfvx_1x2v_ser_p2_mode(int m, int edge, double dv10,
  const double *flux, double* GKYL_RESTRICT out) 
{ 
  const int a = vst_1x2v_ser_p2_ho_prj_v0_pm_a[m]; 
  const int b = vst_1x2v_ser_p2_ho_prj_v0_pm_b[m]; 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 4; ++j) { 
      t += vst_1x2v_ser_p2_ho_prj_v0_Vw[j*3 + b]*flux[0 + i*4 + j]; 
    } 
    g += vst_1x2v_ser_p2_ho_prj_v0_Cw[i*3 + a]*t; 
  } 
  if (edge == -1) { 
    out[m] += dv10*vst_1x2v_ser_p2_ho_prj_v0_pm_cr[m]*g; 
  } else { 
    out[m] += dv10*vst_1x2v_ser_p2_ho_prj_v0_pm_cl[m]*g; 
  } 
} 

GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x2v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  for (int k = 0; k < 8; ++k) { 
  const int a = vst_1x2v_ser_p2_ho_prj_v0_kamap[k]; 
  const int b = vst_1x2v_ser_p2_ho_prj_v0_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 4; ++j) { 
      t += vst_1x2v_ser_p2_ho_prj_v0_Vw[j*3 + b]*flux[0 + i*4 + j]; 
    } 
    g += vst_1x2v_ser_p2_ho_prj_v0_Cw[i*3 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_1x2v_ser_p2_ho_prj_v0_out_off[k]; q < vst_1x2v_ser_p2_ho_prj_v0_out_off[k+1]; ++q) { 
      out[vst_1x2v_ser_p2_ho_prj_v0_out_mode[q]] += dv10*vst_1x2v_ser_p2_ho_prj_v0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x2v_ser_p2_ho_prj_v0_out_off[k]; q < vst_1x2v_ser_p2_ho_prj_v0_out_off[k+1]; ++q) { 
      out[vst_1x2v_ser_p2_ho_prj_v0_out_mode[q]] += dv10*vst_1x2v_ser_p2_ho_prj_v0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
