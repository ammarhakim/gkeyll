#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x3v_ser_p2.h> 
GKYL_CU_DH void vlasov_boundary_ho_surfvz_2x3v_ser_p2_mode(int k, int edge, double dv12,
  const double *flux, double* GKYL_RESTRICT out) 
{ 
  const int a = vst_2x3v_ser_p2_ho_prj_v2_kamap[k]; 
  const int b = vst_2x3v_ser_p2_ho_prj_v2_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 16; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 16; ++j) { 
      t += vst_2x3v_ser_p2_ho_prj_v2_Vw[j*8 + b]*flux[512 + i*16 + j]; 
    } 
    g += vst_2x3v_ser_p2_ho_prj_v2_Cw[i*8 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_2x3v_ser_p2_ho_prj_v2_out_off[k]; q < vst_2x3v_ser_p2_ho_prj_v2_out_off[k+1]; ++q) { 
      out[vst_2x3v_ser_p2_ho_prj_v2_out_mode[q]] += dv12*vst_2x3v_ser_p2_ho_prj_v2_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x3v_ser_p2_ho_prj_v2_out_off[k]; q < vst_2x3v_ser_p2_ho_prj_v2_out_off[k+1]; ++q) { 
      out[vst_2x3v_ser_p2_ho_prj_v2_out_mode[q]] += dv12*vst_2x3v_ser_p2_ho_prj_v2_out_cl[q]*g; 
    } 
  } 
} 

GKYL_CU_DH double vlasov_boundary_ho_surfvz_2x3v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[4]; 
  for (int k = 0; k < 48; ++k) { 
    vlasov_boundary_ho_surfvz_2x3v_ser_p2_mode(k, edge, dv12, flux, out); 
  } 
  return 0.0;
} 
