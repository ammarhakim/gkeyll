#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_2x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_2x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfy_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_y_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_2x1v_ser_p2(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_2x1v_ser_p2(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_2x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfy_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_y_2x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_2x1v_ser_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_2x1v_ser_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_3x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_3x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_3x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_3x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfy_3x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_y_3x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfz_3x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_z_3x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_3x1v_ser_p2(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_3x1v_ser_p2(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_3x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfy_3x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_y_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfz_3x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_z_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_tensor_p2(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_1x1v_tensor_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_1x1v_tensor_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_1x1v_tensor_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_1x1v_tensor_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_2x1v_tensor_p2(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_2x1v_tensor_p2(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1); 
GKYL_CU_DH double vlasov_pkpm_surfx_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_x_2x1v_tensor_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfy_2x1v_tensor_p2(const double *w, const double *dxv, 
    const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_div_ppar_y_2x1v_tensor_p2(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar); 
GKYL_CU_DH double vlasov_pkpm_surfvpar_2x1v_tensor_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_2x1v_tensor_p2(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

EXTERN_C_END 
