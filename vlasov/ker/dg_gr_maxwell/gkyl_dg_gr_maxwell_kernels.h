#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
#include <gkyl_maxwell_kernels.h> 
typedef struct { double chi, gamma, K_phi, K_psi; } gkyl_dg_gr_maxwell_inp; 
EXTERN_C_BEG 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_ser_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_ser_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_ser_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_ser_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_ser_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_ser_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_tensor_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_1x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_tensor_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_tensor_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_tensor_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gr_maxwell_geom_source_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, 
      const double *geom_factor_con_nodal, const double *fields, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 


GKYL_CU_DH void A_dQ_x_calc(const double lapse, const double shift_x, const double shift_y, const double shift_z, 
   const double h_xx, const double h_xy, const double h_xz, const double h_yy, 
   const double h_yz, const double h_zz, const double J_c, const double *lambda, const double *dQ, double* GKYL_RESTRICT A_dQ); 
GKYL_CU_DH void A_dQ_y_calc(const double lapse, const double shift_x, const double shift_y, const double shift_z, 
   const double h_xx, const double h_xy, const double h_xz, const double h_yy, 
   const double h_yz, const double h_zz, const double J_c, const double *lambda, const double *dQ, double* GKYL_RESTRICT A_dQ); 
GKYL_CU_DH void A_dQ_z_calc(const double lapse, const double shift_x, const double shift_y, const double shift_z, 
   const double h_xx, const double h_xy, const double h_xz, const double h_yy, 
   const double h_yz, const double h_zz, const double J_c, const double *lambda, const double *dQ, double* GKYL_RESTRICT A_dQ); 
EXTERN_C_END 
