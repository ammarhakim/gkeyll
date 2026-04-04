#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
#include <gkyl_maxwell_kernels.h> 
EXTERN_C_BEG 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_1x_ser_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_1x_ser_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_2x_ser_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_2x_ser_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_2x_ser_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_2x_ser_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_3x_ser_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_3x_ser_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_z_3x_ser_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_3x_ser_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_3x_ser_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_z_3x_ser_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_1x_tensor_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_1x_tensor_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_2x_tensor_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_2x_tensor_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_2x_tensor_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_2x_tensor_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_3x_tensor_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_3x_tensor_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_z_3x_tensor_p1(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_3x_tensor_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_y_3x_tensor_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_nodal_to_modal_z_3x_tensor_p2(const double *dxv, 
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

EXTERN_C_END 
