#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
#include <gkyl_maxwell_kernels.h> 
EXTERN_C_BEG 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_ser_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_ser_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_ser_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out);
GKYL_CU_DH double outflow_flux_x_2x_ser_p1(const double *dxv, int edge,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
      const double *det_h_nodal, const double *field_con_skin,
      double* GKYL_RESTRICT conf_flux_surf);
GKYL_CU_DH double outflow_flux_x_2x_ser_p1_curved(const double *dxv, int edge,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
      const double *det_h_nodal, const double *field_con_skin,
      double* GKYL_RESTRICT conf_flux_surf);
// Tetrad-frame characteristic-outflow flux: paper eq. 57 wave decomposition
// in tetrad basis, outgoing-only projection, back-transform via M^{-1}.
GKYL_CU_DH double outflow_flux_x_2x_ser_p1_tetrad_roe(const double *dxv, int edge,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
      const double *det_h_nodal, const double *field_con_skin,
      double* GKYL_RESTRICT conf_flux_surf);
GKYL_CU_DH double lax_flux_x_2x_ser_p1_curved(const double *dxv, const int theta_pole,
      const double *h_ij_nodal, const double *J_c,
      const double *flux_l, const double *flux_r, const double *alpha_quad,
      const double *field_con_l, const double *field_con_r,
      double* GKYL_RESTRICT conf_flux_surf);
// Tetrad-frame Roe flux on the radial face: per-quad-node Roe wave decomposition
// in the locally flat tetrad basis (paper eq. 57 eigenvectors), with moving-
// interface correction tilde beta^x/alpha. Self-contained kernel; computes its
// own wave speeds internally (no alpha_quad input needed).
GKYL_CU_DH double lax_flux_x_2x_ser_p1_tetrad_roe(const double *dxv, const int theta_pole,
      const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *det_h,
      const double *flux_l, const double *flux_r,
      const double *field_con_l, const double *field_con_r,
      double* GKYL_RESTRICT conf_flux_surf);
GKYL_CU_DH double lax_flux_y_2x_ser_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out);
GKYL_CU_DH double outflow_flux_y_2x_ser_p1(const double *dxv, int edge,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
      const double *det_h_nodal, const double *field_con_skin,
      double* GKYL_RESTRICT conf_flux_surf);
GKYL_CU_DH double outflow_flux_y_2x_ser_p1_curved(const double *dxv, int edge,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
      const double *det_h_nodal, const double *field_con_skin,
      double* GKYL_RESTRICT conf_flux_surf);
GKYL_CU_DH double lax_flux_y_2x_ser_p1_curved(const double *dxv, const int theta_pole,
      const double *h_ij_nodal, const double *J_c,
      const double *flux_l, const double *flux_r, const double *alpha_quad,
      const double *field_con_l, const double *field_con_r,
      double* GKYL_RESTRICT conf_flux_surf);
// Tetrad-frame Roe flux on the theta-direction face: per-quad-node Roe wave
// decomposition with y-direction flat-Maxwell eigenvectors. theta_pole=1
// triggers a centered-flux-only reflective-like fallback (M^{-1} diverges at
// the coordinate singularity).
GKYL_CU_DH double lax_flux_y_2x_ser_p1_tetrad_roe(const double *dxv, const int theta_pole,
      const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *det_h,
      const double *flux_l, const double *flux_r,
      const double *field_con_l, const double *field_con_r,
      double* GKYL_RESTRICT conf_flux_surf);

// Tetrad-frame characteristic slope limiter (radial / x-direction face).
// In-place edit of the cell's r-slope mode q[1] for every conserved component;
// optionally zeros q[3] (cross slope) when any wave amplitude is clipped.
// Cell-average q[0] is preserved exactly. Reference geometry is the right-
// face nodal arrays of the cell (averaged over 2 GL nodes to a single face
// representative).
GKYL_CU_DH void slope_limiter_x_2x_ser_p1(double limiter_fac, const int theta_pole_lo, const int theta_pole_up,
      const double *lapse_nodal_face, const double *shift_nodal_face,
      const double *h_ij_nodal_face, const double *det_h_face,
      const double *field_con_l, const double *field_con_c, const double *field_con_r,
      double *field_con_out);

// Tetrad-frame characteristic slope limiter (theta / y-direction face).
// Limits q[2] (theta-slope) and optionally zeros q[3]. Theta-pole adjacency
// (theta_pole_lo or theta_pole_up != 0) zeros q[2] and q[3] outright; the
// inverse tetrad transform diverges at the pole and no physical gradient
// information is available across it.
GKYL_CU_DH void slope_limiter_y_2x_ser_p1(double limiter_fac, const int theta_pole_lo, const int theta_pole_up,
      const double *lapse_nodal_face, const double *shift_nodal_face,
      const double *h_ij_nodal_face, const double *det_h_face,
      const double *field_con_l, const double *field_con_c, const double *field_con_r,
      double *field_con_out);

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J);
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_ser_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_ser_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_ser_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_ser_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_ser_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_ser_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_ser_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_ser_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_tensor_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_1x_tensor_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_1x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_1x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_2x_tensor_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_2x_tensor_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_2x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_2x_tensor_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_2x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_2x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_2x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_tensor_p1(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_tensor_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfz_3x_tensor_p1(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J); 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_3x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields); 
GKYL_CU_DH void gr_maxwell_vol_3x_tensor_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
      const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_x_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_x_3x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_3x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfx_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_y_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_y_3x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_y_3x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
GKYL_CU_DH double gr_maxwell_surfy_3x_tensor_p2(const double *w, const double *dxv,
      const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lax_flux_z_3x_tensor_p2(const double *dxv, const int theta_pole,
      const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH double roe_flux_z_3x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) ; 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_z_3x_tensor_p2(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse, const double *shift, const double *h_ij, 
      const double *det_h, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 
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
