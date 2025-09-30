#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_1x1v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_1x1v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_1x1v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_1x1v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_1x1v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_1x1v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_1x1v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_1x1v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_x_1x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node32_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node32_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node33_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node33_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node34_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node34_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node35_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node35_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node36_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node36_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node37_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node37_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node38_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node38_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node39_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node39_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node40_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node40_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node41_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node41_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node42_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node42_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node43_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node43_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node44_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node44_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node45_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node45_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node46_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node46_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node47_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node47_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node48_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node48_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node49_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node49_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node50_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node50_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node51_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node51_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node52_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node52_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node53_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node53_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node54_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node54_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node55_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node55_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node56_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node56_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node57_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node57_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node58_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node58_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node59_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node59_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node60_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node60_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node61_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node61_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node62_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node62_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node63_x_1x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node63_x_1x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_x_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_x_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_y_2x2v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_y_2x2v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node32_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node32_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node33_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node33_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node34_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node34_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node35_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node35_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node36_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node36_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node37_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node37_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node38_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node38_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node39_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node39_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node40_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node40_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node41_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node41_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node42_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node42_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node43_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node43_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node44_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node44_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node45_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node45_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node46_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node46_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node47_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node47_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node48_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node48_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node49_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node49_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node50_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node50_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node51_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node51_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node52_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node52_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node53_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node53_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node54_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node54_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node55_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node55_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node56_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node56_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node57_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node57_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node58_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node58_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node59_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node59_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node60_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node60_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node61_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node61_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node62_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node62_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node63_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node63_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node64_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node64_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node65_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node65_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node66_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node66_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node67_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node67_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node68_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node68_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node69_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node69_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node70_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node70_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node71_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node71_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node72_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node72_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node73_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node73_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node74_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node74_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node75_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node75_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node76_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node76_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node77_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node77_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node78_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node78_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node79_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node79_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node80_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node80_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node81_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node81_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node82_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node82_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node83_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node83_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node84_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node84_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node85_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node85_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node86_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node86_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node87_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node87_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node88_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node88_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node89_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node89_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node90_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node90_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node91_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node91_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node92_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node92_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node93_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node93_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node94_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node94_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node95_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node95_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node96_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node96_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node97_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node97_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node98_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node98_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node99_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node99_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node100_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node100_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node101_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node101_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node102_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node102_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node103_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node103_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node104_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node104_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node105_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node105_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node106_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node106_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node107_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node107_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node108_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node108_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node109_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node109_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node110_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node110_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node111_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node111_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node112_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node112_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node113_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node113_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node114_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node114_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node115_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node115_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node116_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node116_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node117_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node117_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node118_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node118_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node119_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node119_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node120_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node120_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node121_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node121_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node122_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node122_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node123_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node123_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node124_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node124_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node125_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node125_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node126_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node126_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node127_x_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node127_x_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node32_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node32_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node33_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node33_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node34_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node34_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node35_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node35_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node36_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node36_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node37_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node37_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node38_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node38_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node39_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node39_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node40_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node40_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node41_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node41_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node42_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node42_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node43_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node43_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node44_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node44_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node45_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node45_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node46_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node46_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node47_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node47_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node48_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node48_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node49_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node49_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node50_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node50_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node51_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node51_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node52_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node52_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node53_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node53_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node54_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node54_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node55_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node55_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node56_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node56_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node57_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node57_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node58_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node58_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node59_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node59_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node60_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node60_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node61_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node61_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node62_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node62_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node63_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node63_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node64_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node64_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node65_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node65_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node66_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node66_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node67_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node67_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node68_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node68_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node69_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node69_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node70_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node70_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node71_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node71_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node72_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node72_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node73_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node73_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node74_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node74_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node75_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node75_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node76_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node76_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node77_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node77_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node78_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node78_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node79_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node79_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node80_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node80_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node81_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node81_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node82_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node82_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node83_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node83_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node84_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node84_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node85_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node85_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node86_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node86_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node87_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node87_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node88_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node88_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node89_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node89_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node90_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node90_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node91_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node91_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node92_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node92_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node93_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node93_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node94_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node94_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node95_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node95_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node96_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node96_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node97_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node97_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node98_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node98_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node99_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node99_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node100_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node100_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node101_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node101_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node102_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node102_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node103_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node103_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node104_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node104_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node105_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node105_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node106_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node106_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node107_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node107_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node108_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node108_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node109_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node109_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node110_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node110_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node111_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node111_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node112_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node112_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node113_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node113_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node114_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node114_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node115_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node115_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node116_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node116_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node117_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node117_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node118_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node118_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node119_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node119_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node120_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node120_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node121_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node121_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node122_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node122_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node123_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node123_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node124_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node124_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node125_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node125_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node126_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node126_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node127_y_2x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node127_y_2x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node32_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node32_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node33_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node33_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node34_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node34_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node35_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node35_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node36_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node36_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node37_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node37_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node38_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node38_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node39_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node39_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node40_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node40_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node41_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node41_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node42_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node42_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node43_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node43_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node44_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node44_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node45_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node45_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node46_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node46_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node47_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node47_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node48_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node48_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node49_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node49_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node50_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node50_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node51_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node51_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node52_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node52_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node53_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node53_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node54_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node54_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node55_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node55_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node56_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node56_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node57_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node57_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node58_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node58_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node59_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node59_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node60_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node60_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node61_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node61_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node62_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node62_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node63_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node63_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node64_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node64_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node65_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node65_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node66_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node66_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node67_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node67_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node68_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node68_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node69_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node69_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node70_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node70_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node71_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node71_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node72_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node72_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node73_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node73_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node74_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node74_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node75_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node75_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node76_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node76_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node77_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node77_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node78_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node78_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node79_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node79_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node80_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node80_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node81_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node81_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node82_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node82_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node83_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node83_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node84_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node84_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node85_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node85_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node86_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node86_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node87_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node87_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node88_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node88_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node89_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node89_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node90_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node90_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node91_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node91_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node92_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node92_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node93_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node93_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node94_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node94_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node95_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node95_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node96_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node96_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node97_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node97_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node98_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node98_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node99_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node99_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node100_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node100_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node101_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node101_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node102_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node102_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node103_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node103_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node104_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node104_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node105_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node105_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node106_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node106_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node107_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node107_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node108_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node108_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node109_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node109_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node110_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node110_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node111_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node111_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node112_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node112_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node113_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node113_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node114_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node114_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node115_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node115_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node116_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node116_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node117_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node117_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node118_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node118_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node119_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node119_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node120_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node120_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node121_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node121_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node122_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node122_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node123_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node123_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node124_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node124_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node125_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node125_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node126_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node126_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node127_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node127_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node128_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node128_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node129_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node129_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node130_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node130_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node131_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node131_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node132_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node132_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node133_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node133_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node134_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node134_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node135_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node135_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node136_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node136_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node137_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node137_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node138_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node138_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node139_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node139_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node140_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node140_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node141_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node141_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node142_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node142_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node143_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node143_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node144_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node144_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node145_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node145_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node146_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node146_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node147_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node147_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node148_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node148_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node149_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node149_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node150_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node150_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node151_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node151_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node152_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node152_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node153_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node153_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node154_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node154_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node155_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node155_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node156_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node156_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node157_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node157_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node158_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node158_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node159_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node159_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node160_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node160_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node161_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node161_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node162_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node162_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node163_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node163_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node164_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node164_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node165_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node165_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node166_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node166_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node167_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node167_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node168_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node168_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node169_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node169_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node170_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node170_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node171_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node171_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node172_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node172_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node173_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node173_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node174_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node174_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node175_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node175_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node176_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node176_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node177_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node177_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node178_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node178_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node179_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node179_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node180_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node180_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node181_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node181_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node182_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node182_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node183_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node183_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node184_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node184_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node185_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node185_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node186_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node186_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node187_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node187_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node188_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node188_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node189_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node189_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node190_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node190_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node191_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node191_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node192_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node192_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node193_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node193_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node194_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node194_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node195_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node195_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node196_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node196_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node197_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node197_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node198_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node198_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node199_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node199_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node200_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node200_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node201_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node201_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node202_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node202_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node203_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node203_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node204_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node204_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node205_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node205_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node206_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node206_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node207_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node207_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node208_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node208_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node209_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node209_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node210_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node210_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node211_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node211_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node212_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node212_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node213_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node213_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node214_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node214_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node215_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node215_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node216_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node216_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node217_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node217_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node218_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node218_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node219_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node219_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node220_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node220_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node221_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node221_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node222_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node222_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node223_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node223_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node224_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node224_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node225_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node225_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node226_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node226_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node227_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node227_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node228_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node228_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node229_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node229_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node230_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node230_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node231_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node231_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node232_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node232_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node233_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node233_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node234_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node234_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node235_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node235_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node236_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node236_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node237_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node237_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node238_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node238_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node239_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node239_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node240_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node240_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node241_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node241_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node242_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node242_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node243_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node243_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node244_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node244_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node245_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node245_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node246_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node246_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node247_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node247_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node248_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node248_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node249_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node249_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node250_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node250_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node251_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node251_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node252_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node252_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node253_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node253_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node254_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node254_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node255_x_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node255_x_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node32_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node32_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node33_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node33_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node34_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node34_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node35_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node35_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node36_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node36_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node37_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node37_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node38_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node38_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node39_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node39_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node40_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node40_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node41_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node41_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node42_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node42_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node43_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node43_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node44_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node44_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node45_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node45_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node46_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node46_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node47_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node47_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node48_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node48_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node49_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node49_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node50_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node50_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node51_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node51_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node52_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node52_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node53_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node53_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node54_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node54_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node55_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node55_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node56_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node56_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node57_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node57_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node58_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node58_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node59_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node59_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node60_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node60_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node61_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node61_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node62_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node62_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node63_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node63_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node64_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node64_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node65_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node65_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node66_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node66_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node67_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node67_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node68_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node68_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node69_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node69_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node70_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node70_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node71_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node71_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node72_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node72_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node73_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node73_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node74_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node74_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node75_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node75_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node76_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node76_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node77_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node77_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node78_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node78_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node79_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node79_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node80_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node80_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node81_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node81_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node82_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node82_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node83_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node83_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node84_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node84_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node85_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node85_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node86_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node86_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node87_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node87_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node88_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node88_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node89_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node89_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node90_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node90_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node91_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node91_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node92_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node92_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node93_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node93_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node94_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node94_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node95_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node95_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node96_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node96_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node97_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node97_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node98_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node98_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node99_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node99_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node100_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node100_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node101_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node101_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node102_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node102_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node103_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node103_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node104_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node104_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node105_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node105_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node106_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node106_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node107_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node107_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node108_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node108_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node109_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node109_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node110_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node110_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node111_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node111_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node112_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node112_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node113_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node113_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node114_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node114_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node115_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node115_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node116_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node116_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node117_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node117_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node118_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node118_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node119_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node119_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node120_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node120_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node121_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node121_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node122_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node122_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node123_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node123_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node124_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node124_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node125_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node125_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node126_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node126_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node127_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node127_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node128_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node128_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node129_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node129_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node130_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node130_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node131_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node131_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node132_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node132_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node133_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node133_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node134_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node134_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node135_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node135_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node136_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node136_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node137_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node137_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node138_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node138_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node139_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node139_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node140_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node140_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node141_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node141_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node142_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node142_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node143_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node143_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node144_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node144_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node145_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node145_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node146_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node146_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node147_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node147_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node148_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node148_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node149_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node149_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node150_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node150_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node151_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node151_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node152_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node152_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node153_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node153_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node154_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node154_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node155_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node155_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node156_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node156_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node157_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node157_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node158_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node158_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node159_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node159_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node160_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node160_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node161_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node161_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node162_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node162_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node163_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node163_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node164_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node164_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node165_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node165_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node166_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node166_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node167_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node167_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node168_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node168_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node169_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node169_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node170_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node170_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node171_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node171_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node172_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node172_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node173_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node173_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node174_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node174_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node175_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node175_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node176_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node176_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node177_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node177_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node178_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node178_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node179_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node179_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node180_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node180_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node181_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node181_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node182_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node182_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node183_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node183_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node184_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node184_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node185_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node185_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node186_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node186_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node187_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node187_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node188_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node188_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node189_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node189_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node190_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node190_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node191_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node191_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node192_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node192_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node193_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node193_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node194_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node194_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node195_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node195_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node196_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node196_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node197_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node197_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node198_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node198_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node199_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node199_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node200_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node200_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node201_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node201_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node202_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node202_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node203_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node203_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node204_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node204_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node205_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node205_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node206_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node206_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node207_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node207_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node208_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node208_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node209_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node209_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node210_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node210_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node211_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node211_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node212_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node212_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node213_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node213_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node214_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node214_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node215_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node215_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node216_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node216_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node217_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node217_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node218_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node218_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node219_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node219_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node220_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node220_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node221_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node221_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node222_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node222_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node223_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node223_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node224_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node224_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node225_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node225_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node226_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node226_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node227_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node227_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node228_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node228_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node229_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node229_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node230_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node230_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node231_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node231_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node232_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node232_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node233_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node233_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node234_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node234_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node235_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node235_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node236_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node236_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node237_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node237_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node238_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node238_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node239_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node239_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node240_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node240_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node241_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node241_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node242_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node242_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node243_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node243_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node244_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node244_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node245_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node245_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node246_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node246_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node247_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node247_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node248_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node248_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node249_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node249_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node250_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node250_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node251_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node251_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node252_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node252_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node253_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node253_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node254_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node254_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node255_y_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node255_y_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node32_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node32_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node33_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node33_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node34_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node34_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node35_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node35_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node36_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node36_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node37_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node37_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node38_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node38_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node39_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node39_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node40_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node40_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node41_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node41_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node42_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node42_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node43_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node43_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node44_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node44_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node45_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node45_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node46_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node46_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node47_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node47_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node48_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node48_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node49_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node49_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node50_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node50_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node51_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node51_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node52_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node52_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node53_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node53_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node54_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node54_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node55_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node55_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node56_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node56_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node57_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node57_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node58_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node58_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node59_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node59_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node60_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node60_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node61_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node61_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node62_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node62_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node63_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node63_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node64_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node64_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node65_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node65_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node66_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node66_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node67_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node67_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node68_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node68_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node69_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node69_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node70_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node70_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node71_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node71_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node72_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node72_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node73_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node73_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node74_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node74_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node75_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node75_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node76_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node76_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node77_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node77_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node78_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node78_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node79_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node79_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node80_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node80_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node81_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node81_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node82_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node82_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node83_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node83_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node84_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node84_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node85_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node85_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node86_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node86_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node87_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node87_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node88_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node88_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node89_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node89_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node90_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node90_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node91_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node91_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node92_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node92_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node93_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node93_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node94_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node94_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node95_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node95_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node96_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node96_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node97_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node97_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node98_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node98_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node99_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node99_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node100_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node100_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node101_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node101_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node102_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node102_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node103_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node103_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node104_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node104_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node105_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node105_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node106_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node106_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node107_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node107_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node108_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node108_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node109_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node109_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node110_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node110_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node111_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node111_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node112_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node112_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node113_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node113_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node114_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node114_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node115_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node115_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node116_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node116_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node117_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node117_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node118_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node118_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node119_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node119_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node120_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node120_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node121_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node121_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node122_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node122_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node123_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node123_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node124_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node124_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node125_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node125_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node126_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node126_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node127_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node127_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node128_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node128_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node129_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node129_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node130_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node130_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node131_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node131_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node132_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node132_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node133_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node133_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node134_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node134_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node135_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node135_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node136_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node136_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node137_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node137_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node138_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node138_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node139_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node139_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node140_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node140_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node141_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node141_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node142_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node142_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node143_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node143_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node144_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node144_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node145_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node145_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node146_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node146_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node147_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node147_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node148_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node148_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node149_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node149_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node150_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node150_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node151_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node151_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node152_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node152_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node153_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node153_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node154_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node154_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node155_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node155_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node156_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node156_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node157_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node157_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node158_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node158_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node159_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node159_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node160_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node160_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node161_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node161_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node162_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node162_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node163_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node163_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node164_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node164_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node165_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node165_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node166_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node166_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node167_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node167_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node168_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node168_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node169_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node169_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node170_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node170_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node171_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node171_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node172_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node172_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node173_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node173_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node174_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node174_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node175_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node175_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node176_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node176_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node177_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node177_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node178_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node178_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node179_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node179_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node180_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node180_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node181_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node181_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node182_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node182_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node183_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node183_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node184_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node184_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node185_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node185_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node186_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node186_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node187_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node187_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node188_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node188_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node189_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node189_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node190_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node190_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node191_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node191_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node192_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node192_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node193_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node193_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node194_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node194_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node195_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node195_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node196_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node196_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node197_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node197_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node198_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node198_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node199_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node199_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node200_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node200_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node201_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node201_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node202_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node202_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node203_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node203_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node204_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node204_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node205_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node205_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node206_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node206_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node207_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node207_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node208_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node208_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node209_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node209_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node210_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node210_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node211_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node211_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node212_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node212_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node213_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node213_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node214_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node214_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node215_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node215_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node216_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node216_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node217_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node217_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node218_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node218_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node219_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node219_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node220_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node220_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node221_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node221_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node222_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node222_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node223_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node223_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node224_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node224_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node225_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node225_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node226_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node226_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node227_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node227_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node228_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node228_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node229_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node229_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node230_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node230_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node231_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node231_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node232_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node232_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node233_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node233_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node234_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node234_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node235_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node235_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node236_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node236_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node237_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node237_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node238_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node238_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node239_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node239_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node240_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node240_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node241_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node241_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node242_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node242_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node243_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node243_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node244_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node244_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node245_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node245_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node246_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node246_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node247_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node247_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node248_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node248_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node249_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node249_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node250_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node250_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node251_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node251_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node252_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node252_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node253_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node253_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node254_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node254_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node255_z_3x3v_tensor_p1(const double *w, const double *dxv,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
     double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l, double* GKYL_RESTRICT Ghat_r ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node255_z_3x3v_tensor_p1(const double *w, const double *dxv, const int edge,
     const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Ghat ); 

EXTERN_C_END 
