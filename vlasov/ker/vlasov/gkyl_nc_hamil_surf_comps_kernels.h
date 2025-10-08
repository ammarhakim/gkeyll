#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 


GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_y_2x2v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_y_2x2v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 


GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node0_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node0_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node1_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node1_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node2_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node2_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node3_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node3_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node4_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node4_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node5_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node5_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node6_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node6_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node7_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node7_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node8_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node8_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node9_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node9_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node10_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node10_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node11_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node11_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node12_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node12_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node13_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node13_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node14_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node14_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node15_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node15_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node16_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node16_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node17_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node17_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node18_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node18_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node19_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node19_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node20_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node20_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node21_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node21_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node22_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node22_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node23_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node23_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node24_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node24_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node25_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node25_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node26_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node26_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node0_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node0_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node1_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node1_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node2_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node2_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node3_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node3_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node4_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node4_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node5_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node5_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node6_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node6_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node7_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node7_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node8_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node8_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node9_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node9_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node10_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node10_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node11_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node11_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node12_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node12_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node13_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node13_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node14_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node14_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node15_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node15_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node16_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node16_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node17_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node17_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node18_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node18_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node19_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node19_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node20_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node20_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node21_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node21_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node22_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node22_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node23_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node23_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node24_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node24_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node25_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node25_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node26_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node26_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node27_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node27_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node28_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node28_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node29_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node29_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node30_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node30_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node31_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node31_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node32_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node32_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node33_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node33_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node34_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node34_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node35_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node35_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node36_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node36_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node37_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node37_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node38_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node38_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node39_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node39_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node40_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node40_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node41_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node41_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node42_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node42_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node43_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node43_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node44_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node44_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node45_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node45_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node46_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node46_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node47_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node47_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node48_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node48_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node49_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node49_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node50_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node50_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node51_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node51_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node52_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node52_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node53_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node53_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node54_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node54_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node55_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node55_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node56_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node56_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node57_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node57_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node58_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node58_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node59_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node59_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node60_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node60_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node61_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node61_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node62_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node62_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node63_x_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node63_x_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node0_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node0_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node1_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node1_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node2_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node2_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node3_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node3_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node4_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node4_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node5_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node5_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node6_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node6_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node7_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node7_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node8_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node8_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node9_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node9_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node10_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node10_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node11_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node11_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node12_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node12_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node13_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node13_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node14_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node14_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node15_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node15_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node16_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node16_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node17_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node17_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node18_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node18_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node19_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node19_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node20_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node20_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node21_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node21_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node22_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node22_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node23_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node23_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node24_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node24_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node25_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node25_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node26_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node26_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node0_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node0_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node1_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node1_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node2_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node2_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node3_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node3_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node4_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node4_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node5_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node5_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node6_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node6_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node7_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node7_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node8_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node8_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node9_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node9_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node10_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node10_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node11_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node11_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node12_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node12_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node13_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node13_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node14_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node14_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node15_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node15_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node16_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node16_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node17_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node17_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node18_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node18_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node19_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node19_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node20_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node20_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node21_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node21_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node22_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node22_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node23_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node23_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node24_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node24_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node25_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node25_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node26_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node26_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node27_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node27_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node28_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node28_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node29_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node29_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node30_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node30_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node31_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node31_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node32_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node32_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node33_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node33_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node34_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node34_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node35_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node35_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node36_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node36_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node37_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node37_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node38_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node38_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node39_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node39_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node40_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node40_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node41_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node41_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node42_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node42_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node43_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node43_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node44_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node44_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node45_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node45_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node46_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node46_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node47_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node47_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node48_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node48_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node49_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node49_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node50_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node50_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node51_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node51_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node52_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node52_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node53_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node53_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node54_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node54_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node55_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node55_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node56_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node56_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node57_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node57_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node58_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node58_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node59_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node59_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node60_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node60_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node61_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node61_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node62_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node62_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node63_y_2x2v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node63_y_2x2v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_x_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_x_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 


GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_y_2x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_y_2x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 


GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node0_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node0_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node1_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node1_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node2_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node2_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node3_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node3_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node4_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node4_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node5_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node5_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node6_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node6_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node7_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node7_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node8_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node8_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node9_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node9_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node10_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node10_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node11_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node11_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node12_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node12_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node13_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node13_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node14_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node14_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node15_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node15_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node16_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node16_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node17_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node17_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node18_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node18_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node19_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node19_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node20_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node20_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node21_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node21_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node22_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node22_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node23_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node23_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node24_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node24_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node25_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node25_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node26_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node26_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node27_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node27_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node28_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node28_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node29_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node29_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node30_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node30_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node31_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node31_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node32_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node32_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node33_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node33_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node34_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node34_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node35_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node35_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node36_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node36_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node37_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node37_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node38_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node38_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node39_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node39_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node40_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node40_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node41_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node41_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node42_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node42_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node43_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node43_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node44_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node44_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node45_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node45_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node46_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node46_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node47_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node47_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node48_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node48_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node49_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node49_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node50_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node50_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node51_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node51_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node52_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node52_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node53_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node53_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node54_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node54_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node55_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node55_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node56_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node56_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node57_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node57_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node58_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node58_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node59_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node59_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node60_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node60_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node61_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node61_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node62_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node62_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node63_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node63_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node64_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node64_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node65_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node65_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node66_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node66_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node67_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node67_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node68_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node68_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node69_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node69_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node70_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node70_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node71_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node71_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node72_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node72_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node73_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node73_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node74_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node74_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node75_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node75_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node76_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node76_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node77_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node77_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node78_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node78_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node79_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node79_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node80_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node80_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node0_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node0_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node1_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node1_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node2_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node2_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node3_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node3_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node4_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node4_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node5_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node5_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node6_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node6_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node7_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node7_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node8_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node8_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node9_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node9_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node10_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node10_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node11_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node11_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node12_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node12_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node13_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node13_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node14_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node14_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node15_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node15_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node16_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node16_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node17_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node17_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node18_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node18_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node19_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node19_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node20_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node20_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node21_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node21_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node22_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node22_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node23_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node23_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node24_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node24_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node25_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node25_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node26_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node26_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node27_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node27_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node28_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node28_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node29_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node29_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node30_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node30_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node31_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node31_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node32_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node32_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node33_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node33_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node34_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node34_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node35_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node35_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node36_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node36_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node37_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node37_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node38_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node38_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node39_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node39_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node40_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node40_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node41_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node41_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node42_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node42_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node43_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node43_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node44_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node44_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node45_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node45_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node46_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node46_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node47_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node47_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node48_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node48_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node49_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node49_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node50_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node50_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node51_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node51_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node52_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node52_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node53_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node53_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node54_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node54_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node55_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node55_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node56_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node56_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node57_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node57_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node58_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node58_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node59_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node59_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node60_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node60_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node61_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node61_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node62_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node62_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node63_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node63_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node64_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node64_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node65_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node65_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node66_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node66_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node67_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node67_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node68_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node68_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node69_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node69_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node70_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node70_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node71_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node71_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node72_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node72_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node73_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node73_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node74_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node74_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node75_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node75_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node76_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node76_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node77_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node77_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node78_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node78_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node79_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node79_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node80_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node80_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node81_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node81_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node82_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node82_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node83_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node83_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node84_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node84_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node85_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node85_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node86_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node86_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node87_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node87_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node88_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node88_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node89_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node89_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node90_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node90_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node91_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node91_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node92_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node92_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node93_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node93_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node94_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node94_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node95_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node95_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node96_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node96_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node97_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node97_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node98_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node98_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node99_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node99_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node100_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node100_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node101_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node101_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node102_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node102_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node103_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node103_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node104_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node104_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node105_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node105_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node106_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node106_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node107_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node107_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node108_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node108_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node109_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node109_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node110_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node110_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node111_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node111_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node112_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node112_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node113_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node113_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node114_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node114_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node115_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node115_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node116_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node116_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node117_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node117_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node118_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node118_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node119_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node119_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node120_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node120_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node121_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node121_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node122_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node122_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node123_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node123_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node124_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node124_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node125_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node125_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node126_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node126_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node127_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node127_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node128_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node128_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node129_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node129_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node130_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node130_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node131_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node131_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node132_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node132_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node133_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node133_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node134_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node134_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node135_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node135_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node136_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node136_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node137_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node137_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node138_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node138_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node139_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node139_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node140_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node140_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node141_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node141_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node142_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node142_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node143_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node143_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node144_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node144_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node145_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node145_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node146_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node146_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node147_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node147_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node148_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node148_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node149_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node149_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node150_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node150_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node151_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node151_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node152_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node152_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node153_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node153_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node154_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node154_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node155_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node155_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node156_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node156_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node157_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node157_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node158_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node158_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node159_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node159_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node160_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node160_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node161_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node161_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node162_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node162_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node163_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node163_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node164_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node164_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node165_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node165_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node166_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node166_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node167_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node167_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node168_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node168_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node169_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node169_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node170_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node170_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node171_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node171_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node172_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node172_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node173_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node173_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node174_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node174_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node175_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node175_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node176_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node176_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node177_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node177_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node178_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node178_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node179_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node179_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node180_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node180_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node181_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node181_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node182_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node182_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node183_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node183_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node184_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node184_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node185_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node185_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node186_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node186_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node187_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node187_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node188_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node188_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node189_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node189_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node190_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node190_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node191_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node191_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node192_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node192_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node193_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node193_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node194_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node194_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node195_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node195_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node196_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node196_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node197_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node197_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node198_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node198_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node199_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node199_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node200_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node200_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node201_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node201_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node202_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node202_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node203_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node203_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node204_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node204_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node205_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node205_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node206_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node206_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node207_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node207_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node208_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node208_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node209_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node209_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node210_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node210_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node211_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node211_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node212_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node212_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node213_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node213_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node214_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node214_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node215_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node215_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node216_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node216_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node217_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node217_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node218_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node218_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node219_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node219_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node220_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node220_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node221_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node221_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node222_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node222_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node223_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node223_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node224_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node224_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node225_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node225_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node226_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node226_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node227_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node227_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node228_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node228_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node229_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node229_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node230_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node230_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node231_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node231_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node232_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node232_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node233_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node233_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node234_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node234_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node235_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node235_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node236_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node236_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node237_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node237_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node238_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node238_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node239_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node239_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node240_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node240_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node241_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node241_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node242_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node242_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node243_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node243_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node244_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node244_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node245_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node245_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node246_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node246_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node247_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node247_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node248_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node248_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node249_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node249_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node250_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node250_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node251_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node251_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node252_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node252_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node253_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node253_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node254_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node254_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node255_x_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node255_x_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node0_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node0_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node1_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node1_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node2_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node2_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node3_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node3_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node4_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node4_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node5_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node5_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node6_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node6_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node7_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node7_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node8_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node8_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node9_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node9_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node10_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node10_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node11_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node11_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node12_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node12_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node13_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node13_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node14_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node14_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node15_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node15_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node16_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node16_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node17_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node17_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node18_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node18_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node19_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node19_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node20_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node20_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node21_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node21_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node22_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node22_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node23_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node23_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node24_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node24_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node25_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node25_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node26_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node26_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node27_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node27_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node28_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node28_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node29_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node29_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node30_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node30_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node31_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node31_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node32_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node32_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node33_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node33_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node34_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node34_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node35_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node35_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node36_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node36_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node37_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node37_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node38_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node38_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node39_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node39_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node40_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node40_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node41_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node41_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node42_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node42_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node43_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node43_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node44_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node44_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node45_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node45_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node46_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node46_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node47_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node47_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node48_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node48_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node49_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node49_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node50_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node50_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node51_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node51_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node52_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node52_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node53_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node53_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node54_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node54_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node55_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node55_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node56_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node56_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node57_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node57_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node58_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node58_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node59_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node59_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node60_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node60_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node61_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node61_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node62_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node62_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node63_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node63_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node64_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node64_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node65_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node65_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node66_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node66_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node67_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node67_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node68_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node68_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node69_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node69_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node70_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node70_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node71_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node71_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node72_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node72_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node73_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node73_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node74_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node74_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node75_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node75_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node76_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node76_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node77_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node77_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node78_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node78_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node79_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node79_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_lo_node80_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_lo_node80_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node0_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node0_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node1_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node1_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node2_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node2_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node3_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node3_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node4_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node4_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node5_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node5_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node6_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node6_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node7_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node7_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node8_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node8_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node9_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node9_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node10_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node10_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node11_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node11_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node12_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node12_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node13_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node13_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node14_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node14_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node15_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node15_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node16_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node16_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node17_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node17_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node18_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node18_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node19_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node19_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node20_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node20_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node21_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node21_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node22_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node22_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node23_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node23_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node24_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node24_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node25_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node25_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node26_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node26_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node27_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node27_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node28_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node28_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node29_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node29_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node30_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node30_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node31_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node31_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node32_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node32_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node33_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node33_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node34_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node34_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node35_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node35_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node36_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node36_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node37_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node37_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node38_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node38_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node39_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node39_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node40_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node40_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node41_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node41_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node42_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node42_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node43_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node43_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node44_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node44_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node45_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node45_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node46_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node46_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node47_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node47_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node48_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node48_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node49_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node49_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node50_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node50_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node51_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node51_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node52_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node52_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node53_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node53_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node54_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node54_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node55_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node55_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node56_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node56_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node57_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node57_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node58_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node58_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node59_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node59_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node60_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node60_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node61_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node61_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node62_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node62_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node63_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node63_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node64_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node64_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node65_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node65_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node66_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node66_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node67_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node67_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node68_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node68_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node69_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node69_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node70_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node70_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node71_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node71_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node72_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node72_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node73_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node73_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node74_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node74_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node75_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node75_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node76_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node76_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node77_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node77_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node78_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node78_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node79_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node79_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node80_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node80_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node81_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node81_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node82_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node82_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node83_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node83_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node84_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node84_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node85_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node85_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node86_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node86_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node87_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node87_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node88_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node88_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node89_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node89_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node90_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node90_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node91_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node91_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node92_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node92_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node93_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node93_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node94_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node94_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node95_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node95_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node96_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node96_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node97_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node97_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node98_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node98_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node99_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node99_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node100_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node100_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node101_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node101_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node102_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node102_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node103_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node103_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node104_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node104_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node105_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node105_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node106_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node106_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node107_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node107_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node108_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node108_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node109_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node109_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node110_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node110_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node111_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node111_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node112_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node112_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node113_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node113_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node114_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node114_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node115_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node115_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node116_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node116_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node117_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node117_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node118_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node118_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node119_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node119_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node120_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node120_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node121_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node121_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node122_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node122_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node123_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node123_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node124_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node124_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node125_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node125_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node126_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node126_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node127_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node127_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node128_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node128_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node129_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node129_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node130_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node130_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node131_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node131_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node132_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node132_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node133_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node133_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node134_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node134_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node135_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node135_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node136_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node136_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node137_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node137_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node138_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node138_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node139_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node139_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node140_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node140_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node141_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node141_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node142_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node142_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node143_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node143_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node144_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node144_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node145_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node145_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node146_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node146_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node147_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node147_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node148_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node148_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node149_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node149_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node150_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node150_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node151_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node151_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node152_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node152_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node153_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node153_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node154_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node154_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node155_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node155_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node156_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node156_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node157_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node157_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node158_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node158_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node159_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node159_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node160_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node160_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node161_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node161_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node162_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node162_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node163_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node163_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node164_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node164_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node165_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node165_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node166_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node166_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node167_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node167_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node168_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node168_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node169_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node169_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node170_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node170_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node171_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node171_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node172_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node172_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node173_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node173_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node174_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node174_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node175_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node175_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node176_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node176_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node177_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node177_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node178_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node178_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node179_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node179_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node180_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node180_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node181_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node181_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node182_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node182_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node183_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node183_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node184_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node184_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node185_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node185_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node186_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node186_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node187_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node187_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node188_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node188_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node189_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node189_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node190_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node190_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node191_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node191_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node192_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node192_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node193_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node193_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node194_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node194_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node195_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node195_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node196_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node196_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node197_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node197_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node198_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node198_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node199_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node199_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node200_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node200_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node201_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node201_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node202_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node202_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node203_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node203_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node204_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node204_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node205_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node205_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node206_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node206_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node207_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node207_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node208_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node208_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node209_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node209_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node210_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node210_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node211_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node211_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node212_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node212_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node213_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node213_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node214_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node214_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node215_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node215_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node216_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node216_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node217_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node217_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node218_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node218_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node219_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node219_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node220_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node220_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node221_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node221_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node222_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node222_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node223_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node223_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node224_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node224_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node225_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node225_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node226_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node226_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node227_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node227_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node228_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node228_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node229_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node229_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node230_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node230_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node231_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node231_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node232_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node232_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node233_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node233_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node234_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node234_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node235_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node235_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node236_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node236_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node237_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node237_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node238_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node238_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node239_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node239_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node240_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node240_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node241_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node241_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node242_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node242_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node243_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node243_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node244_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node244_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node245_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node245_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node246_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node246_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node247_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node247_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node248_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node248_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node249_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node249_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node250_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node250_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node251_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node251_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node252_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node252_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node253_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node253_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node254_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node254_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_ho_node255_y_2x3v_ser_p2(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_ho_node255_y_2x3v_ser_p2(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 

GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_x_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_x_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 


GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_y_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_y_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 


GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node1_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node1_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node2_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node3_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node3_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node4_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node5_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node5_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node6_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node6_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node7_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node7_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node8_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node8_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node9_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node9_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node10_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node10_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node11_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node11_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node12_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node12_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node13_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node13_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node14_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node14_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node15_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node15_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node16_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node16_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node17_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node17_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node18_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node18_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node19_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node19_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node20_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node20_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node21_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node21_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node22_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node22_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node23_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node23_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node24_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node24_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node25_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node25_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node26_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node26_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node27_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node27_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node28_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node28_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node29_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node29_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node30_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node30_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node31_z_3x3v_ser_p1(const double *w, const double *dxv,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc, const double *fr,
        double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Fhat_l_nodal, double* GKYL_RESTRICT Fhat_r_nodal ); 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node31_z_3x3v_ser_p1(const double *w, const double *dxv, const int edge,
        const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fr, double* GKYL_RESTRICT Fhat_nodal ); 


EXTERN_C_END 
