// Gkyl ------------------------------------------------------------------------
//
// Header file for parallel FEM projection operator that couples cells along x.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once 
#include <gkyl_util.h> 
#include <gkyl_mat.h> 
EXTERN_C_BEG 

// This needs to be inside EXTERN_C
#include <gkyl_mat_triples.h> 

long fem_parproj_couplex_num_nodes_global_1x_ser_p1_periodicx(const int *num_cells);
long fem_parproj_couplex_num_nodes_global_1x_ser_p1_nonperiodicx(const int *num_cells);

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_inx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_inx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_upx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_upx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);

void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_inx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_inx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_inx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_inx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_couplex_sol_stencil_1x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_couplex_num_nodes_global_1x_ser_p2_periodicx(const int *num_cells);
long fem_parproj_couplex_num_nodes_global_1x_ser_p2_nonperiodicx(const int *num_cells);

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_inx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_inx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_upx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_upx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs);

void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_inx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_inx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_inx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_inx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_couplex_sol_stencil_1x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


long fem_parproj_couplex_num_nodes_global_2x_ser_p1_periodicy(const int *num_cells);
long fem_parproj_couplex_num_nodes_global_2x_ser_p1_nonperiodicy(const int *num_cells);

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);

void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_couplex_sol_stencil_2x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_couplex_num_nodes_global_2x_ser_p2_periodicy(const int *num_cells);
long fem_parproj_couplex_num_nodes_global_2x_ser_p2_nonperiodicy(const int *num_cells);

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs);

void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_couplex_sol_stencil_2x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


long fem_parproj_couplex_num_nodes_global_3x_ser_p1_periodicz(const int *num_cells);
long fem_parproj_couplex_num_nodes_global_3x_ser_p1_nonperiodicz(const int *num_cells);

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);

void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_couplex_sol_stencil_3x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_couplex_num_nodes_global_3x_ser_p2_periodicz(const int *num_cells);
long fem_parproj_couplex_num_nodes_global_3x_ser_p2_nonperiodicz(const int *num_cells);

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);
GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs);

void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_couplex_sol_stencil_3x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


EXTERN_C_END 
