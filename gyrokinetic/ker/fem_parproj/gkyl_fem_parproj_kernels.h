// Gkyl ------------------------------------------------------------------------
//
// Header file for parallel FEM projection operator.
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

long fem_parproj_num_nodes_global_1x_ser_p1_periodicx(int numCellsPar);
long fem_parproj_num_nodes_global_1x_ser_p1_nonperiodicx(int numCellsPar);

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_inx_periodicx(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_inx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_upx_periodicx(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_upx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs);

void fem_parproj_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p1_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p1_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p1_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p1_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p1_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p1_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p1_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p1_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p1_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p1_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p1_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p1_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p1_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p1_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p1_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p1_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p1_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p1_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_sol_stencil_1x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_1x_ser_p2_periodicx(int numCellsPar);
long fem_parproj_num_nodes_global_1x_ser_p2_nonperiodicx(int numCellsPar);

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_inx_periodicx(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_inx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_upx_periodicx(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_upx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs);

void fem_parproj_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p2_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p2_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p2_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_1x_ser_p2_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p2_lox_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p2_lox_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p2_upx_nondirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_1x_ser_p2_upx_dirichletx(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p2_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p2_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p2_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p2_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_1x_ser_p2_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p2_inx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p2_lox_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p2_lox_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p2_upx_nondirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_1x_ser_p2_upx_dirichletx(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_sol_stencil_1x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


long fem_parproj_num_nodes_global_2x_ser_p1_periodicy(int numCellsPar);
long fem_parproj_num_nodes_global_2x_ser_p1_nonperiodicy(int numCellsPar);

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_iny_periodicy(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_iny_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_upy_periodicy(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_upy_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs);

void fem_parproj_lhs_stencil_noweight_2x_ser_p1_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p1_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p1_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p1_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p1_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p1_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p1_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p1_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p1_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p1_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p1_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p1_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p1_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p1_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p1_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p1_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p1_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p1_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p1_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p1_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p1_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p1_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p1_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p1_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_sol_stencil_2x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_2x_ser_p2_periodicy(int numCellsPar);
long fem_parproj_num_nodes_global_2x_ser_p2_nonperiodicy(int numCellsPar);

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_iny_periodicy(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_iny_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_upy_periodicy(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_upy_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs);

void fem_parproj_lhs_stencil_noweight_2x_ser_p2_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p2_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p2_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p2_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p2_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_2x_ser_p2_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p2_iny_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p2_iny_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p2_loy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p2_loy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p2_upy_nondirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_2x_ser_p2_upy_dirichlety(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p2_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p2_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p2_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p2_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p2_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_2x_ser_p2_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p2_iny_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p2_iny_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p2_loy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p2_loy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p2_upy_nondirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_2x_ser_p2_upy_dirichlety(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_sol_stencil_2x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


long fem_parproj_num_nodes_global_3x_ser_p1_periodicz(int numCellsPar);
long fem_parproj_num_nodes_global_3x_ser_p1_nonperiodicz(int numCellsPar);

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_inz_periodicz(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_inz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_upz_periodicz(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_upz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs);

void fem_parproj_lhs_stencil_noweight_3x_ser_p1_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p1_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p1_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p1_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p1_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p1_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p1_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p1_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p1_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p1_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p1_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p1_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p1_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p1_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p1_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p1_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p1_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p1_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p1_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p1_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p1_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p1_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p1_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p1_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_sol_stencil_3x_ser_p1(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);

long fem_parproj_num_nodes_global_3x_ser_p2_periodicz(int numCellsPar);
long fem_parproj_num_nodes_global_3x_ser_p2_nonperiodicz(int numCellsPar);

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_inz_periodicz(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_inz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_upz_periodicz(int numCellsPar, int parIdx, long *globalIdxs);
GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_upz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs);

void fem_parproj_lhs_stencil_noweight_3x_ser_p2_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p2_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p2_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p2_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p2_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_noweight_3x_ser_p2_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p2_inz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p2_inz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p2_loz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p2_loz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p2_upz_nondirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);
void fem_parproj_lhs_stencil_weighted_3x_ser_p2_upz_dirichletz(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p2_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p2_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p2_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p2_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p2_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_noweight_3x_ser_p2_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p2_inz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p2_inz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p2_loz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p2_loz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p2_upz_nondirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_parproj_src_stencil_weighted_3x_ser_p2_upz_dirichletz(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_parproj_sol_stencil_3x_ser_p2(const double *sol_nodal_global, long nodeOff, const long *globalIdxs, double *sol_modal_local);


EXTERN_C_END 
