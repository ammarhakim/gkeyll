// Gkyl ------------------------------------------------------------------------
//
// Header file for perpendicular FEM Helmholtz/Poisson solver.
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
 
long fem_poisson_perp_num_nodes_global_2x_ser_p1_periodicx(const int *numCells);
long fem_poisson_perp_num_nodes_global_2x_ser_p1_nonperiodicx(const int *numCells);

GKYL_CU_DH void fem_poisson_perp_local_to_global_2x_ser_p1_inx_periodicx(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_2x_ser_p1_inx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_2x_ser_p1_upx_periodicx(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_2x_ser_p1_upx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs);

void fem_poisson_perp_lhs_stencil_2x_ser_p1_inx_periodicx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_2x_ser_p1_lox_periodicx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_2x_ser_p1_lox_dirichletx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_2x_ser_p1_lox_neumannx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_2x_ser_p1_upx_periodicx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_2x_ser_p1_upx_dirichletx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_2x_ser_p1_upx_neumannx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);

GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_inx_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_lox_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_lox_dirichletx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_lox_neumannx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_upx_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_upx_dirichletx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_upx_neumannx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_poisson_perp_sol_stencil_2x_ser_p1(const double *sol_nodal_global, long perpOff, const long *globalIdxs, double *sol_modal_local);


long fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_periodicy(const int *numCells);
long fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy(const int *numCells);
long fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy(const int *numCells);
long fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells);

GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs);
GKYL_CU_DH void fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs);

void fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_iny_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_iny_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_iny_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_iny_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_loy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_loy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_loy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_loy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_loy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_loy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_upy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_upy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_upy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_upy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_upy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_upy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_loy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_loy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_loy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_loy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_loy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_loy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_upy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_upy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_upy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_upy_periodicy(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_upy_dirichlety(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);
void fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_upy_neumanny(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri);

GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_loy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_loy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_loy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_upy_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_upy_dirichlety(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_upy_neumanny(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc);

GKYL_CU_DH void fem_poisson_perp_sol_stencil_3x_ser_p1(const double *sol_nodal_global, long perpOff, const long *globalIdxs, double *sol_modal_local);




EXTERN_C_END 
