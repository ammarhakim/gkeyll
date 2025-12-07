#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_set_prob_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_get_sol_1x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_pressure_1x_ser_p1(double gas_gamma, const double *moms, const double *udrift, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_set_prob_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_get_sol_1x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_set_prob_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_get_sol_1x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_set_prob_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_get_sol_1x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_pressure_1x_ser_p2(double gas_gamma, const double *moms, const double *udrift, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_set_prob_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_get_sol_2x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_pressure_2x_ser_p1(double gas_gamma, const double *moms, const double *udrift, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_set_prob_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_get_sol_2x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_set_prob_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_get_sol_2x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_set_prob_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_get_sol_2x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_pressure_3x_ser_p1(double gas_gamma, const double *moms, const double *udrift, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms, double gas_gamma, double mass); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *moms); 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, double* GKYL_RESTRICT out); 

EXTERN_C_END 
