#pragma once 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void binop_cross_mul_1d_2d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_2d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_2d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_2d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_2d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_2d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_2d_ser_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_2d_ser_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_2d_ser_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_3d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_3d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_3d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_3d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_3d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_3d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_3d_ser_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_3d_ser_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_3d_ser_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_4d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_4d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_4d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_4d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_4d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_4d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1d_4d_ser_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1d_4d_ser_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1d_4d_ser_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_3d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_3d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_3d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_3d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_3d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_3d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_3d_ser_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_3d_ser_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_3d_ser_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_4d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_4d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_4d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_4d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_4d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_4d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_4d_ser_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_4d_ser_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_4d_ser_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_5d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_5d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_5d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2d_5d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2d_5d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2d_5d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3d_4d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3d_4d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3d_4d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3d_4d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3d_4d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3d_4d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3d_4d_ser_p3(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3d_4d_ser_p3(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3d_4d_ser_p3(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3d_5d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3d_5d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3d_5d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3d_5d_ser_p2(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3d_5d_ser_p2(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3d_5d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3d_6d_ser_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3d_6d_ser_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3d_6d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2); 

EXTERN_C_END 
