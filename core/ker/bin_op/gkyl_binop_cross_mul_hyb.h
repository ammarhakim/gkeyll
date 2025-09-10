#pragma once 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void binop_cross_mul_1x1v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1x1v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1x1v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1x2v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1x2v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1x2v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_1x3v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_1x3v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_1x3v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2x1v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2x1v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2x1v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2x2v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2x2v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2x2v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_2x3v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_2x3v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_2x3v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3x1v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3x1v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3x1v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3x2v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3x2v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3x2v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

GKYL_CU_DH void binop_cross_mul_3x3v_hyb_p1(const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_3x3v_hyb_p1(double a, const double *f, const double *g, double *fg); 
GKYL_CU_DH void binop_cross_mul_accumulate_comp_par_3x3v_hyb_p1(double a, const double *f, const double *g, double *fg, int linc2); 

EXTERN_C_END 
