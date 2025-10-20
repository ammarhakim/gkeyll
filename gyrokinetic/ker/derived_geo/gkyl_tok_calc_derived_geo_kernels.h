#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void tok_derived_geo_3x_Ser_p1(const double *gij, const double *bmag, const double *J, double *Jinv, double *grij, double *bi, double *cmag, double *Jtot, double *Jtotinv, double *gxxJ, double *gxyJ, double *gyyJ, double *gxzJ, double *eps2); 
GKYL_CU_DH void tok_derived_geo_quad_3x_Ser_p1(const double *gij, const double *bmag, const double *J, double *Jinv, double *grij, double *bi, double *cmag, double *Jtot, double *Jtotinv, double *gxxJ, double *gxyJ, double *gyyJ, double *gxzJ, double *eps2); 
GKYL_CU_DH void tok_derived_geo_3x_Ser_p2(const double *gij, const double *bmag, const double *J, double *Jinv, double *grij, double *bi, double *cmag, double *Jtot, double *Jtotinv, double *gxxJ, double *gxyJ, double *gyyJ, double *gxzJ, double *eps2); 
EXTERN_C_END 
