#pragma once 
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH
void gk_neut_hamil_1x3v_tensor_p1(const double *w, const double *dxv, const double *gij, double* GKYL_RESTRICT hamil);

GKYL_CU_DH
void gk_neut_hamil_2x3v_tensor_p1(const double *w, const double *dxv, const double *gij, double* GKYL_RESTRICT hamil);

GKYL_CU_DH
void gk_neut_hamil_3x3v_tensor_p1(const double *w, const double *dxv, const double *gij, double* GKYL_RESTRICT hamil);

EXTERN_C_END
