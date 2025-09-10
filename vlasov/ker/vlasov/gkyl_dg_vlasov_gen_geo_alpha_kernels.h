#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_gen_geo_cot_vec_3x_ser_p1(const double *tvComp, const double *gij, double* GKYL_RESTRICT cot_vec); 

GKYL_CU_DH int vlasov_gen_geo_alpha_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int vlasov_gen_geo_alpha_edge_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 

GKYL_CU_DH int vlasov_gen_geo_alpha_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int vlasov_gen_geo_alpha_edge_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 

GKYL_CU_DH int vlasov_gen_geo_alpha_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int vlasov_gen_geo_alpha_edge_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 

EXTERN_C_END 
