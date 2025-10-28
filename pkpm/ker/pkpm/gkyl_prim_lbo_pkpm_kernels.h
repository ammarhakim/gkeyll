#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 

GKYL_CU_DH void pkpm_self_prim_moments_1x1v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void pkpm_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void pkpm_self_prim_moments_2x1v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void pkpm_self_prim_moments_3x1v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void pkpm_self_prim_moments_1x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu); 

GKYL_CU_DH void pkpm_self_prim_moments_2x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu); 

EXTERN_C_END 
