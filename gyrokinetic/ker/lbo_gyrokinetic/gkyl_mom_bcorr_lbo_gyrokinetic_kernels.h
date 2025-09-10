#pragma once
#include <math.h>
#include <gkyl_eqn_type.h>
#include <gkyl_util.h>

EXTERN_C_BEG 

GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_1x1v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_1x1v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_1x2v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_1x2v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_2x2v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_2x2v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_3x2v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_3x2v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_1x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_1x2v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_2x2v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out); 

EXTERN_C_END 
