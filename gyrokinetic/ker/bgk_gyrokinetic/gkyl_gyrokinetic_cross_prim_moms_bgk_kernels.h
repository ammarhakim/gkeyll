#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_1x1v_ser_p1(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_1x1v_ser_p2(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_1x2v_ser_p1(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_1x2v_ser_p2(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_2x2v_ser_p1(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_2x2v_ser_p2(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_3x2v_ser_p1(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_3x2v_ser_p2(const double beta, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross); 

EXTERN_C_END 
