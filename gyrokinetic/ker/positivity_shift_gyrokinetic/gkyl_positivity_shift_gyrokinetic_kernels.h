#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH bool positivity_shift_gyrokinetic_conf_pos_check_1x_ser_p1(const double *fld);

GKYL_CU_DH bool positivity_shift_gyrokinetic_shift_only_1x1v_ser_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_gyrokinetic_MRS_limiter_1x1v_ser_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_gyrokinetic_shift_only_1x2v_ser_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_gyrokinetic_MRS_limiter_1x2v_ser_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_gyrokinetic_conf_pos_check_2x_ser_p1(const double *fld);

GKYL_CU_DH bool positivity_shift_gyrokinetic_shift_only_2x2v_ser_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_gyrokinetic_MRS_limiter_2x2v_ser_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_gyrokinetic_conf_pos_check_3x_ser_p1(const double *fld);

GKYL_CU_DH bool positivity_shift_gyrokinetic_shift_only_3x2v_ser_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_gyrokinetic_MRS_limiter_3x2v_ser_p1(double ffloor, double *distf);


EXTERN_C_END
