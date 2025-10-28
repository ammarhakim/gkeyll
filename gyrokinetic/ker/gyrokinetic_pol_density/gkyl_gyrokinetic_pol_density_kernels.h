#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void gkyl_gyrokinetic_pol_density_1x_ser_p1(const double *dx, const double *epsilon, const double *phi, double *out);
GKYL_CU_DH void gkyl_gyrokinetic_pol_density_2x_ser_p1(const double *dx, const double *epsilon, const double *phi, double *out);
GKYL_CU_DH void gkyl_gyrokinetic_pol_density_3x_ser_p1(const double *dx, const double *epsilon, const double *phi, double *out);


EXTERN_C_END

