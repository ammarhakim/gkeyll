#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void translate_dim_1x_ser_p1_to_2x_p1(const double *fdo, double *ftar);

GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_dirx_lo(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_dirx_mid(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_dirx_up(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_diry_lo(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_diry_mid(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_diry_up(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_2x_ser_p1_to_3x_p1(const double *fdo, double *ftar);

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirx_lo(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirx_mid(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirx_up(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_diry_lo(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_diry_mid(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_diry_up(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirz_lo(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirz_mid(const double *fdo, double *ftar);
GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirz_up(const double *fdo, double *ftar);

GKYL_CU_DH void translate_dim_gyrokinetic_2x2v_ser_p1_from_1x2v_p1(const double *flow, double *fout);

GKYL_CU_DH void translate_dim_gyrokinetic_3x2v_ser_p1_from_1x2v_p1(const double *flow, double *fout);
GKYL_CU_DH void translate_dim_gyrokinetic_3x2v_ser_p1_from_2x2v_p1(const double *flow, double *fout);


EXTERN_C_END
