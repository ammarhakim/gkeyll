#pragma once

#include <gkyl_util.h>
#include <math.h>

EXTERN_C_BEG

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x_ser_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p2_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p2_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p2_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p2_eval_dirs_2(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_13(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_013(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_2(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_23(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_023(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_3(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_03(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_24(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_024(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0124(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_124(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_34(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_034(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0134(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_134(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_4(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_04(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_014(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_14(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_34(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2345(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_034(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0134(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_134(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_35(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_035(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0135(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01235(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0235(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_135(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1235(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_235(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_4(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_45(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_045(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0145(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01245(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0245(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_145(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1245(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_245(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_04(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_014(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0124(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_024(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_14(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_124(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_24(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_5(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_05(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_015(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0125(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_025(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_15(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_125(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_25(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_2(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_0123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_13(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_013(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_2(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_23(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_023(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_3(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_03(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x1v_gkhyb_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_4(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_34(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_034(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0134(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_134(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_234(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_04(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_014(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0124(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_024(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_14(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_124(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_24(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar);


EXTERN_C_END
