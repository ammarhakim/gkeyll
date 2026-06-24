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

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar);
GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar);


EXTERN_C_END
