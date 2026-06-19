#pragma once
#include <gkyl_util.h>
EXTERN_C_BEG

// 1D serendipity
GKYL_CU_DH void dg_differentiate_1d_ser_p1_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_1d_ser_p2_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_1d_ser_p2_dir0_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_1d_ser_p3_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_1d_ser_p3_dir0_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_1d_ser_p3_dir0_ord3(const double *f, double *df);

// 2D serendipity
GKYL_CU_DH void dg_differentiate_2d_ser_p1_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p1_dir1_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p2_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p2_dir0_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p2_dir1_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p2_dir1_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p3_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p3_dir0_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p3_dir0_ord3(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p3_dir1_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p3_dir1_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_2d_ser_p3_dir1_ord3(const double *f, double *df);

// 3D serendipity
GKYL_CU_DH void dg_differentiate_3d_ser_p1_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p1_dir1_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p1_dir2_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p2_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p2_dir0_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p2_dir1_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p2_dir1_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p2_dir2_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p2_dir2_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir0_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir0_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir0_ord3(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir1_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir1_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir1_ord3(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir2_ord1(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir2_ord2(const double *f, double *df);
GKYL_CU_DH void dg_differentiate_3d_ser_p3_dir2_ord3(const double *f, double *df);

EXTERN_C_END
