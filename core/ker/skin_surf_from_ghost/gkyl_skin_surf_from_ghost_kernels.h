#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void skin_surf_from_ghost_lowerx_1x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_upperx_1x_ser_p1(const double *fghost, double *fskin);

GKYL_CU_DH void skin_surf_from_ghost_lowerx_2x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_upperx_2x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_lowery_2x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_uppery_2x_ser_p1(const double *fghost, double *fskin);

GKYL_CU_DH void skin_surf_from_ghost_lowerx_3x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_upperx_3x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_lowery_3x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_uppery_3x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_lowerz_3x_ser_p1(const double *fghost, double *fskin);
GKYL_CU_DH void skin_surf_from_ghost_upperz_3x_ser_p1(const double *fghost, double *fskin);


EXTERN_C_END
