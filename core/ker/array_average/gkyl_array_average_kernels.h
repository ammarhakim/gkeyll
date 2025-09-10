#pragma once

#include <gkyl_util.h>
EXTERN_C_BEG

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_1x_ser_p2_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgx(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgxy(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgxz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgxyz(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);

EXTERN_C_END