#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double conf_diffusion_vol_1x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfx_1x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfx_1x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

GKYL_CU_DH double conf_diffusion_vol_2x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfx_2x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfx_2x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfy_2x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfy_2x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

GKYL_CU_DH double conf_diffusion_vol_3x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfx_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfx_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfy_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfy_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfz_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfz_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

EXTERN_C_END
