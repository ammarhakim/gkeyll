#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double conf_diffusion_vol_1x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfxx_1x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfxx_1x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

GKYL_CU_DH double conf_diffusion_vol_2x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfxx_2x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfxy_2x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfyx_2x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfyy_2x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfxx_2x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfxy_2x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfyx_2x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfyy_2x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

GKYL_CU_DH double conf_diffusion_vol_3x_ser_p1(const double *w, const double *dx, const double *K, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfxx_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfxy_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfxz_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfyx_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfyy_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfyz_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfzx_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfzy_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_surfzz_3x_ser_p1(const double *wc, const double *dxc, const double *Kl, const double *Kc, const double *Kr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfxx_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfxy_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfxz_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfyx_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfyy_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfyz_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfzx_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfzy_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
GKYL_CU_DH double conf_diffusion_boundary_surfzz_3x_ser_p1(const double *wskin, const double *dxskin, const double *Kedge, const double *Kskin, int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

EXTERN_C_END
