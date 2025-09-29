#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double gk_anomalous_diffusion_vol_2x2v_ser_p1(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out);

GKYL_CU_DH double gk_anomalous_diffusion_surfx_2x2v_ser_p1(const double *wc, const double *dxc, const double *nul, const double *nuc, const double *nur, const double *jacobgeo_invl, const double *jacobgeo_invc, const double *jacobgeo_invr, const double *Jfl, const double *Jfc, const double *Jfr, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_surfx_lower_zero_flux_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_surfx_lower_bound_local_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_diagx_lower_bound_local_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_diagx_lower_bound_recovery_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_surfx_upper_zero_flux_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_surfx_upper_bound_local_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin, int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_diagx_upper_bound_local_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out);
GKYL_CU_DH double gk_anomalous_diffusion_boundary_diagx_upper_bound_recovery_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_invSkin, const double *jacobgeo_invGhost, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out);


EXTERN_C_END
