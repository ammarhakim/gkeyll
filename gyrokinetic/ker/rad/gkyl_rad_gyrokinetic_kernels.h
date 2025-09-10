#pragma once

#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double rad_gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap_prime, 
            const double *nvnu, const double *nvsqnu, 
            const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH void rad_gyrokinetic_drag_nuvpar_1x2v_ser_p1(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_numu_1x2v_ser_p1(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_nI_nu_1x2v_ser_p1(const double *vnu_surf, const double *vnu,
            const double *vsqnu_surf, const double *vsqnu, const double *nI, 
            double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
            double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu); 

GKYL_CU_DH double rad_gyrokinetic_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap_prime, 
            const double *nvnu, const double *nvsqnu, 
            const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH void rad_gyrokinetic_drag_nuvpar_1x2v_ser_p2(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_numu_1x2v_ser_p2(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_nI_nu_1x2v_ser_p2(const double *vnu_surf, const double *vnu,
            const double *vsqnu_surf, const double *vsqnu, const double *nI, 
            double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
            double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu); 

GKYL_CU_DH double rad_gyrokinetic_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap_prime, 
            const double *nvnu, const double *nvsqnu, 
            const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH void rad_gyrokinetic_drag_nuvpar_2x2v_ser_p1(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_numu_2x2v_ser_p1(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_nI_nu_2x2v_ser_p1(const double *vnu_surf, const double *vnu,
            const double *vsqnu_surf, const double *vsqnu, const double *nI, 
            double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
            double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu); 

GKYL_CU_DH double rad_gyrokinetic_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap_prime, 
            const double *nvnu, const double *nvsqnu, 
            const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfmu_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfmu_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH void rad_gyrokinetic_drag_nuvpar_2x2v_ser_p2(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_numu_2x2v_ser_p2(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_nI_nu_2x2v_ser_p2(const double *vnu_surf, const double *vnu,
            const double *vsqnu_surf, const double *vsqnu, const double *nI, 
            double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
            double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu); 

GKYL_CU_DH double rad_gyrokinetic_vol_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap_prime, 
            const double *nvnu, const double *nvsqnu, 
            const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
            const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
            const double *fl, const double *fc, const double *fr, 
            double* GKYL_RESTRICT out); 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
            const double *vmap_prime_edge, const double *vmap_prime_skin, 
            const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
            const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) ; 
GKYL_CU_DH void rad_gyrokinetic_drag_nuvpar_3x2v_ser_p1(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_numu_3x2v_ser_p1(const double *vmap, const double *vmapSq,
            double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
            const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 
GKYL_CU_DH void rad_gyrokinetic_drag_nI_nu_3x2v_ser_p1(const double *vnu_surf, const double *vnu,
            const double *vsqnu_surf, const double *vsqnu, const double *nI, 
            double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
            double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu); 


EXTERN_C_END
