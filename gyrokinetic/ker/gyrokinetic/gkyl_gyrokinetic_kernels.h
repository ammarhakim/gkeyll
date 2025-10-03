#pragma once

#include <math.h>
#include <gkyl_util.h>
#include <gkyl_dg_geom.h>
#include <gkyl_gk_dg_geom.h>

EXTERN_C_BEG

GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x1v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x1v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfx_1x1v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfx_1x1v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_1x1v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfx_1x1v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfvpar_1x1v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfvpar_1x1v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x1v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x1v_ser_p2(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfx_1x1v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfx_1x1v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_1x1v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfx_1x1v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfvpar_1x1v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfvpar_1x1v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x2v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfx_1x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfx_1x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_1x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfx_1x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfvpar_1x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfvpar_1x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x2v_ser_p2(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfx_1x2v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfx_1x2v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_1x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfx_1x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfvpar_1x2v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfvpar_1x2v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_2x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_2x2v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfx_2x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfx_2x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfx_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfy_2x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfy_2x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfy_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfy_2x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfvpar_2x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfvpar_2x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_2x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_2x2v_ser_p2(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfx_2x2v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfx_2x2v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfx_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfy_2x2v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfy_2x2v_ser_p2(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfy_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfy_2x2v_ser_p2(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfvpar_2x2v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfvpar_2x2v_ser_p2(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_3x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_3x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *phi,
            const double *dualcurlbhatoverB, const double *rtg33inv, const double* bioverJB,
            const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_3x2v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfx_3x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfx_3x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfx_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfy_3x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfy_3x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfy_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfy_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfz_3x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfz_3x2v_ser_p1(
              const double *w, const double *dxv,
              
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_edge_surfz_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR, 
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_edge_surfz_3x2v_ser_p1(
                const double *w, const double *dxv,
                
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
                const double *bmag, const double *phi, const double *JfL, const double *JfR,
                double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_flux_surfvpar_3x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_flux_no_by_surfvpar_3x2v_ser_p1(
              const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_r, 
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
              const double *bmag, const double *phi, const double *JfL, const double *JfR, 
              double* GKYL_RESTRICT flux_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *flux_surf_l, const double *flux_surf_r, 
              double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *flux_surf_edge, const double *flux_surf_skin, 
              const int edge, double* GKYL_RESTRICT out); 


EXTERN_C_END
