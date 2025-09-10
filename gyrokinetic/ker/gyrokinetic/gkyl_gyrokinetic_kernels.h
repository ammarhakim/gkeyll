#pragma once

#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x1v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x1v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_1x1v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_1x1v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x1v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x1v_ser_p2(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_1x1v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_1x1v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x2v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_1x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_1x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_1x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_1x2v_ser_p2(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_1x2v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_1x2v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_2x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_2x2v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfy_2x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfy_2x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_2x2v_ser_p2(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_2x2v_ser_p2(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_2x2v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_2x2v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfy_2x2v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfy_2x2v_ser_p2(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_2x2v_ser_p2(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_3x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_vol_3x2v_ser_p1(const double *w, const double *dxv,
            const double *vmap, const double *vmapSq, const double q_, const double m_,
            const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi,
            const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_step2_vol_3x2v_ser_p1(const double *w, const double *dxv,
            const double q_, const double m_, const double *apardot, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_3x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_3x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfx_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfy_3x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfy_3x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfz_3x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfz_3x2v_ser_p1(const double *w, const double *dxv,
                const double *vmap, const double *vmapSq, const double q_, const double m_,
                const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
                const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfz_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_, 
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap, const double *vmapSq, const double q_, const double m_,
              const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i,
              const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double gyrokinetic_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, 
              const double *alpha_surf_l, const double *alpha_surf_r, 
              const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
              const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
              const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
              const double *vmap_prime_edge, const double *vmap_prime_skin, 
              const double *alpha_surf_edge, const double *alpha_surf_skin, 
              const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
              const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
              const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 


EXTERN_C_END
