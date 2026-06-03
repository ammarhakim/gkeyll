#pragma once

#include <gkyl_kann_net.h>
#include <gkyl_bc_sheath_gyrokinetic_gyraze_surrogate.h> 
#include <gkyl_knutils.h>
#include <gkyl_util.h>
#include <math.h>

// approximation for inverse Langevin function 
GKYL_CU_DH
static inline double invL(double x) {
  // from Kroger 
  return (3.*x-x*x*x*(6. + x*x - 2.*x*x*x*x)/5.)/(1.-x*x); 
}

EXTERN_C_BEG

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_2x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_2x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_vparcut_calc_lower_1x2v_ser_p1(const double *vmap, const double *nn_out, const double *temperature, const double *bmag, double *vcut_fact_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vparcut_calc_upper_1x2v_ser_p1(const double *vmap, const double *nn_out, const double *temperature, const double *bmag, double *vcut_fact_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vparcut_calc_lower_2x2v_ser_p1(const double *vmap, const double *nn_out, const double *temperature, const double *bmag, double *vcut_fact_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vparcut_calc_upper_2x2v_ser_p1(const double *vmap, const double *nn_out, const double *temperature, const double *bmag, double *vcut_fact_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vparcut_calc_lower_3x2v_ser_p1(const double *vmap, const double *nn_out, const double *temperature, const double *bmag, double *vcut_fact_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vparcut_calc_upper_3x2v_ser_p1(const double *vmap, const double *nn_out, const double *temperature, const double *bmag, double *vcut_fact_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_infer_lower_1x2v_ser_p1(struct gkyl_kann_net *net, struct gkyl_kn_vec *inp_k, struct gkyl_kn_vec *out_k, const double *vmap, const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, double *out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_infer_upper_1x2v_ser_p1(struct gkyl_kann_net *net, struct gkyl_kn_vec *inp_k, struct gkyl_kn_vec *out_k, const double *vmap, const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, double *out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_infer_lower_2x2v_ser_p1(struct gkyl_kann_net *net, struct gkyl_kn_vec *inp_k, struct gkyl_kn_vec *out_k, const double *vmap, const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, double *out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_infer_upper_2x2v_ser_p1(struct gkyl_kann_net *net, struct gkyl_kn_vec *inp_k, struct gkyl_kn_vec *out_k, const double *vmap, const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, double *out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_infer_lower_3x2v_ser_p1(struct gkyl_kann_net *net, struct gkyl_kn_vec *inp_k, struct gkyl_kn_vec *out_k, const double *vmap, const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, double *out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_infer_upper_3x2v_ser_p1(struct gkyl_kann_net *net, struct gkyl_kn_vec *inp_k, struct gkyl_kn_vec *out_k, const double *vmap, const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, double *out) ; 

EXTERN_C_END

