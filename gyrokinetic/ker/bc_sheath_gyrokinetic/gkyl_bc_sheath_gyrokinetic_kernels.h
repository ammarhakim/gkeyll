#pragma once

#include <gkyl_kann_net.h>
#include <gkyl_bc_sheath_gyrokinetic_gyraze_surrogate.h> 
#include <gkyl_util.h>
#include <math.h>

// approximation for inverse Langevin function 
GKYL_CU_DH
static inline double invL(double x) {
  // from Kroger 
  return (3.*x-x*x*x*(6. + x*x - 2.*x*x*x*x)/5.)/(1.-x*x); 
}

EXTERN_C_BEG

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_2x2v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_2x2v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1(const double *vmap, const double *vcutsq, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_const_lower_1x2v_ser_p1(const double *phi, const double *phi_wall, double q2Dm, double *vcutSq_n) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_const_upper_1x2v_ser_p1(const double *phi, const double *phi_wall, double q2Dm, double *vcutSq_n) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_const_lower_2x2v_ser_p1(const double *phi, const double *phi_wall, double q2Dm, double *vcutSq_n) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_const_upper_2x2v_ser_p1(const double *phi, const double *phi_wall, double q2Dm, double *vcutSq_n) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_const_lower_3x2v_ser_p1(const double *phi, const double *phi_wall, double q2Dm, double *vcutSq_n) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_const_upper_3x2v_ser_p1(const double *phi, const double *phi_wall, double q2Dm, double *vcutSq_n) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_surr_lower_1x2v_ser_p1(const double *vmap, const float *nn_out, int n_out, const double *temperature, const double *bmag, double *vcutsq_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_surr_upper_1x2v_ser_p1(const double *vmap, const float *nn_out, int n_out, const double *temperature, const double *bmag, double *vcutsq_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_surr_lower_2x2v_ser_p1(const double *vmap, const float *nn_out, int n_out, const double *temperature, const double *bmag, double *vcutsq_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_surr_upper_2x2v_ser_p1(const double *vmap, const float *nn_out, int n_out, const double *temperature, const double *bmag, double *vcutsq_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_surr_lower_3x2v_ser_p1(const double *vmap, const float *nn_out, int n_out, const double *temperature, const double *bmag, double *vcutsq_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_vcutsq_surr_upper_3x2v_ser_p1(const double *vmap, const float *nn_out, int n_out, const double *temperature, const double *bmag, double *vcutsq_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_build_input_lower_1x2v_ser_p1(const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, int n_inp, float *nn_inp_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_build_input_upper_1x2v_ser_p1(const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, int n_inp, float *nn_inp_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_build_input_lower_2x2v_ser_p1(const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, int n_inp, float *nn_inp_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_build_input_upper_2x2v_ser_p1(const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, int n_inp, float *nn_inp_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_build_input_lower_3x2v_ser_p1(const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, int n_inp, float *nn_inp_out) ; 
GKYL_CU_DH void bc_sheath_gyrokinetic_build_input_upper_3x2v_ser_p1(const double *phi, const double *phi_wall, const double *density, const double *temperature, const double *bmag, const double *bimpact_angle, int n_inp, float *nn_inp_out) ; 

EXTERN_C_END

