#pragma once

#include <gkyl_const.h>
#include <gkyl_kann_net.h>
#include <gkyl_util.h>
#include <math.h>

// Number of points in the fixed mu-grid that the GYRAZE sheath surrogate
// predicts v_par_cut on (see gkeyll_sheath_ai).
#define SRGRZ_N_MU 20

// approximation for inverse Langevin function
GKYL_CU_DH
static inline double invL(double x) {
  // from Kroger
  return (3.*x-x*x*x*(6. + x*x - 2.*x*x*x*x)/5.)/(1.-x*x);
}

// Linear interpolation of the surrogate output vcut[SRGRZ_N_MU] (defined on the
// fixed mu-grid below) onto mu_new[n], normalising mu by mu_ref (= T/B). Results
// are written into out[n]; clamps at the grid boundaries.
GKYL_CU_DH
static inline void
bc_sheath_gyrokinetic_surr_interpf(const float *vcut, const double *mu_new, int n, double mu_ref, double *out)
{
  const double mu_grid[SRGRZ_N_MU] = {
    0.00, 0.02, 0.08, 0.18, 0.32, 0.50, 0.72, 0.98, 1.28, 1.62,
    2.00, 2.42, 2.88, 3.38, 3.92, 4.50, 5.12, 5.78, 6.48, 7.22};
  const int ng = SRGRZ_N_MU;
  for (int i = 0; i < n; i++) {
    double mu = mu_new[i]/mu_ref;
    if (mu <= mu_grid[0])      { out[i] = vcut[0];      continue; }
    if (mu >= mu_grid[ng - 1]) { out[i] = vcut[ng - 1]; continue; }
    // binary search for the bracketing interval
    int lo = 0, hi = ng - 1;
    while (hi - lo > 1) {
      int mid = (lo + hi) >> 1;
      if (mu_grid[mid] <= mu) lo = mid; else hi = mid;
    }
    double t = (mu - mu_grid[lo]) / (mu_grid[hi] - mu_grid[lo]);
    out[i] = vcut[lo] + t * (vcut[hi] - vcut[lo]);
  }
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

