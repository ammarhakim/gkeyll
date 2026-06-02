/*
 * gkyl_bc_sheath_gyrokinetic_gyraze_surrogate.h  -  GYRAZE surrogate model public API generated from gkeyll_sheath_ai @ 5459dfb
 * Sources:
 *   nn model      : gkeyll_sheath_ai/model/nn_model_conv_MPE.pth
 *   normalization : gkeyll_sheath_ai/model/normalization_conv_MPE.npz
 *   svm model     : gkeyll_sheath_ai/model/svm_model.pkl
 */
#pragma once

#include <gkyl_const.h>
#include <gkyl_util.h>
#include <kann.h>

/* Number of points in the fixed mu-grid. */
#define SRGRZ_N_MU 20

/* --- weight storage type --- **/
typedef struct {
  double W0[225], b0[75];
  double W1[5625], b1[75];
  double W2[5625], b2[75];
  double W3[1500], b3[20];
  double Y_mu[20], Y_sigma[20];
  double MU_GRID[20];
} srgrz_weights_t;

EXTERN_C_BEG

/**
 * Set the active KANN model used by all subsequent srgrz_predict / srgrz_eval_fact calls.
 * Call once at construction with the loaded model, and again with NULL at destruction.
 * The caller retains ownership of the model.
 *
 * @param model Pointer to a loaded kann_t, or NULL to clear.
 */
void bc_sheath_gyrokinetic_srgrz_set_model(kann_t *model);

/**
 * Runs the NN regression; writes SRGRZ_N_MU predicted v_par_cut values into out[].
 *
 * @param alpha: impact angle in degrees
 * @param gamma: normalised plasma density parameter
 * @param phi:   normalised sheath potential drop (e * (phi - phi_wall) / T_e)
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_infer(double alpha, double gamma, double phi, double out[SRGRZ_N_MU]);

/**
 * Linear interpolation of vcut[SRGRZ_N_MU] (on the fixed mu-grid) onto
 * mu_new[n]; results are written into out[n]. Clamps at the grid boundaries.
 *
 * @param vcut:    input array of size SRGRZ_N_MU containing values at the fixed mu-grid
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param mu_ref:  reference mu value for normalisation (e.g. temperature / Bmag)
 * @param out:     output array of size n where interpolated values are written
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_interp(const double *vcut, const double *mu_new, int n, double mu_ref, double *out);

/**
 * Eval the sheath BC surrogate to get vcutsq, normalises by 2 * e * (phi - phi_wall) / T.
 *
 * @param model:   pointer to the loaded KANN model
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param phi:     sheath potential (V)
 * @param phi_wall: wall potential (V)
 * @param density:  electron density (m^-3)
 * @param temperature:  electron temperature (eV)
 * @param q2Dm:     2 x charge-to-mass ratio (C/kg)
 * @param bmag:    magnetic field strength (T)
 * @param impact_angle: magnetic impact angle (radians)
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval(kann_t *model, const double *mu_new, int n, double phi, double phi_wall,
    double density, double temperature, double q2Dm, double bmag, double impact_angle, double *out);

EXTERN_C_END

