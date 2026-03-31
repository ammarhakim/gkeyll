/* gkyl_bc_sheath_gyrokinetic_gyraze_surrogate.h  –  GYRAZE surrogate model public API generated from gkeyll_sheath_ai @ cbd6ef1 */
#pragma once

#include <math.h>
#include <gkyl_const.h>
#include <gkyl_util.h>

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
 * Returns 1 if GYRAZE is predicted to converge, 0 otherwise.
 *
 * @param alpha: impact angle in degrees
 * @param gamma: normalised plasma density parameter
 * @param phi:   normalised sheath potential drop (e * (phi - phi_wall) / T_e)
 */
GKYL_CU_DH int bc_sheath_gyrokinetic_srgrz_converged(double alpha, double gamma, double phi);

/**
 * Runs the NN regression; writes SRGRZ_N_MU predicted v_par_cut values into out[].
 *
 * @param alpha: impact angle in degrees
 * @param gamma: normalised plasma density parameter
 * @param phi:   normalised sheath potential drop (e * (phi - phi_wall) / T_e)
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_predict(double alpha, double gamma, double phi, double out[SRGRZ_N_MU]);

/**
 * Copies the SRGRZ_N_MU-element mu-grid into out[] and returns out.
 *
 * @param out: output array of size SRGRZ_N_MU
 */
GKYL_CU_DH double *bc_sheath_gyrokinetic_srgrz_grid(double *out);

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
 * Returns the prediction of a custom mu grid of size n
 *
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param mu_ref:  reference mu value for normalisation (e.g. temperature / Bmag)
 * @param alpha:   impact angle in degrees
 * @param gamma:   normalised plasma density parameter
 * @param phi:     normalised sheath potential drop (e * (phi - phi_wall) / T_e)
 * @param out:     output array of size n where interpolated values are written
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval(const double *mu_new, int n, double mu_ref, double alpha, double gamma, double phi, double *out);

/**
 * Like bc_sheath_gyrokinetic_srgrz_eval, but projects (alpha, gamma, phi) onto the nearest
 * convergent point in parameter space when GYRAZE is predicted not to converge.
 * The projection minimises svm_score(x)^2 + 1e-3*||x-x0||^2 via gradient
 * descent with Armijo backtracking (mirrors find_nearest() in surrogate_proj.py).
 *
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param mu_ref:  reference mu value for normalisation (e.g. temperature / Bmag)
 * @param alpha:   impact angle in degrees
 * @param gamma:   normalised plasma density parameter
 * @param phi:     normalised sheath potential drop (e * (phi - phi_wall) / T_e)
 * @param out:     output array of size n where interpolated values are written
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_proj_eval(const double *mu_new, int n, double mu_ref,
    double alpha, double gamma, double phi, double *out);

/**
 * Converts from physical parameters and evaluates on a custom mu grid.
 * Conversion formulas:
 *   munorm  = mu*Bmag / temperature
 *   gamma   = (1/Bmag) * sqrt(m_e * density / eps0)
 *   phinorm = e * phi / temperature
 *
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
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical(const double *mu_new, int n, double phi, double phi_wall,
    double density, double temperature, double q2Dm, double bmag, double impact_angle, double *out);

/**
 * Same as srgrz_eval_physical, but normalises output by sqrt(2 * e * (phi - phi_wall) / mass)
 *
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
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical_vcut_fact(const double *mu_new, int n, double phi, double phi_wall,
    double density, double temperature, double q2Dm, double bmag, double impact_angle, double *out);

/**
 * Same as srgrz_eval_physical_vcut_fact, but normalises return 0 if gyraze is not converging.
 *
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
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical_vcut_fact_converged(const double *mu_new, int n, double phi, double phi_wall,
    double density, double temperature, double q2Dm, double bmag, double impact_angle, double *out);

EXTERN_C_END

