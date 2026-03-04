#pragma once

#include <math.h>
#include <gkyl_const.h>
#include <gkyl_util.h>

/* Number of points in the fixed mu-grid. */
#define SRGRZ_N_MU 20

EXTERN_C_BEG

/**
 * Returns 1 if GYRAZE is predicted to converge, 0 otherwise.
 *
 * @param alpha: impact angle in degrees
 * @param gamma: normalised plasma density parameter
 * @param phi:   normalised sheath potential drop (e * (phi - phi_wall) / T_e)
 */
GKYL_CU_DH int  bc_sheath_gyrokinetic_srgrz_converged(double alpha, double gamma, double phi);

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
 * @param out:     output array of size n where interpolated values are written
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_interp(const double *vcut, const double *mu_new, int n, double *out);

/**
 * Returns the prediction of a custom mu grid of size n
 *
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param alpha:   impact angle in degrees
 * @param gamma:   normalised plasma density parameter
 * @param phi:     normalised sheath potential drop (e * (phi - phi_wall) / T_e)
 * @param out:     output array of size n where interpolated values are written
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval(const double *mu_new, int n, double alpha, double gamma, double phi, double *out);

/**
 * Converts from physical parameters and evaluates on a custom mu grid.
 * Conversion formulas:
 *   munorm  = mu*Bmag / temp_e
 *   gamma   = (1/Bmag) * sqrt(m_e * dens_e / eps0)
 *   phinorm = e * phi / temp_e
 *
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param phi:     sheath potential (V)
 * @param phi_wall: wall potential (V)
 * @param dens_e:  electron density (m^-3)
 * @param temp_e:  electron temperature (eV)
 * @param m_e:     electron mass (kg)
 * @param bmag:    magnetic field strength (T)
 * @param bimpact_angle: magnetic impact angle (radians)
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical(const double *mu_new, int n, double phi, double phi_wall,
    double dens_e, double temp_e, double m_e, double bmag, double bimpact_angle, double *out);

/**
 * Same as srgrz_eval_physical, but normalises output by sqrt(2 * e * (phi - phi_wall) / m_e)
 *
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param phi:     sheath potential (V)
 * @param phi_wall: wall potential (V)
 * @param dens_e:  electron density (m^-3)
 * @param temp_e:  electron temperature (eV)
 * @param m_e:     electron mass (kg)
 * @param bmag:    magnetic field strength (T)
 * @param bimpact_angle: magnetic impact angle (radians)
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical_vcut_fact(const double *mu_new, int n, double phi, double phi_wall,
    double dens_e, double temp_e, double m_e, double bmag, double bimpact_angle, double *out);

/**
 * Same as srgrz_eval_physical_vcut_fact, but normalises return 0 if gyraze is not converging.
 *
 * @param mu_new:  input array of size n containing the new mu points
 * @param n:       number of points in mu_new and out
 * @param phi:     sheath potential (V)
 * @param phi_wall: wall potential (V)
 * @param dens_e:  electron density (m^-3)
 * @param temp_e:  electron temperature (eV)
 * @param m_e:     electron mass (kg)
 * @param bmag:    magnetic field strength (T)
 * @param bimpact_angle: magnetic impact angle (radians)
 */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical_vcut_fact_converged(const double *mu_new, int n, double phi, double phi_wall,
    double dens_e, double temp_e, double m_e, double bmag, double bimpact_angle, double *out);

