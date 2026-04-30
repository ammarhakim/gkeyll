#pragma once

#include <stdbool.h>

#define TOV_MAX_POINTS 200000

// Output struct for BL evaluation.
struct tov_ultra_rel_eval_bl {
    double P;
    double e;
    double m;
    double Phi;
};

// Output struct for fluid evaluation (CKS/SKS).
struct tov_ultra_rel_eval_fluid {
    double P;
    double e;
};

// TOV solution object.
struct gkyl_tov_ultra_rel;

/**
 * Create a new ultra-relativistic TOV solution object.
 * Integrates the TOV equations for EOS P = (gas_gamma - 1) e.
 *
 * @param gas_gamma Adiabatic index
 * @param e_c Central energy density
 * @param dr Radial step size
 * @return TOV solution object
 */
struct gkyl_tov_ultra_rel* gkyl_tov_ultra_rel_new(double gas_gamma, double e_c, double dr);


/**
 * Evaluate TOV solution in Boyer-Lindquist coordinates at areal radius r.
 *
 * @param tov TOV solution object
 * @param r Areal radius
 * @param out Output struct
 * @return true on success
 */
bool gkyl_tov_ultra_rel_eval_bl(const struct gkyl_tov_ultra_rel *tov, double r, struct tov_ultra_rel_eval_bl *out);

/**
 * Evaluate fluid variables at position (x,y,z) relative to star center.
 * r = sqrt(x^2+y^2+z^2) is the areal radius shared between Schwarzschild and CKS.
 *
 * @param tov TOV solution object
 * @param x, y, z Position relative to star center
 * @param out Output struct
 * @return true on success
 */
bool gkyl_tov_ultra_rel_eval(const struct gkyl_tov_ultra_rel *tov, double x, double y, double z,
    struct tov_ultra_rel_eval_fluid *out);

/**
 * Get gravitational mass of the star.
 */
double gkyl_tov_ultra_rel_star_mass(const struct gkyl_tov_ultra_rel *tov);

/**
 * Get areal radius of the stellar surface.
 */
double gkyl_tov_ultra_rel_star_radius(const struct gkyl_tov_ultra_rel *tov);

/**
 * Refresh lapse potential Phi and enclosed mass m from the current
 * one-dimensional GR TOV state.
 *
 * The state q is assumed to be contiguous with 8 components per cell.
 */
void gkyl_gr_tov_refresh_geometry_from_state(double gas_gamma, double p_atm, int ncells, double *q);

/**
 * Release TOV solution object and free memory.
 */
void gkyl_tov_ultra_rel_release(struct gkyl_tov_ultra_rel *tov);
