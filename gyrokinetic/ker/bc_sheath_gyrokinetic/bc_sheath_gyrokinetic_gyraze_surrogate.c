/*
 * bc_sheath_gyrokinetic_gyraze_surrogate.c  --  GYRAZE surrogate model, Gkeyll kernel generated from /Users/ahoffman/gkeyll_sheath_ai/src_py
 */
#include "gkyl_bc_sheath_gyrokinetic_gyraze_surrogate.h"

/*   Mu-grid */
/* Fixed evaluation mu-grid (20 points) */
static const double MU_GRID[20] = {0.00000000f, 0.02000000f, 0.08000000f, 0.18000000f, 0.32000000f, 0.50000000f, 0.72000000f, 0.98000000f, 1.28000000f, 1.62000000f, 2.00000000f, 2.42000000f, 2.88000000f, 3.38000000f, 3.92000000f, 4.50000000f, 5.12000000f, 5.78000000f, 6.48000000f, 7.22000000f};

/*  Neural-network forward pass  */
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_predict(double alpha, double gamma, double phi, double out[20])
{
}

static double svm_score(double * input) {
    return 0;
}

GKYL_CU_DH int bc_sheath_gyrokinetic_srgrz_converged(double alpha, double gamma, double phi)
{
    double input[3] = {alpha, gamma, phi};
    return (svm_score(input) >= 0.5) ? 1 : 0;
}

GKYL_CU_DH double *srgrz_grid(double *out)
{
    for (int i = 0; i < SRGRZ_N_MU; i++) {
        out[i] = MU_GRID[i];
    }
    return out;
}

GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_interp(const double *vcut, const double *mu_new, int n, double *out)
{
    int ng = SRGRZ_N_MU;
    for (int i = 0; i < n; i++) {
        double mu = mu_new[i];
        if (mu <= MU_GRID[0])          { out[i] = vcut[0];          continue; }
        if (mu >= MU_GRID[ng - 1])     { out[i] = vcut[ng - 1];     continue; }
        /* binary search for the bracketing interval */
        int lo = 0, hi = ng - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) >> 1;
            if (MU_GRID[mid] <= mu) lo = mid; else hi = mid;
        }
        double t = (mu - MU_GRID[lo]) / (MU_GRID[hi] - MU_GRID[lo]);
        out[i] = vcut[lo] + t * (vcut[hi] - vcut[lo]);
    }
}

GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval(const double *mu_new, int n, double alpha, double gamma, double phi, double *out)
{
    double vcut[SRGRZ_N_MU];
    bc_sheath_gyrokinetic_srgrz_predict(alpha, gamma, phi, vcut);
    bc_sheath_gyrokinetic_srgrz_interp(vcut, mu_new, n, out);
}

GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical(const double *mu_new, int n, double phi, double phi_wall, double dens_e,
    double temp_e, double m_e, double bmag, double bimpact_angle, double *out)
{
    double munorm[n];
    for (int i = 0; i < n; i++) {
        munorm[i] = mu_new[i] * bmag / temp_e;
    }
    double gamma   = (1.0 / bmag) * sqrt(m_e * dens_e / GKYL_EPSILON0);
    double phinorm = (GKYL_ELEMENTARY_CHARGE * phi) / temp_e;
    double alpha = bimpact_angle * 180/M_PI;
    bc_sheath_gyrokinetic_srgrz_eval(munorm, n, alpha, gamma, phinorm, out);
}
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical_vcut_fact(const double *mu_new,  int n, double phi, double phi_wall,
    double dens_e, double temp_e, double m_e, double bmag, double bimpact_angle, double *out)
{
    double vcut_const = sqrt(GKYL_ELEMENTARY_CHARGE * (phi - phi_wall) /temp_e);
    bc_sheath_gyrokinetic_srgrz_eval_physical(mu_new, n, phi, phi_wall, dens_e, temp_e, m_e, bmag, bimpact_angle, out);
    for (int i = 0; i < n; i++) {
        out[i] = out[i]/vcut_const;
    }
}
GKYL_CU_DH void bc_sheath_gyrokinetic_srgrz_eval_physical_vcut_fact_converged(const double *mu_new,  int n, double phi, double phi_wall,
    double dens_e, double temp_e, double m_e, double bmag, double bimpact_angle, double *out)
{
    double gamma   = (1.0 / bmag) * sqrt(m_e * dens_e / GKYL_EPSILON0);
    double phinorm = (GKYL_ELEMENTARY_CHARGE * phi) / temp_e;
    double alpha = bimpact_angle * 180/M_PI;
    int converged = bc_sheath_gyrokinetic_srgrz_converged(alpha, gamma, phinorm);
    if (converged) {
        bc_sheath_gyrokinetic_srgrz_eval_physical(mu_new, n, phi, phi_wall, dens_e, temp_e, m_e, bmag, bimpact_angle, out);
    } else {
        for (int i = 0; i < n; i++) {
            out[i] = 0.0;
        }
    }
}
