#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <gkyl_alloc.h>

#include "tov_solver.h"

struct gkyl_tov {
    int N; // how many radial points

    double K; // polytripic EOS - normalizaiton constant
    double Gamma; // polytripic EOS - adiabatic index

    //radial profiles
    double *r_areal; // areal radius
    double *m; 
    double *P;
    double *rho; // rest mass density
    double *eps; // specific internal energy
    double *e; // total energy density
    double *Phi; // potential (gravitaional potential, metric potential: alpha = exp(Phi))

    // global star properties
    double M_star; // total mass of the star (gravitational mass)
    double R_areal; // areal radius of the star surface
    double rho_c; // density at the center
    double P_c; // pressure at the center
    

    double dr; // radial step

};

// EOS is a single polytrope, so: P = K rho^Gamma

static inline double eos_pressure(double K, double Gamma, double rho)
{
    double P;
    if (rho > 0) {
        P = K * pow(rho, Gamma);
    }
    else P = 0.0;
    return P;
}

// eps = K rho^(gamma-1) / (gamma - 1)
static inline double eos_eps(double K, double Gamma, double rho)
{
    double eps;
    if (rho > 0) {
        eps = K * pow(rho, Gamma - 1.0) / (Gamma - 1.0);
    }
    else eps = 0.0;
    return eps;
}

// e = rho * (1 + eps)
static inline double eos_energy(double K, double Gamma, double rho)
{
    double e;
    if (rho > 0) {
        e = rho + K * pow(rho, Gamma) / (Gamma - 1.0);
    }
    else e = 0.0;
    return e;
}

// rho = (P / K)^(1/Gamma)
static inline double eos_rho_from_pressure(double K, double Gamma, double P)
{
    double rho;
    if (P > 0) {
        rho = pow(P / K, 1.0 / Gamma);
    }
    else rho = 0.0;
    return rho;
}

// TOV integrator
static void tov_rhs(double r, const double y[3], double f[3], double K, double Gamma)
{
    double m_enc = y[0]; // mass enclosed
    double P = y[1]; 

    if (P <= 0.0 || r <= 0.0) {
        f[0] = f[1] = f[2] = 0.0;
        return;
    }

    double rho = eos_rho_from_pressure(K, Gamma, P);
    double e = eos_energy(K, Gamma, rho);

    if (r * (r - 2.0 * m_enc) <= 0.0) { // denominator from dP/dr
        f[0] = f[1] = f[2] = 0.0;  // if (r * (r - 2.0 * m_enc) <= 0.0) then we're at/inside the horizon - unphysical for a stable star
        return;
    }
    
    f[0] = 4.0 * M_PI * r * r * e; // dm/dr RHS
    f[1] = -(e + P) * (m_enc + 4.0 * M_PI * r * r * r * P) / (r * (r - 2.0 * m_enc)); // dP/dr RHS
    f[2] = (m_enc + 4.0 * M_PI * r * r * r * P) / (r * (r - 2.0 * m_enc)); // dPhi/dr 

}

static void implicit_tov_source_update(double r, double y[3], double h, double K, double Gamma) // it's an explicit method, not implicit, it's ssp-rk3
{
    double y_old[3], y_stage1[3], y_stage2[3], y_new[3];
    double f[3];

    for (int i = 0; i < 3; i ++) {
        y_old[i] = y[i]; // y_old = y^n
    }

    // stage 1 — Euler step from y_old
    tov_rhs(r, y_old, f, K, Gamma);
    for (int i = 0; i < 3; i++) {
        y_stage1[i] = y_old[i] + h * f[i];
    }

    // stage 2 — RHS at stage1, then mix with y_old
    tov_rhs(r + h, y_stage1, f, K, Gamma);
    for (int i = 0; i < 3; i++) {
        //y_stage2[i] = 0.75*y_old[i] + 0.25*(y_old[i] + h * f[i]);
        y_stage2[i] = 0.75*y_old[i] + 0.25*(y_stage1[i] + h * f[i]);
    }

    // stage 3 — RHS at stage2, final mix
    tov_rhs(r + 0.5*h, y_stage2, f, K, Gamma);
    for (int i = 0; i < 3; i++) {
        y_new[i] = (1.0/3.0)*y_old[i] + (2.0/3.0)*(y_stage2[i] + h * f[i]);
    }

    for (int i = 0; i < 3; i ++) {
        y[i] = y_new[i];
    }

}

// Interpolate on the Boyler Lindquist 
static double interp(const struct gkyl_tov *tov, const double *val, double r)
{
    if (r <= tov->r_areal[0]) {
        return val[0];
    }
    if (r >= tov->r_areal[tov->N - 1]) {
        return val[tov->N - 1];
    }

    int lo = 0, hi = tov->N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (tov->r_areal[mid] <= r) lo = mid;
        else hi = mid;
    }

    double dr = tov->r_areal[hi] - tov->r_areal[lo];
    double t = 0.0;
    if (dr > 0.0) {
        t = (r - tov->r_areal[lo]) / dr;
    }
    else t = 0.0;

    return (1.0 - t) * val[lo] + t * val[hi];
}

// eval TOV in BL coords at arbitrary radius
bool gkyl_tov_eval_bl(const struct gkyl_tov *tov, double r, struct tov_eval_bl *out)
{
    if (!tov || !out || r < 0.0) {
        return false;
    }

    out->rho = interp(tov, tov->rho, r);
    out->P = interp(tov, tov->P, r);
    out->eps = interp(tov, tov->eps, r);
    out->e = interp(tov, tov->e, r);
    out->m = interp(tov, tov->m, r);
    out->Phi = interp(tov, tov->Phi, r);

    return true;
}


// eval fluid relative to star center - r is the areal radius shared between Schwarzschild and Kerr-Schild
bool gkyl_tov_eval(const struct gkyl_tov *tov, double x, double y, double z,
    struct tov_eval_fluid *out)
{
    if (!tov || !out) return false;
    double r = sqrt(x*x + y*y + z*z);
    
    if (r >= tov->R_areal) {
        out->rho = 0.0;
        out->P   = 0.0;
        out->eps = 0.0;
        out->e   = 0.0;
        return true;
    }
    
    out->rho = interp(tov, tov->rho, r);
    out->P   = interp(tov, tov->P,   r);
    out->eps = interp(tov, tov->eps, r);
    out->e   = interp(tov, tov->e,   r);
    return true;
}

// helpers to acces star mass and radius
double gkyl_tov_star_mass(const struct gkyl_tov *tov)
{
    if(tov) {
        return tov->M_star;
    }
    else {
        return 0.0;
    }
}

double gkyl_tov_star_radius(const struct gkyl_tov *tov)
{
    if(tov) {
        return tov->R_areal;
    }
    else {
        return 0.0;
    }
}

//main
struct gkyl_tov * gkyl_tov_new(double K, double Gamma, double rho_c, double dr)
{
    struct gkyl_tov *tov = (struct gkyl_tov *) gkyl_malloc(sizeof(struct gkyl_tov));

    if (!tov) {
        return NULL;
    }

    tov->K = K;
    tov->Gamma = Gamma;
    tov->rho_c = rho_c;
    tov->dr = dr;

    tov->r_areal = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->m = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->P  = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->rho = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->eps = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->e = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));
    tov->Phi = gkyl_malloc(sizeof(double[TOV_MAX_POINTS]));

    if (!tov->r_areal || !tov->m || !tov->P || !tov->rho || !tov->eps || !tov->e || !tov->Phi) {
        gkyl_tov_solution_release(tov);
        return NULL;
    }

    // First we integrate the TOV eqns in the areal coordinates
    double P_c = eos_pressure(K, Gamma, rho_c);
    double e_c = eos_energy(K, Gamma, rho_c);
    tov->P_c = P_c;

    // center
    tov->r_areal[0] = 0.0;
    tov->m[0] = 0.0;
    tov->P[0] = P_c;
    tov->rho[0] = rho_c;
    tov->eps[0] = eos_eps(K, Gamma, rho_c);
    tov->e[0] = e_c;
    tov->Phi[0] = 0.0;  // "provisional"; shifted below

    // Taylor expansion at r = dr  (since we can't start the integration at 0)
    double r_n = dr;
    double y[3];
    y[0] = (4.0 / 3.0) * M_PI * e_c * r_n * r_n * r_n; //m 
    y[1] = P_c - (2.0 * M_PI / 3.0) * (e_c + P_c) * (e_c / 3.0 + P_c) * r_n * r_n; //P
    y[2] = (2.0 * M_PI / 3.0) * (e_c / 3.0 + P_c) * r_n * r_n; //Phi

    if (y[1] <= 0.0) {
            y[1] = 0.0;
        }

    double rho_n = eos_rho_from_pressure(K, Gamma, y[1]);
    tov->r_areal[1] = r_n;
    tov->m[1]       = y[0];
    tov->P[1]       = y[1];
    tov->rho[1]     = rho_n;
    tov->eps[1]     = eos_eps(K, Gamma, rho_n);
    tov->e[1]       = eos_energy(K, Gamma, rho_n);
    tov->Phi[1]     = y[2];
    int n = 2;     

    // Actual integration call - call to the SSP-RK3 (implicit_tov_source_update)
    double r_max = 1500.0;
    while (r_n < r_max && y[1] > 0.0 && n < TOV_MAX_POINTS) {
        implicit_tov_source_update(r_n, y, dr, K, Gamma);
        r_n += dr;

        if (y[1] <= 0.0) {
            y[1] = 0.0;
        }

        rho_n = eos_rho_from_pressure(K, Gamma, y[1]);

        tov->r_areal[n] = r_n;
        tov->m[n] = y[0];
        tov->P[n]  = y[1];
        tov->rho[n] = rho_n;
        tov->eps[n] = eos_eps(K, Gamma, rho_n);
        tov->e[n] = eos_energy(K, Gamma, rho_n);
        tov->Phi[n] = y[2];
        n++;

        // printf("r=%g m=%g P=%g Phi=%g\n", r_n, y[0], y[1], y[2]);
        if (y[1] <= 0.0) break;
    }

    tov->N = n;
    int i_surf = n - 1; // stellar surface
    tov->M_star  = tov->m[i_surf];
    tov->R_areal = tov->r_areal[i_surf];

    // We fix lapse: Phi(R) = 0.5 ln(1 - 2M/R). Lapse = e^Phi. At the stellar surface, the interior must match exterior Schwarzschild.
    // g_tt = -(1 - 2M/r) = -e^2Phi. => Phi(R) = 0.5 * ln(1 - 2M/R)
    double Phi_exact = 0.5 * log(1.0 - 2.0 * tov->M_star / tov->R_areal);
    //double Phi_exact = 0.0;
    double Phi_shift = Phi_exact - tov->Phi[i_surf]; // We integrated with Phi(0) = Phi_c = 0 before (arbitrary), so we need to shift the profile by a constant
    for (int i = 0; i < n; i++) {
        tov->Phi[i] += Phi_shift;
    }   
    // After the Phi shift -- sanity check that the convention is what we expect
    // double alpha_surf_int = exp(tov->Phi[i_surf]) / sqrt(1.0 + 2.0*tov->M_star/tov->R_areal);
    // double alpha_surf_ext = 1.0 / sqrt(1.0 + 2.0*tov->M_star/tov->R_areal);
    double alpha_surf_int = exp(tov->Phi[i_surf]);
    double alpha_surf_ext = sqrt(1.0 - 2.0*tov->M_star/tov->R_areal);

    assert(fabs(alpha_surf_int - alpha_surf_ext) < 1e-12 && "Lapse doesn't match on the surface"); // this should pass

    // double V_star = 2.0 * tov->M_star / tov->R_areal;
    // double a_surf = exp(tov->Phi[i_surf]);               // BL lapse at R*
    // double X_surf = a_surf*a_surf * V_star*V_star / (1.0 - V_star);
    // double alpha_cks_int = a_surf / sqrt(1.0 - X_surf);  // CKS interior lapse
    // double alpha_cks_ext = 1.0 / sqrt(1.0 + V_star);     // exterior KS lapse
    // assert(fabs(alpha_cks_int - alpha_cks_ext) < 1e-12); // this should pass
    //assert(fabs(alpha_surf_int - 1.0 / sqrt(1.0 + 2.0*tov->M_star/tov->R_areal)) < 1e-12);

    // This passes iff Phi(R_star) = 0, which is the CKS normalization convention.


    // Extend to vacuum exterior out to r_max
    while (r_n < r_max && n < TOV_MAX_POINTS) {
        r_n += dr;
        double V_ext = 2.0 * tov->M_star / r_n;
        tov->r_areal[n] = r_n;
        tov->m[n]       = tov->M_star;
        tov->P[n]       = 0.0;
        tov->rho[n]     = 0.0;
        tov->eps[n]     = 0.0;
        tov->e[n]       = 0.0;

        tov->Phi[n] = (1.0 - V_ext > 0.0)
                  ? 0.5 * log(1.0 - V_ext)
                  : tov->Phi[n-1];
        n++;
    }

    tov->N = n;

    // for (int i=0; i < tov->N && i < 10; ++i) {
    //     printf("%d r=%e P=%e rho=%e m=%e Phi=%e\n",
    //         i, tov->r_areal[i], tov->P[i], tov->rho[i], tov->m[i], tov->Phi[i]);
    // }

    // FILE *fp = fopen("tov_profile.txt", "w");
    // if (fp) {
    //     fprintf(fp, "# r rho P eps e m Phi alpha\n");
    //     for (int i = 0; i < tov->N && tov->r_areal[i] < tov->R_areal + 1.0; i++) {
    //         fprintf(fp, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
    //             tov->r_areal[i], tov->rho[i], tov->P[i],
    //             tov->eps[i], tov->e[i], tov->m[i], tov->Phi[i], exp(tov->Phi[i]));
    //     }
    //     fclose(fp);
    // }

    return tov;
}

void gkyl_tov_solution_release(struct gkyl_tov *tov)
{
    if (!tov) return;
    gkyl_free(tov->r_areal);
    gkyl_free(tov->m);
    gkyl_free(tov->P);
    gkyl_free(tov->rho);
    gkyl_free(tov->eps);
    gkyl_free(tov->e);
    gkyl_free(tov->Phi);
    gkyl_free(tov);
}

