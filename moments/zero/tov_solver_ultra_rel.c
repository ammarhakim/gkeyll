#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <gkyl_alloc.h>

#include "tov_solver_ultra_rel.h"

struct gkyl_tov_ultra_rel {
    int N; // how many radial points

    double Gamma; // Adiabatic index.

    //radial profiles
    double *r_areal; // areal radius
    double *m; 
    double *P;
    double *e; // energy density
    double *Phi; // potential (gravitaional potential, metric potential: alpha = exp(Phi))

    // global star properties
    double M_star; // total mass of the star (gravitational mass)
    double R_areal; // areal radius of the star surface
    double e_c; // central energy density
    double P_c; // pressure at the center
    double P_surface; // pressure cutoff used to define the stellar surface
    double r_max; // maximum tabulated radius
    bool surface_found; // true if the pressure reached P_surface before table end
    bool table_exhausted; // true if the table ended before finding P_surface
    

    double dr; // radial step

};

// Ultra-relativistic EOS: P = (Gamma - 1) e.

static inline double eos_pressure(double Gamma, double e)
{
    if (e > 0) {
        return (Gamma - 1.0) * e;
    }
    return 0.0;
}

static inline double eos_energy(double Gamma, double e)
{
    if (e > 0) {
        return e;
    }
    return 0.0;
}

static inline double eos_e_from_pressure(double Gamma, double P)
{
    if (P > 0) {
        return P / (Gamma - 1.0);
    }
    return 0.0;
}

// TOV integrator
static void tov_rhs(double r, const double y[3], double f[3], double Gamma)
{
    double m_enc = y[0]; // mass enclosed
    double P = y[1]; 

    if (P <= 0.0 || r <= 0.0) {
        f[0] = f[1] = f[2] = 0.0;
        return;
    }

    double e = eos_e_from_pressure(Gamma, P);

    if (r * (r - 2.0 * m_enc) <= 0.0) { // denominator from dP/dr
        f[0] = f[1] = f[2] = 0.0;  // if (r * (r - 2.0 * m_enc) <= 0.0) then we're at/inside the horizon - unphysical for a stable star
        return;
    }
    
    f[0] = 4.0 * M_PI * r * r * e; // dm/dr RHS
    f[1] = -(e + P) * (m_enc + 4.0 * M_PI * r * r * r * P) / (r * (r - 2.0 * m_enc)); // dP/dr RHS
    f[2] = (m_enc + 4.0 * M_PI * r * r * r * P) / (r * (r - 2.0 * m_enc)); // dPhi/dr 

}

static void implicit_tov_source_update(double r, double y[3], double h, double Gamma) // it's an explicit method, not implicit, it's ssp-rk3
{
    double y_old[3], y_stage1[3], y_stage2[3], y_new[3];
    double f[3];

    for (int i = 0; i < 3; i ++) {
        y_old[i] = y[i]; // y_old = y^n
    }

    // stage 1 — Euler step from y_old
    tov_rhs(r, y_old, f, Gamma);
    for (int i = 0; i < 3; i++) {
        y_stage1[i] = y_old[i] + h * f[i];
    }

    // stage 2 — RHS at stage1, then mix with y_old
    tov_rhs(r + h, y_stage1, f, Gamma);
    for (int i = 0; i < 3; i++) {
        //y_stage2[i] = 0.75*y_old[i] + 0.25*(y_old[i] + h * f[i]);
        y_stage2[i] = 0.75*y_old[i] + 0.25*(y_stage1[i] + h * f[i]);
    }

    // stage 3 — RHS at stage2, final mix
    tov_rhs(r + 0.5*h, y_stage2, f, Gamma);
    for (int i = 0; i < 3; i++) {
        y_new[i] = (1.0/3.0)*y_old[i] + (2.0/3.0)*(y_stage2[i] + h * f[i]);
    }

    for (int i = 0; i < 3; i ++) {
        y[i] = y_new[i];
    }

}

// Interpolate on the Boyler Lindquist 
static double interp(const struct gkyl_tov_ultra_rel *tov, const double *val, double r)
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
bool gkyl_tov_ultra_rel_eval_bl(const struct gkyl_tov_ultra_rel *tov, double r, struct tov_ultra_rel_eval_bl *out)
{
    if (!tov || !out || r < 0.0) {
        return false;
    }

    out->P = interp(tov, tov->P, r);
    out->e = interp(tov, tov->e, r);
    out->m = interp(tov, tov->m, r);
    out->Phi = interp(tov, tov->Phi, r);

    return true;
}


// eval fluid relative to star center - r is the areal radius shared between Schwarzschild and Kerr-Schild
bool gkyl_tov_ultra_rel_eval(const struct gkyl_tov_ultra_rel *tov, double x, double y, double z,
    struct tov_ultra_rel_eval_fluid *out)
{
    if (!tov || !out) return false;
    double r = sqrt(x*x + y*y + z*z);
    
    if (r >= tov->R_areal) {
        out->P   = 0.0;
        out->e   = 0.0;
        return true;
    }
    
    out->P   = interp(tov, tov->P,   r);
    out->e   = interp(tov, tov->e,   r);
    return true;
}

// helpers to acces star mass and radius
double gkyl_tov_ultra_rel_star_mass(const struct gkyl_tov_ultra_rel *tov)
{
    if(tov) {
        return tov->M_star;
    }
    else {
        return 0.0;
    }
}

double gkyl_tov_ultra_rel_star_radius(const struct gkyl_tov_ultra_rel *tov)
{
    if(tov) {
        return tov->R_areal;
    }
    else {
        return 0.0;
    }
}

bool gkyl_tov_ultra_rel_surface_found(const struct gkyl_tov_ultra_rel *tov)
{
    if (tov) {
        return tov->surface_found;
    }
    else {
        return false;
    }
}

bool gkyl_tov_ultra_rel_table_exhausted(const struct gkyl_tov_ultra_rel *tov)
{
    if (tov) {
        return tov->table_exhausted;
    }
    else {
        return true;
    }
}

double gkyl_tov_ultra_rel_max_radius(const struct gkyl_tov_ultra_rel *tov)
{
    if (tov) {
        return tov->r_max;
    }
    else {
        return 0.0;
    }
}

double gkyl_tov_ultra_rel_surface_pressure(const struct gkyl_tov_ultra_rel *tov)
{
    if (tov) {
        return tov->P_surface;
    }
    else {
        return 0.0;
    }
}

//main
struct gkyl_tov_ultra_rel * gkyl_tov_ultra_rel_new(double Gamma, double e_c, double dr)
{
    struct gkyl_tov_ultra_rel *tov = (struct gkyl_tov_ultra_rel *) gkyl_malloc(sizeof(struct gkyl_tov_ultra_rel));

    if (!tov) {
        return NULL;
    }

    tov->Gamma = Gamma;
    tov->e_c = e_c;
    tov->dr = dr;

    tov->r_areal = gkyl_malloc(sizeof(double[TOV_ULTRA_REL_MAX_POINTS]));
    tov->m = gkyl_malloc(sizeof(double[TOV_ULTRA_REL_MAX_POINTS]));
    tov->P  = gkyl_malloc(sizeof(double[TOV_ULTRA_REL_MAX_POINTS]));
    tov->e = gkyl_malloc(sizeof(double[TOV_ULTRA_REL_MAX_POINTS]));
    tov->Phi = gkyl_malloc(sizeof(double[TOV_ULTRA_REL_MAX_POINTS]));

    if (!tov->r_areal || !tov->m || !tov->P || !tov->e || !tov->Phi) {
        gkyl_tov_ultra_rel_release(tov);
        return NULL;
    }

    // First we integrate the TOV eqns in the areal coordinates
    double P_c = eos_pressure(Gamma, e_c);
    tov->P_c = P_c;
    // Cutoff-defined surface: the ultra-rel EOS has a long pressure tail,
    // unlike a polytrope with a clean zero-pressure boundary.
    tov->P_surface = eos_pressure(Gamma, TOV_ULTRA_REL_SURFACE_E_FRAC * e_c);
    tov->r_max = (TOV_ULTRA_REL_MAX_POINTS - 1) * dr;
    tov->surface_found = false;
    tov->table_exhausted = false;

    // center
    tov->r_areal[0] = 0.0;
    tov->m[0] = 0.0;
    tov->P[0] = P_c;
    tov->e[0] = e_c;
    tov->Phi[0] = 0.0;  // "provisional"; shifted below

    // Taylor expansion at r = dr  (because the TOV equations are singular-looking at r = 0)
    double r_n = dr;
    double y[3];
    y[0] = (4.0 / 3.0) * M_PI * e_c * r_n * r_n * r_n; //m 
    y[1] = P_c - (2.0 * M_PI / 3.0) * (e_c + P_c) * (e_c / 3.0 + P_c) * r_n * r_n; //P
    y[2] = (2.0 * M_PI / 3.0) * (e_c / 3.0 + P_c) * r_n * r_n; //Phi

    if (y[1] <= 0.0) {
            y[1] = 0.0;
        }

    double e_n = eos_e_from_pressure(Gamma, y[1]);
    tov->r_areal[1] = r_n;
    tov->m[1]       = y[0];
    tov->P[1]       = y[1];
    tov->e[1]     = e_n;
    tov->e[1]       = eos_energy(Gamma, e_n);
    tov->Phi[1]     = y[2];
    int n = 2;     

    // Actual integration call - call to the SSP-RK3 (implicit_tov_source_update)
    double r_max = tov->r_max;
    while (r_n < r_max && y[1] > tov->P_surface && n < TOV_ULTRA_REL_MAX_POINTS) {
        implicit_tov_source_update(r_n, y, dr, Gamma);
        r_n += dr;

        if (y[1] <= tov->P_surface) {
            y[1] = 0.0;
            tov->surface_found = true;
        }

        e_n = eos_e_from_pressure(Gamma, y[1]);

        tov->r_areal[n] = r_n;
        tov->m[n] = y[0];
        tov->P[n]  = y[1];
        tov->e[n] = e_n;
        tov->e[n] = eos_energy(Gamma, e_n);
        tov->Phi[n] = y[2];
        n++;

        if (tov->surface_found) break;
    }

    if (!tov->surface_found) {
        tov->table_exhausted = true;
    }

    tov->N = n;
    int i_surf = n - 1; // stellar surface
    tov->M_star  = tov->m[i_surf];
    tov->R_areal = tov->r_areal[i_surf];

    // We fix lapse: Phi(R) = 0.5 ln(1 - 2M/R). Lapse = e^Phi. At the stellar surface, the interior must match exterior Schwarzschild.
    // g_tt = -(1 - 2M/r) = -e^2Phi. => Phi(R) = 0.5 * ln(1 - 2M/R)
    double Phi_exact = 0.5 * log(1.0 - 2.0 * tov->M_star / tov->R_areal);
    double Phi_shift = Phi_exact - tov->Phi[i_surf]; // We integrated with Phi(0) = Phi_c = 0 before (arbitrary), so we need to shift the profile by a constant
    for (int i = 0; i < n; i++) {
        tov->Phi[i] += Phi_shift;
    }   
    // After the Phi shift -- sanity check that the convention is what we expect
    double alpha_surf_int = exp(tov->Phi[i_surf]);
    double alpha_surf_ext = sqrt(1.0 - 2.0*tov->M_star/tov->R_areal);

    assert(fabs(alpha_surf_int - alpha_surf_ext) < 1e-12 && "Lapse doesn't match on the surface"); // this should pass

    // Extend to vacuum exterior out to r_max
    while (r_n < r_max && n < TOV_ULTRA_REL_MAX_POINTS) {
        r_n += dr;
        double V_ext = 2.0 * tov->M_star / r_n;
        tov->r_areal[n] = r_n;
        tov->m[n]       = tov->M_star;
        tov->P[n]       = 0.0;
        tov->e[n]       = 0.0;

        tov->Phi[n] = (1.0 - V_ext > 0.0)
                  ? 0.5 * log(1.0 - V_ext)
                  : tov->Phi[n-1];
        n++;
    }

    tov->N = n;

    return tov;
}

void
gkyl_gr_tov_ultra_rel_refresh_geometry_from_state(double gas_gamma, double p_atm, int ncells, double *q)
{
  double e_before = 0.0;
  double m_before = 0.0;
  double phi_before = 0.0;

  double r_before = 0.0;
  double dphi_dr_before = 0.0;
  int surface_idx = ncells - 1;

  for (int i = 0; i < ncells; ++i) {
    double *qi = &q[8*i];

    double Etot = qi[1];
    double mom_r = qi[2];
    double r = qi[5];

    double beta = 0.25 * (2.0 - gas_gamma);
    double p_term = (4.0 * beta * beta * Etot * Etot) + ((gas_gamma - 1.0) * ((Etot * Etot) - (mom_r * mom_r)));
    double p = -(2.0 * beta * Etot) + sqrt(fmax(0.0, p_term));
    bool in_atmosphere = false;
    if (p_atm > 0.0 && p <= (1.0 + 1.0e-12) * p_atm) {
        in_atmosphere = true;
    }
    double p_geom = 0.0;
    if (p_atm > 0.0) {
        p_geom = fmax(0.0, p - p_atm);
    }
    else {
        p_geom = p;
    }

    double e_geom = p_geom / (gas_gamma - 1.0);

    double vel = 0.0;
    double denom = Etot + p;
    if (denom > 0.0)
      vel = mom_r / denom;

    double W = 1.0 / sqrt(fmax(1.0e-8, 1.0 - vel * vel));

    if (in_atmosphere && surface_idx == ncells - 1) {
      surface_idx = i > 0 ? i - 1 : 0;
    }

    // double m = qi[4];
    // for refreshing m start:
    double m = 0.0;
    if (!in_atmosphere) {
      double dm_dr = 4.0 * M_PI * r * r * e_geom;
      if (i == 0) {
        m = (4.0 / 3.0) * M_PI * e_geom * r * r * r;
      }
      else {
        double dm_dr_before = 4.0 * M_PI * r_before * r_before * e_before;
        double dr = r - r_before;
        m = m_before + 0.5 * dr * (dm_dr_before + dm_dr);
      }
    }
    else {
      m = m_before;
    }
    qi[4] = m;
    // for refreshing m stop

    double dphi_dr = 0.0;
    if (!in_atmosphere && r > 0.0 && r - 2.0 * m > 0.0) {
      double radial_stress = mom_r * vel + p_geom; // fluid is moving, no longer just pressure
      dphi_dr = (m + 4.0 * M_PI * r * r * r * radial_stress) / (r * (r - 2.0 * m));
    }

    double phi = 0.0;
    if (i == 0) {
      phi = 0.0; // provisional; shift after sweep
    }
    else {
      double dr = r - r_before;
      phi = phi_before + 0.5*dr*(dphi_dr_before + dphi_dr);
    }

    qi[3] = phi;

    r_before = r;
    e_before = e_geom;
    m_before = m;
    phi_before = phi;
    dphi_dr_before = dphi_dr;
  }

  double *qsurf = &q[8*surface_idx];
  double R = qsurf[5];
  double M = qsurf[4];

  if (R > 0.0 && 1.0 - 2.0*M/R > 0.0) {
    double phi_match = 0.5 * log(1.0 - 2.0 * M / R);
    double phi_shift = phi_match - qsurf[3];

    for (int i = 0; i <= surface_idx; ++i) {
      q[8*i + 3] += phi_shift;
    }
  }

  for (int i = surface_idx + 1; i < ncells; ++i) {
    double *qi = &q[8*i];
    double r = qi[5];
    double m = qi[4];
    if (r > 0.0 && 1.0 - 2.0 * m / r > 0.0) {
      qi[3] = 0.5 * log(1.0 - 2.0 * m / r);
    }
  }
}

void gkyl_tov_ultra_rel_release(struct gkyl_tov_ultra_rel *tov)
{
    if (!tov) return;
    gkyl_free(tov->r_areal);
    gkyl_free(tov->m);
    gkyl_free(tov->P);
    gkyl_free(tov->e);
    gkyl_free(tov->Phi);
    gkyl_free(tov);
}
