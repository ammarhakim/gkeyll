// TCV-like inner-wall limited plasma, 3x2v gyrokinetic simulation with adaptive
// sources, using the extended analytical magnetic equilibrium of
//   F. Riva et al., Plasma Phys. Control. Fusion 59, 035001 (2017),
// section 2.3 and appendix B (a solution of the Grad-Shafranov equation in the
// r/R0 -> 0 limit). Unlike the simplified Miller geometry used in
// rt_gk_tcv_iwl_adapt_source_3x2v_p1.c, here the elongation (kappa) and
// triangularity (delta) imposed at the LCFS vary radially through the shaping
// functions S2(r) and S3(r) (eq. (B9)), with a self-consistent Shafranov shift
// Delta(r) (eq. (B10)) and a paramagnetic correction to the toroidal field
// function F(r) (eq. (27)). The equilibrium assumes a parabolic safety factor
// profile q(r) = qaxis + (qLCFS-qaxis)*(r/a)^2 (eq. (29)).
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Number of nodes of the radial tables of the Shafranov shift Delta(r) and of
// dPsidr(r), built once in create_ctx. They avoid nested numerical quadratures
// (Delta appears inside the integrands of dPsidr and alpha) which make the
// geometry initialization orders of magnitude slower.
#define GK_DELTA_TAB_N 257
#define GK_PSI_TAB_N 65

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
  int cdim, vdim;
  // Geometry and magnetic field parameters
  double Z_axis, R_axis, R_geo, R0, a_mid, x_inner, r0, B_axis, B0, kappa, delta,
         qaxis, qLCFS, S2a, S3a, q0, Cy, Bref, x_LCFS;
  // Tabulated Shafranov shift Delta(r) on r in [0, delta_rmax].
  bool delta_tab_set;
  int delta_nr;
  double delta_rmax;
  double delta_tab[GK_DELTA_TAB_N];
  // Tabulated dPsidr(r) on r in [psi_rmin, psi_rmax].
  bool psi_tab_set;
  int psi_nr;
  double psi_rmin, psi_rmax;
  double psi_tab[GK_PSI_TAB_N];
  // Plasma parameters
  int num_species;
  double me, qe, mi, qi, n0, Te0, Ti0;
  // Collision parameters
  double nuFrac, nuElc, nuIon;
  // Source parameters
  int num_sources;
  bool adapt_energy_srcCORE, adapt_particle_srcCORE;
  double center_srcCORE[3], sigma_srcCORE[3];
  double energy_srcCORE, particle_srcCORE;
  double floor_srcCORE;
  bool adapt_energy_srcRECY, adapt_particle_srcRECY;
  double center_srcRECY[3], sigma_srcRECY[3];
  double energy_srcRECY, particle_srcRECY;
  double floor_srcRECY;
  // Grid parameters
  double Lx, Ly, Lz;
  double x_min, x_max, y_min, y_max, z_min, z_max;
  int num_cell_x, num_cell_y, num_cell_z, num_cell_vpar, num_cell_mu;
  int cells[GKYL_MAX_DIM], poly_order;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  // Simulation control parameters
  double final_time, write_phase_freq;
  int num_frames, int_diag_calc_num, num_failures_max;
  double dt_failure_tol;
};

// Geometry related functions
double r_x(double x, double a_mid, double x_inner)
{
  return x+a_mid-x_inner;
}

// Parabolic safety factor profile (Riva et al. eq. (29)).
double qprofile(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double qaxis = app->qaxis;
  double qLCFS = app->qLCFS;
  double a_mid = app->a_mid;
  return qaxis + (qLCFS - qaxis)*pow(r/a_mid,2);
}

// Normalized internal inductance for the parabolic q profile (Riva et al. eq. (B9)).
double l_i(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double qaxis = app->qaxis;
  double q = qprofile(r, ctx);
  double dq = (q - qaxis)/qaxis;
  if (fabs(dq) < 1.e-4)
    return 0.5 + dq/3.; // Series expansion, regular at r=0.
  return pow(q,2)*(log(q/qaxis) + qaxis/q - 1.)/pow(q - qaxis,2);
}

struct integrand_ctx {
  struct gk_app_ctx *app_ctx;
  double r;
};

double shift_integrand(double t, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  struct gk_app_ctx *app = inctx->app_ctx;
  return t*l_i(t, app);
}

double dDeltadr(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double R_geo = app->R_geo;
  return -r*l_i(r,ctx)/(2.*R_geo);
}

// Shafranov shift, with Delta(a)=0 (Riva et al. eq. (B10)).
double Delta(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double R_geo = app->R_geo;

  if (app->delta_tab_set && 0. <= r && r <= app->delta_rmax) {
    // Cubic Hermite interpolation of the tabulated shift, using the exact
    // derivative dDeltadr at the nodes (C^1, O(dr^4) accurate).
    double dr = app->delta_rmax/(app->delta_nr-1);
    int i = (int) (r/dr);
    if (i > app->delta_nr-2) i = app->delta_nr-2;
    double t = (r - i*dr)/dr;
    double y0 = app->delta_tab[i], y1 = app->delta_tab[i+1];
    double m0 = dDeltadr(i*dr, ctx), m1 = dDeltadr((i+1)*dr, ctx);
    double t2 = t*t, t3 = t2*t;
    return (2.*t3-3.*t2+1.)*y0 + (t3-2.*t2+t)*dr*m0
         + (-2.*t3+3.*t2)*y1 + (t3-t2)*dr*m1;
  }

  if (fabs(r - a_mid) < 1.e-14)
    return 0.;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  if (r < a_mid) {
    integral = gkyl_dbl_exp(shift_integrand, &tmp_ctx, r, a_mid, 7, 1e-10);
  } else {
    integral = gkyl_dbl_exp(shift_integrand, &tmp_ctx, a_mid, r, 7, 1e-10);
    integral.res = -integral.res;
  }
  return integral.res/(2.*R_geo);
}

// Shaping functions, S2 (elongation, m=2) and S3 (triangularity, m=3)
// (Riva et al. eq. (B9)), with S2(a) = a*(1-kappa)/(1+kappa) and S3(a) = delta*a/4.
double Sm(double r, int m, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double qaxis = app->qaxis, qLCFS = app->qLCFS, a_mid = app->a_mid;
  double Sm_a = m == 2 ? app->S2a : app->S3a;
  double cm = qaxis*(m+1.)/(m-1.);
  double q = qprofile(r, ctx);
  return Sm_a*pow(r/a_mid,m-1)*((q - qaxis) + cm)/((qLCFS - qaxis) + cm);
}

double dSmdr(double r, int m, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double qaxis = app->qaxis, qLCFS = app->qLCFS, a_mid = app->a_mid;
  double Sm_a = m == 2 ? app->S2a : app->S3a;
  double cm = qaxis*(m+1.)/(m-1.);
  double q = qprofile(r, ctx);
  double dqdr = 2.*(qLCFS - qaxis)*r/pow(a_mid,2);
  return Sm_a/((qLCFS - qaxis) + cm)*( (m-1.)*pow(r/a_mid,m-2)*((q - qaxis) + cm)/a_mid
                                      +pow(r/a_mid,m-1)*dqdr );
}

// Correction factor P(r), ensuring psi'(r) = (r/R0)*F(psi)/q(r) + O((r/R0)^4)
// (Riva et al. eq. (B3)).
double Pfunc(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double R_geo = app->R_geo;
  double S2 = Sm(r,2,ctx), S3 = Sm(r,3,ctx);
  return pow(r,3)/(8.*pow(R_geo,2)) - r*Delta(r,ctx)/(2.*R_geo)
    - pow(S2,2)/(2.*r) - pow(S3,2)/r;
}

double dPfuncdr(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double R_geo = app->R_geo;
  double S2 = Sm(r,2,ctx), S3 = Sm(r,3,ctx);
  double dS2 = dSmdr(r,2,ctx), dS3 = dSmdr(r,3,ctx);
  return 3.*pow(r,2)/(8.*pow(R_geo,2)) - Delta(r,ctx)/(2.*R_geo) - r*dDeltadr(r,ctx)/(2.*R_geo)
    - S2*dS2/r + pow(S2,2)/(2.*pow(r,2)) - 2.*S3*dS3/r + pow(S3,2)/pow(r,2);
}

// Flux surfaces (Riva et al. eqs. (B1)-(B2)).
double R_rtheta(double r, double theta, void *ctx)
{
  // Major radius as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double R_geo = app->R_geo;
  return R_geo + Delta(r,ctx) + (r + Sm(r,2,ctx) - Pfunc(r,ctx))*cos(theta)
    + Sm(r,3,ctx)*cos(2.*theta);
}

double Z_rtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Z_axis = app->Z_axis;
  return Z_axis + (r - Sm(r,2,ctx) - Pfunc(r,ctx))*sin(theta)
    - Sm(r,3,ctx)*sin(2.*theta);
}

double dRdr(double r, double theta, void *ctx)
{
  return dDeltadr(r,ctx) + (1. + dSmdr(r,2,ctx) - dPfuncdr(r,ctx))*cos(theta)
    + dSmdr(r,3,ctx)*cos(2.*theta);
}

double dRdtheta(double r, double theta, void *ctx)
{
  return -(r + Sm(r,2,ctx) - Pfunc(r,ctx))*sin(theta)
    - 2.*Sm(r,3,ctx)*sin(2.*theta);
}

double dZdr(double r, double theta, void *ctx)
{
  return (1. - dSmdr(r,2,ctx) - dPfuncdr(r,ctx))*sin(theta)
    - dSmdr(r,3,ctx)*sin(2.*theta);
}

double dZdtheta(double r, double theta, void *ctx)
{
  return (r - Sm(r,2,ctx) - Pfunc(r,ctx))*cos(theta)
    - 2.*Sm(r,3,ctx)*cos(2.*theta);
}

double Jr(double r, double theta, void *ctx)
{
  return R_rtheta(r,theta,ctx)*( dRdr(r,theta,ctx) * dZdtheta(r,theta,ctx)
		                -dRdtheta(r,theta,ctx) * dZdr(r,theta,ctx) );
}

double integrand(double t, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  double r = inctx->r;
  struct gk_app_ctx *app = inctx->app_ctx;
  return Jr(r,t,app) / pow(R_rtheta(r,t,app),2);
}

// Toroidal field function F(r) = B_axis*R_axis*(1 + F2(r)), with the
// paramagnetic correction F2(r) (Riva et al. eq. (27)).
double Ffunc(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double R_geo = app->R_geo, a_mid = app->a_mid;
  double qaxis = app->qaxis, qLCFS = app->qLCFS;
  double B_axis = app->B_axis, R_axis = app->R_axis;
  double q = qprofile(r, ctx);
  double F2 = -pow(r,2)/(2.*pow(R_geo,2)*pow(q,2))*(2. + (qLCFS - qaxis)/qaxis*pow(r/a_mid,2));
  return B_axis*R_axis*(1. + F2);
}

double Bphi(double r, double R, void *ctx)
{
  // Toroidal magnetic field.
  return Ffunc(r,ctx)/R;
}

double dPsidr(double r, void *ctx)
{
  // Radial derivative of the poloidal flux function (Riva et al. eq. (28)).
  struct gk_app_ctx *app = ctx;

  if (app->psi_tab_set && app->psi_rmin <= r && r <= app->psi_rmax) {
    // Cubic Hermite interpolation with Catmull-Rom slopes.
    int n = app->psi_nr;
    double dr = (app->psi_rmax - app->psi_rmin)/(n-1);
    int i = (int) ((r - app->psi_rmin)/dr);
    if (i > n-2) i = n-2;
    double t = (r - (app->psi_rmin + i*dr))/dr;
    const double *y = app->psi_tab;
    double m0 = i == 0 ? (y[1]-y[0])/dr : (y[i+1]-y[i-1])/(2.*dr);
    double m1 = i == n-2 ? (y[n-1]-y[n-2])/dr : (y[i+2]-y[i])/(2.*dr);
    double t2 = t*t, t3 = t2*t;
    return (2.*t3-3.*t2+1.)*y[i] + (t3-2.*t2+t)*dr*m0
         + (-2.*t3+3.*t2)*y[i+1] + (t3-t2)*dr*m1;
  }

  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);

  return ( Ffunc(r,ctx)/(2.*M_PI*qprofile(r,ctx)) )*integral.res;
}

double alpha(double r, double theta, double phi, void *ctx)
{
  double twrap = theta;
  while (twrap < -M_PI) twrap = twrap+2.*M_PI;
  while (M_PI < twrap) twrap = twrap-2.*M_PI;

  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  if (0. < twrap) {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., twrap, 7, 1e-10);
  } else {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, twrap, 0., 7, 1e-10);
    integral.res = -integral.res;
  }

  return phi - Ffunc(r,ctx)*integral.res/dPsidr(r,ctx);
}

double gradr(double r, double theta, void *ctx)
{
  return (R_rtheta(r,theta,ctx)/Jr(r,theta,ctx))*sqrt(pow(dRdtheta(r,theta,ctx),2) + pow(dZdtheta(r,theta,ctx),2));
}

void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Density initial condition (like TCV exp profile)
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double n0 = 5e19;
  double x0 = -0.03;
  double c1 = 0.5;
  double c2 = 8.0;
  double c3 = 0.005;
  fout[0] = n0*(c1*(1.+tanh(c2*(-10*(x+x0))))+c3);
}

// Electron temperature initial conditions
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double T0 = 200 * GKYL_ELEMENTARY_CHARGE;
  double x0 = -0.03; // position of the transition region
  double c0 = 1.3; // multiplicative factor
  double c1 = 0.5; // control the temperature at the _core
  double c2 = 8.0; // control the width of the transition region
  double c3 = 0.1; // control the temperature at the SOL
  fout[0] = c0*T0*(c1*(1.+tanh(c2*(-10*(x+x0))))+c3);
}

// Ion temperature initial conditions
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double T0 = 200 * GKYL_ELEMENTARY_CHARGE;
  double x0 = -0.04; // position of the transition region
  double c0 = 1.0; // multiplicative factor
  double c1 = 0.5; // control the temperature at the _core
  double c2 = 3.0; // control the width of the transition region
  double c3 = 0.2; // control the temperature at the SOL
  fout[0] = c0*T0*(c1*(1.+tanh(c2*(-10*(x+x0))))+c3);
}

// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];

  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double Cy = app->Cy;

  double r = r_x(x,a_mid,x_inner);

  // Map to cylindrical (R, Z, phi) coordinates.
  double R = R_rtheta(r, z, ctx);
  double Z = Z_rtheta(r, z, ctx);
  double phi = y/Cy + alpha(r, z, 0, ctx);

  // Map to Cartesian (X, Y, Z) coordinates.
  double X = R*cos(phi);
  double Y = R*sin(phi);
  xp[0] = X; xp[1] = Y; xp[2] = Z;
}

// Taken from rt gk d3d 3x2c, is this the non uniform v grid mapping?
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double Cy = app->Cy;
  double x_inner = app->x_inner;
  double r = r_x(x,a_mid,x_inner);
  double R = R_rtheta(r, z, ctx);
  double Bt = Bphi(r, R, ctx);
  double Bp = dPsidr(r,ctx)/R*gradr(r,z,ctx);

  double drdtheta = dRdtheta(r,z,ctx);
  double dzdtheta = dZdtheta(r,z,ctx);
  double den = sqrt(pow(drdtheta,2) + pow(dzdtheta,2));
  double B_r = Bp*drdtheta/den;
  double B_z = Bp*dzdtheta/den;
  double phi = y/Cy + alpha(r, z, 0, ctx);

  // xc are computational coords.
  // Set Cartesian components of magnetic field.
  fout[0] = B_r * cos(phi) + Bt * sin(phi);
  fout[1] = B_r * sin(phi) - Bt * cos(phi);
  fout[2] = B_z;
}

void bc_shift_func_lo(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;

  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double Cy = app->Cy;
  double z_min = app->z_min;
  double z_max = app->z_max;

  double r = r_x(x, a_mid, x_inner);

  fout[0] = Cy*( alpha(r, z_min, 0.0, ctx) - alpha(r, z_max, 0.0, ctx) );
}

void bc_shift_func_up(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;

  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double Cy = app->Cy;
  double z_min = app->z_min;
  double z_max = app->z_max;

  double r = r_x(x, a_mid, x_inner);

  fout[0] = -Cy*( alpha(r, z_min, 0.0, ctx) - alpha(r, z_max, 0.0, ctx) );
}

struct gk_app_ctx create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double Z_axis = 0.1414361745; // Magnetic axis height [m].
  double R_axis = 0.8867856264; // Magnetic axis major radius [m].
  double B_axis = 1.4; // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 1.0870056099999; // Major radius of the LCFS at the outboard midplane [m].
  double x_inner = 0.04; // Radial extent inside LCFS
  double x_outer = 0.08; // Radial extent outside LCFS
  double Rmid_min = R_LCFSmid - x_inner; // Minimum midplane major radius of simulation box [m].
  double Rmid_max = R_LCFSmid + x_outer; // Maximum midplane major radius of simulation box [m].
  double R0 = 0.5*(Rmid_min+Rmid_max); // Major radius of the simulation box [m].

  double kappa = 1.4; // Elongation of the LCFS (=1 for no elongation).
  double delta = -0.38; // Triangularity of the LCFS (=0 for no triangularity).
  // Parabolic safety factor profile parameters (Riva et al. eq. (29)).
  // Note: the simulation domain only covers r/a_mid in ~[0.86, 1.27], so qaxis acts
  // as a profile-shape parameter. It is set above the experimental value because for
  // lower qaxis the radial growth of the shaping functions S2(r), S3(r) (proportional
  // to q(r)-qaxis, eq. (B9)) makes the flux surfaces cross (negative Jacobian) in the
  // far SOL, where the equilibrium (only valid up to r=a_mid) is extrapolated.
  double qaxis = 1.6; // Safety factor at the magnetic axis.
  double qLCFS = 2.49; // Safety factor at the LCFS (TCV NT #65130 q fit value at R_LCFSmid).

  // Solve for R_geo (major radius of the LCFS center, Riva et al.'s R0) and
  // a_mid (LCFS minor radius) such that the magnetic axis sits at R_axis,
  // i.e. R_geo + Delta(0) = R_axis, and the LCFS crosses the outboard midplane
  // at R_LCFSmid, i.e. R_rtheta(a_mid, 0) = R_LCFSmid.
  struct gk_app_ctx geo_ctx = {
    .R_axis = R_axis, .B_axis = B_axis,
    .qaxis = qaxis, .qLCFS = qLCFS,
  };
  double R_geo = R_axis;
  double a_mid = R_LCFSmid - R_axis;
  for (int i=0; i<10; i++) {
    geo_ctx.R_geo = R_geo;
    geo_ctx.a_mid = a_mid;
    geo_ctx.S2a = a_mid*(1.-kappa)/(1.+kappa);
    geo_ctx.S3a = delta*a_mid/4.;
    R_geo = R_axis - Delta(0., &geo_ctx);
    a_mid = R_LCFSmid - R_geo - geo_ctx.S2a - geo_ctx.S3a + Pfunc(a_mid, &geo_ctx);
  }
  double S2a = geo_ctx.S2a;
  double S3a = geo_ctx.S3a;

  // Tabulate the Shafranov shift Delta(r) and dPsidr(r) once (by quadrature,
  // since the tables are not set yet). All subsequent geometry evaluations
  // interpolate these tables, avoiding nested quadratures.
  geo_ctx.delta_nr = GK_DELTA_TAB_N;
  geo_ctx.delta_rmax = a_mid + 3.*x_outer;
  for (int i=0; i<geo_ctx.delta_nr; i++)
    geo_ctx.delta_tab[i] = Delta(i*geo_ctx.delta_rmax/(geo_ctx.delta_nr-1), &geo_ctx);
  geo_ctx.delta_tab_set = true;

  geo_ctx.psi_nr = GK_PSI_TAB_N;
  geo_ctx.psi_rmin = a_mid - 2.*x_inner; // Covers the radial domain plus ghost cells.
  geo_ctx.psi_rmax = a_mid + 2.*x_outer;
  for (int i=0; i<geo_ctx.psi_nr; i++)
    geo_ctx.psi_tab[i] = dPsidr(geo_ctx.psi_rmin
      + i*(geo_ctx.psi_rmax - geo_ctx.psi_rmin)/(geo_ctx.psi_nr-1), &geo_ctx);
  geo_ctx.psi_tab_set = true;

  double r0 = a_mid + 0.5*(x_outer - x_inner); // Minor radius at the center of the simulation box [m].
  double B0 = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  int num_species = 2;
  double AMU = 2.01410177811;
  double mi = mp*AMU; // Deuterium ions.
  double Te0 = 100*eV; // Ion reference temperature [J].
  double Ti0 = 100*eV; // Electron reference temperature [J].
  double n0  = 2.0e19; // Reference density [1/m^3].
  double Bref = 1.129; // Reference magnetic field [T].
  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi);
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Configuration domain parameters
  double Lx = Rmid_max-Rmid_min; // Domain size along x.
  double x_min = 0.;
  double x_max = Lx;
  double x_LCFS = R_LCFSmid - Rmid_min; // Radial location of the last closed flux surface.
  double q0 = qprofile(r0, &geo_ctx); // Safety factor in the center of domain.
  double Cy = r0/q0; // Normalization in binormal coordinate.

  double Ly = 150*rho_s; // Domain size along y.
  // Adjust the domain size along y to have integer toroidal mode number.
  // We need: 2*pi*Cy/Ly = integer.
  Ly = 2.*M_PI*Cy/round(2.*M_PI*Cy/Ly);
  double y_min = -Ly/2.;
  double y_max =  Ly/2.;

  double Lz = 2.*M_PI-1e-10; // Domain size along magnetic field.
  double z_min = -Lz/2.;
  double z_max =  Lz/2.;

  // Collision frequencies
  double nuFrac = 0.1;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4) * n0 /
    (6*sqrt(2.) * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(me) * pow(Te0,3./2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));

  // Source parameters
  int num_sources = 2;
  double P_exp = 0.34e6; // P_sol measured [W]
  double vol_frac = 1.0/(2.*M_PI*Cy/Ly); // Volume fraction of the simulation box (=0.5 here).
  double P_inj = P_exp * vol_frac / num_species; // Injection power normalized to the volume fraction and per species [W]
  // Core source:
  // - Injects energy only in the core region (0.25MW per species).
  // - The particles injection is only the one that are lost through the inner radial boundary.
  bool adapt_energy_srcCORE = true; // The source will compensate the losses in energy according to given boundaries.
  bool adapt_particle_srcCORE = true; // The source will compensate the losses in particle according to given boundaries.
  double energy_srcCORE = P_inj; // What the source must inject in energy [W]
  double particle_srcCORE = 0.0; // What the source must inject in particle [1/s]
  double center_srcCORE[3] = {x_min, 0.0, -Lz/4}; // This is the position of the ion source,
  double sigma_srcCORE[3] = {0.03*Lx, 0.0, Lz/6}; // the electron source will be at +Lz/2.
  double floor_srcCORE = 1e-10;
  // Recycling source:
  // - Reinjects particles that are absorbed by the wall.
  // - Energy is free to leave the system.
  bool adapt_energy_srcRECY = false;
  bool adapt_particle_srcRECY = true;
  double energy_srcRECY = 0.0; // [W]
  double particle_srcRECY = 0.0; // [1/s]
  double center_srcRECY[3] = {0.5*x_LCFS, 0.0, M_PI};
  double sigma_srcRECY[3] = {0.25*x_LCFS, 0.0, 0.05*Lz};
  double floor_srcRECY = 1e-10;

  // Grid parameters (reduced resolution for the regression test, minimal recommended values in comments)
  int num_cell_x = 9; // (24) The LCFS is positionned at 1/3 of the domain -> the resolution must be divisible by 3.
  int num_cell_y = 4; // (16)
  int num_cell_z = 8; // (12)
  int num_cell_vpar = 8; // (12)
  int num_cell_mu = 8; // (8)
  int poly_order = 1;
  // Velocity box dimensions
  double vpar_max_elc = 5.*vte;
  double mu_max_elc   = 1.*me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 5.*vti;
  double mu_max_ion   = 1.*mi*pow(4*vti,2)/(2*B0);
  double final_time = 1.e-7;
  int num_frames = 1;
  double write_phase_freq = 1.0;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .Z_axis = Z_axis,
    .R_axis = R_axis,
    .R_geo = R_geo,
    .R0 = R0,
    .a_mid = a_mid,
    .x_inner = x_inner,
    .r0 = r0,
    .B_axis = B_axis,
    .B0 = B0,
    .kappa = kappa,
    .delta = delta,
    .qaxis = qaxis,
    .qLCFS = qLCFS,
    .S2a = S2a,
    .S3a = S3a,
    .q0 = q0,
    .Cy = Cy,
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .x_min = x_min, .x_max = x_max,
    .y_min = y_min, .y_max = y_max,
    .z_min = z_min, .z_max = z_max,
    .Bref = Bref,
    .x_LCFS = x_LCFS,
    .num_species = num_species,
    .me = me, .qe = qe,
    .mi = mi, .qi = qi,
    .n0 = n0, .Te0 = Te0, .Ti0 = Ti0,
    .nuFrac = nuFrac, .nuElc = nuElc, .nuIon = nuIon,
    .num_sources = num_sources,
    .adapt_energy_srcCORE = adapt_energy_srcCORE,
    .adapt_particle_srcCORE = adapt_particle_srcCORE,
    .center_srcCORE = {center_srcCORE[0], center_srcCORE[1], center_srcCORE[2]},
    .sigma_srcCORE = {sigma_srcCORE[0], sigma_srcCORE[1], sigma_srcCORE[2]},
    .energy_srcCORE = energy_srcCORE, .particle_srcCORE = particle_srcCORE,
    .floor_srcCORE = floor_srcCORE,
    .adapt_energy_srcRECY = adapt_energy_srcRECY,
    .adapt_particle_srcRECY = adapt_particle_srcRECY,
    .center_srcRECY = {center_srcRECY[0], center_srcRECY[1], center_srcRECY[2]},
    .sigma_srcRECY = {sigma_srcRECY[0], sigma_srcRECY[1], sigma_srcRECY[2]},
    .energy_srcRECY = energy_srcRECY,  .particle_srcRECY = particle_srcRECY,
    .floor_srcRECY = floor_srcRECY,
    .num_cell_x = num_cell_x,
    .num_cell_y = num_cell_y,
    .num_cell_z = num_cell_z,
    .num_cell_vpar = num_cell_vpar,
    .num_cell_mu = num_cell_mu,
    .cells = {num_cell_x, num_cell_y, num_cell_z, num_cell_vpar, num_cell_mu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc, .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion, .mu_max_ion = mu_max_ion,
    .write_phase_freq = write_phase_freq,
    .final_time = final_time, .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  // Copy the geometry tables built above.
  ctx.delta_nr = geo_ctx.delta_nr;
  ctx.delta_rmax = geo_ctx.delta_rmax;
  for (int i=0; i<geo_ctx.delta_nr; i++)
    ctx.delta_tab[i] = geo_ctx.delta_tab[i];
  ctx.delta_tab_set = true;
  ctx.psi_nr = geo_ctx.psi_nr;
  ctx.psi_rmin = geo_ctx.psi_rmin;
  ctx.psi_rmax = geo_ctx.psi_rmax;
  for (int i=0; i<geo_ctx.psi_nr; i++)
    ctx.psi_tab[i] = geo_ctx.psi_tab[i];
  ctx.psi_tab_set = true;

  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Definition of the sources projections.
  struct gkyl_gyrokinetic_projection proj_srcCORE_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcCORE[0], ctx.center_srcCORE[1], -ctx.center_srcCORE[2]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1], ctx.sigma_srcCORE[2]},
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0*ctx.Te0,
    .temp_min = 0.1*ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_projection proj_srcCORE_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcCORE[0], ctx.center_srcCORE[1], ctx.center_srcCORE[2]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1], ctx.sigma_srcCORE[2]},
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0*ctx.Ti0,
    .temp_min = 0.1*ctx.Ti0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_projection proj_srcRECY_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcRECY[0], ctx.center_srcRECY[1], ctx.center_srcRECY[2]},
    .gaussian_std_dev = {ctx.sigma_srcRECY[0], ctx.sigma_srcRECY[1], ctx.sigma_srcRECY[2]},
    .total_num_particles = ctx.particle_srcRECY,
    .total_kin_energy = ctx.energy_srcRECY,
    .temp_max = 5.0*ctx.Te0,
    .temp_min = 0.1*ctx.Te0,
    .f_floor = ctx.floor_srcRECY,
  };
  struct gkyl_gyrokinetic_projection proj_srcRECY_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcRECY[0], ctx.center_srcRECY[1], ctx.center_srcRECY[2]},
    .gaussian_std_dev = {ctx.sigma_srcRECY[0], ctx.sigma_srcRECY[1], ctx.sigma_srcRECY[2]},
    .total_num_particles = ctx.particle_srcRECY,
    .total_kin_energy = ctx.energy_srcRECY,
    .temp_max = 5.0*ctx.Ti0,
    .temp_min = 0.1*ctx.Ti0,
    .f_floor = ctx.floor_srcRECY,
  };

  // Definition of the source adaptation structures.
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_e ={
      .adapt_to_species = "elc",
      .adapt_particle = ctx.adapt_particle_srcCORE,
      .adapt_energy = ctx.adapt_energy_srcCORE,
      .num_boundaries = 1, // Only the inner radial boundary.
      .dir = {0},
      .edge = {GKYL_LOWER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_i ={
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1, // Only the inner radial boundary.
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECY_e = {
    .adapt_to_species = "ion", // Adapted to ions to maintain neutrality.
    .adapt_particle = ctx.adapt_particle_srcRECY,
    .adapt_energy = ctx.adapt_energy_srcRECY,
    .num_boundaries = 3, // outer radial boundary + both z boundaries.
    .dir = {0, 2, 2},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECY_i = {
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcRECY,
    .adapt_energy = ctx.adapt_energy_srcRECY,
    .num_boundaries = 3, // outer radial boundary + both z boundaries.
    .dir = {0, 2, 2},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };

  // Species definitions
  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = { 1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    // .init_from_file = {
    //    .type = GKYL_IC_IMPORT_F,
    //    .file_name = "restart-elc.gkyl",
    //    .jacobtot_inv_file_name = "restart-jacobtot_inv.gkyl",
    // },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = zero_func,
      .temp = temp_elc,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "ion"},
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_e,
      .adapt[0] = adapt_srcCORE_e,
      .projection[1] = proj_srcRECY_e,
      .adapt[1] = adapt_srcRECY_e,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN,
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP,
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },

    .time_rate_diagnostics = true,
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = { 1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    // .init_from_file = {
    //    .type = GKYL_IC_IMPORT_F,
    //    .file_name = "restart-ion.gkyl",
    //    .jacobtot_inv_file_name = "restart-jacobtot_inv.gkyl",
    // },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = zero_func,
      .temp = temp_ion,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Ti0, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "elc"},
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_i,
      .adapt[0] = adapt_srcCORE_i,
      .projection[1] = proj_srcRECY_i,
      .adapt[1] = adapt_srcRECY_i,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN,
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP,
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },

    .time_rate_diagnostics = true,
  };

  struct gkyl_poisson_bias_line target_corner_bcs[] = {
    {
     .perp_dirs = {0, 2}, // Directions perpendicular to line.
     .perp_coords = {ctx.x_LCFS, ctx.z_min}, // Coordinates of the line in perpendicular directions.
     .val = 0.0, // Biasing value.
    },
    {
     .perp_dirs = {0, 2}, // Directions perpendicular to line.
     .perp_coords = {ctx.x_LCFS, ctx.z_max}, // Coordinates of the line in perpendicular directions.
     .val = 0.0, // Biasing value.
    },
  };

  struct gkyl_poisson_bias_line_list bias_line_list = {
    .num_bias_line = 2,
    .bl = target_corner_bcs,
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .polarization_bmag = ctx.Bref,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
    .bias_line_list = &bias_line_list,
    .time_rate_diagnostics = true,
  };

  // Geometry
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.},
    .mapc2p = mapc2p, // mapping of cCOREutational to physical space
    .c2p_ctx = &ctx,
    .bfield_func = bfield_func, // magnetic field
    .bfield_ctx = &ctx,
    .has_LCFS = true,
    .x_LCFS = ctx.x_LCFS,
    .parallel_lower_bc_shift_func = bc_shift_func_lo,
    .parallel_upper_bc_shift_func = bc_shift_func_up,
    .parallel_lower_bc_shift_ctx = &ctx,
    .parallel_upper_bc_shift_ctx = &ctx,
  };

  // Parallelism
  struct gkyl_app_parallelism_inp parallelism = {
    .comm = comm,
    .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
    .use_gpu = app_args.use_gpu,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cfl_frac_omegaH = 1.0e9,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim,
    .lower = { ctx.x_min, ctx.y_min, ctx.z_min },
    .upper = { ctx.x_max, ctx.y_max, ctx.z_max },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = geometry,

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = ctx.num_species,
    .species = { elc, ion },

    .field = field,

    .parallelism = parallelism
  };

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping = {
      .t_end = ctx.final_time,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
