// Stellarator-like "rotating ellipse" (stellipse) plasma, 3x2v gyrokinetic
// simulation of a toroidally limited plasma with adaptive sources.
//
// Model magnetic geometry:
// - Nested toroidal surfaces whose elliptical cross section (elongation
//   kappa_ell, semi-axes r*sqrt(kappa_ell) and r/sqrt(kappa_ell)) rotates
//   poloidally by NFP*phi/2 around the torus: the classical l=2 stellarator
//   (see e.g. P. Helander, Rep. Prog. Fusion Phys. 77, 087001 (2014);
//   J. Loizu et al., J. Plasma Phys. 83, 715830601 (2017)).
// - The magnetic field is defined in Clebsch form,
//     B = psi'(r) [ grad(r) x grad(theta) - iota(r) grad(r) x grad(phi) ],
//   with psi'(r) = B0*r. Since the Jacobian of the (r,theta,phi) map is
//   J = r*R exactly (the cross section is area preserving), this reduces to
//     B = (B0/R) * ( iota(r) e_theta + e_phi ),
//   where e_theta, e_phi are the tangent vectors of the map. This field is
//   divergence-free and tangent to the surfaces by construction, and theta
//   is a straight-field-line angle: field lines obey theta = alpha + iota*phi.
// - The rotational transform is set to the on-axis vacuum rotating-ellipse
//   value, iota = (NFP/2)*(a-b)^2/(a^2+b^2) with a,b the ellipse semi-axes
//   (nearly shearless, as the vacuum rotating ellipse).
// This is a self-consistent model field rather than an exact MHD equilibrium:
// |B| carries the helical ripple of the rotating ellipse, but the physical
// toroidal field component is B0 (uniform) instead of ~1/R, so the grad-B
// drive is weaker than in a vacuum field while the curvature drive remains.
//
// Field-aligned computational coordinates (x,y,z):
//   r = x + a_mid - x_inner (radial),  phi = z (toroidal angle),
//   theta = -y/Cy + iota(r)*z,         Cy = r0.
// The tangent d/dz at fixed (x,y) is iota*e_theta + e_phi, parallel to B by
// construction, as required by the gyrokinetic app. The minus sign in the
// y-label makes (x,y,z) right-handed (positive computational Jacobian), as in
// the tokamak IWL tests.
//
// Limited scenario: a toroidal ring limiter at phi = +/- pi spans the SOL
// (x > x_LCFS): sheath BCs at z_min/z_max connect SOL field lines to the
// limiter plates after one toroidal transit, while the core (x < x_LCFS) is
// periodic in z with a twist-shift in y of Cy*iota*Lz, reusing the tokamak
// inner-wall-limited machinery.
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

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
  int cdim, vdim;
  // Geometry and magnetic field parameters
  double Z_axis, R0, a_mid, x_inner, r0, B0, kappa_ell, iota0, Cy, Bref, x_LCFS;
  int NFP;
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

// Rotational transform profile. Set to the on-axis value of the vacuum
// rotating ellipse (nearly shearless), constant in radius.
double iota_func(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  return app->iota0;
}

// Cross-section shape functions: (R,Z) offsets per unit minor radius and
// their theta derivatives, for an ellipse with semi-axes sqrt(kappa_ell) and
// 1/sqrt(kappa_ell) rotated poloidally by NFP*phi/2. The normalization is
// area preserving: u*dv/dtheta - v*du/dtheta = 1 identically, so the Jacobian
// of the (r,theta,phi) -> (X,Y,Z) map is J = r*R exactly.
static inline void ellipse_uv(double theta, double phi, double kappa_ell, int NFP,
  double *u, double *v, double *u_th, double *v_th)
{
  double sa = sqrt(kappa_ell), sb = 1./sqrt(kappa_ell);
  double c = cos(0.5*NFP*phi), s = sin(0.5*NFP*phi);
  double ct = cos(theta), st = sin(theta);
  *u = sa*ct*c - sb*st*s;
  *v = sb*st*c + sa*ct*s;
  *u_th = -sa*st*c - sb*ct*s;
  *v_th =  sb*ct*c - sa*st*s;
}

// Major radius and height on the surface labeled by r, at angles (theta,phi).
double R_rtp(double r, double theta, double phi, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double u, v, u_th, v_th;
  ellipse_uv(theta, phi, app->kappa_ell, app->NFP, &u, &v, &u_th, &v_th);
  return app->R0 + r*u;
}

double Z_rtp(double r, double theta, double phi, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double u, v, u_th, v_th;
  ellipse_uv(theta, phi, app->kappa_ell, app->NFP, &u, &v, &u_th, &v_th);
  return app->Z_axis + r*v;
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

// Geometry evaluation functions for the gk app.
// Computational (x,y,z) -> physical Cartesian (X,Y,Z), with z field aligned:
// following z at fixed (x,y) advances phi by dz and theta by iota*dz.
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];

  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double Cy = app->Cy;

  double r = r_x(x,a_mid,x_inner);
  double phi = z;
  double theta = -y/Cy + iota_func(r,ctx)*z;

  double u, v, u_th, v_th;
  ellipse_uv(theta, phi, app->kappa_ell, app->NFP, &u, &v, &u_th, &v_th);
  double R = app->R0 + r*u;
  double Z = app->Z_axis + r*v;

  xp[0] = R*cos(phi); xp[1] = R*sin(phi); xp[2] = Z;
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

// Magnetic field: B = (B0/R) * ( iota(r) e_theta + e_phi ), the Clebsch field
// B = grad(psi) x grad(alpha) with psi = B0*r^2/2 and alpha = theta - iota*phi.
// Cylindrical components:
//   e_theta = ( r du/dtheta, 0, r dv/dtheta ),
//   e_phi   = ( r du/dphi,   R, r dv/dphi ) = ( -r*(NFP/2)*v, R, r*(NFP/2)*u ),
// in the (Rhat, phihat, Zhat) basis.
void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double Cy = app->Cy;
  double x_inner = app->x_inner;
  double B0 = app->B0;

  double r = r_x(x,a_mid,x_inner);
  double phi = z;
  double iota = iota_func(r,ctx);
  double theta = -y/Cy + iota*z;

  double u, v, u_th, v_th;
  ellipse_uv(theta, phi, app->kappa_ell, app->NFP, &u, &v, &u_th, &v_th);
  double R = app->R0 + r*u;

  double B_R   = (B0/R)*r*( iota*u_th - 0.5*app->NFP*v );
  double B_Z   = (B0/R)*r*( iota*v_th + 0.5*app->NFP*u );
  double B_phi = B0; // (B0/R) * R.

  // Set Cartesian components of the magnetic field.
  fout[0] = B_R*cos(phi) - B_phi*sin(phi);
  fout[1] = B_R*sin(phi) + B_phi*cos(phi);
  fout[2] = B_Z;
}

// Twist-shift in y for the core z-periodicity: two points at z_min and z_max
// coincide physically if their y labels differ by -Cy*iota*(z_max-z_min)
// (theta = -y/Cy + iota*z must match modulo 2*pi).
void bc_shift_func_lo(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;
  double r = r_x(x, app->a_mid, app->x_inner);
  fout[0] = -app->Cy*iota_func(r,ctx)*(app->z_max - app->z_min);
}

void bc_shift_func_up(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;
  double r = r_x(x, app->a_mid, app->x_inner);
  fout[0] = app->Cy*iota_func(r,ctx)*(app->z_max - app->z_min);
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
  double Z_axis = 0.0; // Magnetic axis height [m].
  double R0 = 1.5; // Major radius of the flux-surface centers [m].
  double a_mid = 0.25; // LCFS minor radius [m] (area-equivalent radius of the ellipse).
  double B0 = 1.0; // Toroidal magnetic field magnitude [T] (uniform in this model).
  double kappa_ell = 2.5; // Elongation of the rotating elliptical cross section.
  int NFP = 2; // Number of toroidal field periods.
  double x_inner = 0.05; // Radial extent inside LCFS.
  double x_outer = 0.05; // Radial extent outside LCFS.
  // On-axis rotational transform of the vacuum rotating ellipse,
  // iota = (NFP/2)*(a-b)^2/(a^2+b^2) with semi-axes a,b (Helander 2014).
  double sa = sqrt(kappa_ell), sb = 1./sqrt(kappa_ell);
  double iota0 = 0.5*NFP*pow(sa-sb,2)/(sa*sa+sb*sb);

  double r0 = a_mid + 0.5*(x_outer - x_inner); // Minor radius at the center of the simulation box [m].
  double Cy = r0; // Normalization of the binormal coordinate, y = Cy*(theta - iota*phi).

  // Plasma parameters (model SOL-like plasma, same as the TCV IWL tests).
  int num_species = 2;
  double AMU = 2.01410177811;
  double mi = mp*AMU; // Deuterium ions.
  double Te0 = 100*eV; // Electron reference temperature [J].
  double Ti0 = 100*eV; // Ion reference temperature [J].
  double n0  = 2.0e19; // Reference density [1/m^3].
  double Bref = B0; // Reference magnetic field [T].
  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi);
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Configuration domain parameters
  double Lx = x_inner + x_outer; // Domain size along x.
  double x_min = 0.;
  double x_max = Lx;
  double x_LCFS = x_inner; // Radial location of the last closed flux surface.

  double Ly = 150*rho_s; // Domain size along y.
  // Adjust the domain size along y to be an integer fraction of the poloidal
  // circumference: 2*pi*Cy/Ly = integer.
  Ly = 2.*M_PI*Cy/round(2.*M_PI*Cy/Ly);
  double y_min = -Ly/2.;
  double y_max =  Ly/2.;

  double Lz = 2.*M_PI-1e-10; // One full toroidal transit (limiter to limiter).
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
  double P_exp = 0.34e6; // Model input power [W].
  double vol_frac = 1.0/(2.*M_PI*Cy/Ly); // Volume fraction of the simulation box.
  double P_inj = P_exp * vol_frac / num_species; // Injection power normalized to the volume fraction and per species [W]
  // Core source:
  // - Injects energy only in the core region.
  // - The particles injection is only the one that are lost through the inner radial boundary.
  bool adapt_energy_srcCORE = true; // The source will compensate the losses in energy according to given boundaries.
  bool adapt_particle_srcCORE = true; // The source will compensate the losses in particle according to given boundaries.
  double energy_srcCORE = P_inj; // What the source must inject in energy [W]
  double particle_srcCORE = 0.0; // What the source must inject in particle [1/s]
  double center_srcCORE[3] = {x_min, 0.0, -Lz/4}; // This is the position of the ion source,
  double sigma_srcCORE[3] = {0.03*Lx, 0.0, Lz/6}; // the electron source will be at +Lz/2.
  double floor_srcCORE = 1e-10;
  // Recycling source:
  // - Reinjects particles that are absorbed by the limiter.
  // - Energy is free to leave the system.
  bool adapt_energy_srcRECY = false;
  bool adapt_particle_srcRECY = true;
  double energy_srcRECY = 0.0; // [W]
  double particle_srcRECY = 0.0; // [1/s]
  double center_srcRECY[3] = {0.5*x_LCFS, 0.0, M_PI};
  double sigma_srcRECY[3] = {0.25*x_LCFS, 0.0, 0.05*Lz};
  double floor_srcRECY = 1e-10;

  // Grid parameters (reduced resolution for the regression test, minimal recommended values in comments)
  int num_cell_x = 10; // (24) The LCFS is positionned at 1/2 of the domain -> the resolution must be even.
  int num_cell_y = 4; // (16)
  int num_cell_z = 12; // (24) The cross section rotates with period 2*pi/NFP in z.
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
    .R0 = R0,
    .a_mid = a_mid,
    .x_inner = x_inner,
    .r0 = r0,
    .B0 = B0,
    .kappa_ell = kappa_ell,
    .NFP = NFP,
    .iota0 = iota0,
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
    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = &ctx,
    .bfield_func = bfield_func, // magnetic field
    .bfield_ctx = &ctx,
    .non_axisymmetric = true, // Stellarator geometry: skip the axisymmetry check.
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
