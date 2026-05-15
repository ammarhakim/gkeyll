#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <gkyl_math.h>
#include <rt_arg_parse.h>

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
    char sim_name[128]; // Simulation name.
    int cdim, vdim;
    // Geometry and magnetic field parameters
    double a_shift, Z_axis, R_axis, R0, a_mid, r0, B0, kappa, delta, q0, Cy, qaxis, qlcfs;
    // Plasma parameters
    double me, qe, mi, qi, n0, Te0, Ti0;
    bool static_elc; // Whether to use static electrons.
    // Collision parameters
    double nuFrac;
    // Initial condition parameters
    double Ln, LTe, LTi;
    // Grid parameters
    double Lx, Lz;
    double x_min, x_max, z_min, z_max;
    int num_cell_x, num_cell_z, num_cell_vpar, num_cell_mu;
    int cells[GKYL_MAX_DIM], poly_order;
    double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
    // Simulation control parameters
    double final_time, write_phase_freq;
    int num_frames, int_diag_calc_num, num_failures_max;
    double dt_failure_tol;
    double max_run_time; // Maximum run time in seconds, 0 means no limit.
};

// Geometry related functions 
double r_x(double x, double r0)
{
  return x + r0;
} 

// quadratic q profile
double qprofile(double r, double a_mid, double qaxis, double qlcfs)
{
  // Profile from Grandgirard et al. 2008
  return 1.0 + 2.78 * pow(r/a_mid,2.8);

}

double R_rtheta(double r, double theta, void *ctx)
{
  // Major radius as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double a_shift = app->a_shift;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return R_axis - a_shift*r*r/(2.*R_axis) + r*cos(theta + asin(delta)*sin(theta));
}

double Z_rtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Z_axis = app->Z_axis;
  double kappa = app->kappa;
  return Z_axis + kappa*r*sin(theta);
}

double dRdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double a_shift = app->a_shift;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return -a_shift*r/(R_axis) + cos(theta + asin(delta)*sin(theta));
}

double dRdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double delta = app->delta;
  return -r*sin(theta + asin(delta)*sin(theta))*(1.+asin(delta)*cos(theta));
}

double dZdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*sin(theta);
}

double dZdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*r*cos(theta);
}

double Jr(double r, double theta, void *ctx)
{
  return R_rtheta(r,theta,ctx)*( dRdr(r,theta,ctx) * dZdtheta(r,theta,ctx)
		                -dRdtheta(r,theta,ctx) * dZdr(r,theta,ctx) );
}

struct integrand_ctx {
  struct gk_app_ctx *app_ctx;
  double r;
  double theta;
};

double integrand(double t, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  double r = inctx->r;
  struct gk_app_ctx *app = inctx->app_ctx;
  return Jr(r,t,app) / pow(R_rtheta(r,t,app),2);
}

double Bphi(double R, void *ctx)
{
  // Toroidal magnetic field.
  struct gk_app_ctx *app = ctx;
  double B0 = app->B0;
  double R0 = app->R0;
  return B0*R0/R;
}

double dPsidr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);

  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);
  double R_omp = R_rtheta(r,0.0,ctx);
  return ( R*Bt/(2.*M_PI*qprofile(R_omp, app->a_mid, app->qaxis, app->qlcfs)))*integral.res;
}

double integrant_dpsi(double r, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  struct gk_app_ctx *app = inctx->app_ctx;
  double q = qprofile(r, app->a_mid, app->qaxis, app->qlcfs);
  double dpsi = -r/q; // This is valid in the circular limit. The minus sign comes dBtheta/r < 0
  return dpsi; // I am not sure about this sign, I think it must be from dx/dr = -1.
  // return dPsidr(r, inctx->theta, app); // Exact but way more expensive.
}

double intPsi(double r0, double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .theta = theta};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrant_dpsi, &tmp_ctx, r0, r, 7, 1e-10);
  return integral.res;
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

  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);

  return phi - R*Bt*integral.res/dPsidr(r,theta,ctx);
}

double gradr(double r, double theta, void *ctx)
{
  return (R_rtheta(r,theta,ctx)/Jr(r,theta,ctx))*sqrt(pow(dRdtheta(r,theta,ctx),2) + pow(dZdtheta(r,theta,ctx),2));
}

void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double Cy = app->Cy;
  double r0 = app->r0;
  double r = r_x(x,r0);
  double Bt = Bphi(R_rtheta(r,z,ctx),ctx);
  double Bp = dPsidr(r,z,ctx)/R_rtheta(r,z,ctx)*gradr(r,z,ctx);

  double drdtheta = dRdtheta(r,z,ctx);
  double dzdtheta = dZdtheta(r,z,ctx);
  double den = sqrt(pow(drdtheta,2) + pow(dzdtheta,2));
  double B_r = Bp*drdtheta/den;
  double B_z = Bp*dzdtheta/den;
  double phi = y/Cy + alpha(r, z, 0, ctx);
  double R = R_rtheta(r, z, ctx);

  // xc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = B_r * cos(phi) + Bt * sin(phi);
  fout[1] = B_r * sin(phi) - Bt * cos(phi);
  fout[2] = B_z;
}

double Bmag(double r, double theta, void *ctx)
{
  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);
  double Bp = dPsidr(r,theta,ctx)/R_rtheta(r,theta,ctx)*gradr(r,theta,ctx);

  return sqrt(pow(Bt,2) + pow(Bp,2));
}

double rbar(double m, double q, double r, double theta, double vpar, double mu, void *ctx)
{
  struct gk_app_ctx *app = ctx;

  double e = GKYL_ELEMENTARY_CHARGE;
  double R = R_rtheta(r,theta,ctx);
  double r0 = app->r0;
  double R0 = app->R0;
  double q0 = app->q0;
  double rpsi = q0/r0 * (intPsi(0, r, theta, ctx) - intPsi(0, r0, theta, ctx));
  
  double rmin = r0 - app->x_min;
  double rmax = r0 + app->x_max;
  double zmin = app->z_min;
  double zmax = app->z_max;
  double Bmax = Bmag(rmax, zmax, ctx);
  
  double B0 = app->B0;
  double B = Bmag(r,theta,ctx);
  double E = 0.5 * m * pow(vpar,2) + mu*B;  
  double sgnvpar = vpar>0 ? 1.0 : -1.0;
  double vparbar = E - mu*Bmax > 0 ? 
    sgnvpar*sqrt(2/m)*sqrt(E - mu*Bmax) : 0.0;

  double rvpar = m*q0/(q*B0*r0) * (R*vpar - R0*vparbar);

  double rbar = r0 - rpsi - rvpar;

  return rbar;
}

// Density initial condition (like TCV exp profile)
double density_init(double x, void *ctx)
{
struct gk_app_ctx *app = ctx;
  
  // GYSELA CBC width parameter for density 
  double delta_rho_n = 0.25 * app->Lx; 
  double a_mid = app->a_mid;
  double Ln = app->Ln;
  // The argument for the tanh/cosh functions
  double arg = x / (a_mid * delta_rho_n);
  // Integrating the cosh^-2 gradient yields a tanh profile.
  // The prefactor (a_mid * delta_rho_n / Ln) ensures the peak gradient exactly matches 1/Ln.
  double prof_factor = (a_mid * delta_rho_n) / Ln;
  
  return app->n0 * exp(-prof_factor * tanh(arg));
}

// Electron temperature initial conditions
double temp_init_elc(double x, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  
  // GYSELA CBC width parameter for temperature 
  double delta_rho_T = 0.25 * app->Lx;
  double a_mid = app->a_mid;
  double LTe = app->LTe;
  double arg = x / (a_mid * delta_rho_T);
  double prof_factor = (a_mid * delta_rho_T) / LTe;
  
  return app->Te0 * exp(-prof_factor * tanh(arg));
}

// Ion temperature initial conditions
double temp_init_ion(double x, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  
  // GYSELA CBC width parameter for temperature 
  double delta_rho_T = 0.25 * app->Lx;
  double a_mid = app->a_mid;
  double LTi = app->LTi;

  double arg = x / (a_mid * delta_rho_T);
  double prof_factor = (a_mid * delta_rho_T) / LTi;
  
  return app->Ti0 * exp(-prof_factor * tanh(arg));
}

void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_app_ctx *app = ctx;

  fout[0] = 0.5; // Diffusivity [m^2/s].
}

double maxellian(double m, double n, double T, double E)
{
  return n/pow(2*M_PI*T/m,1.5) * exp(-E/T);
}

void eval_canon_maxwellian_e(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0], z = xn[1], vpar = xn[2], mu = xn[3];
  double m = app->me;

  double r = r_x(x, app->r0);
  double theta = z;

  double rb = rbar(app->me, app->qe, r, theta, vpar, mu, ctx);
  double xb = rb - app->r0;
  double dens = density_init(xb, ctx);
  double temp = temp_init_elc(xb, ctx);

  double energy = 0.5*m*pow(vpar,2) + mu*Bmag(r,theta,ctx);

  fout[0] = maxellian(m, dens, temp, energy);
}

void eval_canon_maxwellian_i(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0], z = xn[1], vpar = xn[2], mu = xn[3];
  double m = app->mi;

  double r = r_x(x, app->r0);
  double theta = z;

  double rb = rbar(app->mi, app->qi, r, theta, vpar, mu, ctx);
  double xb = rb - app->r0;
  double dens = density_init(xb, ctx);
  double temp = temp_init_ion(xb, ctx);

  double energy = 0.5*m*pow(vpar,2) + mu*Bmag(r,theta,ctx);

  fout[0] = maxellian(m, dens, temp, energy);
}

// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];

  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double Cy = app->Cy;

  double r = r_x(x,app->r0);

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

  // Linear map for testing.
  vp[0] = vpar_max_elc*cvpar;
  vp[1] = mu_max_elc*cmu;
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

  // Linear map for testing.
  vp[0] = vpar_max_ion*cvpar;
  vp[1] = mu_max_ion*cmu;
}

struct gk_app_ctx create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field. (Taken from Greenfield et al. 1997 table II)
  double a_shift   = 0.0; // Parameter in Shafranov shift.
  double Z_axis    = 0.0; // Magnetic axis height [m].
  double R_axis    = 1.6714; // Magnetic axis major radius [m].
  double B_axis    = 1.54; // Magnetic field at the magnetic axis [T].
  double a_mid     = 0.604; // Minor radius at OBMP [m].
  double R0        = R_axis + 0.5 * a_mid; // Major radius of the simulation box [m].
  double r0        = 0.5*a_mid; // Minor radius of the simulation box [m].
  double B0        = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double kappa     = 1.0; // cirular geometry for comparison with Dimits et al. 2000.
  double delta     = 0.0; //
  double qaxis     = 1.2; // Safety factor at r=0.
  double qlcfs     = 2.0; // Safety factor at the LCFS.

  // Plasma parameters.
  double AMU = 1.0; // Hydrogen.
  double mi  = mp*AMU;
  double Te0 = 2000*eV;
  double Ti0 = 2000*eV;
  double n0  = 4.5e19; // [1/m^3] according to Fig 5 of Greenfield et al. 1997.
  double nuFrac = 1.0; // Collision factor.

  bool static_elc = false; // Whether to use static electrons.

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;
  double rho_i = vti/omega_ci;
  double q0 = qprofile(r0, a_mid, qaxis, qlcfs);
  double Cy = r0/q0; // Cylindrical coordinate shift for field-alignment.

  // Configuration domain parameters 
  double Lx        = 50*rho_s;   // Domain size along x.
  double x_min     = -Lx/2;
  double x_max     = Lx/2;
  double Lz        = 2.*M_PI-1e-10;       // Domain size along magnetic field.
  double z_min     = -Lz/2.;
  double z_max     =  Lz/2.;

  // Initial conditions and gradients
  // Factor to multiply the gradient because the canonical maxwellian formulation reduce the gradient at HFS.
  double canMaxFactor = 1.7;
  double kTi = 6.92 * canMaxFactor; // R/LTi from Dimits et al. 2000.
  double etai = 3.114; // Lni/LTi from Dimits et al. 2000.
  double LTi = R0/kTi; // Ion temperature gradient scale length [m].
  double Ln = LTi*etai; // Density gradient scale length [m].
  double LTe = LTi; // Electron temperature gradient scale length [m].

  double t_unit = R0/c_s;

  // Check the global shear r/q*dq/dr in the middle of the domain
  double rL = r0 - Lx/2;
  double rR = r0 + Lx/2;
  double qL = qprofile(rL, a_mid, qaxis, qlcfs);
  double qR = qprofile(rR, a_mid, qaxis, qlcfs);
  double dq_dr = (qR-qL)/(rR-rL);
  double s0 = r0/q0*dq_dr;

  // Expected ITG growth rate and time scale
  double gamma_itg_norm = 0.3; // Fig. 1 Dimits et al. 2000 (gamma Lne/vti)
  double ky_itg_norm = 0.3; // Fig. 1 Dimits et al. 2000 (k*rho_i)
  double gamma_itg = gamma_itg_norm*vti/Ln;
  double ky_itg = ky_itg_norm/rho_i;
  double tau_itg = 1./gamma_itg;

  double rhostar = rho_s/a_mid;
  double inv_asp_ratio = a_mid/R0;
  
  printf("-- Simulation parameters and info ---\n");
  printf("vte = %1.2e, vti = %1.2e, c_s = %1.1e [m/s]\n", vte, vti, c_s);
  printf("Lx = %1.2g, rho_s = %1.2g [m]\n", Lx, rho_s);
  printf("x_min = %1.2g, x_max = %1.2g [m]\n", x_min, x_max);
  printf("Ln = %1.2g, LTe = %1.2g, LTi = %1.2g [m]\n", Ln, LTe, LTi);
  printf("q0 = %1.2g, qL = %1.2g, qR = %1.2g, s0 = %1.2g\n", q0, qL, qR, s0);
  printf("ε = %1.2g, ⍴* = 1/%2.2g\n", inv_asp_ratio, 1/rhostar);
  printf("R0/c_s = %1.2e [s]\n", t_unit);

  // Grid parameters
  int num_cell_x = 8;
  int num_cell_z = 8;
  int num_cell_vpar = 8;
  int num_cell_mu = 8;
  int poly_order = 1;
  // Velocity box dimensions
  double vpar_max_elc = 4.*vte;
  double mu_max_elc = 7*Te0/B0;
  double vpar_max_ion = 4.*vti;
  double mu_max_ion = 7*Ti0/B0;
  double final_time = 5*tau_itg;
  int num_frames = 100;
  double write_phase_freq = 0.1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .sim_name = "gk_cbc_2x2v_p1",
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .R_axis = R_axis,
    .R0     = R0    ,
    .a_mid  = a_mid ,
    .r0     = r0    ,
    .B0     = B0    ,
    .kappa  = kappa ,
    .delta  = delta ,
    .q0     = q0    ,
    .Cy     = Cy    ,
    .qaxis  = qaxis ,
    .qlcfs  = qlcfs ,
    .Ln     = Ln   ,
    .LTe    = LTe   ,
    .LTi    = LTi   ,
    .Lx     = Lx    ,
    .Lz     = Lz    ,
    .x_min = x_min,  .x_max = x_max,
    .z_min = z_min,  .z_max = z_max,
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
    .nuFrac = nuFrac,
    .static_elc = static_elc,
    .num_cell_x     = num_cell_x,
    .num_cell_z     = num_cell_z,
    .num_cell_vpar  = num_cell_vpar,
    .num_cell_mu    = num_cell_mu,
    .cells = {num_cell_x, num_cell_z, num_cell_vpar, num_cell_mu},
    .poly_order   = poly_order,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,
    .write_phase_freq = write_phase_freq,
    .final_time = final_time,  .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

int 
main(int argc, char **argv)
{
  struct timespec timer_global = gkyl_wall_clock();
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
  int my_rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    gkyl_comm_get_rank(comm, &my_rank);
#endif

  // electrons
  struct gkyl_gyrokinetic_projection elc_bc = {
    .proj_id = GKYL_PROJ_FUNC,
    .func = eval_canon_maxwellian_e,
    .ctx_func = &ctx,
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = eval_canon_maxwellian_e,
      .ctx_func = &ctx,
    },

    .is_static = ctx.static_elc,

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bc, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bc, },
    },
    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
    .time_rate_diagnostics = true,
  };

  // ions
  struct gkyl_gyrokinetic_projection ion_bc = {
    .proj_id = GKYL_PROJ_FUNC,
    .func = eval_canon_maxwellian_i,
    .ctx_func = &ctx,    
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = eval_canon_maxwellian_i,
      .ctx_func = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bc, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bc, },
    },
    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
    .time_rate_diagnostics = true,
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
    .time_rate_diagnostics = true,
  };

  // Geometry
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.},
    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = &ctx,
    .bfield_func = bfield_func, // magnetic field magnitude
    .bfield_ctx = &ctx,
  };

  // Parallelism
  struct gkyl_app_parallelism_inp parallelism = {
    .comm = comm,
    .cuts = {app_args.cuts[0], app_args.cuts[1]},
    .use_gpu = app_args.use_gpu,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cfl_frac_omegaH = 1.0,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim,
    .lower = { ctx.x_min, ctx.z_min },
    .upper = { ctx.x_max, ctx.z_max },
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = geometry,

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 2,
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
    .print_verbosity = {
      .disable_timings = true,
      .enabled = true,
      .frequency = 100,
    }
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