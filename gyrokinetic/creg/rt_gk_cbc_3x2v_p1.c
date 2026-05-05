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
    double a_shift, Z_axis, R_axis, R0, a_mid, r0, B0, kappa, delta, q0, qaxis, qlcfs, Cy;
    // Plasma parameters
    double me, qe, mi, qi, n0, Te0, Ti0;
    // Collision parameters
    double nuFrac;
    // Initial condition parameters
    double Ln, LTe, LTi;
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
  // This profile ensures:
  // - dq/dr = 0 at r=0
  // - q = qaxis at r=0
  // - q = qlcfs at r=a_mid
  return qaxis + (qlcfs-qaxis)*(r*r)/(a_mid*a_mid);
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
  return - a_shift*r/(R_axis) + cos(theta + asin(delta)*sin(theta));
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

void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Density initial condition (like TCV exp profile)
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct gk_app_ctx *app = ctx;
  fout[0] = app->n0 * (1 - x/app->Ln);
}

// Electron temperature initial conditions
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct gk_app_ctx *app = ctx;
  fout[0] = app->Te0 * (1 - x/app->LTe);
}

// Ion temperature initial conditions
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct gk_app_ctx *app = ctx;
  
  fout[0] = app->Ti0 * (1 - x/app->LTi);
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

void bc_shift_func_lo(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;

  double a_mid = app->a_mid;
  double Cy = app->Cy;
  double z_min = app->z_min;
  double z_max = app->z_max;

  double r = r_x(x, app->r0);

  fout[0] = Cy*( alpha(r, z_min, 0.0, ctx) - alpha(r, z_max, 0.0, ctx) );
}

void bc_shift_func_up(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;

  double a_mid = app->a_mid;
  double Cy = app->Cy;
  double z_min = app->z_min;
  double z_max = app->z_max;

  double r = r_x(x, app->r0);

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

  // Geometry and magnetic field. (Taken from Greenfield et al. 1997 table II)
  double a_shift   = 0.0; // Parameter in Shafranov shift.
  double Z_axis    = 0.1414361745; // Magnetic axis height [m].
  double R_axis    = 1.67; // Magnetic axis major radius [m].
  double B_axis    = 1.54; // Magnetic field at the magnetic axis [T].
  double a_mid     = 0.64; // Minor radius at OBMP [m].
  double R_LCFSmid = R_axis + a_mid; // Major radius of the LCFS at the outboard midplane [m].
  double R0        = R_axis + 0.5 * a_mid; // Major radius of the simulation box [m].
  double r0        = 0.5*a_mid; // Minor radius of the simulation box [m].
  double B0        = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double kappa     = 1.0; // cirular geometry for comparison with Dimits et al. 2000.
  double delta     = 0.0; //
  double qaxis     = 1.2; // Safety factor at r=0.
  double qlcfs     = 3.1; // Safety factor at the LCFS.

  // Plasma parameters.
  double AMU = 1.0; // Hydrogen.
  double mi  = mp*AMU;
  double Te0 = 2000*eV;
  double Ti0 = 2000*eV;
  double n0  = 4.5e19; // [1/m^3] according to Fig 5 of Greenfield et al. 1997.
  double nuFrac = 0.1; // Collision factor.

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;
  double q0 = qprofile(r0, a_mid, qaxis, qlcfs);
  double Cy = r0/q0; // Cylindrical coordinate shift for field-alignment.

  // Configuration domain parameters 
  double Lx        = 100*rho_s;   // Domain size along x.
  double x_min     = -Lx/2;
  double x_max     = Lx/2;
  double Ly        = 100*rho_s;           // Domain size along y.
  // Adjust the domain size along y to have integer toroidal mode number.
  // We need: 2*pi*Cy/Ly = integer (Cy = r0/q0)
  Ly = 2.*M_PI*r0/q0/round(2.*M_PI*r0/q0/Ly); 
  double y_min     = -Ly/2.;
  double y_max     =  Ly/2.;  
  double Lz        = 2.*M_PI-1e-10;       // Domain size along magnetic field.
  double z_min     = -Lz/2.;
  double z_max     =  Lz/2.;
  double vol_frac = 1.0/(2.*M_PI*r0/q0/Ly);

  // Initial conditions and gradients
  double LTi = R0/6.92; // Ion temperature gradient scale length [m].
  double etai = 3.114; // Lni/LTi
  double Ln = LTi*etai; // Density gradient scale length [m].
  double LTe = LTi; // Electron temperature gradient scale length [m].

  double t_unit = R0/c_s;

  printf("R0/c_s = %g\n", t_unit);
  printf("rhos = %g\n", rho_s);
  printf("Lx = %g\n", Lx);
  printf("Ly = %g\n", Ly);

  // Grid parameters
  int num_cell_x = 8;
  int num_cell_y = 8;
  int num_cell_z = 8;
  int num_cell_vpar = 8;
  int num_cell_mu = 6;
  int poly_order = 1;
  // Velocity box dimensions
  double vpar_max_elc = 5.*vte;
  double mu_max_elc   = 1.5*me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 5.*vti;
  double mu_max_ion   = 1.5*mi*pow(4*vti,2)/(2*B0);
  double final_time = 1*t_unit;
  int num_frames = 1;
  double write_phase_freq = 1.0;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .sim_name = "gk_cbc_3x2v_p1",
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
    .Ly     = Ly    ,
    .Lz     = Lz    ,
    .x_min = x_min,  .x_max = x_max,
    .y_min = y_min,  .y_max = y_max,
    .z_min = z_min,  .z_max = z_max,
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
    .nuFrac = nuFrac,
    .num_cell_x     = num_cell_x,
    .num_cell_y     = num_cell_y,
    .num_cell_z     = num_cell_z,
    .num_cell_vpar  = num_cell_vpar,
    .num_cell_mu    = num_cell_mu,
    .cells = {num_cell_x, num_cell_y, num_cell_z, num_cell_vpar, num_cell_mu},
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
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = density_init,
    .ctx_upar = &ctx,
    .upar= zero_func, 
    .ctx_temp = &ctx,
    .temp = temp_elc,      
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

/*
    .init_from_file = {
       .type = GKYL_IC_IMPORT_F,
       .file_name = "restart-elc.gkyl"
    },
//*/  
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
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "ion"},
      .nu_frac = ctx.nuFrac,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bc, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bc, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT, .aux_profile = bc_shift_func_lo, .aux_ctx = &ctx, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT, .aux_profile = bc_shift_func_up, .aux_ctx = &ctx, },
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
  };

  // ions
  struct gkyl_gyrokinetic_projection ion_bc = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = density_init,
    .ctx_upar = &ctx,
    .upar= zero_func, 
    .ctx_temp = &ctx,
    .temp = temp_ion,      
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

/*
    .init_from_file = {
       .type = GKYL_IC_IMPORT_F,
       .file_name = "restart-ion.gkyl"
    },
//*/
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
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "elc"},
      .nu_frac = ctx.nuFrac,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bc, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bc, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT, .aux_profile = bc_shift_func_lo, .aux_ctx = &ctx, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT, .aux_profile = bc_shift_func_up, .aux_ctx = &ctx, },
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
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = 0.0 },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = 0.0 },
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
    .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
    .use_gpu = app_args.use_gpu,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cfl_frac_omegaH = 1.0,
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
      .frequency = 1,
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