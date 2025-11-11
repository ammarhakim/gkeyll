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

// Define the context of the simulation. This is basically all the globals
struct gk_app_ctx {
  int cdim, vdim; // Dimensionality.
  double mu0;

  // Geometry and magnetic field.
  double R; // Major radius of the simulation box [m].
  double a; // Minor radius at outboard midplane [m].
  double r; // Minor radius of the simulation box [m].
  double B; // Magnetic field magnitude in the simulation box [T].
  double z_min;  double z_max;

  double kperp;  // Perpendicular wavenumber.

  // Plasma parameters.
  double me;  double qe;
  double mi;  double qi;
  double n0;  double Te0;  double Ti0; 

  // Grid parameters.
  double kz_min;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double vpar_max_elc;  double mu_max_elc;
  double vpar_max_ion;  double mu_max_ion;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Ion initial conditions.
void density_init_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->n0;
}
void upar_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->Ti0;
}

// Electron initial conditions.
void density_init_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double z = xn[0];
  fout[0] = app->n0 * (1 + 1e-6*cos(app->kz_min*z));
}
void upar_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->Te0;
}

// Slab magnetic equilibrium.
static inline void
mapc2p(double t, const double* xc, double* GKYL_RESTRICT xp, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xc[0], y = xc[1], z = xc[2];
  xp[0] = x; xp[1] = y; xp[2] = z;
}
void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B;
}

struct gk_app_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge
  double me = GKYL_ELECTRON_MASS;
  double mp = GKYL_PROTON_MASS;
  double mu0 = GKYL_MU0;
  // Geometry, magnetic field and gradients.
  double R0 = 1.313; // Major radius of the simulation box [m].
  double a = 0.4701; // Minor radius at outboard midplane [m].
  double B0 = 1.91; // Magnetic field magnitude in the simulation box [T].
  double Ti0 = 2072*eV;
  double eps_n = 0.2;
  // double eta_e = 2.0; // Comes from rt-alfven-p1.lua but unused.
  // double eta_i = 2.5; // Comes from rt-alfven-p1.lua but unused.
  double ky_rhos = 0.01;
  double kz_Ln = 1.0;
  double beta_hat = 10.0; // beta_e/2*m_p/m_e.
  double beta_e = beta_hat*2.0*me/mp;
  double tau = 1.0; // Ti/Te.
  // Derived paramters.
  double r0 = 0.5*a; // Minor radius of the simulation box [m].
  double R = R0 + r0; // Major radius at center of simulation box [m].
  double B = B0 * R0/R; // Magnetic field at center of simulation box [T].
  double Te0 = Ti0/tau; // Electron temperature [J].
  double n0 = beta_e*B*B/(2*mu0*Te0); // Density [m^-3].
  double vte = sqrt(Te0/me); // Electron thermal speed [m/s].
  double vti = sqrt(Ti0/mp); // Ion thermal speed [m/s].
  double c_s = sqrt(Te0/mp); // Sound speed [m/s].
  double omega_ci = fabs(qi*B/mp); // Ion cyclotron frequency [rad/s].
  double rho_s = c_s/omega_ci; // Ion sound gyroradius [m
  double ky_min = ky_rhos/rho_s; // Minimum perpendicular wavenumber [1/m].
  double L_n = R*eps_n; // Density gradient scale length [m].
  double kz_min = kz_Ln/L_n; // Minimum parallel wavenumber [1/m].
  double L_parallel = 2*M_PI/kz_min; // Domain size in z [m].
  // Expected frequency from dispersion relation:
  double omega = kz_min*c_s*sqrt(2/beta_e/(1+beta_hat*pow(ky_rhos,2)));

  // Velocity grid parameters.
  double vpar_max_elc = 6.0*vte; // Max electron parallel velocity [m/s].
  double vpar_min_elc = -vpar_max_elc; // Min electron parallel velocity [m/s].
  double vpar_max_ion = 6.0*vti; // Max ion parallel velocity [m/s].
  double vpar_min_ion = -vpar_max_ion; // Min ion parallel
  double mu_max_elc = me*pow(4*vte,2)/(2*B); // Max electron mu [J/T].
  double mu_max_ion = mp*pow(4*vti,2)/(2*B); // Max ion mu [J/T].

  // Grid parameters
  int Nz = 8;
  int Nvpar = 8;
  int Nmu = 6;
  int poly_order = 1;

  double t_end = 2*M_PI/omega; // End time [s].
  int num_frames = 1;
  double write_phase_freq = 0.1; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 50; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mu0 = mu0,
    .R = R,
    .B = B,
    .z_min = -M_PI/kz_min,  .z_max = M_PI/kz_min,
  
    .kperp = ky_min,
    .kz_min = kz_min,
  
    .me = me,  .qe = qe,
    .mi = mp,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
    
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,

    .t_end = t_end,  .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}


int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = density_init_elc,
      .ctx_upar = &ctx,
      .upar= upar_elc,
      .ctx_temp = &ctx,
      .temp = temp_elc,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_EM,
    },

    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = density_init_ion,
      .ctx_upar = &ctx,
      .upar= upar_ion,
      .ctx_temp = &ctx,
      .temp = temp_ion,      
    },

    .is_static = true,
    
    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .kperpSq = ctx.kperp * ctx.kperp,
    .mu0 = ctx.mu0,
  };

  // GK app
  struct gkyl_gk app_inp = {
    .name = "gk_em_alfven_1x2v_p1",

    .cdim = ctx.cdim,
    .lower = { ctx.z_min },
    .upper = { ctx.z_max },
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0., 0.},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func, // magnetic field
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
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
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  return 0;
}
