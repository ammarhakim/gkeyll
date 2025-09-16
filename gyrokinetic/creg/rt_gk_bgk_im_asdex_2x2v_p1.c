#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_asdex_ctx {
  int cdim, vdim; // Dimensionality

  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double c_s; // sound speed
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double B0; // reference magnetic field
  double n0; // reference density
  // Source parameters
  double lambda_source;
  double x_source;
  // Domain parameters.            
  int Nx, Nz; // Number of cells in x,y,z.
  int Nvpar, Nmu; // Number of cells in vpar,mu.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double Lx; // Box size in x
  double Lz; // Box size in z
  // Physical velocity space limits
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  // Computational velocity space limits
  double vpar_min_elc_c, vpar_max_elc_c;
  double mu_min_elc_c, mu_max_elc_c;
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

double random0to1() 
{
  return (double)rand() / (double)RAND_MAX;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double floor = 0.01;

  //fout[0] = fmax(exp(-128.48398224*x+62.52416972), floor); // The coefficients from fitting the data points on Fig. 7 Carralero 2017 Nuclear Fusion.
  fout[0] = fmax(exp(-269.89663849*x+85.77849169), floor); // High density case.
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double eV = GKYL_ELEMENTARY_CHARGE;

  //fout[0] = (-2000.0*x + 360.0)*eV;   // [16, 40] eV
  fout[0] = (-50*tanh((x-0.15)/0.022) + 40) * eV;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double eV = GKYL_ELEMENTARY_CHARGE;

  //fout[0] = (-3000.0*x + 552.0)*eV;   // [36, 72] eV
  fout[0] = (-90*tanh((x-0.15)/0.022) + 72) * eV;
}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double Lz = app->Lz;
  double z_source = -Lz/4.0;
  double S0 = 24.08e22;
  double source_floor = 0.01*S0;

  if (fabs(z-z_source)<Lz/8.0)
    fout[0] = fmax(S0*exp(-(x-x_source)*(x-x_source)/(2*lambda_source*lambda_source)), source_floor);
  else
    fout[0] = source_floor;
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double Lz = app->Lz;
  double z_source = -Lz/4.0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if ((x < x_source + 3*lambda_source) && (fabs(z-z_source)<Lz/8.0))
    fout[0] = 40.0*eV;
  else
    fout[0] = 1.920292202211762*eV;
}

void
eval_temp_ion_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double Lz = app->Lz;
  double z_source = -Lz/4.0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if ((x < x_source + 3*lambda_source) && (fabs(z-z_source)<Lz/8.0))
    fout[0] = 72.0*eV;
  else
    fout[0] = 3.4565259639811785*eV;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  fout[0] = app->nuIon;
}

// Velocity space mappings.
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
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
  struct gk_asdex_ctx *app = ctx;
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

struct gk_asdex_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 40.0*eV;
  double Ti = 72.0*eV;
  double B0 = 2.57; // Magnetic field magnitude in Tesla
  double n0 = 2.0e19; // Particle density in 1/m^3

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;


  // Collision parameters.
  double nuFrac = 1.0;  
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (psi, theta).
  double Lx = 0.172 - 0.150;
  double Lz = (M_PI-1.0e-14)*2;

  // Source parameters.
  double x_source = 0.1534;   
  double lambda_source = 0.0011;  

  // Physical velocity space limits
  double vpar_max_elc = 6.0*vtElc;
  double mu_max_elc = me*(4.0*vtElc)*(4.0*vtElc)/(2.0*B0);

  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = mi*(4.0*vtIon)*(4.0*vtIon)/(2.0*B0);

  // Computational velocity space limits.
  double vpar_min_ion_c = -1.0/sqrt(2.0);
  double vpar_max_ion_c = 1.0/sqrt(2.0);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  // Computational velocity space limits.
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c = 1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  int Nx = 16; // Number of cells in x.
  int Nz = 16; // Number of cells in z, originally 8.
  int Nvpar = 16; // Number of cells in vpar.
  int Nmu = 8; // Number of cells in mu.

  double t_end = 1.0e-8;
  double num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_asdex_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .c_s = c_s, 
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .B0 = B0, 
    .n0 = n0, 
    .Lx = Lx, 
    .Lz = Lz, 
    .lambda_source = lambda_source,
    .x_source = x_source,
    // Physical velocity space limits
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    // Computational velocity space limits
    .vpar_min_elc_c = vpar_min_elc_c,
    .vpar_max_elc_c = vpar_max_elc_c,
    .mu_min_elc_c = mu_min_elc_c,
    .mu_max_elc_c = mu_max_elc_c,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool is_restart_IC, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    if (!is_restart_IC) {
      gkyl_gyrokinetic_app_write_field_energy(app);
      gkyl_gyrokinetic_app_write_integrated_mom(app);
      gkyl_gyrokinetic_app_write_dt(app);
    }
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_asdex_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // electrons
  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density,
    .ctx_upar = &ctx,
    .upar= eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_elc,      
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = elc_ic,
  
    .correct = {
      .correct_all_moms = true, 
      .use_last_converged = true, 
      .iter_eps = 1e-12,
      .max_iter = 10,
    }, 
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .self_nu = evalNuElc,
      .ctx = &ctx,
      .has_implicit_coll_scheme = true,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .write_diagnostics = true, 
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      //.write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,      
      }, 
    },

    .diffusion = {
      .num_diff_dir = 1,
      .diff_dirs = { 0 },
      .D = { 0.1 },
      .order = 2,
    },
      
    .bcx = {
      .lower={
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = elc_ic,
      },
      .upper={.type = GKYL_SPECIES_ABSORB,}, // Changed from ZERO_FLUX on 07/03/2024.
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },

    .num_diag_moments = 9,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_MAXWELLIAN },
  };

  // ions
  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density,
    .ctx_upar = &ctx,
    .upar= eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_ion,      
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = ion_ic,

    .correct = {
      .correct_all_moms = true, 
      .use_last_converged = true, 
      .iter_eps = 1e-12,
      .max_iter = 10,
    }, 
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .self_nu = evalNuIon,
      .ctx = &ctx,
      .has_implicit_coll_scheme = true,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .write_diagnostics = true, 
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      //.write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
    },

    .diffusion = {
      .num_diff_dir = 1,
      .diff_dirs = { 0 },
      .D = { 0.1 },
      .order = 2,
    },
    
    .bcx = {
      .lower={
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = ion_ic,
      },
      .upper={.type = GKYL_SPECIES_ABSORB,}, // Changed from ZERO_FLUX on 07/03/2024.
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },

    .num_diag_moments = 9,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_MAXWELLIAN },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC}, 
                    .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC}, 
                    .lo_value = {0.0, 0.0}, .up_value = {0.0, 0.0}}, 
  };

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk", // equilibrium to use
    .rz_poly_order = 2,                      // polynomial order for psi(R,Z) used for field line tracing
    .flux_poly_order = 1,                    // polynomial order for fpol(psi)
  };

  struct gkyl_tok_geo_grid_inp grid_inp = {
    .ftype = GKYL_LSN_SOL,                     // type of geometry
    .rclose = 2.5,                             // closest R to region of interest
    .rright = 2.5,                             // Closest R to outboard SOL
    .rleft = 0.7,                              // closest R to inboard SOL
    .rmax = 2.5,                               // largest R in machine
    .rmin = 0.7,                               // smallest R in machine
    .zmin = -1.3,                              // Lower Z boundary 
    .zmax = 1.0,                               // Upper Z boundary
    .zmin_left = -1.0,                         // Z of inboard divertor plate
    .zmin_right = -1.0,                        // Z of outboard divertor plate
  };

  // GK app
  struct gkyl_gk app_inp = {
    .name = "gk_bgk_im_asdex_2x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { 0.15, -ctx.Lz/2.0 },   // Originally 0.16
    .upper = { 0.172, ctx.Lz/2.0 },   // Originally 0.175
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 0.05,

    .enforce_positivity = true,

    .geometry = {
      .geometry_id = GKYL_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = grid_inp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  double t_curr = 0.0, t_end = ctx.t_end; // Initial and final simulation times.
  int frame_curr = 0; // Initialize simulation.

  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(ctx.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, app_args.is_restart, false, -1.0);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, app_args.is_restart, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, t_curr > t_end, status.dt_actual);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_gyrokinetic_app_stat_write(app);

  // Fetch simulation statistics.
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
