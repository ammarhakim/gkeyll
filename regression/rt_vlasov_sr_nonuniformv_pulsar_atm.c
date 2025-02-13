#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct pulsar_atm_ctx 
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Positron mass.
  double charge_ion; // Positron charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double n0; // Reference density.
  double T; // Temperature (units of mc^2).
  double J0; // doubles as curl(B) and gravity
  double grav; // gravitational scale height = g/T
  double noise_amp; // Noise level for perturbation.
  int mode_init; // Initial wave mode to perturb with noise.
  int mode_final; // Final wave mode to perturb with noise.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max; // Domain boundary (velocity space: vx-direction).
  double nonuniform_v_pow; // How nonuniform is velocity space? 
                           // velocity space grid points are spaced as v^n power.
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct pulsar_atm_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Positron mass.
  double charge_ion = 1.0; // Positron charge.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double n0 = 1.0; // Reference density.
  double T = 0.1; // Reference temperature (in units of mc^2).
  double J0 = 2.0; // Reference Goldreich-Julian current (>1 = super-GJ).
  double grav = 1.0; // Reference strength of gravity.

  double noise_amp = 1.0e-6; // Noise level for perturbation.
  int mode_init = 1; // Initial wave mode to perturb with noise.
  int mode_final = 4; // Final wave mode to perturb with noise.

  // Simulation parameters.
  int Nx = 64; // Cell count (configuration space: x-direction).
  int Nvx = 64; // Cell count (velocity space: vx-direction).
  double Lx = 100.0; // Domain size (configuration space: x-direction).
  double vx_max = 512.0; // Domain boundary (velocity space: vx-direction).
  double nonuniform_v_pow = 2.0; // How nonuniform is velocity space: grid points are spaced as v^n power
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 10.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  return (struct pulsar_atm_ctx) {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .n0 = n0, 
    .T = T,
    .J0 = J0,
    .grav = grav,
    .noise_amp = noise_amp,
    .mode_init = mode_init,
    .mode_final = mode_final,
    .Nx = Nx,
    .Nvx = Nvx,
    .Lx = Lx,
    .vx_max = vx_max,
    .nonuniform_v_pow = nonuniform_v_pow, 
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0];

  double n = app->n0;
  double grav = app->grav;

  // Set density (no perturbation).
  fout[0] = 10000.0*exp(-grav*x);
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct pulsar_atm_ctx *app = ctx;

  double T = app->T;

  // Set temperature.
  fout[0] = T;
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct pulsar_atm_ctx *app = ctx;

  // Set drift (four-) velocity.
  fout[0] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0];
  double grav = app->grav;

  pcg64_random_t rng = gkyl_pcg64_init(0); // RNG for use in IC
  double Lx = app->Lx;
  double kx = 2.0*M_PI/Lx;
  double noise_amp = app->noise_amp;

  // linear profile in E_x corresponds to potential drop
  double E_x = x;
  /* // add perturbations to Ex to excite instabilities */
  /* for (int i=mode_init; i<mode_final; ++i)  */
  /*   E_x += noise_amp*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng)); */

  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalAccelFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0];
  double grav = app->grav;
  double T = app->T;
  fout[0] = -grav*T; fout[1] = 0.0, fout[2] = 0.0;
}

void
evalCurrentFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0];
  double J0 = app->J0;
  fout[0] = J0; fout[1] = 0.0, fout[2] = 0.0;
}

static void
mapc2p_vel(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double vx_max = app->vx_max;
  double nonuniform_v_pow = app->nonuniform_v_pow;

  if (vc[0] < 0.0) {
    vp[0] = -vx_max*pow(vc[0], nonuniform_v_pow);
  }
  else {
    vp[0] =  vx_max*pow(vc[0], nonuniform_v_pow);  
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, iot->curr - 1);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, iot->curr - 1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct pulsar_atm_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  // Create global range.
  int ccells[] = { NX };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);
  struct gkyl_range cglobal_r;
  gkyl_create_global_range(cdim, ccells, &cglobal_r);

  // Create decomposition.
  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif  
    
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &cglobal_r);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -1.0 },
    .upper = { 1.0 }, 
    .cells = { NVX },

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
      .use_last_converged = true, 
    },

    .mapc2p_vel[0] = {
      .mapc2p_vel_func = mapc2p_vel, 
      .mapc2p_vel_ctx = &ctx, 
    },

    .app_accel = evalAccelFunc,
    .app_accel_ctx = &ctx,

    .bcx = {
      .lower = { .type = GKYL_SPECIES_FIXED_FUNC, },
      .upper = { .type = GKYL_SPECIES_COPY, },
    },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // positrons
  struct gkyl_vlasov_species pos = {
    .name = "pos",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -1.0 },
    .upper = { 1.0 }, 
    .cells = { NVX },

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
    },

    .mapc2p_vel[0] = {
      .mapc2p_vel_func = mapc2p_vel, 
      .mapc2p_vel_ctx = &ctx, 
    },

    .app_accel = evalAccelFunc,
    .app_accel_ctx = &ctx,

    .bcx = {
      .lower = { .type = GKYL_SPECIES_FIXED_FUNC, },
      .upper = { .type = GKYL_SPECIES_COPY, },
    },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,

    .app_current = evalCurrentFunc,
    .app_current_ctx = &ctx,
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {
    .name = "vlasov_sr_nonuniformv_pulsar_atm",
    
    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, pos },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };
  
  // Create app object.
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_vlasov_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  gkyl_vlasov_app_calc_field_energy(app, t_curr);  
  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
    gkyl_vlasov_app_calc_field_energy(app, t_curr);  
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_vlasov_app_write_integrated_mom(app);
  gkyl_vlasov_app_write_field_energy(app);  
  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // Free resources after simulation completion.
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
