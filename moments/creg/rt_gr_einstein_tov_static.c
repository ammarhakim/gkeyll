// GR Euler Ultra Rel in Spherical Coorindates in 1D for a static TOV
// Euler eqns from: https://arxiv.org/pdf/gr-qc/9904052

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_moment.h>
#include <gkyl_moment_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_tov.h>
#include "tov_solver_ultra_rel.h"

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPIh
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct gr_tov_static_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double kappa; // Stress-energy prefactor in the Einstein field equations.

  double K_poly;
  double rho_c;
  double dr_tov;

  double rho_atm;
  double p_atm;

  struct gkyl_tov_ultra_rel *tov;
  double M_star;
  double R_star;

  // Simulation parameters.
  int Nx; // Cell count (r-direction).
  double Lx; // Domain size (r-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gr_tov_static_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 4.0 / 3.0; // Adiabatic index.
  double kappa = 8.0 * pi; // Stress-energy prefactor in the Einstein field equations.

  double K_poly = 100.0; 
  double rho_c = 5e-4;
  double dr_tov = 0.01;

  double rho_atm = 1e-10 * rho_c;
  double p_atm = (gas_gamma - 1.0) * rho_atm;

  struct gkyl_tov_ultra_rel *tov = gkyl_tov_ultra_rel_new(gas_gamma, rho_c, dr_tov);

  double M_star = gkyl_tov_ultra_rel_star_mass(tov);
  printf("M_star = %e \n", M_star);
  double R_star = gkyl_tov_ultra_rel_star_radius(tov);
  printf("R_star = %e \n", R_star);

  // Simulation parameters.
  int Nx = 4096; // Cell count (r-direction).
  double Lx = 2500.0; // Domain size (r-direction).
  double cfl_frac = 0.8; // CFL coefficient.

  double t_end = 1500; // Final simulation time.
  int num_frames = 1000; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gr_tov_static_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .kappa = kappa,
    .K_poly = K_poly,
    .rho_c = rho_c,
    .dr_tov = dr_tov,
    .rho_atm = rho_atm,
    .p_atm = p_atm,
    .tov = tov,
    .M_star = M_star,
    .R_star = R_star,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalGRTovInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double r = xn[0]; // r_relative
  struct gr_tov_static_ctx *app = ctx;

  struct gkyl_tov_ultra_rel *tov = app->tov;
  
  struct tov_ultra_rel_eval_fluid fluid = {0}; 

  struct tov_ultra_rel_eval_bl bl;
  gkyl_tov_ultra_rel_eval_bl(tov, r, &bl);

  double lapse = exp(bl.Phi);
  double a = 1.0 / sqrt(1.0 - (2.0 * bl.m / r));

  double v_con[3] = {0.0, 0.0, 0.0}; //fluid coordinate velocity

  double v_sq = 0.0; // For a static TOV
  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     v_sq += spatial_metric[i][j] * v_con[i] * v_con[j];
  //   }
  // }

  double W = 1.0 / (sqrt(1.0 - v_sq)); // Lorentz gamma function
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
  }

  double u_con[4] = {0.0}; //fluid four-velocity
  double shift[3] = {0.0, 0.0, 0.0};
  u_con[0] = W / lapse;
  for (int i = 1; i < 4; i ++) {
    u_con[i] = 0.0; // for a static TOV
  }

  double p = fmax(bl.P, app->p_atm);
  double rho = p / (app->gas_gamma - 1.0); // fluid-frame total energy density

  double Etot = ((rho + p) * W * W) - p; // Eulerian-frame (conserved) energy dneisty 
  double mom_r = ((rho + p) * W * W) * v_con[0]; // radial momenutm density

  //Evolved conservative variables: dummy for D, tau, momentum_r, lapse

  fout[0] = 0.0; // Set a dummy for D in the future non ultra-rel case
  fout[1] = Etot;
  fout[2] = mom_r;
  fout[3] = bl.Phi; // Set Phi appearing in dt metric term
  fout[4] = bl.m;
  fout[5] = r;
  fout[6] = 0;
  fout[7] = 0;

  // for (int i = 0; i < 3; i++) {
  //   gkyl_free(spatial_metric[i]);
  // }
  // gkyl_free(spatial_metric);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
  }
}

static void
write_energy_exact_error_norms(const struct gr_tov_static_ctx *ctx, gkyl_moment_app *app, double t_curr, const char *fname)
{
  struct gkyl_array *q = gkyl_moment_app_get_write_array_species(app, 0);
  double dx = app->grid.dx[0];

  double l1_local = 0.0;
  double l2sq_local = 0.0;
  double linf_local = 0.0;
  double l1_int_local = 0.0;
  double l2sq_int_local = 0.0;
  double linf_int_local = 0.0;
  double r_int_cut = 0.8 * ctx->R_star;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qnum = gkyl_array_cfetch(q, loc);

    double xc[GKYL_MAX_DIM] = { 0.0 };
    double qexact[8] = { 0.0 };
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    evalGRTovInit(0.0, xc, qexact, (void*) ctx);

    double etot_num = qnum[1];
    double etot_exact = qexact[1];
    double err = etot_num - etot_exact;
    double abs_err = fabs(err);

    l1_local += dx * abs_err;
    l2sq_local += dx * err * err;
    linf_local = fmax(linf_local, abs_err);
    if (xc[0] <= r_int_cut) {
      l1_int_local += dx * abs_err;
      l2sq_int_local += dx * err * err;
      linf_int_local = fmax(linf_int_local, abs_err);
  }
  }

  double sum_local[2] = { l1_local, l2sq_local };
  double sum_global[2] = { 0.0, 0.0 };
  double linf_global = 0.0;
  double sum_int_local[2] = { l1_int_local, l2sq_int_local };
  double sum_int_global[2] = { 0.0, 0.0 };
  double linf_int_global = 0.0;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2, sum_local, sum_global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &linf_local, &linf_global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2, sum_int_local, sum_int_global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &linf_int_local, &linf_int_global);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    bool write_header = false;
    FILE *check = fopen(fname, "r");
    if (check == NULL) {
      write_header = true;
    }
    else {
      fclose(check);
    }

    FILE *fp = fopen(fname, "a");
    if (fp == NULL) {
      fprintf(stderr, "Could not open %s for convergence diagnostics.\n", fname);
      return;
    }

    if (write_header) {
      fprintf(fp, "# kind nx nx_fine time L1 L2 Linf\n");
    }
    fprintf(fp, "exact %d 0 %.16e %.16e %.16e %.16e\n",
      app->grid.cells[0], t_curr, sum_global[0], sqrt(sum_global[1]), linf_global);
    fprintf(fp, "interior_exact %d 0 %.16e %.16e %.16e %.16e\n",
      app->grid.cells[0], t_curr, sum_int_global[0], sqrt(sum_int_global[1]), linf_int_global);
    fclose(fp);

    printf("\nFinal Etot exact-error norms written to %s:\n", fname);
    printf("  NX = %d, t = %.16e\n", app->grid.cells[0], t_curr);
    printf("  Global:   L1 = %.16e, L2 = %.16e, Linf = %.16e\n",
      sum_global[0], sqrt(sum_global[1]), linf_global);
    printf("  Interior r <= 0.8 R_star = %.16e:\n", r_int_cut);
    printf("            L1 = %.16e, L2 = %.16e, Linf = %.16e\n",
      sum_int_global[0], sqrt(sum_int_global[1]), linf_int_global);
    printf("  L1 = %.16e, L2 = %.16e, Linf = %.16e\n",
      sum_global[0], sqrt(sum_global[1]), linf_global);
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

  struct gr_tov_static_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  double lower = 0.5 * (ctx.Lx / NX);
  double upper = ctx.Lx + 0.5 * (ctx.Lx / NX);

  // Fluid equations.
  struct gkyl_wv_eqn *gr_tov = gkyl_wv_gr_tov_new(ctx.gas_gamma, ctx.kappa, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "gr_tov",
    .equation = gr_tov,
    
    .init = evalGRTovInit,
    .ctx = &ctx,

    .has_gr_tov = true,
    .tov_gas_gamma = ctx.gas_gamma,
    .tov_kappa = ctx.kappa,
    .tov_p_atm = ctx.p_atm,
    .has_dynamic_lapse = false,

    .force_low_order_flux = false,
    .limiter = GKYL_MIN_MOD,

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "gr_einstein_tov_static",

    .ndim = 1,
    .lower = { lower },
    .upper = { upper },
    .cells = { NX },

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  write_energy_exact_error_norms(&ctx, app, t_curr, "gr_tov_static_energy_conv.dat");
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(gr_tov);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
