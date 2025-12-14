#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_elasticity.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct elasticity_fivewave_ctx
{
  // Physical constants (using normalized code units).
  double rho_ref; // Reference density (unstressed configuration).
  double T_ref; // Reference temperature (unstressed configuration).
  double sound_speed; // Sound speed (speed of pressure waves).
  double shear_speed; // Shear speed (speed of shear waves).
  double heat_capacity; // Specific heat capacity at constant volume.

  double alpha_param; // Alpha parameter (nonlinear elasticity).
  double beta_param; // Beta parameter (nonlinear elasticity).
  double gamma_param; // Gamma parameter (nonlinear elasticity).

  double vel_y_l; // Left solid y-velocity.
  double Fxx_l; // Left solid deformation gradient tensor (xx-component).
  double Fxy_l; // Left solid deformation gradient tensor (xy-component).
  double entropy_l; // Left solid entropy.

  double vel_y_r; // Right solid y-velocity.
  double Fxx_r; // Right solid deformation gradient tensor (xx-component).
  double Fxy_r; // Right solid deformation gradient tensor (xy-component).
  double entropy_r; // Right solid entropy.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct elasticity_fivewave_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double rho_ref = 8.93; // Reference density (unstressed configuration).
  double T_ref = 300.0; // Reference temperature (unstressed configuration).
  double sound_speed = 4.6; // Sound speed (speed of pressure waves).
  double shear_speed = 2.1; // Shear speed (speed of shear waves).
  double heat_capacity = 3.9 * pow(10.0, -4.0); // Specific heat capacity at constant volume.

  double alpha_param = 1.0; // Alpha parameter (nonlinear elasticity).
  double beta_param = 3.0; // Beta parameter (nonlinear elasticity).
  double gamma_param = 2.0; // Gamma parameter (nonlinear elasticity).

  double vel_y_l = 1.0; // Left solid y-velocity.
  double Fxx_l = 0.95; // Left solid deformation gradient tensor (xx-component).
  double Fxy_l = 0.05; // Left solid deformation gradient tensor (xy-component).
  double entropy_l = pow(10.0, -3.0); // Left solid entropy.

  double vel_y_r = 0.0; // Right solid y-velocity.
  double Fxx_r = 1.0; // Right solid deformation gradient tensor (xx-component).
  double Fxy_r = 0.0; // Right solid deformation gradient tensor (xy-component).
  double entropy_r = 0.0; // Right solid entropy.

  // Simulation parameters.
  int Nx = 4096; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double cfl_frac = 0.5; // CFL coefficient.

  double t_end = 0.06; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct elasticity_fivewave_ctx ctx = {
    .rho_ref = rho_ref,
    .T_ref = T_ref,
    .sound_speed = sound_speed,
    .shear_speed = shear_speed,
    .heat_capacity = heat_capacity,
    .alpha_param = alpha_param,
    .beta_param = beta_param,
    .gamma_param = gamma_param,
    .vel_y_l = vel_y_l,
    .Fxx_l = Fxx_l,
    .Fxy_l = Fxy_l,
    .entropy_l = entropy_l,
    .vel_y_r = vel_y_r,
    .Fxx_r = Fxx_r,
    .Fxy_r = Fxy_r,
    .entropy_r = entropy_r,
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
evalElasticityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct elasticity_fivewave_ctx *app = ctx;

  double rho_ref = app->rho_ref;
  double T_ref = app->T_ref;
  double sound_speed = app->sound_speed;
  double shear_speed = app->shear_speed;
  double heat_capacity = app->heat_capacity;

  double alpha_param = app->alpha_param;
  double beta_param = app->beta_param;
  double gamma_param = app->gamma_param;

  double vel_y_l = app->vel_y_l;
  double Fxx_l = app->Fxx_l;
  double Fxy_l = app->Fxy_l;
  double entropy_l = app->entropy_l;

  double vel_y_r = app->vel_y_r;
  double Fxx_r = app->Fxx_r;
  double Fxy_r = app->Fxy_r;
  double entropy_r = app->entropy_r;

  double vel_y = 0.0;
  double Fxx = 0.0;
  double Fxy = 0.0;
  double entropy = 0.0;

  if (x < 0.5) {
    vel_y = vel_y_l; // Solid y-velocity (left).
    Fxx = Fxx_l; // Solid deformation gradient tensor (xx-component, left).
    Fxy = Fxy_l; // Solid deformation gradient tensor (xy-component, left).
    entropy = entropy_l; // Solid entropy (left).
  }
  else {
    vel_y = vel_y_r; // Solid y-velocity (right).
    Fxx = Fxx_r; // Solid deformation gradient tensor (xx-component, right).
    Fxy = Fxy_r; // Solid deformation gradient tensor (xy-component, right).
    entropy = entropy_r; // Solid entropy (right).
  }

  double deformation_gradient[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        deformation_gradient[i][j] = 1.0;
      }
      else {
        deformation_gradient[i][j] = 0.0;
      }
    }
  }

  deformation_gradient[0][0] = Fxx;
  deformation_gradient[0][1] = Fxy;

  double deformation_gradient_det = (deformation_gradient[0][0] * ((deformation_gradient[1][1] * deformation_gradient[2][2]) - (deformation_gradient[2][1] * deformation_gradient[1][2]))) -
    (deformation_gradient[0][1] * ((deformation_gradient[1][0] * deformation_gradient[2][2]) - (deformation_gradient[1][2] * deformation_gradient[2][0]))) +
    (deformation_gradient[0][2] * ((deformation_gradient[1][0] * deformation_gradient[2][1]) - (deformation_gradient[1][1] * deformation_gradient[2][0])));

  double bulk_modulus = (sound_speed * sound_speed) - ((4.0 / 3.0) * (shear_speed * shear_speed));
  double shear_modulus = shear_speed * shear_speed;

  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += deformation_gradient[i][i];
  }

  double deformation_gradient_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        deformation_gradient_sq[i][j] += deformation_gradient[i][k] * deformation_gradient[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += deformation_gradient_sq[i][i];
  }

  double identity_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        identity_tensor[i][j] = 1.0;
      }
      else {
        identity_tensor[i][j] = 0.0;
      }
    }
  }
  
  double inv_deformation_gradient[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient[i][j] = (1.0 / deformation_gradient_det) *
        ((0.5 * ((trace * trace) - sq_trace) * identity_tensor[i][j]) - (trace * deformation_gradient[i][j]) + deformation_gradient_sq[i][j]);
    }
  }

  double inv_deformation_gradient_transpose[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_transpose[i][j] = inv_deformation_gradient[j][i];
    }
  }

  double strain_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor[i][j] += inv_deformation_gradient_transpose[i][k] * inv_deformation_gradient[k][j];
      }
    }
  }

  double strain_tensor_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_sq[i][j] += strain_tensor[i][k] * strain_tensor[k][j];
      }
    }
  }

  double strain_invariant1 = 0.0;
  for (int i = 0; i < 3; i++) {
    strain_invariant1 += strain_tensor[i][i];
  }

  double sq_strain_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace += strain_tensor_sq[i][i];
  }

  double strain_invariant2 = 0.5 * ((strain_invariant1 * strain_invariant1) - sq_strain_trace);

  double strain_invariant3 = (strain_tensor[0][0] * ((strain_tensor[1][1] * strain_tensor[2][2]) - (strain_tensor[2][1] * strain_tensor[1][2]))) -
    (strain_tensor[0][1] * ((strain_tensor[1][0] * strain_tensor[2][2]) - (strain_tensor[1][2] * strain_tensor[2][0]))) +
    (strain_tensor[0][2] * ((strain_tensor[1][0] * strain_tensor[2][1]) - (strain_tensor[1][1] * strain_tensor[2][0])));

  double internal_energy = ((bulk_modulus / (2.0 * (alpha_param * alpha_param))) * (pow(strain_invariant3, 0.5 * alpha_param) - 1.0) * (pow(strain_invariant3, 0.5 * alpha_param) - 1.0)) +
    (heat_capacity * T_ref * pow(strain_invariant3, 0.5 * gamma_param) * (exp(entropy / heat_capacity) - 1.0)) + (0.5 * shear_modulus * pow(strain_invariant3, 0.5 * beta_param) *
    (((strain_invariant1 * strain_invariant1) / 3.0) - strain_invariant2));

  double rho = rho_ref / deformation_gradient_det; // Solid mass density.
  double mom_x = 0.0; // Solid momentum density (x-direction).
  double mom_y = rho * vel_y; // Solid momentum density (y-direction).
  double mom_z = 0.0; // Solid momentum density (z-direction).

  double deformation_gradient_cons[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_cons[i][j] = rho * deformation_gradient[i][j]; // Solid deformation gradient tensor.
    }
  }

  double E_tot = rho * (internal_energy + (0.5 * (vel_y * vel_y))); // Solid total energy density.

  // Set solid mass density.
  fout[0] = rho;
  // Set solid momentum density.
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;

  // Set solid deformation gradient tensor
  fout[4] = deformation_gradient_cons[0][0]; fout[5] = deformation_gradient_cons[0][1]; fout[6] = deformation_gradient_cons[0][2];
  fout[7] = deformation_gradient_cons[1][0]; fout[8] = deformation_gradient_cons[1][1]; fout[9] = deformation_gradient_cons[1][2];
  fout[10] = deformation_gradient_cons[2][0]; fout[11] = deformation_gradient_cons[2][1]; fout[12] = deformation_gradient_cons[2][2];

  // Set solid total energy density.
  fout[13] = E_tot;

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

  struct elasticity_fivewave_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Fluid equations.
  struct gkyl_wv_eqn *elasticity = gkyl_wv_elasticity_new(ctx.T_ref, ctx.sound_speed, ctx.shear_speed, ctx.heat_capacity, ctx.alpha_param, ctx.beta_param,
    ctx.gamma_param, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "elasticity",
    .equation = elasticity,
    
    .init = evalElasticityInit,
    .ctx = &ctx,
    .force_low_order_flux = true,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
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
    .name = "elasticity_fivewave",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
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
  gkyl_wv_eqn_release(elasticity);
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
