#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

void
calc_integrated_diagnostics_singleb(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    if ( !(dt < 0.0) ) {
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
    }
  }
}

void
write_data_singleb(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
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

void gyrokinetic_run_singleb_simulation(struct gkyl_gyrokinetic_run_inp inp)
{
  struct gkyl_gyrokinetic_timings timing = inp.timing;

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(inp.app_inp);
 
  double t_curr = 0.0, t_end = timing.t_end; // Initial and final simulation times.
  int frame_curr = 0; // Initialize simulation.

  if (timing.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, timing.restart_frame);

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
  int num_frames = timing.num_frames, num_int_diag_calc = timing.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(timing.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
                                                .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics_singleb(&trig_calc_intdiag, app, t_curr, timing.is_restart, false, -1.0);
  write_data_singleb(&trig_write_conf, &trig_write_phase, app, t_curr, timing.is_restart, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = timing.dt_failure_tol;
  int num_failures = 0, num_failures_max = timing.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= timing.num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics_singleb(&trig_calc_intdiag, app, t_curr, false, t_curr > t_end, status.dt_actual);
    write_data_singleb(&trig_write_conf, &trig_write_phase, app, t_curr, false, t_curr > t_end);

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
        calc_integrated_diagnostics_singleb(&trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
        write_data_singleb(&trig_write_conf, &trig_write_phase, app, t_curr, false, true);
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
}

void
calc_integrated_diagnostics_multib(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_multib_app_save_dt(app, t_curr, dt);
  }
}

void
write_data_multib(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_multib_app* app, double t_curr, bool is_restart_IC, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_multib_app_write_conf(app, t_curr, frame);
    if (!is_restart_IC) {
      gkyl_gyrokinetic_multib_app_write_field_energy(app);
      gkyl_gyrokinetic_multib_app_write_integrated_mom(app);
      gkyl_gyrokinetic_multib_app_write_dt(app);
    }
  }
  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_multib_app_write_phase(app, t_curr, frame);
  }
}

void gyrokinetic_run_multib_simulation(struct gkyl_gyrokinetic_run_inp inp)
{
  struct gkyl_gyrokinetic_timings timing = inp.timing;

  // Create app object.
  gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(inp.multib_app_inp);

  double t_curr = 0.0, t_end = timing.t_end; // Initial and final simulation times.
  int frame_curr = 0; // Initialize simulation.

  if (timing.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_multib_app_read_from_frame(app, timing.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_multib_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_multib_app_apply_ic(app, t_curr);
  }

  // Create triggers for IO.
  int num_frames = timing.num_frames, num_int_diag_calc = timing.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(timing.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics_multib(&trig_calc_intdiag, app, t_curr, timing.is_restart, false, -1.0);
  write_data_multib(&trig_write_conf, &trig_write_phase, app, t_curr, timing.is_restart, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = timing.dt_failure_tol;
  int num_failures = 0, num_failures_max = timing.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= timing.num_steps)) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_multib_update(app, dt);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_multib_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics_multib(&trig_calc_intdiag, app, t_curr, false, t_curr > t_end, status.dt_actual);
    write_data_multib(&trig_write_conf, &trig_write_phase, app, t_curr, false, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_multib_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics_multib(&trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
        write_data_multib(&trig_write_conf, &trig_write_phase, app, t_curr, false, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }
    step += 1;
  }

  gkyl_gyrokinetic_multib_app_stat_write(app);

  // Fetch simulation statistics.
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_multib_app_stat(app);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of write calls %ld.\n", stat.n_io);
  gkyl_gyrokinetic_multib_app_print_timings(app, stdout);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release(app);
}

void gkyl_gyrokinetic_run_simulation(struct gkyl_gyrokinetic_run_inp inp)
{
  if (inp.app_inp != 0) {
    gyrokinetic_run_singleb_simulation(inp);
  }
  else if (inp.multib_app_inp != 0) {
    gyrokinetic_run_multib_simulation(inp);
  }
  else {
    fprintf(stderr, "Error: No valid application input provided to run simulation.\n");
    exit(EXIT_FAILURE);
  }

}
