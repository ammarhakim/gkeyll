#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

// This file contains standard routines to run gyrokinetic simulations.
// If a user would like to customize the time-stepping loop, they can copy the relevant methods
// from this file into their input file and modify as needed. Make sure not to #include this file
// in that case.

// Step message context.
struct message_trigs {
  int log_count; // Number of times logging called.
  int tenth, p1c; // 10% and 1% counters.
  struct gkyl_tm_trigger log_trig; // 10% trigger.
  struct gkyl_tm_trigger log_trig_1p; // 1% trigger.
  long io_period; // How often to print verbose messages.
  double t_end; // End time of simulation.
  bool estimate_completion_time; // Whether to estimate completion time.
  struct timespec tm_loop_start; // Wall clock time at start of timing window.
  double t_loop_start; // Simulation time at start of timing window.
};

//
// ............. Single block simulations ............... //
// 

// Function pointer types for output handlers.
typedef void (*write_message_pre_update_singleb_t)(const struct gkyl_gyrokinetic_app*, long, double,
  struct message_trigs*);
typedef void (*write_message_post_update_singleb_t)(const struct gkyl_gyrokinetic_app*, long, double,
  double, struct message_trigs*);

// Write log message to console.
// Verbose output functions for single-block
static void
write_message_pre_update_verbose_singleb(const struct gkyl_gyrokinetic_app *app, long step, double t_curr, 
  struct message_trigs *trigs)
{
  if (((step % trigs->io_period) == 0) || step == 1) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %#11.8g ...", step, t_curr);
  }
}

static void
write_message_post_update_verbose_singleb(const struct gkyl_gyrokinetic_app *app, long step, 
  double t_curr, double dt_actual, struct message_trigs *trigs)
{
  if (step == 1) {
    gkyl_gyrokinetic_app_cout(app, stdout, "\tdt = %.6e\n", dt_actual);
  }
  else if ((step % trigs->io_period == 0)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "\tdt = %.6e ", dt_actual);
    double pct_complete = 100.0 * (t_curr - trigs->t_loop_start) / (trigs->t_end - trigs->t_loop_start);
    
    if (trigs->estimate_completion_time) {
      double wall_time_elapsed = gkyl_time_diff_now_sec(trigs->tm_loop_start);
      double sim_time_progressed = t_curr - trigs->t_loop_start;
      double sim_time_remaining = trigs->t_end - t_curr;
      double wall_time_per_sim_time = wall_time_elapsed / sim_time_progressed;
      double wall_time_remaining = wall_time_per_sim_time * sim_time_remaining;

      int hours = (int)(wall_time_remaining / 3600.0);
      int minutes = (int)((wall_time_remaining - hours*3600.0) / 60.0);
      int seconds = (int)(wall_time_remaining - hours*3600.0 - minutes*60.0);

      gkyl_gyrokinetic_app_cout(app, stdout, "\t(%.1f%% complete, est. %dh %dm %ds remaining)\n", 
        pct_complete, hours, minutes, seconds);
      
      // Reset timing window for next batch
      trigs->tm_loop_start = gkyl_wall_clock();
      trigs->t_loop_start = t_curr;
    }
    else {
      gkyl_gyrokinetic_app_cout(app, stdout, "\t(%.1f%% complete)\n", pct_complete);
    }
  }
}

static void
write_message_pre_update_nonverbose_singleb(const struct gkyl_gyrokinetic_app *app, long step, double t_curr,
  struct message_trigs *trigs)
{
  // No output before update in non-verbose mode
}

static void
write_message_post_update_nonverbose_singleb(const struct gkyl_gyrokinetic_app *app, long step,
  double t_curr, double dt_actual, struct message_trigs *trigs)
{
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig, t_curr)) {
    if (trigs->log_count > 0) {
      gkyl_gyrokinetic_app_cout(app, stdout, "\tStep %ld at time %#11.8g.  Time-step  %.6e.  Completed %g%s\n",
        step, t_curr, dt_actual, trigs->tenth * 10.0, "%");
    }
    else {
      trigs->log_count += 1;
    }
    
    trigs->tenth += 1;
  }
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig_1p, t_curr)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "%d", trigs->p1c);
    trigs->p1c = (trigs->p1c+1) % 10;
    fflush(stdout);
  }
}

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

void
gyrokinetic_run_singleb_simulation(struct gkyl_gyrokinetic_run_inp* inp)
{
  struct gkyl_gyrokinetic_time_stepping_inp time_stepping = inp->time_stepping;
  struct gkyl_gyrokinetic_run_verbosity_inp verbose = inp->print_verbosity;
  struct timespec tm_init = gkyl_wall_clock();

  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&inp->app_inp);

  if (verbose.enabled) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Gyrokinetic simulation initialized...\n");
  }
 
  double t_curr = 0.0, t_end = time_stepping.t_end; // Initial and final simulation times.
  int frame_curr = 0; // Initialize simulation.

  if (time_stepping.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, time_stepping.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
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
  int num_frames = time_stepping.num_frames, num_int_diag_calc = time_stepping.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = 
    { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = 
    { .dt = t_end/(time_stepping.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = 
    { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc), .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics_singleb(&trig_calc_intdiag, app, t_curr, time_stepping.is_restart, false, -1.0);
  write_data_singleb(&trig_write_conf, &trig_write_phase, app, t_curr, time_stepping.is_restart, false);

  if (verbose.enabled) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Initialization completed in %g sec\n\n", gkyl_time_diff_now_sec(tm_init));
  }

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = time_stepping.dt_failure_tol;
  int num_failures = 0, num_failures_max = time_stepping.num_failures_max;

  // Set up function pointers based on verbosity mode
  struct message_trigs m_trig = {
    .log_count = 0,
    .tenth = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 10.0) : 0.0,
    .p1c = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 100.0) % 10 : 0.0,
    .log_trig = { .dt = t_end / 10.0, .tcurr = t_curr },
    .log_trig_1p = { .dt = t_end / 100.0, .tcurr = t_curr },
    .io_period = (long)(1/verbose.frequency),
    .t_end = t_end,
    .estimate_completion_time = verbose.estimate_completion_time,
    .tm_loop_start = gkyl_wall_clock(),
    .t_loop_start = t_curr,
  };

  write_message_pre_update_singleb_t write_message_pre_update;
  write_message_post_update_singleb_t write_message_post_update;

  if (verbose.enabled) {
    write_message_pre_update = write_message_pre_update_verbose_singleb;
    write_message_post_update = write_message_post_update_verbose_singleb;
  }
  else {
    write_message_pre_update = write_message_pre_update_nonverbose_singleb;
    write_message_post_update = write_message_post_update_nonverbose_singleb;
  }

  long step = 1;
  while ((t_curr < t_end) && (step <= time_stepping.num_steps)) {
    write_message_pre_update(app, step, t_curr, &m_trig);

    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    write_message_post_update(app, step, t_curr, status.dt_actual, &m_trig);

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
  gkyl_gyrokinetic_app_release(app);
}

//
// ............. Single block multirate simulations ............... //
//


void reset_io_triggers(struct gkyl_gyrokinetic_time_stepping_inp *time_stepping, struct multirate_time_frame_state *tfs,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag)
{
  // Reset I/O triggers:
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;
  int frame_curr = tfs->frame_curr;
  int num_frames = tfs->num_frames;
  int num_int_diag_calc = time_stepping->int_diag_calc_num;

  // Prevent division by zero when frame_curr equals num_frames
  int frames_remaining = GKYL_MAX2(1, num_frames - frame_curr);
  double time_remaining = GKYL_MAX2(1e-12, t_end - t_curr);

  trig_write_conf->dt = time_remaining / frames_remaining;
  trig_write_conf->tcurr = t_curr;
  trig_write_conf->curr = frame_curr;

  trig_write_phase->dt = time_remaining / (time_stepping->write_phase_freq * frames_remaining);
  trig_write_phase->tcurr = t_curr;
  trig_write_phase->curr = frame_curr;

  int diag_frames = GKYL_MAX2(frames_remaining, (num_int_diag_calc/GKYL_MAX2(1, num_frames)) * frames_remaining);
  trig_calc_intdiag->dt = time_remaining / diag_frames;
  trig_calc_intdiag->tcurr = t_curr;
  trig_calc_intdiag->curr = frame_curr;
}

void run_phase(gkyl_gyrokinetic_app* app, struct gkyl_gyrokinetic_time_stepping_inp *time_stepping, struct gkyl_gyrokinetic_run_verbosity_inp *verbose,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag, struct multirate_time_frame_state *tfs,
  struct gkyl_gyrokinetic_multirate_phase_params *pparams, int phase_idx)
{
  tfs->t_end = tfs->t_curr + pparams->duration;
  tfs->num_frames = tfs->frame_curr + pparams->num_frames;

  // Run an OAP or FDP.
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;
  
  // Reset I/O triggers:
  reset_io_triggers(time_stepping, tfs, trig_write_conf, trig_write_phase, trig_calc_intdiag);

  // Reset simulation parameters and function pointers.
  struct gkyl_gyrokinetic_collisionless collisionless_inp = {
    .type = GKYL_GK_COLLISIONLESS_ES,
    .scale_factor = pparams->alpha,
  };
  struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp = {
    .type = pparams->fdot_mult_type,
    .cellwise_const = true,
    .write_diagnostics = true,
  };
  struct gkyl_gyrokinetic_field field_inp = {
    .is_static = pparams->is_static_field,
  };

  gkyl_gyrokinetic_app_reset_species_collisionless(app, t_curr, "ion", collisionless_inp);
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, t_curr, "ion", fdot_mult_inp);
  gkyl_gyrokinetic_app_reset_species_positivity(app, t_curr, "ion", pparams->is_positivity_enabled);
  gkyl_gyrokinetic_app_reset_field(app, t_curr, field_inp);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = time_stepping->dt_failure_tol;
  int num_failures = 0, num_failures_max = time_stepping->num_failures_max;

  // Set up function pointers based on verbosity mode
  struct message_trigs m_trig = {
    .log_count = 0,
    .tenth = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 10.0) : 0.0,
    .p1c = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 100.0) % 10 : 0.0,
    .log_trig = { .dt = t_end / 10.0, .tcurr = t_curr },
    .log_trig_1p = { .dt = t_end / 100.0, .tcurr = t_curr },
    .io_period = (long)(1/verbose->frequency),
    .t_end = t_end,
    .estimate_completion_time = verbose->estimate_completion_time,
    .tm_loop_start = gkyl_wall_clock(),
    .t_loop_start = t_curr,
  };

  write_message_pre_update_singleb_t write_message_pre_update;
  write_message_post_update_singleb_t write_message_post_update;

  if (verbose->enabled) {
    write_message_pre_update = write_message_pre_update_verbose_singleb;
    write_message_post_update = write_message_post_update_verbose_singleb;
  }
  else {
    write_message_pre_update = write_message_pre_update_nonverbose_singleb;
    write_message_post_update = write_message_post_update_nonverbose_singleb;
  }

  if (verbose->multirate_phase_start_info) {
    gkyl_gyrokinetic_app_cout(app, stdout, "\n========== Phase %d ==========\n", phase_idx);
    gkyl_gyrokinetic_app_cout(app, stdout, "Duration: %g | Frames: %d | α: %g | Field: %s | Positivity: %s | fdot_mult: %d\n", 
      pparams->duration, pparams->num_frames, pparams->alpha,
      pparams->is_static_field ? "static" : "evolve",
      pparams->is_positivity_enabled ? "ON" : "OFF",
      pparams->fdot_mult_type);
    gkyl_gyrokinetic_app_cout(app, stdout, "==============================\n");
  }


  long step = 1;
  while ((t_curr < t_end) && (step <= time_stepping->num_steps))
  {
    write_message_pre_update(app, step, t_curr, &m_trig);

    dt = t_end - t_curr; // Ensure we don't step beyond t_end.
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    write_message_post_update(app, step, t_curr, status.dt_actual, &m_trig);

    if (!status.success){
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics_singleb(trig_calc_intdiag, app, t_curr, false, t_curr >= t_end, status.dt_actual);
    write_data_singleb(trig_write_conf, trig_write_phase, app, t_curr, false, t_curr >= t_end);

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
        calc_integrated_diagnostics_singleb(trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
        write_data_singleb(trig_write_conf, trig_write_phase, app, t_curr, false, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  tfs->t_curr = t_curr;
  tfs->frame_curr = tfs->frame_curr+pparams->num_frames;
}

void
gyrokinetic_run_singleb_multirate_simulation(struct gkyl_gyrokinetic_run_inp* inp)
{
  struct gkyl_gyrokinetic_time_stepping_inp time_stepping = inp->time_stepping;
  struct gkyl_gyrokinetic_run_verbosity_inp verbose = inp->print_verbosity;
  struct gkyl_gyrokinetic_multirate multirate = inp->multirate;
  struct timespec tm_init = gkyl_wall_clock();

  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&inp->app_inp);

  if (verbose.enabled) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Gyrokinetic simulation initialized...\n");
  }

  // Triggers for IO.
  struct gkyl_tm_trigger trig_write_conf, trig_write_phase, trig_calc_intdiag;

  struct multirate_time_frame_state tfs = {
    .t_curr = 0.0, // Initial simulation time.
    .frame_curr = 0, // Initial frame.
    .t_end = multirate.phases[0].duration, // Final time of 1st phase.
    .num_frames = multirate.phases[0].num_frames, // Number of frames in 1st phase.
  };

  int phase_idx_init = 0, phase_idx_end = multirate.num_phases; // Initial and final phase index.
  if (time_stepping.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, time_stepping.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    tfs.frame_curr = status.frame;
    tfs.t_curr = status.stime;

    // Find out what phase we are in.
    double time_count = 0.0;
    int frame_count = 0;
    int pit_curr = 0;
    for (int pit=0; pit<multirate.num_phases; pit++) {
      time_count += multirate.phases[pit].duration;
      frame_count += multirate.phases[pit].num_frames;
      if ((tfs.t_curr <= time_count) && (tfs.frame_curr <= frame_count)) {
        pit_curr = pit;
        break;
      }
    };
    phase_idx_init = pit_curr;

    // Change the duration and number frames so this phase reaches the expected
    // time and number of frames and not beyond.
    struct gkyl_gyrokinetic_multirate_phase_params *pparams = &multirate.phases[phase_idx_init];
    pparams->num_frames = frame_count - tfs.frame_curr;
    pparams->duration = time_count - tfs.t_curr;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", tfs.frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", tfs.t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, tfs.t_curr);

    // Write out ICs.
    reset_io_triggers(&time_stepping, &tfs, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag);

    calc_integrated_diagnostics_singleb(&trig_calc_intdiag, app, tfs.t_curr, false, true, -1.0);
    write_data_singleb(&trig_write_conf, &trig_write_phase, app, tfs.t_curr, false, true);
  }
  
  if (verbose.enabled) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Initialization completed in %g sec\n\n", gkyl_time_diff_now_sec(tm_init));
  }

  if (time_stepping.num_steps != INT_MAX)
    phase_idx_end = 1;

  // Loop over number of number of phases;
  for (int pit=phase_idx_init; pit<phase_idx_end; pit++) {
    struct gkyl_gyrokinetic_multirate_phase_params *phase_params = &multirate.phases[pit];
    run_phase(app, &time_stepping, &verbose, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag, &tfs, phase_params, pit);
  }

  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app); // fetch simulation statistics
  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0)
  {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

  freeresources:
  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
}

//
// ............. Multi-block simulations ............... //
// 

typedef void (*write_message_pre_update_multib_t)(const struct gkyl_gyrokinetic_multib_app*, 
  long, double, struct message_trigs*);
typedef void (*write_message_post_update_multib_t)(const struct gkyl_gyrokinetic_multib_app*, 
  long, double, double, struct message_trigs*);

// Verbose output functions for multi-block
static void
write_message_pre_update_verbose_multib(const struct gkyl_gyrokinetic_multib_app *app, long step,
  double t_curr, struct message_trigs *trigs)
{
  if (((step % trigs->io_period) == 0) || step == 1) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Taking time-step %ld at t = %#11.8g ...", step, t_curr);
  }
}

static void
write_message_post_update_verbose_multib(const struct gkyl_gyrokinetic_multib_app *app, long step,
  double t_curr, double dt_actual, struct message_trigs *trigs)
{
  if (step == 1) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "\tdt = %.6e\n", dt_actual);
  }
  else if ((step % trigs->io_period == 0)) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "\tdt = %.6e ", dt_actual);
    double pct_complete = 100.0 * t_curr / trigs->t_end;
    
    if (trigs->estimate_completion_time) {
      double wall_time_elapsed = gkyl_time_diff_now_sec(trigs->tm_loop_start);
      double sim_time_progressed = t_curr - trigs->t_loop_start;
      double sim_time_remaining = trigs->t_end - t_curr;
      double wall_time_per_sim_time = wall_time_elapsed / sim_time_progressed;
      double wall_time_remaining = wall_time_per_sim_time * sim_time_remaining;

      int hours = (int)(wall_time_remaining / 3600.0);
      int minutes = (int)((wall_time_remaining - hours*3600.0) / 60.0);
      int seconds = (int)(wall_time_remaining - hours*3600.0 - minutes*60.0);

      gkyl_gyrokinetic_multib_app_cout(app, stdout, "\t(%.1f%% complete, est. %dh %dm %ds remaining)\n",
        pct_complete, hours, minutes, seconds);
      
      // Reset timing window for next batch
      trigs->tm_loop_start = gkyl_wall_clock();
      trigs->t_loop_start = t_curr;
    }
    else {
      gkyl_gyrokinetic_multib_app_cout(app, stdout, "\t(%.1f%% complete)\n", pct_complete);
    }
  }
}

static void
write_message_pre_update_nonverbose_multib(const struct gkyl_gyrokinetic_multib_app *app, long step, 
  double t_curr, struct message_trigs *trigs)
{
  // No output before update in non-verbose mode
}

static void
write_message_post_update_nonverbose_multib(const struct gkyl_gyrokinetic_multib_app *app, long step,
  double t_curr, double dt_actual, struct message_trigs *trigs)
{
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig, t_curr)) {
    if (trigs->log_count > 0) {
      gkyl_gyrokinetic_multib_app_cout(app, stdout, "\tStep %ld at time %#11.8g.  Time-step  %.6e.  Completed %g%s\n",
        step, t_curr, dt_actual, trigs->tenth * 10.0, "%");
    }
    else {
      trigs->log_count += 1;
    }
    
    trigs->tenth += 1;
  }
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig_1p, t_curr)) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "%d", trigs->p1c);
    trigs->p1c = (trigs->p1c+1) % 10;
    fflush(stdout);
  }
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

void 
gyrokinetic_run_multib_simulation(struct gkyl_gyrokinetic_run_inp* inp)
{
  struct gkyl_gyrokinetic_time_stepping_inp time_stepping = inp->time_stepping;
  struct gkyl_gyrokinetic_run_verbosity_inp verbose = inp->print_verbosity;
  struct timespec tm_init = gkyl_wall_clock();

  gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(&inp->multib_app_inp);

  if (verbose.enabled) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Gyrokinetic simulation initialized...\n");
  }

  double t_curr = 0.0, t_end = time_stepping.t_end; // Initial and final simulation times.
  int frame_curr = 0; // Initialize simulation.

  if (time_stepping.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_multib_app_read_from_frame(app, time_stepping.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_multib_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
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
  int num_frames = time_stepping.num_frames, num_int_diag_calc = time_stepping.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = 
    { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = 
    { .dt = t_end/(time_stepping.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = 
    { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc), .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics_multib(&trig_calc_intdiag, app, t_curr, time_stepping.is_restart, false, -1.0);
  write_data_multib(&trig_write_conf, &trig_write_phase, app, t_curr, time_stepping.is_restart, false);

  if (verbose.enabled) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Initialization completed in %g sec\n\n",
      gkyl_time_diff_now_sec(tm_init));
  }

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = time_stepping.dt_failure_tol;
  int num_failures = 0, num_failures_max = time_stepping.num_failures_max;

  // Set up function pointers based on verbosity mode
  struct message_trigs m_trig = {
    .log_count = 0,
    .tenth = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 10.0) : 0.0,
    .p1c = t_curr > 0.0 ?  (int) floor(t_curr / t_end * 100.0) % 10 : 0.0,
    .log_trig = { .dt = t_end / 10.0, .tcurr = t_curr },
    .log_trig_1p = { .dt = t_end / 100.0, .tcurr = t_curr },
    .io_period = (long)(1/verbose.frequency),
    .t_end = t_end,
    .estimate_completion_time = verbose.estimate_completion_time,
    .tm_loop_start = gkyl_wall_clock(),
    .t_loop_start = t_curr
  };

  write_message_pre_update_multib_t write_message_pre_update;
  write_message_post_update_multib_t write_message_post_update;
  
  if (verbose.enabled) {
    write_message_pre_update = write_message_pre_update_verbose_multib;
    write_message_post_update = write_message_post_update_verbose_multib;
  }
  else {
    write_message_pre_update = write_message_pre_update_nonverbose_multib;
    write_message_post_update = write_message_post_update_nonverbose_multib;
  }

  long step = 1;
  while ((t_curr < t_end) && (step <= time_stepping.num_steps)) {
    write_message_pre_update(app, step, t_curr, &m_trig);

    struct gkyl_update_status status = gkyl_gyrokinetic_multib_update(app, dt);

    write_message_post_update(app, step, t_curr, status.dt_actual, &m_trig);

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
  gkyl_gyrokinetic_multib_app_release(app);
}

void
gkyl_gyrokinetic_run_simulation(struct gkyl_gyrokinetic_run_inp* inp)
{
  if (inp->print_verbosity.frequency == 0.0)
    inp->print_verbosity.frequency = 0.1; // Default to logging every 10 steps.

  if (inp->app_type == GKYL_GK_SINGLEB) {
    gyrokinetic_run_singleb_simulation(inp);
  }
  else if (inp->app_type == GKYL_GK_MULTIB) {
    gyrokinetic_run_multib_simulation(inp);
  }
  else if (inp->app_type == GKYL_GK_SINGLEB_MULTIRATE) {
    gyrokinetic_run_singleb_multirate_simulation(inp);
  }
  else {
    gkyl_exit("gyrokinetic_run: No valid application input provided to run simulation.");
  }
}