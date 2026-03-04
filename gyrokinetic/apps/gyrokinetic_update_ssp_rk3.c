#include <gkyl_gyrokinetic_priv.h>

static void
gyrokinetic_forward_euler(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, struct gkyl_update_status *st)
{
  struct timespec wst_fe = gkyl_wall_clock();
  // Take a forward Euler step with the suggested time-step dt. This may
  // not be the actual time-step taken. However, the function will never
  // take a time-step larger than dt even if it is allowed by
  // stability. The actual time-step and dt_suggested are returned in
  // the status object.
  app->stat.nfeuler += 1;

  // Compute the time rate of change of the distributions, df/dt.
  gyrokinetic_rhs(app, tcurr, dt, fdot_args, st);

  // Complete update of distribution functions.
  struct timespec wst = gkyl_wall_clock();
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;

  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  double dta = st->dt_actual;

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_step_f(gks, fout[i], dta, fin[i]);
    gk_species_bflux_step_f(app, &gks->bflux, bflux_out[i], dta, bflux_in[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_step_f(gkns, fout_neut[i], dta, fin_neut[i]);
    gk_neut_species_bflux_step_f(app, &gkns->bflux, bflux_out_neut[i], dta, bflux_in_neut[i]);
  }
  app->stat.fwd_euler_step_f_tm += gkyl_time_diff_now_sec(wst);

  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst_fe);
}

static void
gyrokinetic_bflux_scale(gkyl_gyrokinetic_app* app, struct gkyl_array ***bflux, struct gkyl_array ***bflux_neut, double fac)
{
  // Scale boundary fluxes by a factor.
  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_scale(app, &gks->bflux, bflux[i], fac);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_bflux_scale(app, &gkns->bflux, bflux_neut[i], fac);
  }
}

void
gyrokinetic_pre_process_step_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of a SSP-RK step.
  // NOTE: do not use dt here because its size is not yet determined.
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Adapt sources.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_source_adapt(app, gks, &gks->src, gks->lte.f_lte, bflux_in, tcurr);
  }

  // Copy boundary fluxes to a buffer in case a later stage fails (used to
  // recompute e.g. Boltzmann elc fields) .
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_copy(app, &gks->bflux, gks->bflux.f_copy, bflux_in[i]);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_species_bflux_copy(app, &gkns->bflux, gkns->bflux.f_copy, bflux_in_neut[i]);
  }
}

static void
gyrokinetic_pre_process_rk_stage_initial_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of the first SSP-RK stage.
  // NOTE: do not use dt here because its size is not yet determined.
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Clear input boundary fluxes so we compute fluxes through the boundary in one step.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_clear(app, &gks->bflux, bflux_in[i], 0.0);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_bflux_clear(app, &gkns->bflux, bflux_in_neut[i], 0.0);
  }
}

void
gyrokinetic_pre_process_rk_stage_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an SSP-RK stage.
  if (stage_idx == 0) {
    // First RK stage.
    gyrokinetic_pre_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);
  }
}

static void
gyrokinetic_post_process_rk_stage_initial_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the first SSP-RK stage.
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Divide new bflux by 1/dt (since forward Euler multiplied by dt).
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, 1.0/dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    // Compute moment of f_old to later compute moment of df/dt.
    // Do it before the fields are updated, but after dt is calculated.
    gk_species_calc_int_mom_dt(app, gks, fin[i], dt, gks->fdot_mom_old);
  }

  // Compute field energy divided by dt for energy balance diagnostics.
  gk_field_calc_energy_dt(app, app->field, dt, app->field->em_energy_red_old);

  // Compute the fields and apply BCs.
  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, bflux_out, fout_neut);

  // Multiply new bflux by dt (so forward Euler multiplied can add dt*bflux).
  wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);
}

static void
gyrokinetic_post_process_rk_stage_intermediate_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of intermediate SSP-RK stages.
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  // Divide new bflux by 1/dt (since forward Euler multiplied by dt).
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, 1.0/dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  // Update fields and apply BCs.
  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, bflux_out, fout_neut);

  // Multiply new bflux by dt (so forward Euler multiplied can add dt*bflux).
  wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);
}

static void
gyrokinetic_post_process_rk_stage_final_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the last SSP-RK stage.
}

void
gyrokinetic_post_process_rk_stage_ssprk(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an SSP-RK stage.
  if (stage_idx == 0) {
    // First RK stage.
    gyrokinetic_post_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);
  }
  else if (stage_idx == num_stages-1) {
    // Last RK stage.
//    printf("      post_process_rk_stage_final dt=%.7e\n",dt);
    gyrokinetic_post_process_rk_stage_final_ssprk(app, tcurr, dt, fdot_args);
  }
  else {
    // Other RK stages.
    gyrokinetic_post_process_rk_stage_intermediate_ssprk(app, tcurr, dt, fdot_args);
  }
}

void
gyrokinetic_post_process_step_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the SSP-RK step.
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Divide new bflux by 1/dt (since forward Euler multiplied by dt).
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_bflux_scale(app, bflux_out, bflux_out_neut, 1.0/dt);
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  // Compute volume- and time-integrated boundary fluxes.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_calc_voltime_integrated_mom(app, gks, &gks->bflux, tcurr);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_bflux_calc_voltime_integrated_mom(app, gkns, &gkns->bflux, tcurr);
  }

  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    // Scale species according to balance between recycling and reactions.
    gk_neut_species_recycle_react_scale_apply(app, gkns, &gkns->rrs, fout_neut[i], bflux_out);
  }

  // Apply positivity shift if requested.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_positivity_apply(app, gks, &gks->positivity, gks->fnew, fout[i]);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_positivity_apply(app, gkns, &gkns->positivity, gkns->fnew, fout_neut[i]);
  }

  // Enforce quasineutrality of the positivity shifts.
  gyrokinetic_post_positivity_quasineut(app, fout);

  // Compute the fields and apply BCs
  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, bflux_out, fout_neut);

  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    // Compute moment of f_new to compute moment of df/dt.
    // Need to do it after the fields are updated.
    gk_species_calc_int_mom_dt(app, gks, fout[i], dt, gks->fdot_mom_new);
  }

  // Compute field energy divided by dt for energy balance diagnostics.
  gk_field_calc_energy_dt(app, app->field, dt, app->field->em_energy_red_new);
}

void
gyrokinetic_post_process_failed_step_ssprk(gkyl_gyrokinetic_app* app, double tcurr,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations performed after a failed SSP-RK step.
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Copy boundary fluxes from the beginning of the step back into bflux_in.
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_copy(app, &gks->bflux, bflux_in[i], gks->bflux.f_copy);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_species_bflux_copy(app, &gkns->bflux, bflux_in_neut[i], gkns->bflux.f_copy);
  }
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  // Recalculate the field.
  gyrokinetic_calc_field(app, tcurr, fin, bflux_in);
}

struct gkyl_update_status
gyrokinetic_update_ssp_rk3(gkyl_gyrokinetic_app* app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  struct gkyl_gyrokinetic_fdot_args *fdot_args = &app->fdot_args;
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  struct gkyl_update_status st = { .success = true };
 
  int num_rk_stages = 3; // Number of RK stages.
  enum { RK_STAGE_1 = 0, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1; // Time-stepper state.

  double tcurr = app->tcurr, dt = dt0;

  // Pre-process step.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    fin[i] = gks->f;
    bflux_in[i] = gks->bflux.f;
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fin_neut[i] = gkns->f;
    bflux_in_neut[i] = gkns->bflux.f;
  }
  gyrokinetic_pre_process_step_ssprk(app, tcurr, dt, fdot_args);

  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<ns_charged; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f;
          fout[i] = gks->f1;
          bflux_in[i] = gks->bflux.f;
          bflux_out[i] = gks->bflux.f1;
        }
        for (int i=0; i<ns_neut; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f;
          fout_neut[i] = gkns->f1;
          bflux_in_neut[i] = gkns->bflux.f;
          bflux_out_neut[i] = gkns->bflux.f1;
        }

        // Pre-process stage.
        gyrokinetic_pre_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);

        // Step solution.
        gyrokinetic_forward_euler(app, tcurr, dt, fdot_args, &st);
        dt = st.dt_actual;

        // Post process stage.
        gyrokinetic_post_process_rk_stage_initial_ssprk(app, tcurr, dt, fdot_args);

        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<ns_charged; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          bflux_in[i] = gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<ns_neut; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          bflux_in_neut[i] = gkns->bflux.f1;
          bflux_out_neut[i] = gkns->bflux.fnew;
        }

        // Step solution.
        gyrokinetic_forward_euler(app, tcurr+dt, dt, fdot_args, &st);

        if (st.dt_actual < dt) {

          // Post process failed stage.
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            fin[i] = gks->f;
            bflux_in[i] = gks->bflux.f;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            fin_neut[i] = gkns->f;
            bflux_in_neut[i] = gkns->bflux.f;
          }
          gyrokinetic_post_process_failed_step_ssprk(app, tcurr, fdot_args);

          // Collect stats.
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // Restart from stage 1.

        } 
        else {
          struct timespec wst = gkyl_wall_clock();
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            gk_species_combine(gks, gks->f1, 3.0/4.0, gks->f, 1.0/4.0, gks->fnew, &gks->local_ext);
            gk_species_bflux_combine(app, &gks->bflux, gks->bflux.f1, 3.0/4.0, gks->bflux.f, 1.0/4.0, gks->bflux.fnew);
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_combine(gkns, gkns->f1, 3.0/4.0, gkns->f, 1.0/4.0, gkns->fnew, &gkns->local_ext);
            gk_neut_species_bflux_combine(app, &gkns->bflux, gkns->bflux.f1, 3.0/4.0, gkns->bflux.f, 1.0/4.0, gkns->bflux.fnew);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Post process stage.
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            fout[i] = gks->f1;
            bflux_out[i] = gks->bflux.f1;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            fout_neut[i] = gkns->f1;
            bflux_out_neut[i] = gkns->bflux.f1;
          }
          gyrokinetic_post_process_rk_stage_intermediate_ssprk(app, tcurr, dt, fdot_args);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<ns_charged; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          bflux_in[i] = gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<ns_neut; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          bflux_in_neut[i] = gkns->bflux.f1;
          bflux_out_neut[i] = gkns->bflux.fnew;
        }

        // Step solution.
        gyrokinetic_forward_euler(app, tcurr+dt/2, dt, fdot_args, &st);

        if (st.dt_actual < dt) {
          // Post process failed stage.
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            fin[i] = gks->f;
            bflux_in[i] = gks->bflux.f;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            fin_neut[i] = gkns->f;
            bflux_in_neut[i] = gkns->bflux.f;
          }
          gyrokinetic_post_process_failed_step_ssprk(app, tcurr, fdot_args);

          // Collect stats.
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // Restart from stage 1.

          app->stat.nstage_2_fail += 1;
        }
        else {
          struct timespec wst = gkyl_wall_clock();
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            // Step f.
            gk_species_combine(gks, gks->f, 1.0/3.0, gks->f, 2.0/3.0, gks->fnew, &gks->local_ext);
            // Step boundary fluxes.
            gk_species_bflux_combine(app, &gks->bflux, gks->bflux.f, 1.0/3.0, gks->bflux.f, 2.0/3.0, gks->bflux.fnew);
          }

          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            // Step f.
            gk_neut_species_combine(gkns, gkns->f, 1.0/3.0, gkns->f, 2.0/3.0, gkns->fnew, &gkns->local_ext);
            // Step boundary fluxes.
            gk_neut_species_bflux_combine(app, &gkns->bflux, gkns->bflux.f, 1.0/3.0, gkns->bflux.f, 2.0/3.0, gkns->bflux.fnew);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Post process stage.
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &app->species[i];
            fout[i] = gks->f;
            bflux_out[i] = gks->bflux.f;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            fout_neut[i] = gkns->f;
            bflux_out_neut[i] = gkns->bflux.f;
          }
          gyrokinetic_post_process_rk_stage_final_ssprk(app, tcurr, dt, fdot_args);

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  // Post-process step.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    fout[i] = gks->f;
    bflux_out[i] = gks->bflux.f;
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fout_neut[i] = gkns->f;
    bflux_out_neut[i] = gkns->bflux.f;
  }
  gyrokinetic_post_process_step_ssprk(app, tcurr, dt, fdot_args);

  return st;
}

void
gyrokinetic_pre_process_step_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of an STS step.
}

void
gyrokinetic_pre_process_rk_stage_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an RK stage.
}

void
gyrokinetic_post_process_rk_stage_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an STS stage.
}

void
gyrokinetic_post_process_step_sts(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the STS step.
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  // Compute the fields and apply BCs
  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, bflux_out, fout_neut);
}

void
gyrokinetic_post_process_failed_step_sts(gkyl_gyrokinetic_app* app, double tcurr,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations performed after a failed STS step.
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Copy boundary fluxes from the beginning of the step back into bflux_in.
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_copy(app, &gks->bflux, bflux_in[i], gks->bflux.f_copy);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_species_bflux_copy(app, &gkns->bflux, bflux_in_neut[i], gkns->bflux.f_copy);
  }
  app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

  // Recalculate the field.
  gyrokinetic_calc_field(app, tcurr, fin, bflux_in);
}

void
gyrokinetic_pre_process_step_opsplit(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of an operator split step.
}

void
gyrokinetic_post_process_step_opsplit(gkyl_gyrokinetic_app* app, double tcurr,
  double dt, struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the operator split step.
}

