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
    gk_species_bflux_accumulate(app, &gks->bflux, bflux_out[i], 1.0, bflux_in[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_step_f(gkns, fout_neut[i], dta, fin_neut[i]);
    gk_neut_species_bflux_accumulate(app, &gkns->bflux, bflux_out_neut[i], 1.0, bflux_in_neut[i]);
  }
  app->stat.fwd_euler_step_f_tm += gkyl_time_diff_now_sec(wst);

  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst_fe);
}

static void
gyrokinetic_pre_process_rk_stage_initial(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the beginning of the first RK stage.
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;

  // Adapt sources.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_source_adapt(app, gks, &gks->src, gks->lte.f_lte, bflux_in, tcurr);
  }
}

void
gyrokinetic_pre_process_rk_stage(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an RK stage.
  if (stage_idx == 0) {
    // First RK stage.
    gyrokinetic_pre_process_rk_stage_initial(app, tcurr, dt, fdot_args);
  }
}

static void
gyrokinetic_post_process_rk_stage_initial(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the first RK stage.
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Subtract boundary flux f from f1 so that we only step boundary
  // fluxes during a given time step, not over all time. And so that the
  // boundary flux in f is kept in case a later RK stage fails.
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bflux_accumulate(app, &gks->bflux, bflux_out[i], -1.0, bflux_in[i]);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_bflux_accumulate(app, &gkns->bflux, bflux_out_neut[i], -1.0, bflux_in_neut[i]);
  }
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
}

static void
gyrokinetic_post_process_rk_stage_intermediate(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of intermediate RK stages.
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;

  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, bflux_out, fout_neut);
}

static void
gyrokinetic_post_process_rk_stage_final(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args)
{
  // Operations at the end of the last RK stage.
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;

  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;

  // Apply positivity shift if requested.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_positivity_apply(app, gks, &gks->positivity, gks->fnew, fout[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_positivity_apply(app, gkns, &gkns->positivity, gkns->fnew, fout_neut[i]);
  }

  // Enforce quasineutrality of the positivity shifts.
  gyrokinetic_post_positivity_quasineut(app, fout);

  // Compute the fields and apply BCs
  gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, bflux_out, fout_neut);

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    // Compute moment of f_new to compute moment of df/dt.
    // Need to do it after the fields are updated.
    gk_species_calc_int_mom_dt(app, gks, fout[i], dt, gks->fdot_mom_new);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    // Scale species according to balance between recycling and reactions.
    gk_neut_species_recycle_react_scale_apply(app, gkns, &gkns->rrs, fout_neut[i], bflux_out);
  }

  // Compute field energy divided by dt for energy balance diagnostics.
  gk_field_calc_energy_dt(app, app->field, dt, app->field->em_energy_red_new);
}

void
gyrokinetic_post_process_rk_stage(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations at the end of an RK stage.
  if (stage_idx == 0) {
    // First RK stage.
    gyrokinetic_post_process_rk_stage_initial(app, tcurr, dt, fdot_args);
  }
  else if (stage_idx < num_stages-1) {
    // Last RK stage.
    gyrokinetic_post_process_rk_stage_final(app, tcurr, dt, fdot_args);
  }
  else {
    // Other RK stages.
    gyrokinetic_post_process_rk_stage_intermediate(app, tcurr, dt, fdot_args);
  }
}

void
gyrokinetic_post_process_failed_rk_stage(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_gyrokinetic_fdot_args *fdot_args, int stage_idx, int num_stages)
{
  // Operations performed after a failed stage.
  const struct gkyl_array **fin = fdot_args->fin;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;

  // Recalculate the field.
  gyrokinetic_calc_field(app, tcurr, fin, bflux_in);
}

struct gkyl_update_status
gyrokinetic_update_ssp_rk3(gkyl_gyrokinetic_app* app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.

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
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f;
          fout[i] = gks->f1;
          bflux_in[i] = gks->bflux.f;
          bflux_out[i] = gks->bflux.f1;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f;
          fout_neut[i] = gkns->f1;
          bflux_in_neut[i] = gkns->bflux.f;
          bflux_out_neut[i] = gkns->bflux.f1;
        }

        // Pre-process stage.
        gyrokinetic_pre_process_rk_stage_initial(app, tcurr, dt, fdot_args);

        // Step solution.
        gyrokinetic_forward_euler(app, tcurr, dt, fdot_args, &st);
        dt = st.dt_actual;

        // Post process stage.
        gyrokinetic_post_process_rk_stage_initial(app, tcurr, dt, fdot_args);

        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          bflux_in[i] = gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
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
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fin[i] = gks->f;
            bflux_in[i] = gks->bflux.f;
          }
          gyrokinetic_post_process_failed_rk_stage(app, tcurr, dt, fdot_args, RK_STAGE_2, num_rk_stages);

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
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            gk_species_combine(gks, gks->f1, 3.0/4.0, gks->f, 1.0/4.0, gks->fnew, &gks->local_ext);
            gk_species_bflux_set(app, &gks->bflux, gks->bflux.f1, 1.0/4.0, gks->bflux.fnew);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_combine(gkns, gkns->f1, 3.0/4.0, gkns->f, 1.0/4.0, gkns->fnew, &gkns->local_ext);
            gk_neut_species_bflux_set(app, &gkns->bflux, gkns->bflux.f1, 1.0/4.0, gkns->bflux.fnew);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Post process stage.
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fout[i] = gks->f1;
            bflux_out[i] = gks->bflux.f1;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f1;
          }
          gyrokinetic_post_process_rk_stage_intermediate(app, tcurr, dt, fdot_args);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          bflux_in[i] = gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
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
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fin[i] = gks->f;
            bflux_in[i] = gks->bflux.f;
          }
          gyrokinetic_post_process_failed_rk_stage(app, tcurr, dt, fdot_args, RK_STAGE_3, num_rk_stages);

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
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            // Step f.
            gk_species_combine(gks, gks->f1, 1.0/3.0, gks->f, 2.0/3.0, gks->fnew, &gks->local_ext);
            gk_species_copy_range(gks, gks->f, gks->f1, &gks->local_ext);
            // Step boundary fluxes.
            gk_species_bflux_set(app, &gks->bflux, gks->bflux.f, 2.0/3.0, gks->bflux.fnew);
            gk_species_bflux_calc_voltime_integrated_mom(app, gks, &gks->bflux, tcurr);
          }

          for (int i=0; i<app->num_neut_species; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            // Step f.
            gk_neut_species_combine(gkns, gkns->f1, 1.0/3.0, gkns->f, 2.0/3.0, gkns->fnew, &gkns->local_ext);
            gk_neut_species_copy_range(gkns, gkns->f, gkns->f1, &gkns->local_ext);
            // Step boundary fluxes.
            gk_neut_species_bflux_set(app, &gkns->bflux, gkns->bflux.f, 2.0/3.0, gkns->bflux.fnew);
            gk_neut_species_bflux_calc_voltime_integrated_mom(app, gkns, &gkns->bflux, tcurr);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Post process stage.
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            fout[i] = gks->f;
            bflux_out[i] = gks->bflux.f;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f;
          }
          gyrokinetic_post_process_rk_stage_final(app, tcurr, dt, fdot_args);

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}
