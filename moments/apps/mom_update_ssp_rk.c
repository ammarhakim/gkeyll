#include <gkyl_moment_priv.h>

// Remap a vacuum-Einstein state vector into the metric slots of a GR-Euler fluid state vector. 
static void
sync_fluid_metric_from_einstein(const double *qe, double *qf, double excision_threshold)
{
  double lapse = qe[9];

  qf[5] = lapse;                                       // lapse
  qf[6] = qe[52]; qf[7] = qe[53]; qf[8] = qe[54];      // shift^i
  for (int n = 0; n < 9; n++) qf[9 + n] = qe[0 + n];   // gamma_ij
  for (int n = 0; n < 9; n++) qf[18 + n] = qe[10 + n]; // K_ij
  qf[27] = (lapse < excision_threshold) ? -1.0 : 1.0;  // excision flag (lapse-based)
  for (int n = 0; n < 3; n++) qf[28 + n] = lapse * qe[46 + n];  // d_i alpha = alpha * A_i
  for (int n = 0; n < 9; n++) qf[31 + n] = 2.0 * qe[55 + n];    // d_i beta^j = 2 * (1/2 d_i beta^j)
  for (int n = 0; n < 27; n++) qf[40 + n] = 2.0 * qe[19 + n];   // d_k gamma_ij = 2 * D_kij
}

// Wire up the GR-Euler <-> vacuum-Einstein coupling for the current RK stage
static void
couple_gr_fluid_spacetime(gkyl_moment_app* app, const struct gkyl_array *fin[])
{
  int fluid_idx = -1, einstein_idx = -1;
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].has_gr_euler) fluid_idx = i;
    if (app->species[i].has_vacuum_einstein) einstein_idx = i;
    app->species[i].coupling_partner_fin = NULL;
  }

  if (fluid_idx < 0 || einstein_idx < 0) {
    return; // uncoupled: leave everything untouched.
  }

  app->species[fluid_idx].coupling_partner_fin = fin[einstein_idx];
  app->species[einstein_idx].coupling_partner_fin = fin[fluid_idx];

  double excision_threshold = app->species[einstein_idx].vacuum_einstein_excision_threshold;
  const struct gkyl_array *einstein_f = fin[einstein_idx];
  struct gkyl_array *fluid_f = (struct gkyl_array *) fin[fluid_idx]; // metric slots are mutable scratch.

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qe = gkyl_array_cfetch(einstein_f, loc);
    double *qf = gkyl_array_fetch(fluid_f, loc);
    sync_fluid_metric_from_einstein(qe, qf, excision_threshold);
  }
}

static void
apply_vacuum_einstein_excision_in_step(gkyl_moment_app* app, const struct moment_species *sp,
  struct gkyl_array *f)
{
  if (!sp->has_vacuum_einstein && !sp->has_vacuum_einstein_conformal) {
    return;
  }

  double excision_threshold = 0.0;
  if (sp->has_vacuum_einstein) {
    excision_threshold = sp->vacuum_einstein_excision_threshold;
  }
  else {
    excision_threshold = sp->vacuum_einstein_conformal_excision_threshold;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double *q = gkyl_array_fetch(f, loc);

    if (q[9] < excision_threshold) {
      for (int m = 0; m < sp->num_equations; m++) {
        q[m] = 0.0;
      }
    }
  }
}

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(gkyl_moment_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *emout,
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  double dtmin = DBL_MAX;

  // Couple GR fluid and dynamical spacetime for this stage: sync the live metric into the fluid state
  // and expose partner stage states (no-op for uncoupled runs).
  couple_gr_fluid_spacetime(app, fin);

  // compute RHS of fluid equations
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = moment_species_rhs(app, &app->species[i], fin[i], fout[i]);
    dtmin = fmin(dtmin, dt1);
  }
  // compute RHS of Maxwell equations
  if (app->has_field) {
    double dt1 = moment_field_rhs(app, &app->field, emin, emout);
    dtmin = fmin(dtmin, dt1);
  }

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of species
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate_range(gkyl_array_scale_range(fout[i], dta, &(app->local)),
      1.0, fin[i], &(app->local));
    apply_vacuum_einstein_excision_in_step(app, &app->species[i], fout[i]);
    moment_species_apply_bc(app, tcurr, &app->species[i], fout[i]);
  }
  if (app->has_field) {
    // complete update of field (even when field is static, it is
    // safest to do this accumulate as it ensure emout = emin)
    gkyl_array_accumulate_range(gkyl_array_scale_range(emout, dta, &(app->local)),
      1.0, emin, &(app->local));
    
    moment_field_apply_bc(app, tcurr, &app->field, emout);
  }
}

// internal function that takes a single time-step using a single-step
// Strang-split scheme
struct gkyl_update_status
moment_update_ssp_rk3(gkyl_moment_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f0;
          fout[i] = app->species[i].f1;
        }
        forward_euler(app, tcurr, dt, fin, app->has_field ? app->field.f0 : 0,
          fout, app->has_field ? app->field.f1 : 0,
          &st);
	dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        forward_euler(app, tcurr+dt, dt, fin, app->has_field ? app->field.f1 : 0,
          fout, app->has_field ? app->field.fnew : 0,
          &st);
	if (st.dt_actual < dt) {

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;
            
          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        } else {
          for (int i=0; i<app->num_species; ++i)
            array_combine(app->species[i].f1,
              3.0/4.0, app->species[i].f0, 1.0/4.0, app->species[i].fnew, &app->local_ext);
          if (app->has_field)
            array_combine(app->field.f1,
              3.0/4.0, app->field.f0, 1.0/4.0, app->field.fnew, &app->local_ext);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        forward_euler(app, tcurr+dt/2, dt, fin, app->has_field ? app->field.f1 : 0,
          fout, app->has_field ? app->field.fnew : 0,
          &st);
        if (st.dt_actual < dt) {
          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;
            
          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1
            
          app->stat.nstage_2_fail += 1;
        }
        else {
          for (int i=0; i<app->num_species; ++i) {
            array_combine(app->species[i].f1,
              1.0/3.0, app->species[i].f0, 2.0/3.0, app->species[i].fnew, &app->local_ext);
            gkyl_array_copy_range(app->species[i].f0, app->species[i].f1, &app->local_ext);
          }
          if (app->has_field) {
            array_combine(app->field.f1,
              1.0/3.0, app->field.f0, 2.0/3.0, app->field.fnew, &app->local_ext);
            gkyl_array_copy_range(app->field.f0, app->field.f1, &app->local_ext);
          }
          
          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }
  
  return st;
}
