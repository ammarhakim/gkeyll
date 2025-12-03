#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <string.h>

void
gk_field_energy_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Allocate energy reduction arrays.
  if (app->use_gpu) {
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    f->em_energy_red = gkyl_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_malloc(sizeof(double[1]));
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_write_call = true;

  f->calc_energy_func = gk_field_calc_energy_enabled;
  f->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  
  if (f->info.time_rate_diagnostics) {
    gk_field_time_rate_diags_new(app, f);
  }


  // Factors for ES energy.
  f->es_energy_fac = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  f->es_energy_fac_1d = 0.0;
}

void
gk_field_time_rate_diags_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->calc_energy_dt_func = gk_field_calc_energy_dt_active;
  if (app->use_gpu)
  {
    f->em_energy_red_new = gkyl_cu_malloc(sizeof(double[1]));
    f->em_energy_red_old = gkyl_cu_malloc(sizeof(double[1]));
    gkyl_cu_memset(f->em_energy_red_new, 0, sizeof(double[1]));
    gkyl_cu_memset(f->em_energy_red_old, 0, sizeof(double[1]));
  }
  else
  {
    f->em_energy_red_new = gkyl_malloc(sizeof(double[1]));
    f->em_energy_red_old = gkyl_malloc(sizeof(double[1]));
    memset(f->em_energy_red_new, 0, sizeof(double[1]));
    memset(f->em_energy_red_old, 0, sizeof(double[1]));
  }
  f->integ_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_dot_write_call = true;
}

void
gk_field_calc_energy_dt_active(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
    1.0/dt, field->es_energy_fac, &app->local, &app->local, energy_reduced);
  app->stat.phidot_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy_dt_none(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
}

void
gk_field_calc_energy_enabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm)
{
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth,
                               1.0, field->es_energy_fac, &app->local, &app->local, field->em_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

  double energy_global[1] = {0.0};
  if (app->use_gpu)
    gkyl_cu_memcpy(energy_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  else
    energy_global[0] = field->em_energy_red_global[0];

  if (app->cdim == 1)
    energy_global[0] *= field->es_energy_fac_1d;

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  if (field->info.time_rate_diagnostics)
  {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_old, field->em_energy_red_global);
    double energy_dot_global_old[1] = {0.0};
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_old, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_old[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_old[0] *= field->es_energy_fac_1d;

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_new, field->em_energy_red_global);
    double energy_dot_global_new[1] = {0.0};
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_new, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_new[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_new[0] *= field->es_energy_fac_1d;

    double energy_dot_global[1] = {0.0};
    energy_dot_global[0] = energy_dot_global_new[0] - energy_dot_global_old[0];

    gkyl_dynvec_append(field->integ_energy_dot, tm, energy_dot_global);
  }
}

void
gk_field_calc_energy_disabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm)
{
  // Do nothing.
}

void
gk_field_energy_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_dynvec_release(f->integ_energy);

  if (app->use_gpu) {
    gkyl_cu_free(f->em_energy_red);
    gkyl_cu_free(f->em_energy_red_global);
  } else {
    gkyl_free(f->em_energy_red);
    gkyl_free(f->em_energy_red_global);
  }

  if (f->info.time_rate_diagnostics) {
    gk_field_time_rate_diags_release(app, f);
  }

  gkyl_array_release(f->es_energy_fac);
}

void
gk_field_time_rate_diags_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  if (app->use_gpu)
  {
    gkyl_cu_free(f->em_energy_red_new);
    gkyl_cu_free(f->em_energy_red_old);
  }
  else
  {
    gkyl_free(f->em_energy_red_new);
    gkyl_free(f->em_energy_red_old);
  }
  gkyl_dynvec_release(f->integ_energy_dot);
}
