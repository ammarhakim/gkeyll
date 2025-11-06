#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_scaling_cross_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  // Do nothing.
}

static void
gk_species_scaling_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // Do nothing.
}

static void
gk_species_scaling_apply_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  struct gk_species *gks_ref = &app->species[sca->ref_species_idx];

  // Reciprocal of the density of this species.
  struct gkyl_array *m0_inv = gks->lte.moms.marr;
  gk_species_moment_calc(&gks->m0, gks->local, app->local, fin);
  gkyl_dg_inv_op_range(app->basis, 0, m0_inv, 0, gks->m0.marr, &app->local);

  // Density of the reference species.
  gk_species_moment_calc(&gks_ref->m0, gks_ref->local, app->local, gks_ref->f);

  // Ratio of the densities.
  struct gkyl_array *m0ratio = gks->m0.marr;
  gkyl_dg_mul_op_range(app->basis, 0, m0ratio, 0, m0_inv, 0, gks_ref->m0.marr, &app->local);

  // Multiply this species by ratio of densities times fixed_fraction.
  gkyl_array_scale_range(m0ratio, sca->fixed_fraction, &app->local);

  gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, fin, m0ratio, fin, &app->local, &gks->local);
}

static void
gk_species_scaling_apply_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  // Do nothing.
}

static void
gk_species_scaling_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  // React diagnostics usually written from gk_species.
  // In the case of static gk_species, write_diagnostics flag
  // can be used to check reaction rates from gk_species.
  struct timespec wst = gkyl_wall_clock();
//  // Compute reaction rate
//  const struct gkyl_array *fin[app->num_species];
//  const struct gkyl_array *fin_neut[app->num_species];
//  for (int i=0; i<app->num_species; ++i) {
//    fin[i] = app->species[i].f;
//  }
//  for (int i=0; i<app->num_species; ++i) {
//    fin_neut[i] = app->species[i].f;
//  }
//  gk_species_scaling_cross_moms(app, gks, gkr, fin, fin_neut);
  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  
  struct timespec wtm = gkyl_wall_clock();
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gk_species_scaling_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  // Do nothing
}

void 
gk_species_scaling_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_scaling *sca)
{
  struct gkyl_gyrokinetic_scaling_inp *sca_inp = &gks->info.scaling;

  sca->type = 0;
  sca->cross_moms_func = gk_species_scaling_cross_moms_disabled;
  sca->rhs_func = gk_species_scaling_rhs_disabled;
  sca->apply_func = gk_species_scaling_apply_disabled;
  sca->write_func = gk_species_scaling_write_disabled;

  if (sca_inp->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
    sca->type = sca_inp->type;
    sca->write_diagnostics = sca_inp->write_diagnostics;

    sca->fixed_fraction = sca_inp->fixed_fraction;

    sca->apply_func = gk_species_scaling_apply_enabled;
    if (sca->write_diagnostics)
      sca->write_func = gk_species_scaling_write_enabled;
    else
      sca->write_func = gk_species_scaling_write_disabled;
  }
}

void 
gk_species_scaling_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
    struct gkyl_gyrokinetic_scaling_inp *sca_inp = &gks->info.scaling;

    // Fetch index of species for indexing arrays.
    sca->ref_species_idx = gk_find_species_idx(app, sca_inp->ref_species_name);
  }
}

void 
gk_species_scaling_apply_ic_cross(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
  }
}

void
gk_species_scaling_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  sca->cross_moms_func(app, gks, sca, fin, fin_neut);
}

void
gk_species_scaling_rhs(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  sca->rhs_func(app, gks, sca, fin, rhs);
}

void
gk_species_scaling_apply(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  sca->apply_func(app, gks, sca, fin, bflux);
}

void
gk_species_scaling_write(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  sca->write_func(app, gks, sca, ridx, tm, frame);
}

void 
gk_species_scaling_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
  } 
}
