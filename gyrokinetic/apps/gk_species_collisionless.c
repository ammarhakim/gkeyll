#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_collisionless_flux_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  // do nothing.
}

static void
gk_species_collisionless_flux_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  // Compute the surface expansion of the phase space flux
  // Note: Each cell stores the *lower* surface expansions of the 
  // phase space flux, so local_ext range needed to index the output
  // values of flux_surf even though we only loop over local ranges
  // to avoid evaluating quantities such as geometry in ghost cells
  // where they are not defined.
  gkyl_gk_collisionless_flux_surf(gkcls->surf_flux_op, 
    &app->local, &species->local, &app->local_ext, &species->local_ext, 
    species->gyro_phi, fin, gkcls->flux_surf, species->cflrate);
}

static void
gk_species_collisionless_add_em_flux_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  gkyl_gk_collisionless_flux_surf(gkcls->surf_flux_op_add_em, 
    &app->local, &species->local, &app->local_ext, &species->local_ext, 
    species->gyro_apardot, fin, gkcls->flux_surf, species->cflrate);
}

static void
gk_species_collisionless_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_species_collisionless_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkcls->flux_func(app, species, gkcls, fin);

  gkyl_dg_updater_gyrokinetic_advance(gkcls->slvr, &species->local, 
    fin, species->cflrate, rhs);

  app->stat.species_collisionless_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_collisionless_add_em_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkcls->add_em_flux(app, species, gkcls, fin);

  gkyl_dg_updater_gyrokinetic_advance(gkcls->slvr_add_em, &species->local, 
    fin, species->cflrate, rhs);

  app->stat.species_collisionless_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_collisionless_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
}

static void
gk_species_collisionless_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_species_collisionless_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_collisionless *gkcls)
{
  gkcls->collisionless_id = gks->info.collisionless.type;
  gkcls->no_by = gks->info.collisionless.no_by;
  gkcls->is_em = app->field->is_em;
  gkcls->write_diagnostics = gks->info.collisionless.write_diagnostics;

  gkcls->write_diags_func = gk_species_collisionless_write_diags_disabled;
  gkcls->flux_func = gk_species_collisionless_flux_disabled;
  gkcls->rhs_func = gk_species_collisionless_rhs_disabled;
  gkcls->add_em_flux = gk_species_collisionless_flux_disabled;
  gkcls->add_em_rhs = gk_species_collisionless_rhs_disabled;

  if (gkcls->collisionless_id) {

    int cdim = app->cdim, vdim = app->vdim;
    int pdim = cdim+vdim;

    // Determine which directions are zero-flux. By default
    // we do not have zero-flux boundary conditions in any direction.
    bool is_zero_flux[2*GKYL_MAX_DIM] = {false};
    for (int dir=0; dir<app->cdim; ++dir) {
      if (gks->lower_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        is_zero_flux[dir] = true;
      if (gks->upper_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        is_zero_flux[dir+pdim] = true;
    }

    // Need to figure out size of flux_surf by finding size of surface basis set 
    struct gkyl_basis surf_basis;
    struct gkyl_basis surf_vpar_basis;
    // Define surface bases
    if (app->poly_order > 1) {
      gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
    }
    else {
      gkyl_cart_modal_serendip(&surf_vpar_basis, pdim-1, app->poly_order);
      if (vdim>1) {
        gkyl_cart_modal_gkhybrid(&surf_basis, cdim-1, vdim); // p=2 in vparallel
      }
      else {
        gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2); // p=2 in vparallel
      }
    }
    int flux_surf_sz = (cdim)*surf_basis.num_basis + surf_vpar_basis.num_basis;

    // Allocate arrays to store surface phase space flux.
    gkcls->flux_surf = mkarr(app->use_gpu, flux_surf_sz, gks->local_ext.volume);

    if (gkcls->is_em) {
      // Parallel component of magnetic vector potential.
      gkcls->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      gkcls->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);    
    }
    else {
      gkcls->apar    = gkyl_array_acquire(app->field->phi_smooth); // Not used.
      gkcls->apardot = gkyl_array_acquire(app->field->phi_smooth); // Not used.
    }

    enum gkyl_gyrokinetic_bc_type bctype_conf[2*GKYL_MAX_CDIM];
    for (int d=0; d<app->cdim; d++) {
      bctype_conf[d] = gks->lower_bc[d].type;
      bctype_conf[GKYL_MAX_CDIM+d] = gks->upper_bc[d].type;
    }

    gkcls->surf_flux_op = gkyl_gk_collisionless_flux_new(&gks->grid, &app->basis, &gks->basis, 
      gks->info.charge, gks->info.mass, gkcls->no_by, false, app->gk_geom, 
      app->dg_geom, app->gk_dg_geom, gks->vel_map, bctype_conf, app->use_gpu);

    struct gkyl_dg_gyrokinetic_auxfields aux_inp = { .flux_surf = gkcls->flux_surf, 
      .phi = gks->gyro_phi, .apar = gkcls->apar, .apardot = gkcls->apardot };
    // Create solver.
    gkcls->slvr = gkyl_dg_updater_gyrokinetic_new(&gks->grid, &app->basis, &gks->basis, 
      &app->local, &gks->local, is_zero_flux, gks->info.charge, gks->info.mass,
      gks->info.skip_cell_threshold, gkcls->no_by, false, app->gk_geom, gks->vel_map, 
      &aux_inp, app->use_gpu);

    // Methods chosen at runtime.
    gkcls->flux_func = gk_species_collisionless_flux_enabled;
    gkcls->rhs_func = gk_species_collisionless_rhs_enabled;
    if (gkcls->is_em) {
      gkcls->add_em_flux = gk_species_collisionless_add_em_flux_enabled;
      gkcls->add_em_rhs = gk_species_collisionless_add_em_rhs_enabled;
    }
    if (gkcls->write_diagnostics) {
      gkcls->write_diags_func = gk_species_collisionless_write_diags_enabled;
    }
  }
}

void
gk_species_collisionless_flux(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  gkcls->flux_func(app, species, gkcls, fin);
}

void
gk_species_collisionless_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkcls->rhs_func(app, species, gkcls, fin, rhs);
}

void
gk_species_collisionless_add_em_flux(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  gkcls->add_em_flux(app, species, gkcls, fin);
}

void
gk_species_collisionless_add_em_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkcls->add_em_rhs(app, species, gkcls, fin, rhs);
}

void
gk_species_collisionless_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  gkcls->write_diags_func(app, gks, gkcls, tm, frame);
}

void
gk_species_collisionless_release(const struct gkyl_gyrokinetic_app *app, const struct gk_collisionless *gkcls)
{
  if (gkcls->collisionless_id) {

    gkyl_array_release(gkcls->flux_surf);
    gkyl_array_release(gkcls->apar);
    gkyl_array_release(gkcls->apardot);
  
    gkyl_gk_collisionless_flux_release(gkcls->surf_flux_op);
    gkyl_dg_updater_gyrokinetic_release(gkcls->slvr);

    if (gkcls->write_diagnostics) {
    }
  }
}
