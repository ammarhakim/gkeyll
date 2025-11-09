#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_collisionless_flux_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  // Do nothing.
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

  gkcls->fdot_scaling(app, species, gkcls, rhs, species->cflrate, &species->local);

  app->stat.species_collisionless_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_collisionless_fdot_scaling_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_range *rng)
{
  // Do nothing.
}

static void
gk_species_collisionless_fdot_scaling_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_range *rng)
{
  // First scale both rhs and cflrate by scale_fac
  gkyl_array_scale_range(rhs, gkcls->scale_fac, rng);
  gkyl_array_scale_range(cflrate, gkcls->scale_fac, rng);

  // Limit the timestep based on omega_H or a user specified minimum timestep dt.
  // This procedure is called time dilation. https://arxiv.org/html/2510.09756
  if (gkcls->cfl_dt_min_omegaH || (gkcls->cfl_dt_min_value > 0.0))
  {
    // Compute the maximum omega from either dt_omegaH or cfl_dt_min_value
    // WARNING: dt_omegaH is DBL_MAX for boltzmann and adiabatic fields! This takes infinite time steps
    double omega_max = DBL_MAX; // A ceiling on omega is a floor on dt;
    if (gkcls->cfl_dt_min_omegaH) {
      omega_max = (gks->dt_omegaH > 1e-30) ? 1.0 / gks->dt_omegaH : DBL_MAX;
    }
    if (gkcls->cfl_dt_min_value > 0.0) {
      double omega_from_user = 1.0 / gkcls->cfl_dt_min_value;
      omega_max = fmin(omega_max, omega_from_user);// Take the largest timestep from either condition
    }

    // Compute scale_fac_array = min(1.0, omega_H / omega_cfl)
    gkyl_dg_inv_op(*gkcls->cfl_basis, 0, gkcls->scale_fac_array, 0, cflrate); // 1/omega_cfl
    gkyl_array_scale(gkcls->scale_fac_array, omega_max); // omega_max / omega_cfl
    gkyl_array_min_by_cell(gkcls->scale_fac_array, 1.0); // min(1.0, omega_max / omega_cfl)
    
    // Apply cell-wise scaling to both rhs and cflrate
    gkyl_array_scale_by_cell(rhs, gkcls->scale_fac_array);
    gkyl_array_scale_by_cell(cflrate, gkcls->scale_fac_array);
  }
}


static void
gk_species_collisionless_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  // Do nothing.
}

static void
gk_species_collisionless_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  // Write scale factor array to file.
  char fileNm[256];

  if ((gkcls->cfl_dt_min_omegaH || (gkcls->cfl_dt_min_value > 0.0))) {
    snprintf(fileNm, sizeof fileNm, "%s-%s_collisionless_scale_fac_%d.gkyl", 
      app->name, gks->info.name, frame);

    struct gkyl_msgpack_data *mt = gk_array_meta_new(
      (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = gkcls->cfl_basis->poly_order,
        .basis_type = gkcls->cfl_basis->id,
      }, GKYL_GK_META_NONE, 0
    );

    if (app->use_gpu) {
      struct gkyl_array *scale_fac_host = mkarr(false, 1, gks->local_ext.volume);
      gkyl_array_copy(scale_fac_host, gkcls->scale_fac_array);
      gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, 
        scale_fac_host, fileNm);
      gkyl_array_release(scale_fac_host);
    }
    else {
      gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, 
        gkcls->scale_fac_array, fileNm);
    }

    gk_array_meta_release(mt);
  }

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}


void 
gk_species_collisionless_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_collisionless *gkcls)
{
  gkcls->collisionless_id = gks->info.collisionless.type;
  gkcls->write_diagnostics = gks->info.collisionless.write_diagnostics;

  gkcls->write_diags_func = gk_species_collisionless_write_diags_disabled;
  gkcls->flux_func = gk_species_collisionless_flux_disabled;
  gkcls->rhs_func = gk_species_collisionless_rhs_disabled;

  if (gkcls->collisionless_id) {

    int cdim = app->cdim, vdim = gks->info.vdim;
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

    if (gkcls->collisionless_id == GKYL_GK_COLLISIONLESS_EM_BPERP) {
      // Parallel component of magnetic vector potential.
      gkcls->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      gkcls->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);    
    }
    else {
      gkcls->apar    = gkyl_array_acquire(app->field->phi_smooth); // Not used.
      gkcls->apardot = gkyl_array_acquire(app->field->phi_smooth); // Not used.
    }

    enum gkyl_gyrokinetic_bc_type bctype_conf[2*GKYL_MAX_CDIM];
    for (int d=0; d<cdim; d++) {
      bctype_conf[d] = gks->lower_bc[d].type;
      bctype_conf[GKYL_MAX_CDIM+d] = gks->upper_bc[d].type;
    }

    gkcls->surf_flux_op = gkyl_gk_collisionless_flux_new(&gks->grid, &app->basis, &gks->basis, 
      gks->info.charge, gks->info.mass, gkcls->collisionless_id, app->gk_geom, 
      app->dg_geom, app->gk_dg_geom, gks->vel_map, bctype_conf, app->use_gpu);

    struct gkyl_dg_gyrokinetic_auxfields aux_inp = { .flux_surf = gkcls->flux_surf, 
      .phi = gks->gyro_phi, .apar = gkcls->apar, .apardot = gkcls->apardot };
    // Create solver.
    gkcls->slvr = gkyl_dg_updater_gyrokinetic_new(&gks->grid, &app->basis, &gks->basis, 
      &app->local, &gks->local, is_zero_flux, gks->info.charge, gks->info.mass,
      gks->info.skip_cell_threshold, gkcls->collisionless_id, app->gk_geom, gks->vel_map, 
      &aux_inp, app->use_gpu);

    gkcls->scale_fac = 1.0; // Not used if scale_factor in input file is not given.
    gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_disabled;
    if (1.0e-16 < fabs(gks->info.collisionless.scale_factor) || 
        gks->info.collisionless.cfl_dt_min_omegaH || 
        (gks->info.collisionless.cfl_dt_min_value > 0.0))
    {
      gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_enabled;
      // We want to set scale_fac to 1 if not provided 
      if (1.0e-16 < fabs(gks->info.collisionless.scale_factor)) {
        gkcls->scale_fac = gks->info.collisionless.scale_factor;
      } else {
        gkcls->scale_fac = 1.0;
      }

      // Allocate array to hold cell-wise scale factors for omega_cfl screening.
      // This array will store min(1.0, omega_H/omega_cfl) for each cell.
      gkcls->cfl_dt_min_omegaH = gks->info.collisionless.cfl_dt_min_omegaH;
      gkcls->cfl_dt_min_value = gks->info.collisionless.cfl_dt_min_value;
      gkcls->scale_fac_array = mkarr(app->use_gpu, 1, gks->local_ext.volume);
      gkyl_array_clear(gkcls->scale_fac_array, 1.0); // Initialize to 1.0.

      // Create P0 basis for cflrate calculations. 1D p=0, because each cell
      // just has a single number
      gkcls->cfl_basis = gkyl_cart_modal_serendip_new(1, 0);
    }

    // Other methods chosen at runtime.
    gkcls->flux_func = gk_species_collisionless_flux_enabled;
    gkcls->rhs_func = gk_species_collisionless_rhs_enabled;
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

    // Release scale_fac_array if it was allocated for omega_cfl screening
    if (gkcls->scale_fac_array) {
      gkyl_array_release(gkcls->scale_fac_array);
      gkyl_free(gkcls->cfl_basis);
    }

    if (gkcls->write_diagnostics) {
    }
  }
}

void
gk_species_collisionless_reset(gkyl_gyrokinetic_app* app, double tm, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_gyrokinetic_collisionless gkcls_inp)
{
  // Update input struct
  gks->info.collisionless.scale_factor = gkcls_inp.scale_factor;
  gks->info.collisionless.cfl_dt_min_omegaH = gkcls_inp.cfl_dt_min_omegaH;
  gks->info.collisionless.cfl_dt_min_value = gkcls_inp.cfl_dt_min_value;

  // Update local struct
  gkcls->scale_fac = gkcls_inp.scale_factor;
  gkcls->cfl_dt_min_omegaH = gkcls_inp.cfl_dt_min_omegaH;
  gkcls->cfl_dt_min_value = gkcls_inp.cfl_dt_min_value;

  // Choose appropriate scaling method
  gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_disabled;
  if (1.0e-16 < fabs(gkcls_inp.scale_factor) ||
            gkcls->cfl_dt_min_omegaH || 
            (gkcls->cfl_dt_min_value > 0.0)) {              
    gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_enabled;
    gkyl_array_clear(gkcls->scale_fac_array, 1.0);
  }
}
