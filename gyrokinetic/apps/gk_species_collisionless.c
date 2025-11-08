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
  gkyl_array_scale_range(rhs, gkcls->scale_fac, rng);
  gkyl_array_scale_range(cflrate, gkcls->scale_fac, rng);
}

static void
gk_species_collisionless_fdot_scaling_omegaH(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_range *rng)
{
  // Omega-based CFL screening: for each cell, compute scale_factor based on the
  // ratio of omega_cfl (from previous step) to omega_H.
  //   - If omega_cfl < omega_H: scale_factor = 1.0 (no screening)
  //   - If omega_cfl > omega_H: scale_factor = omega_H/omega_cfl (slow down to omega_H)
  //
  // omega_H is computed from the previous step as: omega_H = cfl_omegaH / dt_omegaH
  // where dt_omegaH is stored in gks->dt_omegaH

  // Can also be expressed as scale_fac = min(1.0, omega_H / omega_CFL)
  
  // Compute omega_H from previous step (global value)
  double omega_H = (gks->dt_omegaH > 1e-30) ?
    1.0 / gks->dt_omegaH : DBL_MAX;
  
  // This loop makes it clear that this should be inside an updater
  // or use well known array methods
  // The array methods do not have dividing operations, for good reason, 
  
  // To do this with array methods, first compute omega_H / omega_cell
  // We have gkyl_dg_div_op, which divides 2 arrays, but we would need to make a basis
  // and allocate dt_omegaH to an array.
  // We would need to implement a min operation in array_ops, called gkyl_array_min which
  // computes min( a, arr ). I'm not sure how to do this in a general way > p0

  // I think really this is a dg_div_op, then we need to implement some kind of dg_min_op, which takes the minimum of element with a dg array. We could just implement this dg_bin_op_min for p0 only and make the kernels simple. Put a warning in there that it's only implemtned for p0. I'm not sure that this operation even makes sense in p>0 because of weak equality. Do you take the min operation at quadrature nodes? I'm really terrible at maxima kernel generation, so it would take me a long time to implement this p>0. I think I could hardcode p0 without kernels, but it's not really within the framework of dg_bin_ops.

  // I'm also thinking that this logic with storing the omegacfl from previous step will lead to issues with resetting. It's kind of recursive how we have it set up. Instead, I think we should floor the CFL for just the collisionless module explicitly. Since this is the first one called in the gyrokinetic_rhs, we can floor CFL right here, just based on omegaH and input cflrate. We should also make a user option to specity GKYL_GK_COLLISIONLESS_ES_CFL_DT_FLOOR, which floors CFL timestep in general. We could use a boolean whether we want the floor based on omegaH or some user input. Maybe .omegaH_floor = true, or .floor_val = 8 to limit dt to 8 s here. I think when this option is specified without any parameters, it shouldn't do anything. You need to specity one of these parameters, at leasst.

  // I don't think I should use an enum to control this behavior. It doesn't have to do with the equation solver. Maybe it's just a boolean flag like. .enable_cfl_dt_floor = true. However, this boolean doesn't make sense if there is no follow up. So maybe the logic should just be controlled by .cfl_dt_floor_min_omegaH = true, or .cfl_dt_floor = 8. The code should default these options to off, but if something is specified, then use it.

  gkyl_dg_inv_op(cfl_basis, 0, gkcls->inv_cflrate_array, 0, cflrate);
  gkyl_array_scale(gkcls->inv_cflrate_array, omega_H);
  gkyl_array_min(gkcls->scale_fac_array, 1.0);
  
  // Apply cell-wise scaling to both rhs and cflrate
  gkyl_array_scale_by_cell(rhs, gkcls->scale_fac_array);
  gkyl_array_scale_by_cell(cflrate, gkcls->scale_fac_array);

  // Finish by scaling the entire rhs and cflrate by the scale_fac
  gk_species_collisionless_fdot_scaling_enabled(app, gks, gkcls, rhs, cflrate, rng);
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
    for (int d=0; d<app->cdim; d++) {
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
    if (1.0e-16 < fabs(gks->info.collisionless.scale_factor)) {
      gkcls->scale_fac = gks->info.collisionless.scale_factor;
      gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_enabled;
    }

    if (gkcls->collisionless_id == GKYL_GK_COLLISIONLESS_ES_OMEGA_CFL_SCREENING) {
      // Allocate array to hold cell-wise scale factors for omega_cfl screening.
      // This array will store min(1.0, omega_H/omega_cfl) for each cell.
      gkcls->scale_fac_array = mkarr(app->use_gpu, 1, gks->local_ext.volume);
      gkyl_array_clear(gkcls->scale_fac_array, 1.0); // Initialize to 1.0.
      gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_omegaH;
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
    }

    if (gkcls->write_diagnostics) {
    }
  }
}

void
gk_species_collisionless_reset(gkyl_gyrokinetic_app* app, double tm, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_gyrokinetic_collisionless gkcls_inp)
{
  gkcls->scale_fac = 0.0;
  gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_disabled;
  if (1.0e-16 < fabs(gkcls_inp.scale_factor)) {
    gks->info.collisionless.scale_factor = gkcls_inp.scale_factor;

    gkcls->scale_fac = gkcls_inp.scale_factor;
    gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_enabled;
  }
}
