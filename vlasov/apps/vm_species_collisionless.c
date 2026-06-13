#include <assert.h>
#include <gkyl_vlasov_priv.h>

static void
vm_species_collisionless_rhs_enabled(gkyl_vlasov_app *app, struct vm_species *vms,
  struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Set values of q/m*EM and the total potentials based on field and external forces configuration. 
  gkyl_array_clear(cls->qmem, 0.0); 
  if (cls->has_app_accel) {
    gkyl_array_accumulate_range(cls->qmem, 1.0, cls->app_accel, &app->local);
  }

  if (app->has_field) {
    if (app->field->has_ext_em) {
      gkyl_array_accumulate_range(cls->qmem, cls->qbym, app->field->ext_em, &app->local);
    }

    if (vms->field_id == GKYL_FIELD_E_B || vms->field_id == GKYL_FIELD_GR_D_B) {
      gkyl_array_accumulate_range(cls->qmem, cls->qbym, em, &app->local);
    }
    else if (vms->field_id == GKYL_FIELD_PHI) {
      gkyl_array_set_offset(cls->pot_tot, cls->qbym, app->field->phi, 0);
      if (app->field->has_ext_pot) {
        gkyl_array_accumulate_offset(cls->pot_tot, cls->qbym, app->field->ext_pot, 0);
      }
    }
  }

  // Divide out velocity-space Jacobian.
  gkyl_vlasov_velocity_map_divide_jacobvel(vms->vel_map, &app->basis, &vms->basis,
    &vms->local, fin, vms->f_no_J);

  // Compute the surface expansion of the phase space flux in configuration space. 
  if (vms->model_id == GKYL_MODEL_TRIAD || vms->model_id == GKYL_MODEL_TRIAD_GR) {
    gkyl_dg_vlasov_conf_flux_surf_advance(cls->calc_conf_flux, &app->local, &vms->local, &vms->local_ext, 
      vms->conf_poisson_tensor, vms->hamil, fin, vms->cflrate, cls->conf_flux_surf);
  }

  // Compute the surface expansion of the phase space flux in velocity space.
  gkyl_dg_vlasov_vel_flux_surf_advance(cls->calc_vel_flux, &app->local, &vms->local,
    vms->conf_poisson_tensor, vms->hamil, cls->qmem, cls->pot_tot, vms->rad,
    vms->f_no_J, vms->cflrate, cls->vel_flux_surf);

  gkyl_hyper_dg_advance(cls->slvr, &vms->local, fin, vms->cflrate, rhs);

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_collisionless_rhs_disabled(gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs)
{
  // Do nothing.
}

void 
vm_species_collisionless_init(struct gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;
  enum gkyl_basis_type b_type = app->basis.b_type;

  // Allocate array to store q/m*(E,B) or potentials (q/m*phi + m*phi_g, q/m*A) depending on equation system. 
  // Note: the potentials are the total potentials and thus can include both (or either) gravitational
  // or electrostatic interactions. 
  cls->qbym = vms->info.charge/vms->info.mass;
  cls->qmem = mkarr(app->use_gpu, 8*app->basis.num_basis, app->local_ext.volume);
  cls->pot_tot = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);

  // Initialize applied acceleration for use in force update. 
  cls->app_accel = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(cls->app_accel, 0.0);
  cls->has_app_accel = false;
  cls->app_accel_evolve = false;
  // setup applied acceleration
  if (vms->info.app_accel) {
    cls->has_app_accel = true;
    if (vms->info.app_accel_evolve) {
      cls->app_accel_evolve = vms->info.app_accel_evolve;
    }

    // Host-side distribution function for projection on GPUs.
    cls->app_accel_host = app->use_gpu ? mkarr(false, cls->app_accel->ncomp, cls->app_accel->size)
                                       : gkyl_array_acquire(cls->app_accel);      
    cls->app_accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->basis, app->basis.poly_order+1,
      3, vms->info.app_accel, vms->info.app_accel_ctx);
  }

  // Determine which forces we need based on combination of field ID and presence 
  // of applied accelerations and external fields/potentials. 
  cls->has_E = false; 
  cls->has_B = false; 
  cls->has_phi = false; 
  if (app->has_field) {
    if (vms->field_id == GKYL_FIELD_E_B || app->field->has_ext_em || vms->field_id == GKYL_FIELD_GR_D_B) {
      cls->has_E = true; 
      cls->has_B = true; 
    } 
    if (vms->field_id == GKYL_FIELD_PHI) {
      cls->has_phi = true; 
    }
  }
  else if (cls->has_app_accel) {
    cls->has_E = true; 
  }
  cls->use_lo = false; 
  if (vms->info.use_lo == true) {
    cls->use_lo = true; 
  }
  cls->use_vierbein = false; 
  if (vms->info.use_vierbein == true) {
    cls->use_vierbein = true; 
  }
  cls->use_extended_hamil_def = false; 
  if (vms->info.use_extended_hamil_def == true) {
    cls->use_extended_hamil_def = true; 
  }
  cls->use_preset_geom = false;
  if (vms->geom && vms->geom->use_preset_geom) {
    cls->use_preset_geom = true;
  }

  // Select the number of nodes, with case for hybrid-tensor.
  int highorder = cls->use_lo ? 0 : 1;

  // p + 1 is equivalent to p + 2 for ser p1
  if ((app->poly_order == 1) && (b_type == GKYL_BASIS_MODAL_SERENDIPITY)) {
    highorder = 0;
  }
  cls->num_surf_vel_nodes = pow(app->poly_order+1+highorder, pdim - 1);
  if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (app->poly_order == 1)) {
    cls->num_surf_vel_nodes = (int) ( pow(app->poly_order+1+highorder,vdim - 1) + pow(app->poly_order + 1,cdim) );
  }

  // Allocate nodal surface expansion of velocity space flux array (conf). 
  if (vms->model_id == GKYL_MODEL_TRIAD || vms->model_id == GKYL_MODEL_TRIAD_GR) {

    // Compute the number of configuration space nodes, with case for hybrid-tensor.
    cls->num_surf_conf_nodes = pow(app->poly_order+1+highorder,pdim - 1);
    if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (app->poly_order == 1)) {
      cls->num_surf_conf_nodes = (int) ( pow(app->poly_order+1+highorder,vdim) + pow(app->poly_order + 1,cdim- 1) );
    }

    cls->conf_flux_surf = mkarr(app->use_gpu, cdim*cls->num_surf_conf_nodes, vms->local_ext.volume);
    struct gkyl_dg_vlasov_conf_flux_surf_inp inp_conf_flux = {
      .phase_grid = &vms->grid, 
      .conf_basis = &app->basis,
      .phase_basis = &vms->basis,
      .vel_range = &vms->local_vel,
      .hamil_range = &vms->hamil_range,
      .skip_cell_thresh = vms->info.skip_cell_thresh > 0.0 ? vms->info.skip_cell_thresh : 0.0, 
      .model_id = vms->model_id,
      .use_lo = cls->use_lo,
      .use_gpu = app->use_gpu,
    }; 
    cls->calc_conf_flux = gkyl_dg_vlasov_conf_flux_surf_inew(&inp_conf_flux); 
  }

  // Allocate nodal surface expansion of velocity space flux array (vel).
  cls->vel_flux_surf = mkarr(app->use_gpu, vdim*cls->num_surf_vel_nodes, vms->local_ext.volume);
  struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
    .phase_grid = &vms->grid,
    .conf_basis = &app->basis,
    .phase_basis = &vms->basis,
    .vel_map = vms->vel_map,
    .hamil_range = &vms->hamil_range,
    .skip_cell_thresh = vms->info.skip_cell_thresh > 0.0 ? vms->info.skip_cell_thresh : 0.0, 
    .model_id = vms->model_id,
    .has_E = cls->has_E, 
    .has_phi = cls->has_phi, 
    .has_B = cls->has_B, 
    .has_rad = vms->has_rad, 
    .use_lo = cls->use_lo,
    .use_gpu = app->use_gpu,
  }; 
  cls->calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux);

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &app->basis,
    .phase_basis = &vms->basis,
    .conf_range =  &app->local,
    .hamil_range = &vms->hamil_range,
    .phase_range = &vms->local,
    .vel_map = vms->vel_map,
    .skip_cell_thresh = vms->info.skip_cell_thresh > 0.0 ? vms->info.skip_cell_thresh : 0.0,
    .model_id = vms->model_id,
    .has_E = cls->has_E, 
    .has_phi = cls->has_phi, 
    .has_B = cls->has_B, 
    .has_rad = vms->has_rad, 
    .poisson_tensor_conf = vms->conf_poisson_tensor,
    .hamil = vms->hamil,
    .qmem = cls->qmem, 
    .pot_tot = cls->pot_tot, 
    .conf_flux_surf = cls->conf_flux_surf,
    .vel_flux_surf = cls->vel_flux_surf, 
    .f_no_J = vms->f_no_J, 
    .rad = vms->rad, 
    .use_lo = cls->use_lo,
    .use_gpu = app->use_gpu,
  };  
  // Construct Vlasov equation and Hyper DG object for updating equation. 
  cls->eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  int up_dirs[GKYL_MAX_DIM];
  for (int d=0; d<pdim; ++d) {
    up_dirs[d] = d;
  }
  int num_up_dirs = pdim;

  // Default: no zero-flux BCs in any configuration-space direction.
  int zero_flux_flags[2*GKYL_MAX_DIM] = {false}; 
  for (int dir=0; dir<cdim; ++dir) {
    if (vms->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
      zero_flux_flags[dir] = true;
    }
    if (vms->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
      zero_flux_flags[dir+pdim] = true;
    }
  }
  // Default: zero-flux BCs in velocity space
  for (int dir=cdim; dir<pdim; ++dir) {
    zero_flux_flags[dir] = zero_flux_flags[dir+pdim] = 1;
  }

  cls->slvr = gkyl_hyper_dg_new(&vms->grid, &vms->basis, cls->eqn, 
    num_up_dirs, up_dirs, zero_flux_flags, 1, app->use_gpu);

  cls->rhs_func = vm_species_collisionless_rhs_enabled;
  if (vms->info.no_collisionless_terms) {
    cls->rhs_func = vm_species_collisionless_rhs_disabled;
  }
}

void
vm_species_collisionless_app_accel(gkyl_vlasov_app *app, struct vm_collisionless *cls, double tm)
{
  if (cls->has_app_accel) {
    gkyl_proj_on_basis_advance(cls->app_accel_proj, tm, &app->local_ext, cls->app_accel_host);
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(cls->app_accel, cls->app_accel_host);
    }
  }
}

void
vm_species_collisionless_rhs(gkyl_vlasov_app *app, struct vm_species *vms, 
  struct vm_collisionless *cls, const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs)
{
  cls->rhs_func(app, vms, cls, fin, em, rhs);
}

void
vm_species_collisionless_release(const struct gkyl_vlasov_app *app, 
  const struct vm_species *vms, const struct vm_collisionless *cls)
{
  if (vms->model_id == GKYL_MODEL_TRIAD || vms->model_id == GKYL_MODEL_TRIAD_GR) {
    gkyl_dg_vlasov_conf_flux_surf_release(cls->calc_conf_flux);
    gkyl_array_release(cls->conf_flux_surf);
  }
  gkyl_dg_vlasov_vel_flux_surf_release(cls->calc_vel_flux);
  gkyl_array_release(cls->qmem); 
  gkyl_array_release(cls->pot_tot); 
  gkyl_array_release(cls->app_accel);
  if (cls->has_app_accel) {
    gkyl_array_release(cls->app_accel_host);
    gkyl_proj_on_basis_release(cls->app_accel_proj);
  } 
  gkyl_array_release(cls->vel_flux_surf);  
  gkyl_hyper_dg_release(cls->slvr);
  gkyl_dg_eqn_release(cls->eqn);
}
