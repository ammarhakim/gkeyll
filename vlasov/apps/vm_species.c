#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <time.h>

static void 
vm_species_new_hamil(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  int vdim = app->vdim;  
  if (vms->model_id == GKYL_MODEL_DEFAULT || vms->model_id == GKYL_MODEL_SR) {
    // Hamiltonain is only a function of velocity space. 
    vms->hamil_range = vms->local_vel; 
    vms->hamil = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
    gkyl_dg_vlasov_calc_hamil(&vms->grid_vel, &vms->basis_vel, &vms->local_vel, 
      vms->model_id, vms->hamil, app->use_gpu); 

    // If relativistic, allocate additional updater for computing derived relativistic moments,
    // such as the spatial component of the four-velocity and pressure. 
    if (vms->model_id == GKYL_MODEL_SR) {
      vms->gamma_inv = mkarr(app->use_gpu, vms->basis_vel.num_basis, vms->local_vel.volume);
      vms->sr_vars = gkyl_dg_calc_sr_vars_new(&vms->grid, &vms->grid_vel,
        &app->basis,  &vms->basis_vel, &app->local, &vms->local_vel, app->use_gpu);
      // Project inverse of the Hamiltonian. 
      gkyl_calc_sr_vars_init_p_vars(vms->sr_vars, vms->gamma_inv);      
    }
  }
  else {
    // Hamiltonian is a full phase-space array. 
    vms->hamil_range = vms->local; 
    vms->hamil = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
    vms->hamil_host = vms->hamil;
    if (app->use_gpu){
      vms->hamil_host = mkarr(false, vms->basis.num_basis, vms->local_ext.volume);
    }

    // Allocate arrays for specified metric inverse
    vms->h_ij = mkarr(app->use_gpu, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    vms->h_ij_host = vms->h_ij;
    if (app->use_gpu){
      vms->h_ij_host = mkarr(false, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    }

    // Allocate arrays for specified metric inverse
    vms->h_ij_inv = mkarr(app->use_gpu, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    vms->h_ij_inv_host = vms->h_ij_inv;
    if (app->use_gpu){
      vms->h_ij_inv_host = mkarr(false, app->basis.num_basis*vdim*(vdim+1)/2, app->local_ext.volume);
    }

    // Allocate arrays for specified metric determinant
    vms->det_h = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    vms->det_h_host = vms->det_h;
    if (app->use_gpu){
      vms->det_h_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    }

    // Evaluate specified hamiltonian function at nodes to insure continuity of hamiltoniam
    struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&vms->grid, &vms->basis, 1, vms->info.hamil, vms->info.hamil_ctx);
    gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &vms->local_ext, vms->hamil_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->hamil, vms->hamil_host);
    }
    gkyl_eval_on_nodes_release(hamil_proj);

    // Evaluate specified metric function at nodes to insure continuity
    struct gkyl_eval_on_nodes* h_ij_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, vdim*(vdim+1)/2, vms->info.h_ij, vms->info.h_ij_ctx);
    gkyl_eval_on_nodes_advance(h_ij_proj, 0.0, &app->local, vms->h_ij_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->h_ij, vms->h_ij_host);
    }
    gkyl_eval_on_nodes_release(h_ij_proj);

    // Evaluate specified inverse metric function at nodes to insure continuity of the inverse 
    struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, vdim*(vdim+1)/2, vms->info.h_ij_inv, vms->info.h_ij_inv_ctx);
    gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &app->local, vms->h_ij_inv_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->h_ij_inv, vms->h_ij_inv_host);
    }
    gkyl_eval_on_nodes_release(h_ij_inv_proj);

    // Evaluate specified determinant metric function at nodes to insure continuity of the determinant
    struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 1, vms->info.det_h, vms->info.det_h_ctx);
    gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &app->local, vms->det_h_host);
    if (app->use_gpu){
      gkyl_array_copy(vms->det_h, vms->det_h_host);
    }
    gkyl_eval_on_nodes_release(det_h_proj);    
  }
}

static void
vm_species_apply_ic_dynamic(gkyl_vlasov_app *app, struct vm_species *vms, double tm)
{
  if (vms->num_init > 1) {
    gkyl_array_clear(vms->f, 0.0);  
    for (int k=0; k<vms->num_init; k++) {
      vm_species_projection_calc(app, vms, &vms->proj_init[k], vms->f1, tm);
      gkyl_array_accumulate(vms->f, 1.0, vms->f1);
    }
  }
  else {
    vm_species_projection_calc(app, vms, &vms->proj_init[0], vms->f, tm);
  }

  // Pre-compute applied acceleration in case it's time-independent
  vm_species_calc_app_accel(app, vms, tm);

  // we are pre-computing source for now as it is time-independent
  vm_species_source_calc(app, vms, &vms->src, tm);

  if (vms->calc_bflux) {
    vm_species_bflux_rhs(app, vms, &vms->bflux, vms->f, vms->f1);
  }

  // Optional runtime configuration to use BGK collisions but with fixed input 
  // temperature relaxation based on the initial temperature value. 
  if (vms->bgk.fixed_temp_relax) {
    vm_species_bgk_moms_fixed_temp(app, vms, &vms->bgk, vms->f);
  }
  
  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  if (vms->lower_bc[0].type == GKYL_SPECIES_FIXED_FUNC) {
    gkyl_bc_basic_buffer_fixed_func(vms->bc_lo[0], vms->bc_buffer_lo_fixed, vms->f);
  }
  if (vms->upper_bc[0].type == GKYL_SPECIES_FIXED_FUNC) {
    gkyl_bc_basic_buffer_fixed_func(vms->bc_up[0], vms->bc_buffer_up_fixed, vms->f);  
  }
}

static void
vm_species_apply_ic_static(gkyl_vlasov_app *app, struct vm_species *vms, double tm)
{
  if (vms->num_init > 1) {
    gkyl_array_clear(vms->f, 0.0);  
    struct gkyl_array *f_tmp = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
    for (int k=0; k<vms->num_init; k++) {
      vm_species_projection_calc(app, vms, &vms->proj_init[k], f_tmp, tm);
      gkyl_array_accumulate(vms->f, 1.0, f_tmp);
    }
    gkyl_array_release(f_tmp);
  }
  else {
    vm_species_projection_calc(app, vms, &vms->proj_init[0], vms->f, tm);
  }  
}

static void
vm_species_collisionless_rhs_included(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Set values of q/m*EM and the total potentials based on field and external forces configuration. 
  gkyl_array_clear(vms->qmem, 0.0); 
  if (vms->has_app_accel) {
    gkyl_array_accumulate_range(vms->qmem, 1.0, vms->app_accel, &app->local);
  }

  if (app->has_field) {
    if (app->field->has_ext_em) {
      gkyl_array_accumulate_range(vms->qmem, vms->qbym, app->field->ext_em, &app->local);
    }

    if (vms->field_id == GKYL_FIELD_E_B) {
      gkyl_array_accumulate_range(vms->qmem, vms->qbym, em, &app->local);
    }
    else if (vms->field_id == GKYL_FIELD_PHI) {
      gkyl_array_set_offset(vms->pot_tot, vms->qbym, app->field->phi, 0);
      if (app->field->has_ext_pot) {
        gkyl_array_accumulate_offset(vms->pot_tot, vms->qbym, app->field->ext_pot, 0);
      }
    }
  }

  // Compute the surface expansion of the phase space flux in velocity space. 
  gkyl_dg_vlasov_vel_flux_surf_advance(vms->calc_vel_flux, &app->local, &vms->local, 
    vms->hamil, vms->qmem, vms->pot_tot, fin, vms->cflrate, vms->vel_flux_surf);

  gkyl_hyper_dg_advance(vms->slvr, &vms->local, fin, vms->cflrate, rhs);

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_collisionless_rhs_empty(gkyl_vlasov_app *app, struct vm_species *vms, 
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
}

static void
vm_species_collisionless_rhs(gkyl_vlasov_app *app, struct vm_species *vms, 
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  vms->collisionless_rhs_func(app, vms, fin, em, rhs);
}

static double
vm_species_rhs_dynamic(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  gkyl_array_clear(vms->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  vm_species_collisionless_rhs(app, vms, fin, em, rhs);

  if (vms->lbo.collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_rhs(app, vms, &vms->lbo, fin, rhs);
  }
  else if (vms->bgk.collision_id == GKYL_BGK_COLLISIONS && !app->has_implicit_coll_scheme) {
    vms->bgk.implicit_step = false;
    vm_species_bgk_rhs(app, vms, &vms->bgk, fin, rhs);
  }

  if (vms->calc_bflux) {
    vm_species_bflux_rhs(app, vms, &vms->bflux, fin, rhs);
  }

  if (vms->rad.radiation_id == GKYL_VM_COMPTON_RADIATION) {
    vm_species_radiation_rhs(app, vms, &vms->rad, fin, rhs);
  }  

  // Reduce the CFL frequency anc compute stable dt needed by this species.
  app->stat.n_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(vms->omega_cfl, vms->cflrate, GKYL_MAX, &vms->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, vms->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omega_cfl_ho[0] = vms->omega_cfl[0];
  }
  double dt_out = app->cfl/omega_cfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return dt_out;
}

double
vm_species_rhs_implicit_dynamic(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  gkyl_array_clear(vms->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute implicit update and update rhs to new time step.   
  if (vms->bgk.collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_rhs(app, vms, &vms->bgk, fin, rhs);
  }
  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);

  if (vms->calc_bflux) {
    vm_species_bflux_rhs(app, vms, &vms->bflux, fin, rhs);
  }

  // Implicit step does not affect the stable time step. 
  double dt_out = DBL_MAX;  
  return dt_out;
}

static double
vm_species_rhs_static(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  double dt_out = DBL_MAX;  
  return dt_out;
}

static double
vm_species_rhs_implicit_static(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  double dt_out = DBL_MAX;  
  return dt_out;
}

static void
vm_species_apply_bc_dynamic(gkyl_vlasov_app *app, const struct vm_species *vms, struct gkyl_array *f,
  double tcurr)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(vms->comm, &vms->local, &vms->local_ext,
    num_periodic_dir, app->periodic_dirs, f); 
  
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (vms->lower_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(app, vms, &vms->bc_emission_lo, f, tcurr);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(vms->bc_lo[d], vms->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(vms->bc_lo[d], vms->bc_buffer_lo_fixed, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (vms->upper_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(app, vms, &vms->bc_emission_up, f, tcurr);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(vms->bc_up[d], vms->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(vms->bc_up[d], vms->bc_buffer_up_fixed, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(vms->comm, &vms->local, &vms->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_apply_bc_static(gkyl_vlasov_app *app, const struct vm_species *vms, struct gkyl_array *f, 
  double tcurr)
{
  // do nothing
}

static void
vm_species_step_f_dynamic(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

static void
vm_species_step_f_static(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // do nothing
}

static void
vm_species_combine_dynamic(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

static void
vm_species_combine_static(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // do nothing
}

static void
vm_species_copy_range_dynamic(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  gkyl_array_copy_range(out, inp, range);
}

static void
vm_species_copy_range_static(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // do nothing
}

static void
vm_species_write_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  vms->write_cfl_func(app, vms, tm, frame);
  vms->write_lte_func(app, vms, tm, frame);

  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = vms->basis.id
    }
  );
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, frame);
  
  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(vms->f_host, vms->f);
  }
  gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, vms->f_host, fileNm);
    
  vlasov_array_meta_release(mt);  

  app->stat.species_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
vm_species_write_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_write_cfl_enabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = 0,
      .basis_type = vms->basis.id,
    }
  );

  const char *fmt = "%s-%s-cflrate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, frame);
  gkyl_array_copy(vms->cflrate_host, vms->cflrate);
  gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt,
    vms->cflrate_host, fileNm);

  vlasov_array_meta_release(mt);  

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
vm_species_write_cfl_disabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_write_lte_enabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = vms->basis.id
    }
  );

  const char *fmt = "%s-%s_%d_lte.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, frame);

  vm_species_lte(app, vms, &vms->lte, vms->f);
  
  // copy data from device to host before writing it out
  // Just re-use f_host host-side array to avoid allocating 
  // more distribution function-size arrays. 
  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(vms->f_host, vms->lte.f_lte);
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, 
      vms->f_host, fileNm);
  }
  else {
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, 
      vms->lte.f_lte, fileNm);
  }

  vlasov_array_meta_release(mt);  

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

static void
vm_species_write_lte_disabled(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_write_mom_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  for (int m=0; m<vms->info.num_diag_moments; ++m) {
    struct timespec wtm = gkyl_wall_clock();
    vm_species_moment_calc(&vms->moms[m], vms->local, app->local, vms->f);
    app->stat.n_mom += 1;
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wtm);
      
    struct timespec wst = gkyl_wall_clock();
    if (app->use_gpu) {
      gkyl_array_copy(vms->moms[m].marr_host, vms->moms[m].marr);
    }

    const char *fmt = "%s-%s_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name,
      gkyl_distribution_moments_strs[vms->info.diag_moments[m]], frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name,
      gkyl_distribution_moments_strs[vms->info.diag_moments[m]], frame);
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      vms->moms[m].marr_host, fileNm);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_diag_io += 1;
  }
  
  vlasov_array_meta_release(mt); 

  app->stat.n_diag += 1;
}

static void
vm_species_write_mom_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  // do nothing
}

static void
vm_species_calc_integrated_mom_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  int vdim = app->vdim;
  double avals[2+vdim], avals_global[2+vdim];
  
  vm_species_moment_calc(&vms->integ_moms, vms->local, app->local, vms->f);
  app->stat.n_mom += 1;

  // reduce to compute sum over whole domain, append to diagnostics
  if (app->use_gpu) {
    gkyl_array_reduce_range(vms->red_integ_diag, vms->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(avals, vms->red_integ_diag, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    gkyl_array_reduce_range(avals, vms->integ_moms.marr_host, GKYL_SUM, &app->local);
  }

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, avals, avals_global);
  gkyl_dynvec_append(vms->integ_diag, tm, avals_global);  

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_calc_integrated_mom_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  // do nothing
}

static void
vm_species_write_integrated_mom_dynamic(gkyl_vlasov_app *app, struct vm_species *vms)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, "imom");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, "imom");
    
    if (vms->is_first_integ_write_call) {
      gkyl_dynvec_write(vms->integ_diag, fileNm);
      vms->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(vms->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(vms->integ_diag);
  app->stat.n_diag_io += 1;
  
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
vm_species_write_integrated_mom_static(gkyl_vlasov_app* app, struct vm_species *vms)
{
  // do nothing
}

static void
vm_species_calc_L2_dynamic(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_dg_calc_l2_range(vms->basis, 0, vms->L2_f, 0, vms->f, vms->local);
  gkyl_array_scale_range(vms->L2_f, vms->grid.cellVolume, &vms->local);
  
  double L2[1] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_reduce_range(vms->red_L2_f, vms->L2_f, GKYL_SUM, &vms->local);
    gkyl_cu_memcpy(L2, vms->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(L2, vms->L2_f, GKYL_SUM, &vms->local);
  }
  double L2_global[1] = { 0.0 };
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);
  
  gkyl_dynvec_append(vms->integ_L2_f, tm, L2_global);  

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

static void
vm_species_calc_L2_static(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  // do nothing
}

static void
vm_species_write_L2_dynamic(gkyl_vlasov_app* app, struct vm_species *vms)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write the L2 norm.
    const char *fmt = "%s-%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, vms->info.name, "L2");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, vms->info.name, "L2");

    if (vms->is_first_integ_L2_write_call) {
      gkyl_dynvec_write(vms->integ_L2_f, fileNm);
      vms->is_first_integ_L2_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(vms->integ_L2_f, fileNm);
    }
  }
  gkyl_dynvec_clear(vms->integ_L2_f);

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

static void
vm_species_write_L2_static(gkyl_vlasov_app* app, struct vm_species *vms)
{
  // do nothing
}

static void
vm_species_release_dynamic(const gkyl_vlasov_app* app, const struct vm_species *vms)
{
  // Release various arrays used in time stepping and boundary conditions. 
  gkyl_array_release(vms->f1);
  gkyl_array_release(vms->fnew);
  gkyl_array_release(vms->cflrate);
  gkyl_array_release(vms->bc_buffer);
  gkyl_array_release(vms->bc_buffer_lo_fixed);
  gkyl_array_release(vms->bc_buffer_up_fixed);

  vm_species_moment_release(app, &vms->integ_moms); 

  gkyl_array_release(vms->L2_f);
  gkyl_dynvec_release(vms->integ_L2_f);
  gkyl_dynvec_release(vms->integ_diag);

  if (vms->source_id) {
    vm_species_source_release(app, &vms->src);
  }

  if (vms->lbo.collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_release(app, &vms->lbo);
  }
  else if (vms->bgk.collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_release(app, &vms->bgk);
  }

  if (vms->rad.radiation_id == GKYL_VM_COMPTON_RADIATION) {
    vm_species_radiation_release(app, &vms->rad);
  }

  if (vms->calc_bflux) {
    vm_species_bflux_release(app, &vms->bflux);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (vms->lower_bc[d].type == GKYL_SPECIES_EMISSION) {
      vm_species_emission_release(&vms->bc_emission_lo);
    }
    else {
      gkyl_bc_basic_release(vms->bc_lo[d]);
    }
    
    if (vms->upper_bc[d].type == GKYL_SPECIES_EMISSION) {
      vm_species_emission_release(&vms->bc_emission_up);
    }
    else {
      gkyl_bc_basic_release(vms->bc_up[d]);
    }
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(vms->omega_cfl);
    gkyl_cu_free(vms->red_L2_f);
    gkyl_cu_free(vms->red_integ_diag);
  }
  else {
    gkyl_free(vms->omega_cfl);
  }  
}

static void
vm_species_release_static(const gkyl_vlasov_app* app, const struct vm_species *vms)
{
}

// Initialize dynamic species object.
static void
vm_species_new_dynamic(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  // Allocate remaining distribution function arrays for RK time stepping. 
  vms->f1 = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
  vms->fnew = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);

  // Allocate cflrate (scalar array) and maximum frequency for stable time step.
  vms->cflrate = mkarr(app->use_gpu, 1, vms->local_ext.volume);
  // Host-side cflrate for I/O on GPUs.
  vms->cflrate_host = app->use_gpu ? mkarr(false, vms->cflrate->ncomp, vms->cflrate->size)
                                   : gkyl_array_acquire(vms->cflrate);    
  if (app->use_gpu) {
    vms->omega_cfl = gkyl_cu_malloc(sizeof(double));
  }
  else {
    vms->omega_cfl = gkyl_malloc(sizeof(double));
  }

  // Allocate data for integrated moments. 
  vm_species_moment_init(app, vms, &vms->integ_moms, GKYL_F_MOMENT_M0M1M2, true);  
  // Array for storing f^2 in each cell.
  vms->L2_f = mkarr(app->use_gpu, 1, vms->local_ext.volume);
  if (app->use_gpu) {
    vms->red_L2_f = gkyl_cu_malloc(sizeof(double));
    vms->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments and f^2.
  vms->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  vms->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  vms->is_first_integ_L2_write_call = true;
  vms->is_first_integ_write_call = true; 

  if (vms->info.collisions.collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_init(app, vms, &vms->lbo);
  }
  else if (vms->info.collisions.collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_init(app, vms, &vms->bgk);
  }  

  if (vms->info.radiation.radiation_id == GKYL_VM_COMPTON_RADIATION) {
    vm_species_radiation_init(app, vms, &vms->rad);
  }
  
  // Allocate buffer for applying BCs.
  long buff_sz = 0;
  for (int dir=0; dir<cdim; ++dir) {
    long vol = GKYL_MAX2(vms->lower_skin[dir].volume, vms->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  vms->bc_buffer = mkarr(app->use_gpu, vms->basis.num_basis, buff_sz);
  // Buffers for fixed function BCs on distribution function.
  vms->bc_buffer_lo_fixed = mkarr(app->use_gpu, vms->basis.num_basis, buff_sz);
  vms->bc_buffer_up_fixed = mkarr(app->use_gpu, vms->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (vms->lower_bc[d].type == GKYL_SPECIES_EMISSION) {
      vms->emit_lo = true;
      vms->calc_bflux = true;
      vm_species_emission_init(app, &vms->bc_emission_lo, d, GKYL_LOWER_EDGE, vms->lower_bc[d].aux_ctx);
    }
    else {
      if (vms->lower_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (vms->lower_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (vms->lower_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_DISTF_REFLECT;
      else if (vms->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      vms->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, vms->basis_on_dev,
        &vms->lower_skin[d], &vms->lower_ghost[d], vms->f->ncomp, cdim, app->use_gpu);
    }

    // Upper BC updater. Copy BCs by default.
    if (vms->upper_bc[d].type == GKYL_SPECIES_EMISSION) {
      vms->emit_up = true;
      vms->calc_bflux = true;
      vm_species_emission_init(app, &vms->bc_emission_up, d, GKYL_UPPER_EDGE, vms->upper_bc[d].aux_ctx);
    }
    else {
      if (vms->upper_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (vms->upper_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (vms->upper_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_DISTF_REFLECT;
      else if (vms->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      vms->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, vms->basis_on_dev,
        &vms->upper_skin[d], &vms->upper_ghost[d], vms->f->ncomp, cdim, app->use_gpu);
    }
  }

  // Intitalize boundary flux updater if we need boundary fluxes.
  if (vms->calc_bflux) {
    vm_species_bflux_init(app, vms, &vms->bflux);
  }

  // Set function pointers.
  vms->apply_ic_func = vm_species_apply_ic_dynamic; 
  vms->rhs_func = vm_species_rhs_dynamic;
  vms->rhs_implicit_func = vm_species_rhs_implicit_dynamic;
  vms->bc_func = vm_species_apply_bc_dynamic;
  vms->release_func = vm_species_release_dynamic;
  vms->step_f_func = vm_species_step_f_dynamic;
  vms->combine_func = vm_species_combine_dynamic;
  vms->copy_func = vm_species_copy_range_dynamic;
  vms->write_func = vm_species_write_dynamic;
  if (vms->info.correct.output_f_lte) {
    vms->write_lte_func = vm_species_write_lte_enabled;
  }
  else {
    vms->write_lte_func = vm_species_write_lte_disabled;
  }
  if (vms->info.write_omega_cfl) {
    vms->write_cfl_func = vm_species_write_cfl_enabled;
  }
  else {
    vms->write_cfl_func = vm_species_write_cfl_disabled;
  }
  vms->write_mom_func = vm_species_write_mom_dynamic;
  vms->calc_integrated_mom_func = vm_species_calc_integrated_mom_dynamic;
  vms->write_integrated_mom_func = vm_species_write_integrated_mom_dynamic;
  vms->calc_L2_func = vm_species_calc_L2_dynamic;
  vms->write_L2_func = vm_species_write_L2_dynamic;
}

// Initialize static species object.
static void
vm_species_new_static(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  // When using static species, only need one allocation.
  vms->f1 = vms->f;
  vms->fnew = vms->f;
  
  // Set function pointers.
  vms->apply_ic_func = vm_species_apply_ic_static; 
  vms->rhs_func = vm_species_rhs_static;
  vms->rhs_implicit_func = vm_species_rhs_implicit_static;
  vms->bc_func = vm_species_apply_bc_static;
  vms->release_func = vm_species_release_static;
  vms->step_f_func = vm_species_step_f_static;
  vms->combine_func = vm_species_combine_static;
  vms->copy_func = vm_species_copy_range_static;
  vms->write_func = vm_species_write_static;
  vms->write_lte_func = vm_species_write_lte_disabled;
  vms->write_cfl_func = vm_species_write_cfl_disabled;
  vms->write_mom_func = vm_species_write_mom_static;
  vms->calc_integrated_mom_func = vm_species_calc_integrated_mom_static;
  vms->write_integrated_mom_func = vm_species_write_integrated_mom_static;
  vms->calc_L2_func = vm_species_calc_L2_static;
  vms->write_L2_func = vm_species_write_L2_static;
}

// End static function definitions.

// Initialize species object.
void
vm_species_init(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_species *vms)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm_app_inp->cells[d];
    lower[d] = vm_app_inp->lower[d];
    upper[d] = vm_app_inp->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // Full phase-space grid.
    cells[cdim+d] = vms->info.cells[d];
    lower[cdim+d] = vms->info.lower[d];
    upper[cdim+d] = vms->info.upper[d];
    ghost[cdim+d] = 0; // No ghost-cells in velocity space.

    // Only velocity space.
    cells_vel[d] = vms->info.cells[d];
    lower_vel[d] = vms->info.lower[d];
    upper_vel[d] = vms->info.upper[d];
    ghost_vel[d] = 0; // No ghost-cells in velocity space.
  }

  // Allocate device basis if we are using GPUs.
  if (app->use_gpu) {
    vms->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    vms->basis_on_dev = &vms->basis;
  }  

  // Determine basis type from configuration-space basis and create phase-space basis. 
  enum gkyl_basis_type b_type = app->basis.b_type;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&vms->basis, pdim, app->poly_order);
      gkyl_cart_modal_serendip(&vms->basis_surf, pdim-1, app->poly_order);
      gkyl_cart_modal_serendip(&vms->basis_vel, vdim, app->poly_order);
      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(vms->basis_on_dev, pdim, app->poly_order);
      }
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&vms->basis, pdim, app->poly_order); 
      gkyl_cart_modal_tensor(&vms->basis_surf, pdim-1, app->poly_order); 
      gkyl_cart_modal_tensor(&vms->basis_vel, vdim, app->poly_order); 
      if (app->use_gpu) {
        gkyl_cart_modal_tensor_cu_dev(vms->basis_on_dev, pdim, app->poly_order);
      }
      break;    
    default:
      assert(false);
      break;
  }  

  // Full phase-space grid.
  gkyl_rect_grid_init(&vms->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vms->grid, ghost, &vms->global_ext, &vms->global);
  
  // Velocity-space grid.
  gkyl_rect_grid_init(&vms->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&vms->grid_vel, ghost_vel, &vms->local_ext_vel, &vms->local_vel);

  // Phase-space communicator.
  vms->comm = gkyl_comm_extend_comm(app->comm, &vms->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // local = conf-local X local_vel
  gkyl_range_ten_prod(&local, &app->local, &vms->local_vel);
  gkyl_create_ranges(&local, ghost, &vms->local_ext, &vms->local);

  // Determine which directions are not periodic.
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  // Default: no zero-flux BCs in any configuration-space direction.
  int zero_flux_flags[2*GKYL_MAX_DIM] = {false}; 
  for (int dir=0; dir<cdim; ++dir) {
    vms->lower_bc[dir].type = vms->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const struct gkyl_vlasov_bcs *bc;
      if (dir == 0) {
        bc = &vms->info.bcx;
      }
      else if (dir == 1) {
        bc = &vms->info.bcy;
      }
      else {
        bc = &vms->info.bcz;
      }

      vms->lower_bc[dir] = bc->lower;
      vms->upper_bc[dir] = bc->upper;
      if (vms->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        zero_flux_flags[dir] = true;
      }
      if (vms->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        zero_flux_flags[dir+pdim] = true;
      }
    }
  }
  // Default: zero-flux BCs in velocity space
  for (int dir=cdim; dir<pdim; ++dir) {
    zero_flux_flags[dir] = zero_flux_flags[dir+pdim] = 1;
  }

  // Store model type.
  vms->model_id = vms->info.model_id; 
  // Store field type.
  vms->field_id = app->has_field ? app->field->field_id : GKYL_FIELD_NULL; 
  
  // Allocate distribution function array.
  vms->f = mkarr(app->use_gpu, vms->basis.num_basis, vms->local_ext.volume);
  // Host-side distribution function for I/O on GPUs.
  vms->f_host = app->use_gpu ? mkarr(false, vms->f->ncomp, vms->f->size)
                             : gkyl_array_acquire(vms->f);  

  // Allocate array to store q/m*(E,B) or potentials (q/m*phi + m*phi_g, q/m*A) depending on equation system. 
  // Note: the potentials are the total potentials and thus can include both (or either) gravitational
  // or electrostatic interactions. 
  vms->qbym = vms->info.charge/vms->info.mass;
  vms->qmem = mkarr(app->use_gpu, 8*app->basis.num_basis, app->local_ext.volume);
  vms->pot_tot = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);

  // Initialize applied acceleration for use in force update. 
  vms->app_accel = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(vms->app_accel, 0.0);
  vms->has_app_accel = false;
  vms->app_accel_evolve = false;
  // setup applied acceleration
  if (vms->info.app_accel) {
    vms->has_app_accel = true;
    if (vms->info.app_accel_evolve) {
      vms->app_accel_evolve = vms->info.app_accel_evolve;
    }

    // Host-side distribution function for projection on GPUs.
    vms->app_accel_host = app->use_gpu ? mkarr(false, vms->app_accel->ncomp, vms->app_accel->size)
                                       : gkyl_array_acquire(vms->app_accel);      
    vms->app_accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->basis, app->basis.poly_order+1,
      3, vms->info.app_accel, vms->info.app_accel_ctx);
  }

  // Determine which forces we need based on combination of field ID and presence 
  // of applied accelerations and external fields/potentials. 
  vms->has_qmem = false; 
  vms->has_phi = false; 
  if (app->has_field) {
    if (vms->field_id == GKYL_FIELD_E_B || app->field->has_ext_em) {
      vms->has_qmem = true; 
    } 
    if (vms->field_id == GKYL_FIELD_PHI) {
      vms->has_phi = true; 
    }
  }
  else if (vms->has_app_accel) {
    vms->has_qmem = true; 
  }

  // Construct Hamiltonian. 
  vm_species_new_hamil(vm_app_inp, app, vms); 

  // Allocate modal surface expansion of velocity space flux array. 
  vms->vel_flux_surf = mkarr(app->use_gpu, vdim*vms->basis_surf.num_basis, vms->local_ext.volume);
  struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
    .phase_grid = &vms->grid, 
    .conf_basis = &app->basis,
    .phase_basis = &vms->basis,
    .hamil_range = &vms->hamil_range,
    .model_id = vms->model_id,
    .has_qmem = vms->has_qmem, 
    .has_phi = vms->has_phi, 
    .use_gpu = app->use_gpu,
  }; 
  vms->calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux); 

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &app->basis,
    .phase_basis = &vms->basis,
    .conf_range =  &app->local,
    .hamil_range = &vms->hamil_range,
    .phase_range = &vms->local,
    .model_id = vms->model_id,
    .has_qmem = vms->has_qmem, 
    .has_phi = vms->has_phi, 
    .hamil = vms->hamil,
    .qmem = vms->qmem, 
    .pot_tot = vms->pot_tot, 
    .vel_flux_surf = vms->vel_flux_surf, 
    .use_gpu = app->use_gpu,
  };  
  // Construct Vlasov equation and Hyper DG object for updating equation. 
  vms->eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  int up_dirs[GKYL_MAX_DIM];
  for (int d=0; d<pdim; ++d) {
    up_dirs[d] = d;
  }
  int num_up_dirs = pdim;
  vms->slvr = gkyl_hyper_dg_new(&vms->grid, &vms->basis, vms->eqn, 
    num_up_dirs, up_dirs, zero_flux_flags, 1, app->use_gpu);

  vms->collisionless_rhs_func = vm_species_collisionless_rhs_included;
  if (vms->info.no_collisionless_terms) {
    vms->collisionless_rhs_func = vm_species_collisionless_rhs_empty;
  }

  // Allocate data for momentum (for use in current accumulation).
  vm_species_moment_init(app, vms, &vms->m1i, GKYL_F_MOMENT_M1_FROM_H, false);
  // Allocate date for density (for use in charge density accumulation and weak division for V_drift).
  vm_species_moment_init(app, vms, &vms->m0, GKYL_F_MOMENT_M0, false);

  // Allocate data for diagnostic moments.
  int ndm = vms->info.num_diag_moments;
  vms->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m) {
    vm_species_moment_init(app, vms, &vms->moms[m], vms->info.diag_moments[m], false); 
  }

  // Create skin/ghost ranges for applying BCs.
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&vms->lower_skin[dir], &vms->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &vms->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&vms->upper_skin[dir], &vms->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &vms->local_ext, ghost);
  }

  // Global skin and ghost ranges, only valid (i.e. volume>0) in ranges abutting boundaries.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&vms->global_lower_skin[dir], &vms->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE, &vms->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&vms->global_upper_skin[dir], &vms->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE, &vms->global_ext, ghost);

    gkyl_sub_range_intersect(&vms->global_lower_skin[dir], &vms->local_ext, &vms->global_lower_skin[dir]);
    gkyl_sub_range_intersect(&vms->global_upper_skin[dir], &vms->local_ext, &vms->global_upper_skin[dir]);

    gkyl_sub_range_intersect(&vms->global_lower_ghost[dir], &vms->local_ext, &vms->global_lower_ghost[dir]);
    gkyl_sub_range_intersect(&vms->global_upper_ghost[dir], &vms->local_ext, &vms->global_upper_ghost[dir]);
  }

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine.
  // Projection routine optionally corrects all the Maxwellian/LTE moments.
  // This routine is utilized by BGK collisions.
  vms->lte = (struct vm_lte) { };
  struct correct_all_moms_inp corr_inp = { 
    .correct_all_moms = vms->info.correct.correct_all_moms, 
    .max_iter = vms->info.correct.max_iter > 0 ? vms->info.correct.max_iter : 50, 
    .iter_eps = vms->info.correct.iter_eps > 0 ? vms->info.correct.iter_eps : 1e-10, 
    .use_last_converged = vms->info.correct.use_last_converged, 
  };
  vm_species_lte_init(app, vms, &vms->lte, corr_inp);

  // Initialize projection routine for initial conditions. 
  vms->num_init = vms->info.num_init;
  for (int k=0; k<vms->num_init; k++) {
    vm_species_projection_init(app, vms, vms->info.projection[k], &vms->proj_init[k]);
  }

  // Initialize empty structs. New methods will fill them if specified.
  vms->src = (struct vm_source) { };
  vms->lbo = (struct vm_lbo_collisions) { };
  vms->bgk = (struct vm_bgk_collisions) { };
  vms->rad = (struct vm_rad_drag) { };
  vms->bflux = (struct vm_boundary_fluxes) { };

  // Set species source id. 
  vms->source_id = vms->info.source.source_id;
  if (vms->source_id == GKYL_BFLUX_SOURCE) {
    vms->calc_bflux = true;
  }

  if (!vms->info.is_static) {
    vm_species_new_dynamic(vm_app_inp, app, vms);
  }
  else {
    vm_species_new_static(vm_app_inp, app, vms);
  }
}

void
vm_species_calc_app_accel(gkyl_vlasov_app *app, struct vm_species *vms, double tm)
{
  if (vms->has_app_accel) {
    gkyl_proj_on_basis_advance(vms->app_accel_proj, tm, &app->local_ext, vms->app_accel_host);
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(vms->app_accel, vms->app_accel_host);
    }
  }
}

void
vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *vms, double tm)
{
  return vms->apply_ic_func(app, vms, tm); 
}

// Compute the RHS for species update, returning maximum stable time-step.
double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  return vms->rhs_func(app, vms, fin, em, rhs);
}

// Compute the implicit RHS for species update; returns DBL_MAX because implicit
// update does not affect the stable time step.
double
vm_species_rhs_implicit(gkyl_vlasov_app *app, struct vm_species *vms,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  return vms->rhs_implicit_func(app, vms, fin, rhs, dt);
}

// Accummulate function for forward Euler method.
void
vm_species_step_f(struct vm_species *vms, struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  vms->step_f_func(out, a, inp);
}

// Combine function for SSP-RK3 updates.
void
vm_species_combine(struct vm_species *vms, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  vms->combine_func(out, c1, arr1, c2, arr2, rng);
}

// Copy function for SSP-RK3 updates.
void
vm_species_copy_range(struct vm_species *vms, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  vms->copy_func(out, inp, range);
}

// Apply boundary conditions to the distribution function.
void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *vms, struct gkyl_array *f, 
  double tcurr)
{
  vms->bc_func(app, vms, f, tcurr);
}

// Write distribution function. 
void
vm_species_write(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  if (frame == 0) {
    vm_species_write_dynamic(app, vms, tm, frame);
  }
  else {
    vms->write_func(app, vms, tm, frame);
  }
}

// Write distribution function moments. 
void
vm_species_write_mom(gkyl_vlasov_app* app, struct vm_species *vms, double tm, int frame)
{
  if (frame == 0) {
    vm_species_write_mom_dynamic(app, vms, tm, frame);
  }
  else {
    vms->write_mom_func(app, vms, tm, frame);
  }
}

// Calculate integrated moments. 
void
vm_species_calc_integrated_mom(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  vms->calc_integrated_mom_func(app, vms, tm);
}

// Write integrated moments. 
void
vm_species_write_integrated_mom(gkyl_vlasov_app* app, struct vm_species *vms)
{
  vms->write_integrated_mom_func(app, vms);
}

// Calculate integrated f^2. 
void
vm_species_calc_L2(gkyl_vlasov_app* app, struct vm_species *vms, double tm)
{
  vms->calc_L2_func(app, vms, tm);
}

// Write integrated f^2.
void
vm_species_write_L2(gkyl_vlasov_app* app, struct vm_species *vms)
{
  vms->write_L2_func(app, vms);
}

// Write the number of correction iterations to output statistics struct. 
void
vm_species_n_iter_corr(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    app->stat.num_corr[i] = app->species[i].lte.num_corr;
    app->stat.n_iter_corr[i] = app->species[i].lte.n_iter;
  }
}

void
vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *vms)
{
  // Release resources for Vlasov species.
  gkyl_array_release(vms->f);
  gkyl_comm_release(vms->comm);
  if (app->use_gpu) {
    gkyl_array_release(vms->f_host);
    gkyl_cu_free(vms->basis_on_dev);
  }

  for (int k=0; k<vms->num_init; k++) {
    vm_species_projection_release(app, &vms->proj_init[k]);
  }

  gkyl_array_release(vms->qmem); 
  gkyl_array_release(vms->pot_tot); 
  gkyl_array_release(vms->app_accel);
  if (vms->has_app_accel) {
    gkyl_array_release(vms->app_accel_host);
    gkyl_proj_on_basis_release(vms->app_accel_proj);
  } 
  gkyl_array_release(vms->vel_flux_surf); 

  // Release arrays for different types of Vlasov equations.
  if (vms->model_id  == GKYL_MODEL_DEFAULT || vms->model_id  == GKYL_MODEL_SR) {
    if (vms->model_id == GKYL_MODEL_SR) {
      gkyl_dg_calc_sr_vars_release(vms->sr_vars);
      gkyl_array_release(vms->gamma_inv);   
    }
  }
  else  {
    gkyl_array_release(vms->h_ij);
    gkyl_array_release(vms->h_ij_inv);
    gkyl_array_release(vms->det_h);
    if (app->use_gpu){
      gkyl_array_release(vms->hamil_host);
      gkyl_array_release(vms->h_ij_host);
      gkyl_array_release(vms->h_ij_inv_host);
      gkyl_array_release(vms->det_h_host);
    }
  }
  gkyl_array_release(vms->hamil);  

  // Release equation object and solver.
  gkyl_dg_eqn_release(vms->eqn);
  gkyl_hyper_dg_release(vms->slvr);

  // Release moment data.
  vm_species_moment_release(app, &vms->m1i);
  vm_species_moment_release(app, &vms->m0);
  for (int i=0; i<vms->info.num_diag_moments; ++i) {
    vm_species_moment_release(app, &vms->moms[i]);
  }
  gkyl_free(vms->moms);
  
  vm_species_lte_release(app, &vms->lte);

  vms->release_func(app, vms);
}