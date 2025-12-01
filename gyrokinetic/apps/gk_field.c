#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>


void
gk_field_calc_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  field->calc_energy_dt_func(app, field, dt, energy_reduced);
}

void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm)
{
  if (field->has_phi_wall_lo && field->phi_wall_lo_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_lo_proj, tm, &app->local_ext, field->phi_wall_lo_host);
    if (app->use_gpu) { // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(field->phi_wall_lo, field->phi_wall_lo_host);
    }
  }
  if (field->has_phi_wall_up && field->phi_wall_up_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_up_proj, tm, &app->local_ext, field->phi_wall_up_host);
    if (app->use_gpu) { // note: phi_wall_up_host is same as phi_wall_up when not on GPUs
      gkyl_array_copy(field->phi_wall_up, field->phi_wall_up_host);
    }
  }
}

void
gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(field->rho_c, 0.0);
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    gk_species_moment_calc(&s->m0, s->local, app->local, fin[i]);
    field->accumulate_rhoc_func(app, field, s, bflux[i]);
  } 
  app->stat.field_phi_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x).

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}


// Initialize field object.
struct gk_field* 
gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app)
{
  struct gk_field *f = gkyl_malloc(sizeof(struct gk_field));

  f->info = gk->field;

  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD_ES;

  f->calc_init_field = !f->info.zero_init_field;
  f->update_field = !f->info.is_static;
  // The combination update_field=true, calc_init_field=false is not allowed.
  assert(!(f->update_field && (!f->calc_init_field)));

  // allocate arrays for charge density
  f->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->rho_c_global_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // allocate arrays for electrostatic potential
  // global phi (only used in 1x simulations)
  f->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  // local phi (assuming domain decomposition is *only* in z right now)
  f->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    f->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->apardot_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  f->init_phi_pol = false;
  if (f->info.polarization_potential) {
    gk_field_polarization_potential_new(f, app);
  }

  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  f->epsilon = 0;
  f->kSq = 0;  // not currently used by fem_perp_poisson

  // Factors for ES energy. 
  f->es_energy_fac = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  f->es_energy_fac_1d = 0.0;

  // Initialize the individual field solvers.
  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    gk_field_fem_new_boltzmann(app, f);
  } else {
    if (app->cdim == 1) {
      gk_field_fem_new_1x(app, f);
    } else if (app->cdim > 1) {
      gk_field_fem_new_2x3x(app, f);
    }
  }

  f->phi_host = f->phi_smooth;  
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
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

  // Set up biased lower wall (same size as electrostatic potential), by default is 0.0
  gk_field_biassed_wall_new(app, f);

  return f;
}

void
gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Compute the electrostatic potential.
  struct timespec wst = gkyl_wall_clock();
  field->rhs_phi_func(app, field); // set inside the initialization function
  app->stat.field_phi_solve_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_file_import_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file.
  struct gkyl_app_restart_status rstat = header_from_file(app, inp.file_name);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    struct gkyl_app_restart_status rstat;
    rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, inp.file_name);
    gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
  }
  else {
    assert(false);
  }
}

void
gk_field_project_init(struct gkyl_gyrokinetic_app *app)
{
  // Project the initial field.
  struct gkyl_eval_on_nodes *phi_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
    1, app->field->info.init_field_profile, app->field->info.init_field_profile_ctx);
  gkyl_eval_on_nodes_advance(phi_proj, 0.0, &app->local, app->field->phi_host);
  gkyl_eval_on_nodes_release(phi_proj);
  gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();
  field->calc_energy_func(app, field, tm);
  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

// release resources for field
void
gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->apardot_fem);
  }

  if (f->init_phi_pol) {
    gk_field_polarization_potential_release(f);
  }

  gkyl_array_release(f->es_energy_fac);
  gkyl_dynvec_release(f->integ_energy);

  // Release field-solver-specific resources
  f->solver_release_func(app, f);

  // Release common energy diagnostics
  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
    gkyl_cu_free(f->em_energy_red);
    gkyl_cu_free(f->em_energy_red_global);
  } else {
    gkyl_free(f->em_energy_red);
    gkyl_free(f->em_energy_red_global);
  }

  if (f->info.time_rate_diagnostics) {
    gk_field_time_rate_diags_release(app, f);
  }

  gk_field_biassed_wall_release(app, f);

  gkyl_free(f);
}
