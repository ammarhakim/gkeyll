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

// Functions related to setting the potential by adjusting the polarization density

static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}

static void
gk_field_polarization_potential_new(struct gk_field *f, struct gkyl_gyrokinetic_app *app)
{
  // Project the initial potential onto a p+1 tensor basis and compute the polarization
  // density to use use by species in calculating the initial ion density.
  f->init_phi_pol = true;
  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order + 1);

  f->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
  struct gkyl_array *phi_pol_ho = app->use_gpu ? mkarr(false, f->phi_pol->ncomp, f->phi_pol->size)
                                               : gkyl_array_acquire(f->phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &phi_pol_basis,
      .num_ret_vals = 1,
      .eval = f->info.polarization_potential,
      .ctx = f->info.polarization_potential_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
  });

  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &app->local, phi_pol_ho);
  gkyl_array_copy(f->phi_pol, phi_pol_ho);

  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_release(phi_pol_ho);
}

static void
gk_field_polarization_potential_release(struct gk_field *gkf)
{
  gkyl_array_release(gkf->phi_pol);
}

// Functions related to the field energy allocations, diagnostics, and release

static void
gk_field_calc_energy_dt_enabled(gkyl_gyrokinetic_app *app, const struct gk_field *gkf,
  double dt, struct gkyl_array *fields[], double *energy_reduced)
{
  struct timespec wst = gkyl_wall_clock();
  // Calculate ES energy.
  struct gkyl_array *phi_curr = gk_field_get_phi_from_fields(gkf, fields);
  gkyl_array_integrate_advance(gkf->calc_em_energy, phi_curr,
    1.0/dt, gkf->es_energy_fac, &app->local, &app->local, energy_reduced);
  app->stat.phidot_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_calc_energy_dt_disabled(gkyl_gyrokinetic_app *app, const struct gk_field *gkf,
  double dt, struct gkyl_array *fields[], double *energy_reduced)
{
  // Do nothing.
}

static void
gk_field_calc_energy_enabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *gkf, double tm)
{
  gkyl_array_integrate_advance(gkf->calc_em_energy, gkf->phi,
                               1.0, gkf->es_energy_fac, &app->local, &app->local, gkf->em_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, gkf->em_energy_red, gkf->em_energy_red_global);

  double energy_global[1] = {0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(energy_global, gkf->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  } else {
    energy_global[0] = gkf->em_energy_red_global[0];
  }

  if (app->cdim == 1) {
    energy_global[0] *= gkf->es_energy_fac_1d;
  }

  gkyl_dynvec_append(gkf->integ_energy, tm, energy_global);

  if (gkf->info.time_rate_diagnostics) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, gkf->em_energy_red_old, gkf->em_energy_red_global);
    double energy_dot_global_old[1] = {0.0};
    if (app->use_gpu) {
      gkyl_cu_memcpy(energy_dot_global_old, gkf->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    } else {
      energy_dot_global_old[0] = gkf->em_energy_red_global[0];
    }
    if (app->cdim == 1) {
      energy_dot_global_old[0] *= gkf->es_energy_fac_1d;
    }

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, gkf->em_energy_red_new, gkf->em_energy_red_global);
    double energy_dot_global_new[1] = {0.0};
    if (app->use_gpu) {
      gkyl_cu_memcpy(energy_dot_global_new, gkf->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    } else {
      energy_dot_global_new[0] = gkf->em_energy_red_global[0];
    }
    if (app->cdim == 1) {
      energy_dot_global_new[0] *= gkf->es_energy_fac_1d;
    }

    double energy_dot_global[1] = {0.0};
    energy_dot_global[0] = energy_dot_global_new[0] - energy_dot_global_old[0];

    gkyl_dynvec_append(gkf->integ_energy_dot, tm, energy_dot_global);
  }
}

static void
gk_field_calc_energy_disabled(struct gkyl_gyrokinetic_app *app, const struct gk_field *gkf, double tm)
{
  // Do nothing.
}

static void
gk_field_time_rate_diags_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkf->calc_energy_dt_func = gk_field_calc_energy_dt_enabled;
  if (app->use_gpu) {
    gkf->em_energy_red_new = gkyl_cu_malloc(sizeof(double[1]));
    gkf->em_energy_red_old = gkyl_cu_malloc(sizeof(double[1]));
    gkyl_cu_memset(gkf->em_energy_red_new, 0, sizeof(double[1]));
    gkyl_cu_memset(gkf->em_energy_red_old, 0, sizeof(double[1]));
  } else {
    gkf->em_energy_red_new = gkyl_malloc(sizeof(double[1]));
    gkf->em_energy_red_old = gkyl_malloc(sizeof(double[1]));
    memset(gkf->em_energy_red_new, 0, sizeof(double[1]));
    memset(gkf->em_energy_red_old, 0, sizeof(double[1]));
  }
  gkf->integ_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  gkf->is_first_energy_dot_write_call = true;
}

static void
gk_field_energy_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // Allocate energy reduction arrays.
  if (app->use_gpu) {
    gkf->em_energy_red = gkyl_cu_malloc(sizeof(double[1]));
    gkf->em_energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    gkf->em_energy_red = gkyl_malloc(sizeof(double[1]));
    gkf->em_energy_red_global = gkyl_malloc(sizeof(double[1]));
  }

  gkf->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  gkf->is_first_energy_write_call = true;

  gkf->calc_energy_func = gk_field_calc_energy_enabled;
  gkf->calc_energy_dt_func = gk_field_calc_energy_dt_disabled;
  
  if (gkf->info.time_rate_diagnostics) {
    gk_field_time_rate_diags_new(app, gkf);
  }

  // Factors for ES energy.
  gkf->es_energy_fac = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  gkf->es_energy_fac_1d = 0.0;
}

static void
gk_field_time_rate_diags_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkf->calc_energy_dt_func = gk_field_calc_energy_dt_disabled;
  if (app->use_gpu) {
    gkyl_cu_free(gkf->em_energy_red_new);
    gkyl_cu_free(gkf->em_energy_red_old);
  } else {
    gkyl_free(gkf->em_energy_red_new);
    gkyl_free(gkf->em_energy_red_old);
  }
  gkyl_dynvec_release(gkf->integ_energy_dot);
}

static void
gk_field_energy_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_dynvec_release(gkf->integ_energy);

  if (app->use_gpu) {
    gkyl_cu_free(gkf->em_energy_red);
    gkyl_cu_free(gkf->em_energy_red_global);
  } else {
    gkyl_free(gkf->em_energy_red);
    gkyl_free(gkf->em_energy_red_global);
  }

  if (gkf->info.time_rate_diagnostics) {
    gk_field_time_rate_diags_release(app, gkf);
  }

  gkyl_array_release(gkf->es_energy_fac);
}

// Related to enforcing parallel Vlasov boundary conditions
void
gk_field_enforce_parallel_bc_disabled(const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *finout)
{
  // Do nothing.
}

// Initialize field object.
struct gk_field* 
gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app)
{
  struct gk_field *gkf = gkyl_malloc(sizeof(struct gk_field));

  gkf->info = gk->field;

  gkf->gkfield_id = gkf->info.gkfield_id ? gkf->info.gkfield_id : GKYL_GK_FIELD_ES;
  
  gkf->calc_init_field = !gkf->info.zero_init_field;
  gkf->update_field = !gkf->info.is_static;
  // The combination update_field=true, calc_init_field=false is not allowed.
  assert(!(gkf->update_field && (!gkf->calc_init_field)));

  gkf->num_fields = 1;
  if (gkf->gkfield_id == GKYL_GK_FIELD_EM_APAR)
    gkf->num_fields = 2;

  // Allocate array of field pointers.
  gkf->f = gkyl_malloc(gkf->num_fields * sizeof(struct gkyl_array *));
  gkf->f1 = gkyl_malloc(gkf->num_fields * sizeof(struct gkyl_array *));
  gkf->fnew = gkyl_malloc(gkf->num_fields * sizeof(struct gkyl_array *));

  // Initialize polarization potential if needed.
  gkf->init_phi_pol = false;
  if (gkf->info.polarization_potential) {
    gk_field_polarization_potential_new(gkf, app);
  }
  
  // Initialize energy diagnostics.
  gk_field_energy_new(app, gkf);

  // Initialize the field solver
  if (gkf->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    gk_field_fem_new_boltzmann(app, gkf);
  } else {
    if (app->cdim == 1) {
      gk_field_fem_new_1x(app, gkf);
    } else if (app->cdim > 1) {
      gk_field_fem_new_2x3x(app, gkf);
    }
  }

  // Initialize biased walls.
  gk_field_biased_wall_new(app, gkf);

  return gkf;
}

void
gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *gkf, struct gkyl_array *emfields[])
{
  struct timespec wst = gkyl_wall_clock();

  // Compute the electrostatic potential.
  struct gkyl_array *phi_out = gk_field_get_phi_from_fields(gkf, emfields);
  gkf->rhs_phi_func(app, gkf, phi_out);

  app->stat.field_phi_solve_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *gkf)
{
  struct timespec wst = gkyl_wall_clock();
  gkf->calc_energy_func(app, gkf, tm);
  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *gkf,
  double dt, struct gkyl_array *fields[], double *energy_reduced)
{
  gkf->calc_energy_dt_func(app, gkf, dt, fields, energy_reduced);
}

void gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *gkf,
  struct gk_species *s, struct gkyl_array **bflux)
{
  // Gyroaverage the density if needed.
  s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);
  gkyl_array_accumulate_range(gkf->rho_c, s->info.charge, s->m0_gyroavg, &app->local);
  // Add the background (electron) charge density.
  double n_s0 = gkf->info.electron_density;
  double q_s = gkf->info.electron_charge;
  double dg_norm = pow(sqrt(2), app->basis.ndim);
  gkyl_array_shiftc_range(gkf->rho_c, q_s * n_s0 * dg_norm, 0, &app->local);
}

void gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app, struct gk_field *gkf,
  struct gk_species *s, struct gkyl_array **bflux)
{
  // Gyroaverage the density if needed.
  s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);
  gkyl_array_accumulate_range(gkf->rho_c, s->info.charge, s->m0_gyroavg, &app->local);
}

void
gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *gkf, 
  struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(gkf->rho_c, 0.0);
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    gk_species_moment_calc(app, &s->m0, &s->local, &app->local, 0, 0, 0, fin[i]);
    gkf->accumulate_rhoc_func(app, gkf, s, bflux[i]);
  } 
  app->stat.field_phi_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *gkf,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x).

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, gkf->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(gkf->fem_parproj, gkf->rho_c_global_dg, gkf->rho_c_global_dg);
  gkyl_fem_parproj_solve(gkf->fem_parproj, gkf->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, gkf->phi_fem, &app->local, &gkf->global_sub_range);
}

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = GKYL_ARRAY_RIO_FOPEN_FAILED };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gkyl_msgpack_map_elem elem_list[] = {
      { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
      { .key = "time", .elem_type = GKYL_MP_DOUBLE, .cval = 0 },
    };
    int elem_list_len = sizeof(elem_list)/sizeof(elem_list[0]);
    gkyl_msgpack_to_map_elem_list(&(struct gkyl_msgpack_data) {
        .meta = hdr.meta,
        .meta_sz = hdr.meta_size
      }, elem_list_len, elem_list);

    rstat.frame = gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "frame");
    rstat.stime = gkyl_msgpack_map_elem_get_double(elem_list_len, elem_list, "time");

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}

void
gk_field_file_import_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file.
  struct gkyl_app_restart_status rstat = header_from_file(app, inp.file_name);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    struct gkyl_app_restart_status rstat;
    rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, inp.file_name);
    gkyl_array_copy(app->field->phi, app->field->phi_host);
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
  gkyl_array_copy(app->field->phi, app->field->phi_host);
}


// Release resources for field.
void
gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *gkf)
{
  // Release polarization potential if used.
  if (gkf->init_phi_pol) {
    gk_field_polarization_potential_release(gkf);
  }

  // Release solver-specific resources.
  gkf->solver_release_func(app, gkf);

  // Release energy diagnostics.
  gk_field_energy_release(app, gkf);

  // Release biased walls.
  gk_field_biased_wall_release(app, gkf);

  gkyl_free(gkf->f);
  gkyl_free(gkf->f1);
  gkyl_free(gkf->fnew);

  gkyl_free(gkf);
}
