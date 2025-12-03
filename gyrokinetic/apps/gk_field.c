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

  // Initialize polarization potential if needed.
  f->init_phi_pol = false;
  if (f->info.polarization_potential) {
    gk_field_polarization_potential_new(f, app);
  }
  
  // Initialize energy diagnostics.
  gk_field_energy_new(app, f);

  // Initialize the field solver
  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    gk_field_fem_new_boltzmann(app, f);
  } else {
    if (app->cdim == 1) {
      gk_field_fem_new_1x(app, f);
    } else if (app->cdim > 1) {
      gk_field_fem_new_2x3x(app, f);
    }
  }

  // Initialize biased walls.
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
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();
  field->calc_energy_func(app, field, tm);
  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  field->calc_energy_dt_func(app, field, dt, energy_reduced);
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

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = 2 };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gyrokinetic_output_meta meta =
      gk_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }, 
        GKYL_GK_META_NONE,
        0
      );

    rstat.frame = meta.frame;
    rstat.stime = meta.stime;

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


// Release resources for field.
void
gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f)
{

  // Release polarization potential if used.
  if (f->init_phi_pol) {
    gk_field_polarization_potential_release(f);
  }

  // Release solver-specific resources.
  f->solver_release_func(app, f);

  // Release energy diagnostics.
  gk_field_energy_release(app, f);

  // Release biased walls.
  gk_field_biassed_wall_release(app, f);

  gkyl_free(f);
}
