#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_translate_dim.h>

#include <assert.h>
#include <time.h>

// initialize species grid and arrays
// initialize species object

static double
gk_neut_species_rhs_dynamic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  double omega_cfl = 1/DBL_MAX;   
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);
  
  // Collisionless terms.
  struct timespec wst = gkyl_wall_clock();
  gkyl_dg_updater_vlasov_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);
  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);

  if (species->bgk.collision_id == GKYL_BGK_COLLISIONS && !app->has_implicit_coll_scheme) {
    gk_neut_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  }

  if (species->react_neut.num_react) {
    gk_neut_species_react_rhs(app, species, &species->react_neut, fin, rhs);
  }
  
  // Compute and store (in the ghost cell of rhs) the boundary fluxes.
  gk_neut_species_bflux_rhs(app, &species->bflux, fin, rhs);

  // Compute diagnostic moments of the boundary fluxes.
  gk_neut_species_bflux_calc_moms(app, &species->bflux, rhs, bflux_moms);
  
  app->stat.n_neut_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);
  
  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];
  
  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl/omega_cfl;
}

static double
gk_neut_species_rhs_implicit_dynamic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{ 
  double omega_cfl = 1/DBL_MAX;
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute implicit update and update rhs to new time step
  if (species->bgk.collision_id == GKYL_BGK_COLLISIONS) {
    gk_neut_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  }
  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);
  
  app->stat.n_neut_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);
  
  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];
  
  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl/omega_cfl;
}

static double
gk_neut_species_rhs_static(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  double omega_cfl = 1/DBL_MAX;
  return app->cfl/omega_cfl;
}

static double
gk_neut_species_rhs_implicit_static(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  double omega_cfl = 1/DBL_MAX;
  return app->cfl/omega_cfl;
}

static void
gk_neut_species_apply_bc_dynamic(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, app->periodic_dirs, f); 
  
  for (int d=0; d<cdim; ++d) {
    if (species->bc_is_np[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_RECYCLE:
          gk_neut_species_recycle_apply_bc(app, &species->bc_recycle_lo, species, f);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_SPECIES_RECYCLE:
          gk_neut_species_recycle_apply_bc(app, &species->bc_recycle_up, species, f);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_apply_bc_static(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  // empty function
}

static void
gk_neut_species_step_f_dynamic(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

static void
gk_neut_species_step_f_static(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // do nothing
}

static void
gk_neut_species_combine_dynamic(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

static void
gk_neut_species_combine_static(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // do nothing
}

static void
gk_neut_species_copy_range_dynamic(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  gkyl_array_copy_range(out, inp, range);
}

static void
gk_neut_species_copy_range_static(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // do nothing
}

static void
gk_neut_species_apply_pos_shift_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns)
{
  struct timespec wtm = gkyl_wall_clock();
  // Copy f so we can calculate the moments of the change later. 
  gkyl_array_set(gkns->fnew, -1.0, gkns->f);

  // Shift each species.
  gkyl_positivity_shift_vlasov_advance(gkns->pos_shift_op, &app->local, &gkns->local,
    gkns->f, gkns->m0.marr, gkns->ps_delta_m0);
  app->stat.neut_species_pos_shift_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gk_neut_species_apply_pos_shift_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns)
{
}

void
gk_neut_species_apply_pos_shift(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns)
{
  gkns->apply_pos_shift_func(app, gkns);
}

static void
gk_neut_species_write_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = gkns->basis.id
    }, GKYL_GK_META_NONE, 0
  );

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, frame);
  
  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(gkns->f_host, gkns->f);
  }
  
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(gkns->comm, &gkns->grid, &gkns->local, mt, gkns->f_host, fileNm);
  app->stat.neut_species_io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.n_neut_io += 1;
  
  gk_array_meta_release(mt);  
}

static void
gk_neut_species_write_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  // do nothing
}

static void
gk_neut_species_write_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{

  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }, GKYL_GK_META_NONE, 0
  );

  for (int m=0; m<gkns->info.num_diag_moments; ++m) {
    struct timespec wst = gkyl_wall_clock();
    gk_neut_species_moment_calc(&gkns->moms[m], gkns->local, app->local, gkns->f);
    app->stat.n_neut_mom += 1;

    // Rescale moment by inverse of Jacobian. 
    // For LTE (Maxwellian) moments, we only need to re-scale
    // the density (the 0th component).
    gkyl_dg_div_op_range(gkns->moms[m].mem_geo, app->basis, 
      0, gkns->moms[m].marr, 0, gkns->moms[m].marr, 0, 
      app->gk_geom->geo_int.jacobgeo, &app->local);      
    app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

    struct timespec wtm = gkyl_wall_clock();
    if (app->use_gpu)
      gkyl_array_copy(gkns->moms[m].marr_host, gkns->moms[m].marr);

    const char *fmt = "%s-%s_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, 
      gkyl_distribution_moments_strs[gkns->info.diag_moments[m]], frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name,
      gkyl_distribution_moments_strs[gkns->info.diag_moments[m]], frame);

    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gkns->moms[m].marr_host, fileNm);
    app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_neut_diag_io += 1;
  }

  if (gkns->enforce_positivity) {
    struct timespec wst = gkyl_wall_clock();
    // We placed the change in f from the positivity shift in fnew.
    gk_neut_species_moment_calc(&gkns->ps_moms, gkns->local, app->local, gkns->fnew);
    app->stat.n_mom += 1;

    // Rescale moment by inverse of Jacobian.
    gkyl_dg_div_op_range(gkns->ps_moms.mem_geo, app->basis,
      0, gkns->ps_moms.marr, 0, gkns->ps_moms.marr, 0,
      app->gk_geom->geo_int.jacobgeo, &app->local);
    app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

    struct timespec wtm = gkyl_wall_clock();
    if (app->use_gpu)
      gkyl_array_copy(gkns->ps_moms.marr_host, gkns->ps_moms.marr);

    const char *fmt = "%s-%s_positivity_shift_M0_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, frame);

    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gkns->ps_moms.marr_host, fileNm);
    app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_neut_diag_io += 1;
  }
  gk_array_meta_release(mt); 
}

static void
gk_neut_species_write_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  // do nothing
}

static void
gk_neut_species_calc_integrated_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  int vdim = app->vdim+1; // Neutrals are always 3V
  int num_mom = gkns->integ_moms.num_mom;
  double avals_global[num_mom];
  
  gk_neut_species_moment_calc(&gkns->integ_moms, gkns->local, app->local, gkns->f); 
  app->stat.n_neut_mom += 1;

  // reduce to compute sum over whole domain, append to diagnostics
  gkyl_array_reduce_range(gkns->red_integ_diag, gkns->integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
    gkns->red_integ_diag, gkns->red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, gkns->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, gkns->red_integ_diag_global, sizeof(double[num_mom]));
  }
  gkyl_dynvec_append(gkns->integ_diag, tm, avals_global);

  if (gkns->enforce_positivity) {
    // The change in f from the positivity shift is in fnew.
    gk_neut_species_moment_calc(&gkns->integ_moms, gkns->local, app->local, gkns->fnew);
    app->stat.n_mom += 1;

    // Reduce (sum) over whole domain, append to diagnostics.
    gkyl_array_reduce_range(gkns->red_integ_diag, gkns->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom,
      gkns->red_integ_diag, gkns->red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gkns->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gkns->red_integ_diag_global, sizeof(double[num_mom]));
    }
    gkyl_dynvec_append(gkns->ps_integ_diag, tm, avals_global);
  }

  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_neut_diag += 1;
}

static void
gk_neut_species_calc_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm)
{
  // do nothing
}

static void
gk_neut_species_write_integrated_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns)
{
  struct timespec wst = gkyl_wall_clock();
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, "integrated_moms");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, "integrated_moms");
    
    if (gkns->is_first_integ_write_call) {
      gkyl_dynvec_write(gkns->integ_diag, fileNm);
      gkns->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(gkns->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(gkns->integ_diag);
  app->stat.n_neut_diag_io += 1;

  if (gkns->enforce_positivity) {
    if (rank == 0) {
      // Write integrated diagnostic moments.
      const char *fmt = "%s-%s_positivity_shift_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, "integrated_moms");

      if (gkns->is_first_ps_integ_write_call) {
        gkyl_dynvec_write(gkns->ps_integ_diag, fileNm);
        gkns->is_first_ps_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gkns->ps_integ_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(gkns->ps_integ_diag);
    app->stat.n_diag_io += 1;
  }

  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_write_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns)
{
  // do nothing
}

// release all resources for dynamic species
static void
gk_neut_species_release_dynamic(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s)
{
  // release various arrays
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omega_cfl);
  }
  else {
    gkyl_free(s->omega_cfl);
  }

  gk_neut_species_source_release(app, &s->src);

  if (s->bgk.collision_id == GKYL_BGK_COLLISIONS) {
    gk_neut_species_bgk_release(app, &s->bgk);
  }
  if (s->react_neut.num_react) {
    gk_neut_species_react_release(app, &s->react_neut);
  }

  // Release integrated mom data.
  gk_neut_species_moment_release(app, &s->integ_moms); 

  // Release integrated mom diag data.
  gkyl_dynvec_release(s->integ_diag);
  
  if (app->use_gpu) {
    gkyl_cu_free(s->red_integ_diag);
    gkyl_cu_free(s->red_integ_diag_global);
  }
  else {
    gkyl_free(s->red_integ_diag);
    gkyl_free(s->red_integ_diag_global);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d].type == GKYL_SPECIES_RECYCLE) {
      gk_neut_species_recycle_release(app, &s->bc_recycle_lo);
    }
    else { 
      gkyl_bc_basic_release(s->bc_lo[d]);
    }
    if (s->upper_bc[d].type == GKYL_SPECIES_RECYCLE) {
      gk_neut_species_recycle_release(app, &s->bc_recycle_up);
    }
    else {
      gkyl_bc_basic_release(s->bc_up[d]);
    }
  }

  if (s->enforce_positivity) {
    gkyl_array_release(s->ps_delta_m0);
    gkyl_positivity_shift_vlasov_release(s->pos_shift_op);
    gk_species_moment_release(app, &s->ps_moms);
    gkyl_dynvec_release(s->ps_integ_diag);
  }

}

// empty function for static species
static void
gk_neut_species_release_static(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s)
{ 
}

static void
gk_neut_species_new_dynamic(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  int cdim = app->cdim, vdim = app->vdim+1; // neutral species are 3v
  int pdim = cdim+vdim;
  
  // allocate additional distribution function arrays for time stepping
  s->f1 = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
  
  // Allocate cflrate (scalar array).
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);

  if (app->use_gpu) {
    s->omega_cfl = gkyl_cu_malloc(sizeof(double));
  }
  else {
    s->omega_cfl = gkyl_malloc(sizeof(double));
  }

  // Allocate data for integrated moments.
  gk_neut_species_moment_init(app, s, &s->integ_moms, GKYL_F_MOMENT_M0M1M2, true);

  // Allocate data for integrated diagnostics.
  if (app->use_gpu) {
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[s->integ_moms.num_mom]));
    s->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[s->integ_moms.num_mom]));
  } else {
    s->red_integ_diag = gkyl_malloc(sizeof(double[s->integ_moms.num_mom]));
    s->red_integ_diag_global = gkyl_malloc(sizeof(double[s->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, s->integ_moms.num_mom);
  s->is_first_integ_write_call = true;

  // Determine collision type and initialize it.
  if (s->info.collisions.collision_id == GKYL_BGK_COLLISIONS) {
    gk_neut_species_bgk_init(app, s, &s->bgk);
  }
  // Determine reaction type(s) and initialize them. 
  if (s->info.react_neut.num_react) {
    gk_neut_species_react_init(app, s, s->info.react_neut, &s->react_neut);
  }

  // Allocate buffer needed for BCs.
  long buff_sz = 0;
  for (int dir=0; dir<cdim; ++dir) {
    long vol = GKYL_MAX2(s->lower_skin[dir].volume, s->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->use_gpu, s->basis.num_basis, buff_sz);
  // Buffer arrays for fixed function boundary conditions on distribution function.
  s->bc_buffer_lo_fixed = mkarr(app->use_gpu, s->basis.num_basis, buff_sz);
  s->bc_buffer_up_fixed = mkarr(app->use_gpu, s->basis.num_basis, buff_sz);
  
  s->bc_recycle_lo = (struct gk_recycle_wall) {};
  s->bc_recycle_up = (struct gk_recycle_wall) {};
  for (int d=0; d<cdim; ++d) {
    // Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (s->lower_bc[d].type == GKYL_SPECIES_RECYCLE) {
      // Initialize recycling object. 
      gk_neut_species_recycle_init(app, &s->bc_recycle_lo, d, GKYL_LOWER_EDGE, &s->lower_bc[d].emission, s, app->use_gpu);
    }
    else { 
      if (s->lower_bc[d].type == GKYL_SPECIES_COPY) {
        bctype = GKYL_BC_COPY;
      }
      else if (s->lower_bc[d].type == GKYL_SPECIES_ABSORB) {
        bctype = GKYL_BC_ABSORB;
      }
      else if (s->lower_bc[d].type == GKYL_SPECIES_REFLECT) {
        bctype = GKYL_BC_DISTF_REFLECT;
      }
      else if (s->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        bctype = GKYL_BC_FIXED_FUNC;
      }

      s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, s->basis_on_dev,
        &s->lower_skin[d], &s->lower_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);

      if (s->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_lo;
        gk_neut_species_projection_init(app, s, s->lower_bc[d].projection, &gk_proj_bc_lo);
        gk_neut_species_projection_calc(app, s, &gk_proj_bc_lo, s->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(s->bc_lo[d], s->bc_buffer_lo_fixed, s->f1);
        gkyl_array_clear(s->f1, 0.0);
        gk_neut_species_projection_release(app, &gk_proj_bc_lo);
      }
    }

    if (s->upper_bc[d].type == GKYL_SPECIES_RECYCLE) {
      // Initialize recycling object. 
      gk_neut_species_recycle_init(app, &s->bc_recycle_up, d, GKYL_UPPER_EDGE, &s->upper_bc[d].emission, s, app->use_gpu);
    }
    else {
      // Upper BC updater. Copy BCs by default.
      if (s->upper_bc[d].type == GKYL_SPECIES_COPY) {
        bctype = GKYL_BC_COPY;
      }
      else if (s->upper_bc[d].type == GKYL_SPECIES_ABSORB) {
        bctype = GKYL_BC_ABSORB;
      }
      else if (s->upper_bc[d].type == GKYL_SPECIES_REFLECT) {
        bctype = GKYL_BC_DISTF_REFLECT;
      }
      else if (s->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        bctype = GKYL_BC_FIXED_FUNC;
      }

      s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, s->basis_on_dev,
        &s->upper_skin[d], &s->upper_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);

      if (s->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_up;
        gk_neut_species_projection_init(app, s, s->upper_bc[d].projection, &gk_proj_bc_up);
        gk_neut_species_projection_calc(app, s, &gk_proj_bc_up, s->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(s->bc_up[d], s->bc_buffer_up_fixed, s->f1);
        gkyl_array_clear(s->f1, 0.0);
        gk_neut_species_projection_release(app, &gk_proj_bc_up);
      }
    }
  }

  if (app->enforce_positivity || s->info.enforce_positivity) {
    s->enforce_positivity = true;

    // Positivity enforcing by shifting f (ps=positivity shift).
    s->ps_delta_m0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    s->pos_shift_op = gkyl_positivity_shift_vlasov_new(app->basis, s->basis,
      s->grid, &app->local_ext, app->use_gpu);

    // Allocate data for diagnostic moments
    gk_neut_species_moment_init(app, s, &s->ps_moms, GKYL_F_MOMENT_M0, false);

    s->ps_integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, s->integ_moms.num_mom);
    s->is_first_ps_integ_write_call = true;
  }

  // Set function pointers
  s->rhs_func = gk_neut_species_rhs_dynamic;
  s->rhs_implicit_func = gk_neut_species_rhs_implicit_dynamic;
  s->bc_func = gk_neut_species_apply_bc_dynamic;
  s->release_func = gk_neut_species_release_dynamic;
  s->step_f_func = gk_neut_species_step_f_dynamic;
  s->combine_func = gk_neut_species_combine_dynamic;
  s->copy_func = gk_neut_species_copy_range_dynamic;
  if (s->enforce_positivity)
    s->apply_pos_shift_func = gk_neut_species_apply_pos_shift_enabled;
  else
    s->apply_pos_shift_func = gk_neut_species_apply_pos_shift_disabled;
  s->write_func = gk_neut_species_write_dynamic;
  s->write_mom_func = gk_neut_species_write_mom_dynamic;
  s->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_dynamic;
  s->write_integrated_mom_func = gk_neut_species_write_integrated_mom_dynamic;
}

static void
gk_neut_species_new_static(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  // Set pointers for RK methods.
  s->f1 = s->f;
  s->fnew = s->f;

  // Set function pointers
  s->rhs_func = gk_neut_species_rhs_static;
  s->rhs_implicit_func = gk_neut_species_rhs_implicit_static;
  s->bc_func = gk_neut_species_apply_bc_static;
  s->release_func = gk_neut_species_release_static;
  s->step_f_func = gk_neut_species_step_f_static;
  s->combine_func = gk_neut_species_combine_static;
  s->copy_func = gk_neut_species_copy_range_static;
  s->apply_pos_shift_func = gk_neut_species_apply_pos_shift_disabled;
  s->write_func = gk_neut_species_write_static;
  s->write_mom_func = gk_neut_species_write_mom_static;
  s->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_static;
  s->write_integrated_mom_func = gk_neut_species_write_integrated_mom_static;
}

void
gk_neut_species_file_import_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file. Intended options include importing:
  //   1) ICs with same grid.
  //   2) ICs one dimensionality lower (e.g. 2x2v for 3x2v sim).
  //   3) ICs with same grid extents but different resolution (NYI).

  struct gkyl_rect_grid grid = s->grid;
  int pdim = grid.ndim;
  int cdim = app->cdim;
  int vdim = pdim - cdim;
  int poly_order = app->poly_order;

  struct gkyl_rect_grid grid_do; // Donor grid.
  struct gkyl_array_header_info hdr;
  int pdim_do, vdim_do, cdim_do;
  bool same_res = true;

  // Read the header of the input file, extract needed info an create a grid
  // and other things needed.
  FILE *fp;
  with_file(fp, inp.file_name, "r") {

    int status = gkyl_grid_sub_array_header_read_fp(&grid_do, &hdr, fp);

    pdim_do = grid_do.ndim;
    vdim_do = vdim; // Assume velocity space dimensionality is the same.
    cdim_do = pdim_do - vdim_do;

    // Perform some basic checks.
    if (pdim_do == pdim) {
      // Check if the grid resolution is the same.
      for (int d=0; d<pdim; d++)
        same_res = same_res && (grid_do.cells[d] == grid.cells[d]);
    }
    else {
      // Assume the loaded file has one lower conf-space dimension.
      // Primarily meant for loading:
      //   - 1x3v for a 2x3v sim.
      //   - 2x3v for a 3x3v sim.
      assert(pdim_do == pdim-1);
      double tol_eq = 1e-10;
      for (int d=0; d<cdim_do-1; d++) {
        assert(gkyl_compare_double(grid_do.lower[d], grid.lower[d], tol_eq));
        assert(gkyl_compare_double(grid_do.upper[d], grid.upper[d], tol_eq));
        assert(gkyl_compare_double(grid_do.dx[d]   , grid.dx[d]   , tol_eq));
        assert(grid_do.cells[d] == grid.cells[d]);
      }
      assert(gkyl_compare_double(grid_do.lower[cdim_do-1], grid.lower[cdim-1], tol_eq));
      assert(gkyl_compare_double(grid_do.upper[cdim_do-1], grid.upper[cdim-1], tol_eq));
      assert(gkyl_compare_double(grid_do.dx[cdim_do-1]   , grid.dx[cdim-1]   , tol_eq));
      assert(grid_do.cells[cdim_do-1] == grid.cells[cdim-1]);
      for (int d=0; d<vdim; d++) {
        assert(gkyl_compare_double(grid_do.lower[cdim_do+d], grid.lower[cdim+d], tol_eq));
        assert(gkyl_compare_double(grid_do.upper[cdim_do+d], grid.upper[cdim+d], tol_eq));
        assert(gkyl_compare_double(grid_do.dx[cdim_do+d]   , grid.dx[cdim+d]   , tol_eq));
        assert(grid_do.cells[cdim_do+d] == grid.cells[cdim+d]);
      }
    }

    struct gyrokinetic_output_meta meta =
      gk_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }, GKYL_GK_META_NONE, 0
      );
    assert(strcmp(s->basis.id, meta.basis_type_nm) == 0);
    assert(poly_order == meta.poly_order);
    gkyl_grid_sub_array_header_release(&hdr);
  }

  // Donor basis.
  struct gkyl_basis basis_do;
  // Basis is tensor for p=1 and ser for p>1
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis_do, pdim_do, poly_order);
  }
  else if (poly_order == 1) {
    gkyl_cart_modal_tensor(&basis_do, pdim_do, poly_order); // for canonical PB
  }

  // Donor global range.
  int ghost_do[pdim_do];
  for (int d=0; d<cdim_do; d++) ghost_do[d] = 1;
  for (int d=0; d<vdim_do; d++) ghost_do[cdim_do+d] = 0;
  struct gkyl_range global_ext_do, global_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &global_ext_do, &global_do);

  // Create a donor communicator.
  int cuts_tar[GKYL_MAX_DIM] = {-1}, cuts_do[GKYL_MAX_CDIM] = {-1};
  gkyl_rect_decomp_get_cuts(app->decomp, cuts_tar);
  if (cdim_do == cdim-1) {
    for (int d=0; d<cdim_do-1; d++) {
      cuts_do[d] = cuts_tar[d];
    }
    cuts_do[cdim_do-1] = cuts_tar[cdim-1];
  }
  else {
    for (int d=0; d<cdim; d++) {
      cuts_do[d] = cuts_tar[d];
    }
  }
  // Set velocity space cuts to 1 as we do not use MPI in vel-space.
  for (int d=0; d<vdim; d++) {
    cuts_do[cdim_do+d] = 1;
  }

  struct gkyl_rect_decomp *decomp_do = gkyl_rect_decomp_new_from_cuts(pdim_do, cuts_do, &global_do);
  struct gkyl_comm* comm_do = gkyl_comm_split_comm(s->comm, 0, decomp_do);

  // Donor local range.
  int my_rank = 0;
  gkyl_comm_get_rank(comm_do, &my_rank);

  struct gkyl_range local_ext_do, local_do;
  gkyl_create_ranges(&decomp_do->ranges[my_rank], ghost_do, &local_ext_do, &local_do);

  // Donor array.
  struct gkyl_array *fdo = mkarr(app->use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_host = app->use_gpu? mkarr(false, basis_do.num_basis, local_ext_do.volume)
                                            : gkyl_array_acquire(fdo);

  // Read donor distribution function and Jacobian inverse.
  struct gkyl_app_restart_status rstat;
  rstat.io_status = gkyl_comm_array_read(comm_do, &grid_do, &local_do, fdo_host, inp.file_name);

  bool scale_by_jacobgeo = false;
  with_file(fp, inp.jacobtot_inv_file_name, "r") {
    // Set up configuration space donor grid and basis
    struct gkyl_rect_grid conf_grid_do;
    gkyl_rect_grid_init(&conf_grid_do, cdim_do, grid_do.lower, grid_do.upper, grid_do.cells);
    // Create configuration space global ranges
    struct gkyl_range conf_local_ext_do, conf_local_do, conf_global_ext_do, conf_global_do;
    gkyl_create_grid_ranges(&conf_grid_do, ghost_do, &conf_global_ext_do, &conf_global_do);
    // Create configuration space local ranges
    struct gkyl_rect_decomp *conf_decomp_do = gkyl_rect_decomp_new_from_cuts(cdim_do, cuts_do, &conf_global_do);
    gkyl_create_ranges(&conf_decomp_do->ranges[my_rank], ghost_do, &conf_local_ext_do, &conf_local_do);
    // Create a configuration space basis.
    struct gkyl_basis conf_basis_do;
    gkyl_cart_modal_serendip(&conf_basis_do, cdim_do, poly_order);
    // Array for Jacobian inverse
    struct gkyl_array *jacobgeo_inv_do_host = mkarr(false, conf_basis_do.num_basis, conf_local_ext_do.volume);
    rstat.io_status = gkyl_comm_array_read(comm_do, &conf_grid_do, &conf_local_do, jacobgeo_inv_do_host, inp.jacobtot_inv_file_name);
    gkyl_dg_mul_conf_phase_op_range(&conf_basis_do, &basis_do, fdo_host, jacobgeo_inv_do_host, fdo_host, &conf_local_ext_do, &local_ext_do);
    gkyl_array_release(jacobgeo_inv_do_host);
    gkyl_rect_decomp_release(conf_decomp_do);
    scale_by_jacobgeo = true;
  }

  if (app->use_gpu) {
    gkyl_array_copy(fdo, fdo_host);
  }

  if (pdim_do == pdim-1) {
    struct gkyl_translate_dim* transdim = gkyl_translate_dim_new(cdim_do,
      basis_do, cdim, s->basis, -1, GKYL_NO_EDGE, app->use_gpu);
    gkyl_translate_dim_advance(transdim, &local_do, &s->local, fdo, 1, s->f);
    gkyl_translate_dim_release(transdim);
  }
  else {
    if (same_res) {
      gkyl_array_copy(s->f, fdo);
    }
    else {
      // Interpolate the donor distribution to the target grid.
      struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(app->cdim, &s->basis,
        &grid_do, &grid, &local_do, &s->local, ghost_do, app->use_gpu);
      gkyl_dg_interpolate_advance(interp, fdo, s->f);
      gkyl_dg_interpolate_release(interp);
    }
  }

  if (inp.type == GKYL_IC_IMPORT_AF || inp.type == GKYL_IC_IMPORT_AF_B) {
    // Scale f by a conf-space factor.
    gkyl_proj_on_basis *proj_conf_scale = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      poly_order+1, 1, inp.conf_scale, inp.conf_scale_ctx);
    struct gkyl_array *xfac = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *xfac_ho = app->use_gpu? mkarr(false, app->basis.num_basis, app->local_ext.volume)
                                             : gkyl_array_acquire(xfac_ho);
    gkyl_proj_on_basis_advance(proj_conf_scale, 0.0, &app->local, xfac_ho);
    gkyl_array_copy(xfac, xfac_ho);
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &s->basis, s->f, xfac, s->f, &app->local, &s->local);
    gkyl_proj_on_basis_release(proj_conf_scale);
    gkyl_array_release(xfac_ho);
    gkyl_array_release(xfac);
  }
  if (inp.type == GKYL_IC_IMPORT_F_B || inp.type == GKYL_IC_IMPORT_AF_B) {
    // Add a phase factor to f.
    struct gk_proj proj_phase_add;
    gk_neut_species_projection_init(app, s, inp.phase_add, &proj_phase_add);
    gk_neut_species_projection_calc(app, s, &proj_phase_add, s->fnew, 0.0);
    gkyl_array_accumulate_range(s->f, 1.0, s->fnew, &s->local);
    gk_neut_species_projection_release(app, &proj_phase_add);
  }

  // Multiply f by the Jacobian.
  if (scale_by_jacobgeo)
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &s->basis, s->f, app->gk_geom->geo_int.jacobgeo, s->f, &app->local, &s->local);

  gkyl_rect_decomp_release(decomp_do);
  gkyl_comm_release(comm_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_host);
}

static void
gkyl_array_move_comp(struct gkyl_array *out, int cout, struct gkyl_array *in, int cin, struct gkyl_array *tmp)
{
  // Move the 'cin' scalar field from the multi-component (multi-scalar field) array
  // 'in' to 'cin' component of the the multi-component 'out'. This requires a
  // temporary array 'tmp'. This is really a work around for the fact that some
  // of our array offset ops don't work as we'd like for multi-component arrays
  // with the same number of components.
  gkyl_array_set_offset(tmp, 1.0,  in, cin);
  gkyl_array_set_offset(out, 1.0, tmp, cout);
}

void
gk_neut_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  int cdim = app->cdim, vdim = app->vdim+1; // neutral species are 3v
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = gk->cells[d];
    lower[d] = gk->lower[d];
    upper[d] = gk->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = s->info.cells[d];
    lower_vel[d] = s->info.lower[d];
    upper_vel[d] = s->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }

  // Basis is tensor for p=1 and ser for p>1
  if (app->poly_order > 1) {
    gkyl_cart_modal_serendip(&s->basis, pdim, app->poly_order);
  }
  else if (app->poly_order == 1) {
    gkyl_cart_modal_tensor(&s->basis, pdim, app->poly_order); // for canonical PB
  }

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    s->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (app->poly_order > 1) {
      gkyl_cart_modal_serendip_cu_dev(s->basis_on_dev, pdim, app->poly_order);
    }
    else if (app->poly_order == 1) {
      gkyl_cart_modal_tensor_cu_dev(s->basis_on_dev, pdim, app->poly_order); // for canonical PB
    }
  }
  else {
    s->basis_on_dev = &s->basis;
  }
  
  // Full phase space grid.
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->global_ext, &s->global);
  
  // Velocity space grid.
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

  // Phase-space communicator.
  s->comm = gkyl_comm_extend_comm(app->comm, &s->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // local = conf-local X local_vel.
  gkyl_range_ten_prod(&local, &app->local, &s->local_vel);
  gkyl_create_ranges(&local, ghost, &s->local_ext, &s->local);

  s->model_id = GKYL_MODEL_CANONICAL_PB;

  // Velocity space mapping.
  assert(s->info.mapc2p.mapping == 0); // mapped v-space not implemented for neutrals yet.
  s->vel_map = gkyl_velocity_map_new(s->info.mapc2p, s->grid, s->grid_vel,
    s->local, s->local_ext, s->local_vel, s->local_ext_vel, app->use_gpu);

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // modify it in GK_IWL BCs without modifying the app's.
  s->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<s->num_periodic_dir; ++d)
    s->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) s->bc_is_np[d] = true;
  for (int d=0; d<s->num_periodic_dir; ++d)
    s->bc_is_np[s->periodic_dirs[d]] = false;

  // Store the BCs from the input file.
  for (int dir=0; dir<app->cdim; ++dir) {
    s->lower_bc[dir].type = s->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (s->bc_is_np[dir]) {
      const struct gkyl_gyrokinetic_bcs *bc;
      if (dir == 0)
        bc = &s->info.bcx;
      else if (dir == 1)
        bc = &s->info.bcy;
      else
        bc = &s->info.bcz;

      s->lower_bc[dir] = bc->lower;
      s->upper_bc[dir] = bc->upper;
    }
  }

  // Determine which directions are zero-flux. By default
  // we do not have zero-flux boundary conditions in any direction.
  bool is_zero_flux[2*GKYL_MAX_DIM] = {false};
  for (int dir=0; dir<app->cdim; ++dir) {
    if (s->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX)
      is_zero_flux[dir] = true;
    if (s->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX)
      is_zero_flux[dir+pdim] = true;
  }

  // Allocate distribution function array for initialization and I/O.
  s->f = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);

  s->f_host = s->f;
  if (app->use_gpu) {
    s->f_host = mkarr(false, s->basis.num_basis, s->local_ext.volume);
  }

  if (app->cdim < 3) {
    // Reorganize g_ij and gij as done in calculation of Hamiltonian to
    // compute momentum and temperature.
    s->gij = mkarr(app->use_gpu, app->gk_geom->geo_int.gij_neut->ncomp, app->gk_geom->geo_int.gij_neut->size);
    s->g_ij = mkarr(app->use_gpu, app->gk_geom->geo_int.g_ij_neut->ncomp, app->gk_geom->geo_int.g_ij_neut->size);

    // Reorganize the metric tensor so ignorable coordinates are last.
    int metric_reorg_idxs_1x[] = {5, 2, 4, 0, 1, 3};
    int metric_reorg_idxs_2x[] = {0, 2, 1, 5, 4, 3};
    int *metric_reorg_idxs = app->cdim == 1? metric_reorg_idxs_1x : metric_reorg_idxs_2x;
    int num_basis_conf = app->basis.num_basis;
    struct gkyl_array *tmp_arr = mkarr(app->use_gpu, num_basis_conf, app->gk_geom->geo_int.gij_neut->size);
    for (int i=0; i<6; i++) {
      gkyl_array_move_comp(s->gij, i*num_basis_conf, app->gk_geom->geo_int.gij_neut, metric_reorg_idxs[i]*num_basis_conf, tmp_arr);
      gkyl_array_move_comp(s->g_ij, i*num_basis_conf, app->gk_geom->geo_int.g_ij_neut, metric_reorg_idxs[i]*num_basis_conf, tmp_arr);
    }
    gkyl_array_release(tmp_arr);
  }
  else {
    s->gij = gkyl_array_acquire(app->gk_geom->geo_int.gij);
    s->g_ij = gkyl_array_acquire(app->gk_geom->geo_int.g_ij);
  }

  // Allocate array for the Hamiltonian.
  s->hamil = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    
  // Call updater to evaluate hamiltonian
  struct gkyl_dg_calc_gk_neut_hamil* hamil_calc = gkyl_dg_calc_gk_neut_hamil_new(&s->grid, &s->basis, app->cdim, app->use_gpu);
  gkyl_dg_calc_gk_neut_hamil_calc(hamil_calc, &app->local, &s->local, s->gij, s->hamil);
  gkyl_dg_calc_gk_neut_hamil_release(hamil_calc);
    
  s->hamil_host = s->hamil;
  if (app->use_gpu) {
    s->hamil_host = mkarr(false, s->basis.num_basis, s->local_ext.volume);
    gkyl_array_copy(s->hamil_host, s->hamil);
  }

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
  gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, app->poly_order);

  int alpha_surf_sz = (cdim+vdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (cdim+vdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  s->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, s->local_ext.volume);
  s->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, s->local_ext.volume);
  s->const_sgn_alpha = mk_int_arr(app->use_gpu, cdim+vdim, s->local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&s->grid, 
    &app->basis, &s->basis, app->use_gpu);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &app->local, &s->local, &s->local_ext, s->hamil,
    s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);

  struct gkyl_dg_canonical_pb_auxfields aux_inp = {.hamil = s->hamil, .alpha_surf = s->alpha_surf, 
    .sgn_alpha_surf = s->sgn_alpha_surf, .const_sgn_alpha = s->const_sgn_alpha};

  // Set field type and model id for neutral species in GK system and create solver
  s->field_id = GKYL_FIELD_NULL;
  s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->basis, &s->basis, 
    &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);

  // Acquire equation object.
  s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);
  
  // Create skin/ghost ranges fir applying BCs. Only used for dynamic neutrals but included here to avoid
  // code duplication since the "ghost" array is needed.
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->lower_skin[dir], &s->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &s->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->upper_skin[dir], &s->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);
  }

  // Global skin and ghost ranges, only valid (i.e. volume>0) in ranges
  // abutting boundaries.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&s->global_lower_skin[dir], &s->global_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &s->global_ext, ghost); 
    gkyl_sub_range_intersect(&s->global_lower_skin[dir], &s->local_ext, &s->global_lower_skin[dir]);
    gkyl_sub_range_intersect(&s->global_lower_ghost[dir], &s->local_ext, &s->global_lower_ghost[dir]);
    gkyl_skin_ghost_ranges(&s->global_upper_skin[dir], &s->global_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);
    gkyl_sub_range_intersect(&s->global_upper_skin[dir], &s->local_ext, &s->global_upper_skin[dir]);
    gkyl_sub_range_intersect(&s->global_upper_ghost[dir], &s->local_ext, &s->global_upper_ghost[dir]);
  }

  if (s->info.init_from_file.type == 0) {
    // Initialize projection routine for initial conditions.
    gk_neut_species_projection_init(app, s, s->info.projection, &s->proj_init);
  }
  else {
    // Read initial condition from file.
    gk_neut_species_file_import_init(app, s, s->info.init_from_file);
  }
  
  // Allocate object for computing number .density 
  gk_neut_species_moment_init(app, s, &s->m0, GKYL_F_MOMENT_M0, false);

  // Allocate objects for computing diagnostic moments.
  int ndm = s->info.num_diag_moments;
  s->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_neut_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m], false);

  // Initialize boundary fluxes.
  s->bflux = (struct gk_boundary_fluxes) { };
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp) { };
  // Set the operation type for the bflux app.
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  if (s->info.boundary_flux_diagnostics.num_diag_moments > 0 ||
      s->info.boundary_flux_diagnostics.num_integrated_diag_moments > 0) {
    bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS;
  }
  else {
    // Set bflux_type to 
    //   - GK_SPECIES_BFLUX_CALC_FLUX to only put bfluxes in ghost cells of rhs.
    //   - GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS to calc bfluxes and step its moments.
    // The latter also requires that you place the moment you desire in add_bflux_moms_inp below.
  }
  // Introduce new moments into moms_inp if needed.
  gk_neut_species_bflux_init(app, s, &s->bflux, bflux_type, add_bflux_moms_inp);
  
  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  s->lte = (struct gk_lte) { };
  bool correct_all_moms = s->info.correct.correct_all_moms;
  int max_iter = s->info.correct.max_iter > 0 ? s->info.correct.max_iter : 50;
  double iter_eps = s->info.correct.iter_eps > 0 ? s->info.correct.iter_eps  : 1e-10;
  bool use_last_converged = s->info.correct.use_last_converged;
  struct correct_all_moms_inp corr_inp = { .correct_all_moms = correct_all_moms, 
    .max_iter = max_iter, .iter_eps = iter_eps, 
    .use_last_converged = use_last_converged };
  gk_neut_species_lte_init(app, s, &s->lte, corr_inp);

  s->enforce_positivity = false;
  
  // Initialize empty structs. New methods will fill them if specified.
  s->src = (struct gk_source) { };
  s->bgk = (struct gk_bgk_collisions) { };
  s->react_neut = (struct gk_react) { };
  if (!s->info.is_static) {
    gk_neut_species_new_dynamic(gk, app, s);
  }
  else {
    gk_neut_species_new_static(gk, app, s);
  }
}

void
gk_neut_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double t0)
{
  if (species->info.init_from_file.type == 0)
    gk_neut_species_projection_calc(app, species, &species->proj_init, species->f, t0);

  // we are pre-computing source for now as it is time-independent
  gk_neut_species_source_calc(app, species, &species->src, species->lte.f_lte, t0);
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
gk_neut_species_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms) 
{
  return species->rhs_func(app, species, fin, rhs, bflux_moms);
}

// Compute the implicit RHS for species update, returning maximum stable
// time-step.
double
gk_neut_species_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt)
{
  return species->rhs_implicit_func(app, species, fin, rhs, dt);
}

// Accummulate function for forward euler method.
void
gk_neut_species_step_f(struct gk_neut_species *species, struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  species->step_f_func(out, a, inp);
}

// Combine function for rk3 updates.
void
gk_neut_species_combine(struct gk_neut_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  species->combine_func(out, c1, arr1, c2, arr2, rng);
}

// Copy function for rk3 updates.
void
gk_neut_species_copy_range(struct gk_neut_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  species->copy_func(out, inp, range);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
gk_neut_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  species->bc_func(app, species, f);
}

// write functions
void
gk_neut_species_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  if (frame == 0) {
    gk_neut_species_write_dynamic(app, gkns, tm, frame);
  }
  else
    gkns->write_func(app, gkns, tm, frame);
}

void
gk_neut_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  if (frame == 0) {
    gk_neut_species_write_mom_dynamic(app, gkns, tm, frame);
  }
  else
    gkns->write_mom_func(app, gkns, tm, frame);
}

void
gk_neut_species_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm)
{
  gkns->calc_integrated_mom_func(app, gkns, tm);
}

void
gk_neut_species_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns)
{
  gkns->write_integrated_mom_func(app, gkns);
}

void
gk_neut_species_n_iter_corr(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_neut_species; ++i) {
    app->stat.neut_num_corr[i] = app->neut_species[i].lte.num_corr;
    app->stat.neut_n_iter_corr[i] = app->neut_species[i].lte.n_iter;
  }
}

void
gk_neut_species_release(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s)
{
  // Release resources for neutral species.

  gkyl_array_release(s->f);
  if (s->info.init_from_file.type == 0) {
    gk_neut_species_projection_release(app, &s->proj_init);
  }
  gkyl_comm_release(s->comm);

  if (app->use_gpu) {
    gkyl_array_release(s->f_host);
    gkyl_cu_free(s->basis_on_dev);
  }

  gkyl_velocity_map_release(s->vel_map);

  gkyl_array_release(s->alpha_surf);
  gkyl_array_release(s->sgn_alpha_surf);
  gkyl_array_release(s->const_sgn_alpha);

  // Release equation object and solver.
  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_release(s->slvr);

  // Release moment data.
  gk_neut_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    gk_neut_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);

  // Free boundary flux memory.
  gk_neut_species_bflux_release(app, s, &s->bflux);
  
  gk_neut_species_lte_release(app, &s->lte);

  gkyl_array_release(s->gij);
  gkyl_array_release(s->g_ij);

  gkyl_array_release(s->hamil);
  if (app->use_gpu)
    gkyl_array_release(s->hamil_host);

  s->release_func(app, s);
}
