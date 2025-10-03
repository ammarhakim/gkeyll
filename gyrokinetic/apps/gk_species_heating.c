#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static double
gk_heating_volume_integrate(gkyl_gyrokinetic_app *app, struct gk_heating *src, const struct gkyl_array *arrin)
{
  // Compute the volume integral of arrin.
  gkyl_array_integrate_advance(src->vol_integ_op, arrin, 1.0,
    0, &app->local, 0, src->volint_local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, src->volint_local, src->volint_global);
  double volint_global = 0.0;
  if (app->use_gpu)
    gkyl_cu_memcpy(&volint_global, src->volint_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&volint_global, src->volint_global, sizeof(double));
  return volint_global;
}

static void
gk_species_heating_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_heating *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_species_heating_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_heating *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute Maxwellian moments (n, u_par, T/m).
  gk_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, app->basis, 0, species->lte.moms.marr, 
    0, species->lte.moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local);  

  // Volume integrate Jrate times the thermal M2.
  gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_mom, 0, src->Jrate, 0, species->lte.moms.marr, &app->local);
  gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_mom, 0, src->Jrate_mom, 2, species->lte.moms.marr, &app->local);
  double Jrate_M2thermal_int = GKYL_MAX2(0.0, gk_heating_volume_integrate(app, src, src->Jrate_mom));

  // Volume integrate Jrate times the vtsq_shape time M0.
  gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_mom, 0, src->Jrate_vtsq_shape, 0, species->lte.moms.marr, &app->local);
  double Jrate_vtsq_shape_M0_int = GKYL_MAX2(0.0, gk_heating_volume_integrate(app, src, src->Jrate_mom));

  // Thermal speed squared of the Maxwellian.
  src->vtsq_amplitude = (src->norm_power + Jrate_M2thermal_int)/Jrate_vtsq_shape_M0_int;
  gkyl_array_set_offset_range(species->lte.moms.marr, src->vtsq_amplitude, src->vtsq_shape, 2*app->basis.num_basis, &app->local);

  // Compute the Maxwellian.
  gk_species_lte_from_moms(app, species, &species->lte, species->lte.moms.marr);

  // Multiply the Maxwellian by Jrate.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, src->Jrate_fmax, 
    src->Jrate, species->lte.f_lte, &app->local, &species->local);

  // Assemble the BGK-like term and add it to rhs.
  gkyl_bgk_collisions_advance(src->bgk_op, &app->local, &species->local, 
    src->rate, src->Jrate_fmax, fin, src->implicit_step, src->dt_implicit, rhs, species->cflrate);

  app->stat.species_heat_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_heating_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_heating *src, double tm, int frame)
{
}

static void
gk_species_heating_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_heating *src, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // Write the Maxwellian square thermal speed amplitude.
  gkyl_dynvec_append(src->vtsq_amp_diag, tm, &src->vtsq_amplitude);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    const char *fmt = "%s-%s_heating_vtsq_amplitude.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name);
    
    if (src->is_first_diag_dynvec_write_call) {
      gkyl_dynvec_write(src->vtsq_amp_diag, fileNm);
      src->is_first_diag_dynvec_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(src->vtsq_amp_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(src->vtsq_amp_diag);
  app->stat.n_diag_io += 1;
  
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_heating_write_conf_array(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_heating *src, int frame, double stime, char* file_suffix,
  struct gkyl_array *arrout, struct gkyl_array *arrout_host)
{
  // Write out a conf-space array.
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = stime,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }, GKYL_GK_META_NONE, 0
  );
  // Construct the file handles for collision frequency and primitive moments.
  const char *fmt = "%s-%s_%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, file_suffix, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, file_suffix, frame);

  struct gkyl_array *arr_ho;
  if (app->use_gpu) {  
    if (arrout_host)
      arr_ho = gkyl_array_acquire(arrout_host);
    else {
      arr_ho = mkarr(false, arrout->ncomp, arrout->size);
    }
    // Copy data from device to host before writing it out.
    gkyl_array_copy(arr_ho, arrout);
  }
  else {
    if (arrout_host)
      arr_ho = gkyl_array_acquire(arrout_host);
    else
      arr_ho = gkyl_array_acquire(arrout);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, arr_ho, fileNm);
  gk_array_meta_release(mt); 
  gkyl_array_release(arr_ho);
}

void 
gk_species_heating_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_heating *src)
{
  src->heating_id = gks->info.heating.heating_id;
  src->write_diagnostics = gks->info.heating.write_diagnostics;

  src->write_diags_func = gk_species_heating_write_diags_disabled;
  src->rhs_func = gk_species_heating_rhs_disabled;

  if (src->heating_id) {
    src->norm_power = 2.0*gks->info.heating.power/((app->vdim == 1? 1.0 : 3.0)*gks->info.mass);

    // Heating rate.
    src->rate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *rate_host = app->use_gpu? mkarr(false, src->rate->ncomp, src->rate->size)
                                               : gkyl_array_acquire(src->rate);
    gkyl_proj_on_basis *proj_rate = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      app->poly_order+1, 1, gks->info.heating.rate_profile, gks->info.heating.rate_profile_ctx);
    gkyl_proj_on_basis_advance(proj_rate, 0.0, &app->local, rate_host);
    gkyl_array_copy(src->rate, rate_host);
    gkyl_proj_on_basis_release(proj_rate);
    gkyl_array_release(rate_host);
    // Multiply the rate by the conf-space Jacobian.
    src->Jrate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_dg_mul_op_range(app->basis, 0, src->Jrate, 0, app->gk_geom->geo_int.jacobgeo, 0, src->rate, &app->local);

    // Heating rate.
    src->vtsq_shape = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *vtsq_shape_host = app->use_gpu? mkarr(false, src->vtsq_shape->ncomp, src->vtsq_shape->size)
                                                     : gkyl_array_acquire(src->vtsq_shape);
    gkyl_proj_on_basis *proj_vtsq_shape = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      app->poly_order+1, 1, gks->info.heating.temp_shape, gks->info.heating.temp_shape_ctx);
    gkyl_proj_on_basis_advance(proj_vtsq_shape, 0.0, &app->local, vtsq_shape_host);
    gkyl_array_copy(src->vtsq_shape, vtsq_shape_host);
    gkyl_proj_on_basis_release(proj_vtsq_shape);
    gkyl_array_release(vtsq_shape_host);

    // Multiply Jrate by the shape of v_t^2.
    src->Jrate_vtsq_shape = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_dg_mul_op_range(app->basis, 0, src->Jrate_vtsq_shape, 0, src->Jrate, 0, src->vtsq_shape, &app->local);

    // Rate times the Maxwellian.
    src->Jrate_fmax = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
    // Rate times a velocity moment.
    src->Jrate_mom = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    // Volume integrator.
    src->vol_integ_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    if (app->use_gpu) {
      src->volint_local = gkyl_cu_malloc(sizeof(double));
      src->volint_global = gkyl_cu_malloc(sizeof(double));
    } 
    else {
      src->volint_local = gkyl_malloc(sizeof(double));
      src->volint_global = gkyl_malloc(sizeof(double));
    }

    // BGK operator.
    src->bgk_op = gkyl_bgk_collisions_new(&app->basis, &gks->basis, app->use_gpu);
    src->implicit_step = false;
    src->dt_implicit = 1e9;

    if (src->write_diagnostics) {
      src->vtsq_amp_diag = gkyl_dynvec_new(GKYL_DOUBLE, 1);
      // Write out the heating rate and vtsq shape.
      gk_heating_write_conf_array(app, gks, src, 0, 0.0, "heating_rate", src->rate, 0);
      gk_heating_write_conf_array(app, gks, src, 0, 0.0, "heating_temp_shape", src->vtsq_shape, 0);
    }

    // Methods chosen at runtime.
    src->rhs_func = gk_species_heating_rhs_enabled;
    if (src->write_diagnostics) {
      src->write_diags_func = gk_species_heating_write_diags_enabled;
    }
  }
}

void
gk_species_heating_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_heating *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  src->rhs_func(app, species, src, fin, rhs);
}

void
gk_species_heating_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_heating *src, double tm, int frame)
{
  src->write_diags_func(app, gks, src, tm, frame);
}

void
gk_species_heating_release(const struct gkyl_gyrokinetic_app *app, const struct gk_heating *src)
{
  if (src->heating_id) {
    gkyl_array_release(src->rate);
    gkyl_array_release(src->Jrate);
    gkyl_array_release(src->vtsq_shape);
    gkyl_array_release(src->Jrate_vtsq_shape);
    gkyl_array_release(src->Jrate_fmax);
    gkyl_array_release(src->Jrate_mom);

    gkyl_array_integrate_release(src->vol_integ_op);
    if (app->use_gpu) {
      gkyl_cu_free(src->volint_local);
      gkyl_cu_free(src->volint_global);
    }
    else {
      gkyl_free(src->volint_local);
      gkyl_free(src->volint_global);
    }

    gkyl_bgk_collisions_release(src->bgk_op);

    if (src->write_diagnostics) {
      gkyl_dynvec_release(src->vtsq_amp_diag);
    }
  }
}
