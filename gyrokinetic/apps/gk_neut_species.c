#include <gkyl_gk_neut_species_priv.h>

void
gk_neut_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns)
{
  if (ns->info.cells[0] == 0) {
    // Fluid neutrals.
    gk_neut_species_fluid_init(gk, app, ns);
  }
  else {
    // Kinetic neutrals.
    gk_neut_species_kinetic_init(gk, app, ns);
  }
}

void
gk_neut_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double t0)
{
  if (species->info.init_from_file.type == 0)
    gk_neut_species_projection_calc(app, species, &species->proj_init, species->f, t0);

  // We are pre-computing source for now as it is time-independent.
  gk_neut_species_source_calc(app, species, &species->src, species->lte.f_lte, t0);
}

void
gk_neut_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns_self, double t0)
{
  // Store initial density in recycle_react_scale.
  gk_neut_species_recycle_react_scale_apply_ic_cross(app, gkns_self, &gkns_self->rrs);
}

double
gk_neut_species_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms) 
{
  // Compute the RHS for species update, returning maximum stable time-step.
  return species->rhs_func(app, species, fin, rhs, bflux_moms);
}

double
gk_neut_species_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt)
{
  // Compute the implicit RHS for species update, returning maximum stable time-step.
  return species->rhs_implicit_func(app, species, fin, rhs, bflux_moms, dt);
}

void
gk_neut_species_step_f(struct gk_neut_species *species, struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  // Accummulate function for forward euler method.
  species->step_f_func(out, a, inp);
}

void
gk_neut_species_combine(struct gk_neut_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // Combine function for rk3 updates.
  species->combine_func(out, c1, arr1, c2, arr2, rng);
}

void
gk_neut_species_copy_range(struct gk_neut_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // Copy function for rk3 updates.
  species->copy_func(out, inp, range);
}

void
gk_neut_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  species->bc_func(app, species, f);
}

void
gk_neut_species_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  gkns->write_func(app, gkns, tm, frame);
}

void
gk_neut_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
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
gk_neut_species_n_iter_corr(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int sidx)
{
  gkns->report_n_iter_corr_func(app, gkns, sidx);
}

void
gk_neut_species_release(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *gkns)
{
  gkns->release_func(app, gkns);
}

/**
 *
 * Methods common to fluid and kinetic neutrals.
 *
 */

double
gk_neut_species_rhs_static(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  double omega_cfl = 1/DBL_MAX;
  return app->cfl/omega_cfl;
}

double
gk_neut_species_rhs_implicit_static(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt)
{
  double omega_cfl = 1/DBL_MAX;
  return app->cfl/omega_cfl;
}

void
gk_neut_species_apply_bc_static(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  // empty function
}

void
gk_neut_species_step_f_dynamic(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

void
gk_neut_species_step_f_static(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // Do nothing.
}

void
gk_neut_species_combine_dynamic(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

void
gk_neut_species_combine_static(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // Do nothing.
}

void
gk_neut_species_copy_range_dynamic(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  gkyl_array_copy_range(out, inp, range);
}

void
gk_neut_species_copy_range_static(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // Do nothing.
}

void
gk_neut_species_write_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  int io_meta_len[] = {app->io_meta_basic_len, gkns->io_meta_len, app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, gkns->io_meta, app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, frame);
  
  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(gkns->f_host, gkns->f);
  }
  
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(gkns->comm, &gkns->grid, &gkns->local, mt, gkns->f_host, fileNm);
  app->stat.neut_species_io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.n_neut_io += 1;
  
  gk_array_meta_release(mt);  
}

void
gk_neut_species_write_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{

  // Package metadata.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  int io_meta_len[] = {app->io_meta_basic_len, app->io_meta_len, app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, app->io_meta, app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

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
  gk_array_meta_release(mt); 
}

void
gk_neut_species_calc_integrated_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  int num_mom = gkns->integ_moms.num_mom;
  double avals_global[num_mom];
  
  gk_neut_species_moment_calc(&gkns->integ_moms, gkns->local, app->local, gkns->f); 
  app->stat.n_neut_mom += 1;

  // Reduce to compute sum over whole domain, append to diagnostics.
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

  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_neut_diag += 1;
}

void
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

  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_neut_species_write_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  // Do nothing.
}

void
gk_neut_species_write_init_only(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  gk_neut_species_write_dynamic(app, gkns, tm, frame);
  gkns->write_func = gk_neut_species_write_static;
}

void
gk_neut_species_write_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  // Do nothing.
}

void
gk_neut_species_write_mom_init_only(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame)
{
  gk_neut_species_write_mom_dynamic(app, gkns, tm, frame);
  gkns->write_mom_func = gk_neut_species_write_mom_static;
}

void
gk_neut_species_calc_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm)
{
  // Do nothing.
}

void
gk_neut_species_write_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns)
{
  // Do nothing.
}

void
gk_neut_species_n_iter_corr_enabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int sidx)
{
  app->stat.neut_num_corr[sidx] = gkns->lte.num_corr;
  app->stat.neut_n_iter_corr[sidx] = gkns->lte.n_iter;
}

void
gk_neut_species_n_iter_corr_disabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int sidx)
{
}

void
gk_neut_species_release_static(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s)
{ 
  // Do nothing.
}
