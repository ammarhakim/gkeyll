#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic_pol_density.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_translate_dim_gyrokinetic.h>
#include <gkyl_proj_on_basis.h>

#include <assert.h>
#include <time.h>

// Begin static function definitions.
static double
gk_species_omegaH_dt(gkyl_gyrokinetic_app *app, struct gk_species *gks, const struct gkyl_array *fin)
{
  // Compute the time step under which the omega_H mode is stable, dt_omegaH =
  // omegaH_CFL/omega_H.
  // Each species computes its own omega_H as:
  //   omega_H = q_e*sqrt(jacobgeo*n_{s0}/m_s) * omegaH_gf
  // where
  //   - n_{s0} is either a reference, average or max density.
  //   - omegaH_gf = (cmag/(jacobgeo*B^_\parallel))*kpar_max / 
  //                 min(sqrt(k_x^2*eps_xx+k_x*k_y*eps_xy+k_y^2*eps_yy+)).
  // and k_x,k_y,k_par are wavenumbers in computational space, and eps_ij is
  // the polarization weight in our field equation.
 
  if (!(app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN || app->field->gkfield_id == GKYL_GK_FIELD_ADIABATIC)) {
    // Obtain the maximum density (using cell centers).
    gk_species_moment_calc(&gks->m0, gks->local, app->local, fin);
    gkyl_array_reduce_range(gks->m0_max, gks->m0.marr, GKYL_MAX, &app->local);
  
    double m0_max[1];
    if (app->use_gpu)
      gkyl_cu_memcpy(m0_max, gks->m0_max, sizeof(double), GKYL_CU_MEMCPY_D2H);
    else
      m0_max[0] = gks->m0_max[0];
    m0_max[0] *= 1.0/pow(sqrt(2.0),app->cdim);
  
    double omegaH = fabs(gks->info.charge)*sqrt(GKYL_MAX2(0.0,m0_max[0])/gks->info.mass)*app->omegaH_gf;
  
    return omegaH > 1e-20? app->cfl_omegaH/omegaH : DBL_MAX;
  }
  else
    return DBL_MAX;
}

static double
gk_species_rhs_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_set(species->phi, 1.0, app->field->phi_smooth);

  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute the surface expansion of the phase space flux
  // Note: Each cell stores the *lower* surface expansions of the 
  // phase space flux, so local_ext range needed to index the output
  // values of alpha_surf even though we only loop over local ranges
  // to avoid evaluating quantities such as geometry in ghost cells
  // where they are not defined.
  
  gkyl_dg_calc_gyrokinetic_vars_alpha_surf(species->calc_gk_vars, 
    &app->local, &species->local, &species->local_ext, 
    species->phi, species->alpha_surf, species->sgn_alpha_surf, species->const_sgn_alpha);

  gkyl_dg_updater_gyrokinetic_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);

  if (species->lbo.collision_id == GKYL_LBO_COLLISIONS)
    gk_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  if (species->bgk.collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  
  if (species->has_diffusion)
    gkyl_dg_updater_diffusion_gyrokinetic_advance(species->diff_slvr, &species->local, 
      species->diffD, app->gk_geom->jacobgeo_inv, fin, species->cflrate, rhs);

  if (species->rad.radiation_id == GKYL_GK_RADIATION)
    gk_species_radiation_rhs(app, species, &species->rad, fin, rhs);

  if (species->react.num_react)
    gk_species_react_rhs(app, species, &species->react, fin, rhs);
  if (species->react_neut.num_react)
    gk_species_react_rhs(app, species, &species->react_neut, fin, rhs);
  
  app->stat.nspecies_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omega_cfl_ho[0] = species->omega_cfl[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);

  double dt_out = app->cfl/omega_cfl_ho[0];
  
  // Enforce the omega_H constraint on dt.
  double dt_omegaH = gk_species_omegaH_dt(app, species, fin);
  dt_out = fmin(dt_out, dt_omegaH);

  return dt_out;
}

static double
gk_species_rhs_static(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  double omega_cfl = 1/DBL_MAX;
  return app->cfl/omega_cfl;
}

static void
gk_species_apply_bc_dynamic(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = species->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, species->periodic_dirs, f); 
  
  for (int d=0; d<cdim; ++d) {
    if (species->bc_is_np[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_GK_SHEATH:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_lo, app->field->phi_smooth, 
            app->field->phi_wall_lo, f, &app->local);
          break;
        case GKYL_SPECIES_GK_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_lo, app->field->phi_smooth, 
            app->field->phi_wall_lo, f, &app->local);
	  if (cdim == 3)
            gkyl_bc_twistshift_advance(species->bc_ts_lo, f, f);
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
        case GKYL_SPECIES_GK_SHEATH:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_up, app->field->phi_smooth, 
            app->field->phi_wall_up, f, &app->local);
          break;
        case GKYL_SPECIES_GK_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_up, app->field->phi_smooth, 
            app->field->phi_wall_up, f, &app->local);
	  if (cdim == 3)
            gkyl_bc_twistshift_advance(species->bc_ts_up, f, f);
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

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_apply_bc_static(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  // do nothing
}

static void
gk_species_step_f_dynamic(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

static void
gk_species_step_f_static(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp)
{
  // do nothing
}

static void
gk_species_combine_dynamic(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

static void
gk_species_combine_static(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  // do nothing
}

static void
gk_species_copy_range_dynamic(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  gkyl_array_copy_range(out, inp, range);
}

static void
gk_species_copy_range_static(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  // do nothing
}

static void
gk_species_write_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
  
  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(gks->f_host, gks->f);
  }
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->f_host, fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 1;
    
  gk_array_meta_release(mt);  
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

static void
gk_species_write_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // do nothing
}

static void
gk_species_write_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );
  for (int m=0; m<gks->info.num_diag_moments; ++m) {
    gk_species_moment_calc(&gks->moms[m], gks->local, app->local, gks->f);
    app->stat.nmom += 1;
    
    const char *fmt = "%s-%s_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
      gks->info.diag_moments[m], frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
      gks->info.diag_moments[m], frame);
    
    // Rescale moment by inverse of Jacobian if not already re-scaled 
    if (!gks->moms[m].is_bimaxwellian_moms && !gks->moms[m].is_maxwellian_moms) {
      gkyl_dg_div_op_range(gks->moms[m].mem_geo, app->confBasis, 
        0, gks->moms[m].marr, 0, gks->moms[m].marr, 0, 
        app->gk_geom->jacobgeo, &app->local);  
    }    
      
    if (app->use_gpu) {
      gkyl_array_copy(gks->moms[m].marr_host, gks->moms[m].marr);
    }
    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gks->moms[m].marr_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;
  }
  if (app->enforce_positivity) {
    // We placed the change in f from the positivity shift in fnew.
    gk_species_moment_calc(&gks->ps_moms, gks->local, app->local, gks->fnew);
    app->stat.nmom += 1;
    const char *fmt = "%s-%s_positivity_shift_FourMoments_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
    
    // Rescale moment by inverse of Jacobian.
    gkyl_dg_div_op_range(gks->ps_moms.mem_geo, app->confBasis, 
      0, gks->ps_moms.marr, 0, gks->ps_moms.marr, 0, 
      app->gk_geom->jacobgeo, &app->local);  
    if (app->use_gpu) {
      gkyl_array_copy(gks->ps_moms.marr_host, gks->ps_moms.marr);
    }
    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gks->ps_moms.marr_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;
  }
  gk_array_meta_release(mt);   
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

static void
gk_species_write_mom_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // do nothing
}

static void
gk_species_calc_integrated_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  int vdim = app->vdim;
  double avals_global[2+vdim];
  
  gk_species_moment_calc(&gks->integ_moms, gks->local, app->local, gks->f); 
  // Reduce (sum) over whole domain, append to diagnostics.
  gkyl_array_reduce_range(gks->red_integ_diag, gks->integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
    gks->red_integ_diag, gks->red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[2+vdim]));
  }
  gkyl_dynvec_append(gks->integ_diag, tm, avals_global);
  
  if (app->enforce_positivity) {
    // The change in f from the positivity shift is in fnew.
    gk_species_moment_calc(&gks->integ_moms, gks->local, app->local, gks->fnew); 
    // Reduce (sum) over whole domain, append to diagnostics.
    gkyl_array_reduce_range(gks->red_integ_diag, gks->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
      gks->red_integ_diag, gks->red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[2+vdim]));
      }
    gkyl_dynvec_append(gks->ps_integ_diag, tm, avals_global);
  }
  
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

static void
gk_species_calc_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  // do nothing
}

static void
gk_species_write_integrated_mom_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  
  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");
    
    struct timespec wtm = gkyl_wall_clock();
    if (gks->is_first_integ_write_call) {
      gkyl_dynvec_write(gks->integ_diag, fileNm);
      gks->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(gks->integ_diag, fileNm);
    }
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;
  }
  gkyl_dynvec_clear(gks->integ_diag);
  
  if (app->enforce_positivity) {
    if (rank == 0) {
      // Write integrated diagnostic moments.
      const char *fmt = "%s-%s_positivity_shift_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");
      
      struct timespec wtm = gkyl_wall_clock();
      if (gks->is_first_ps_integ_write_call) {
	gkyl_dynvec_write(gks->ps_integ_diag, fileNm);
	gks->is_first_ps_integ_write_call = false;
      }
      else {
	gkyl_dynvec_awrite(gks->ps_integ_diag, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gks->ps_integ_diag);
  }
  
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

static void
gk_species_write_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  // do nothing
}

static void
gk_species_calc_L2norm_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_integrate_advance(gks->integ_wfsq_op, gks->f, (2.0*M_PI)/gks->info.mass,
    app->jacobtot_inv_weak, &gks->local, &app->local, gks->L2norm_local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, gks->L2norm_local, gks->L2norm_global);

  double L2norm_global[] = {0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(L2norm_global, gks->L2norm_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    memcpy(L2norm_global, gks->L2norm_global, sizeof(double));
  }
  
  gkyl_dynvec_append(gks->L2norm, tm, L2norm_global);  

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

static void
gk_species_calc_L2norm_static(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  // do nothing
}

static void
gk_species_write_L2norm_dynamic(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write the L2 norm.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "L2norm");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "L2norm");

    struct timespec wtm = gkyl_wall_clock();
    if (gks->is_first_L2norm_write_call) {
      gkyl_dynvec_write(gks->L2norm, fileNm);
      gks->is_first_L2norm_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(gks->L2norm, fileNm);
    }
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;
  }
  gkyl_dynvec_clear(gks->L2norm);
}

static void
gk_species_write_L2norm_static(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  // do nothing
}

static void
gk_species_release_dynamic(const gkyl_gyrokinetic_app* app, const struct gk_species *s)
{
    // release various arrays and species objects
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  if (s->lbo.collision_id == GKYL_LBO_COLLISIONS)
    gk_species_lbo_release(app, &s->lbo);
  if (s->bgk.collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_release(app, &s->bgk);

  if (s->react.num_react)
    gk_species_react_release(app, &s->react);
  if (s->react_neut.num_react)
    gk_species_react_release(app, &s->react_neut);  

  if (s->rad.radiation_id == GKYL_GK_RADIATION) 
    gk_species_radiation_release(app, &s->rad);

  if (s->has_diffusion) {
    gkyl_array_release(s->diffD);
    gkyl_dg_updater_diffusion_gyrokinetic_release(s->diff_slvr);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d].type == GKYL_SPECIES_GK_SHEATH)
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_lo);
    else if (s->lower_bc[d].type == GKYL_SPECIES_GK_IWL) { 
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_lo);
      if (app->cdim == 3)
        gkyl_bc_twistshift_release(s->bc_ts_lo);
    }
    else 
      gkyl_bc_basic_release(s->bc_lo[d]);
    
    if (s->upper_bc[d].type == GKYL_SPECIES_GK_SHEATH)
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_up);
    else if (s->upper_bc[d].type == GKYL_SPECIES_GK_IWL) {
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_up);
      if (app->cdim == 3)
        gkyl_bc_twistshift_release(s->bc_ts_up);
    }
    else 
      gkyl_bc_basic_release(s->bc_up[d]);
  }
  
  if (app->enforce_positivity) {
    gkyl_array_release(s->ps_delta_m0);
    gkyl_array_release(s->ps_delta_m0s_tot);
    gkyl_array_release(s->ps_delta_m0r_tot);
    gkyl_positivity_shift_gyrokinetic_release(s->pos_shift_op);
    gk_species_moment_release(app, &s->ps_moms);
    gkyl_dynvec_release(s->ps_integ_diag);
  }

  if (app->use_gpu) {
    gkyl_cu_free(s->omega_cfl);
    gkyl_cu_free(s->m0_max);
  }
  else {
    gkyl_free(s->omega_cfl);
    gkyl_free(s->m0_max);
  }

  // Release L2 norm memory.
  gkyl_array_integrate_release(s->integ_wfsq_op);
  gkyl_dynvec_release(s->L2norm);
  if (app->use_gpu) {
    gkyl_cu_free(s->L2norm_local);
    gkyl_cu_free(s->L2norm_global);
  }
  else {
    gkyl_free(s->L2norm_local);
    gkyl_free(s->L2norm_global);
  }
}

static void
gk_species_release_static(const gkyl_gyrokinetic_app* app, const struct gk_species *s)
{
}

// initialize species object
static void
gk_species_new_dynamic(struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int ghost[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    ghost[cdim+d] = 0; // No ghost-cells in velocity space.
  }

  // Allocate distribution function arrays.
  gks->f1 = mkarr(app->use_gpu, app->basis.num_basis, gks->local_ext.volume);
  gks->fnew = mkarr(app->use_gpu, app->basis.num_basis, gks->local_ext.volume);

  // Allocate cflrate (scalar array).
  gks->cflrate = mkarr(app->use_gpu, 1, gks->local_ext.volume);

  if (app->use_gpu) {
    gks->omega_cfl = gkyl_cu_malloc(sizeof(double));
    gks->m0_max = gkyl_cu_malloc(app->confBasis.num_basis*sizeof(double));
  }
  else {
    gks->omega_cfl = gkyl_malloc(sizeof(double));
    gks->m0_max = gkyl_malloc(app->confBasis.num_basis*sizeof(double));
  }

  for (int dir=0; dir<app->cdim; ++dir) {
    if (gks->lower_bc[dir].type == GKYL_SPECIES_GK_IWL || gks->upper_bc[dir].type == GKYL_SPECIES_GK_IWL) {
      // Make the parallel direction periodic so that we sync the core before
      // applying sheath BCs in the SOL.
      gks->periodic_dirs[gks->num_periodic_dir] = app->cdim-1; // The last direction is the parallel one.
      gks->num_periodic_dir += 1;
      // Check that the LCFS is the same on both BCs and that it's on a cell boundary within our grid.
      double xLCFS = gks->lower_bc[dir].aux_parameter;
      assert(fabs(xLCFS-gks->upper_bc[dir].aux_parameter) < 1e-14);
      // Check the split happens within the domain and at a cell boundary.
      assert((app->grid.lower[0]<xLCFS) && (xLCFS<app->grid.upper[0]));
      double needint = (xLCFS-app->grid.lower[0])/app->grid.dx[0];
      assert(floor(fabs(needint-floor(needint))) < 1.);
    }
  }
  
  // Determine collision type and initialize it.
  gks->lbo = (struct gk_lbo_collisions) { };
  gks->bgk = (struct gk_bgk_collisions) { };
  if (gks->info.collisions.collision_id == GKYL_LBO_COLLISIONS) {
    gk_species_lbo_init(app, gks, &gks->lbo);
  }
  if (gks->info.collisions.collision_id == GKYL_BGK_COLLISIONS) {
    gk_species_bgk_init(app, gks, &gks->bgk);
  }

  gks->react = (struct gk_react) { };
  gks->react_neut = (struct gk_react) { };
  if (gks->info.react.num_react) {
    gk_species_react_init(app, gks, gks->info.react, &gks->react, true);
  }
  if (gks->info.react_neut.num_react) {
    gk_species_react_init(app, gks, gks->info.react_neut, &gks->react_neut, false);
  }

  // Determine radiation type to use in gyrokinetic update.
  gks->rad = (struct gk_rad_drag) { };

  // Initialize diffusion if present.
  gks->has_diffusion = false;  
  if (gks->info.diffusion.num_diff_dir) {
    gks->has_diffusion = true;
    int diffusion_order = gks->info.diffusion.order ? gks->info.diffusion.order : 2;

    int szD = cdim*app->confBasis.num_basis;
    gks->diffD = mkarr(app->use_gpu, szD, app->local_ext.volume);
    bool diff_dir[GKYL_MAX_CDIM] = {false};

    int num_diff_dir = gks->info.diffusion.num_diff_dir ? gks->info.diffusion.num_diff_dir : app->cdim;
    // Assuming diffusion along x only for now.
    assert(num_diff_dir == 1);
    assert(gks->info.diffusion.diff_dirs[0] == 0);
    for (int d=0; d<num_diff_dir; ++d) {
      int dir = gks->info.diffusion.diff_dirs[d]; 
      diff_dir[dir] = 1; 
      gkyl_array_shiftc(gks->diffD, gks->info.diffusion.D[d]*pow(sqrt(2),app->cdim), dir);
    }
    // Multiply diffD by g^xx*jacobgeo.
    gkyl_dg_mul_op(app->confBasis, 0, gks->diffD, 0, app->gk_geom->gxxj, 0, gks->diffD);

    // By default, we do not have zero-flux boundary conditions in any direction
    // Determine which directions are zero-flux.
    bool is_zero_flux[2*GKYL_MAX_DIM] = {false};
    for (int dir=0; dir<app->cdim; ++dir) {
      if (gks->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir] = true;
      }
      if (gks->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir+pdim] = true;
      }
    }

    gks->diff_slvr = gkyl_dg_updater_diffusion_gyrokinetic_new(&gks->grid, &app->basis, &app->confBasis, 
      false, diff_dir, diffusion_order, &app->local, is_zero_flux, app->use_gpu);
  }
  
  // Allocate buffer needed for BCs.
  long buff_sz = 0;
  for (int dir=0; dir<cdim; ++dir) {
    long vol = GKYL_MAX2(gks->lower_skin[dir].volume, gks->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  gks->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  // buffer arrays for fixed function boundary conditions on distribution function
  gks->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  gks->bc_buffer_up_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;

    // Lower BC.
    if (gks->lower_bc[d].type == GKYL_SPECIES_GK_SHEATH) {
      gks->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, app->basis_on_dev.basis, 
        &gks->lower_skin[d], &gks->lower_ghost[d], gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else if (gks->lower_bc[d].type == GKYL_SPECIES_GK_IWL) {
      double xLCFS = gks->lower_bc[d].aux_parameter;
      // Index of the cell that abuts the xLCFS from below.
      int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
      gkyl_range_shorten_from_below(&gks->lower_skin_par_sol, &gks->lower_skin[d], 0, app->grid.cells[0]-idxLCFS_m+1);
      gkyl_range_shorten_from_below(&gks->lower_ghost_par_sol, &gks->lower_ghost[d], 0, app->grid.cells[0]-idxLCFS_m+1);

      gks->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, app->basis_on_dev.basis, 
        &gks->lower_skin_par_sol, &gks->lower_ghost_par_sol, gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);

      if (cdim == 3) {
        // For 3x2v we need a twistshift BC in the core.
        // Create a core local range, extended in the BC dir.
        int ndim = cdim+vdim;
        int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
        for (int i=0; i<ndim; i++) {
          lower_bcdir_ext[i] = gks->local.lower[i];
          upper_bcdir_ext[i] = gks->local.upper[i];
        }
        upper_bcdir_ext[0] = idxLCFS_m;
        lower_bcdir_ext[d] = gks->local_ext.lower[d];
        upper_bcdir_ext[d] = gks->local_ext.upper[d];
        gkyl_sub_range_init(&gks->local_par_ext_core, &gks->local_ext, lower_bcdir_ext, upper_bcdir_ext);

        struct gkyl_bc_twistshift_inp tsinp = {
          .bc_dir = d,
          .shift_dir = 1, // y shift.
          .shear_dir = 0, // shift varies with x.
          .edge = GKYL_LOWER_EDGE,
          .cdim = cdim,
          .bcdir_ext_update_r = gks->local_par_ext_core,
          .num_ghost = ghost,
          .basis = app->basis,
          .grid = gks->grid,
          .shift_func = gks->lower_bc[d].aux_profile,
          .shift_func_ctx = gks->lower_bc[d].aux_ctx,
          .use_gpu = app->use_gpu,
        };

        gks->bc_ts_lo = gkyl_bc_twistshift_new(&tsinp);
      }
    }
    else { 
      if (gks->lower_bc[d].type == GKYL_SPECIES_COPY) 
        bctype = GKYL_BC_COPY;
      else if (gks->lower_bc[d].type == GKYL_SPECIES_ABSORB) 
        bctype = GKYL_BC_ABSORB;
      else if (gks->lower_bc[d].type == GKYL_SPECIES_REFLECT) 
        bctype = GKYL_BC_REFLECT;
      else if (gks->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      gks->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.basis,
        &gks->lower_skin[d], &gks->lower_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu);

      if (gks->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_lo;
        gk_species_projection_init(app, gks, gks->lower_bc[d].projection, &gk_proj_bc_lo);
        gk_species_projection_calc(app, gks, &gk_proj_bc_lo, gks->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(gks->bc_lo[d], gks->bc_buffer_lo_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_lo);
      }
    }

    // Upper BC.
    if (gks->upper_bc[d].type == GKYL_SPECIES_GK_SHEATH) {
      gks->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, app->basis_on_dev.basis, 
        &gks->upper_skin[d], &gks->upper_ghost[d], gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else if (gks->upper_bc[d].type == GKYL_SPECIES_GK_IWL) {
      double xLCFS = gks->upper_bc[d].aux_parameter;
      // Index of the cell that abuts the xLCFS from below.
      int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
      gkyl_range_shorten_from_below(&gks->upper_skin_par_sol, &gks->upper_skin[d], 0, app->grid.cells[0]-idxLCFS_m+1);
      gkyl_range_shorten_from_below(&gks->upper_ghost_par_sol, &gks->upper_ghost[d], 0, app->grid.cells[0]-idxLCFS_m+1);

      gks->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, app->basis_on_dev.basis, 
        &gks->upper_skin_par_sol, &gks->upper_ghost_par_sol, gks->vel_map,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);

      if (cdim == 3) {
        // For 3x2v we need a twistshift BC in the core.
        // Create a core local range, extended in the BC dir.
        struct gkyl_bc_twistshift_inp tsinp = {
          .bc_dir = d,
          .shift_dir = 1, // y shift.
          .shear_dir = 0, // shift varies with x.
          .edge = GKYL_UPPER_EDGE,
          .cdim = cdim,
          .bcdir_ext_update_r = gks->local_par_ext_core,
          .num_ghost = ghost,
          .basis = app->basis,
          .grid = gks->grid,
          .shift_func = gks->upper_bc[d].aux_profile,
          .shift_func_ctx = gks->upper_bc[d].aux_ctx,
          .use_gpu = app->use_gpu,
        };

        gks->bc_ts_up = gkyl_bc_twistshift_new(&tsinp);
      }
    }
    else {
      if (gks->upper_bc[d].type == GKYL_SPECIES_COPY) 
        bctype = GKYL_BC_COPY;
      else if (gks->upper_bc[d].type == GKYL_SPECIES_ABSORB) 
        bctype = GKYL_BC_ABSORB;
      else if (gks->upper_bc[d].type == GKYL_SPECIES_REFLECT) 
        bctype = GKYL_BC_REFLECT;
      else if (gks->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      gks->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.basis,
        &gks->upper_skin[d], &gks->upper_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu);

      if (gks->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_up;
        gk_species_projection_init(app, gks, gks->upper_bc[d].projection, &gk_proj_bc_up);
        gk_species_projection_calc(app, gks, &gk_proj_bc_up, gks->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(gks->bc_up[d], gks->bc_buffer_up_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_up);
      }
    }
  }

  if (app->enforce_positivity) {
    // Positivity enforcing by shifting f (ps=positivity shift).
    gks->ps_delta_m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    // Set pointers to total ion/electron Delta m0.
    if (gks->info.charge > 0.0) {
      gks->ps_delta_m0s_tot = gkyl_array_acquire(app->ps_delta_m0_ions);
      gks->ps_delta_m0r_tot = gkyl_array_acquire(app->ps_delta_m0_elcs);
    }
    else {
      gks->ps_delta_m0s_tot = gkyl_array_acquire(app->ps_delta_m0_elcs);
      gks->ps_delta_m0r_tot = gkyl_array_acquire(app->ps_delta_m0_ions);
    }

    gks->pos_shift_op = gkyl_positivity_shift_gyrokinetic_new(app->confBasis, app->basis,
      gks->grid, gks->info.mass, app->gk_geom, gks->vel_map, &app->local_ext, app->use_gpu);

    // Allocate data for diagnostic moments
    gk_species_moment_init(app, gks, &gks->ps_moms, "FourMoments");

    gks->ps_integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, gks->ps_moms.num_mom);
    gks->is_first_ps_integ_write_call = true;
  }

  // set function pointers
  gks->rhs_func = gk_species_rhs_dynamic;
  gks->bc_func = gk_species_apply_bc_dynamic;
  gks->release_func = gk_species_release_dynamic;
  gks->step_f_func = gk_species_step_f_dynamic;
  gks->combine_func = gk_species_combine_dynamic;
  gks->copy_func = gk_species_copy_range_dynamic;
  gks->write_func = gk_species_write_dynamic;
  gks->write_mom_func = gk_species_write_mom_dynamic;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_dynamic;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_dynamic;
  gks->calc_L2norm_func = gk_species_calc_L2norm_dynamic;
  gks->write_L2norm_func = gk_species_write_L2norm_dynamic;
}

// initialize static species object
static void
gk_species_new_static(struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Allocate distribution function arrays.
  gks->f1 = gks->f;
  gks->fnew = gks->f;

  // set function pointers
  gks->rhs_func = gk_species_rhs_static;
  gks->bc_func = gk_species_apply_bc_static;
  gks->release_func = gk_species_release_static;
  gks->step_f_func = gk_species_step_f_static;
  gks->combine_func = gk_species_combine_static;
  gks->copy_func = gk_species_copy_range_static;
  gks->write_func = gk_species_write_static;
  gks->write_mom_func = gk_species_write_mom_static;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_static;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_static;
  gks->calc_L2norm_func = gk_species_calc_L2norm_static;
  gks->write_L2norm_func = gk_species_write_L2norm_static;
}

// End static function definitions.

void
gk_species_file_import_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file. Intended options include importing:
  //   1) ICs with same grid.
  //   2) ICs one dimensionality lower (e.g. 2x2v for 3x2v sim).
  //   3) ICs with same grid extents but different resolution (NYI).

  struct gkyl_rect_grid grid = gks->grid;
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
      for (int d=0; d<pdim; d++) {
        assert(grid_do.lower[d] == grid.lower[d]);
        assert(grid_do.upper[d] == grid.upper[d]);
      }
      // Check if the grid resolution is the same.
      for (int d=0; d<pdim; d++)
        same_res = same_res && (grid_do.cells[d] == grid.cells[d]);
    }
    else {
      // Assume the loaded file has one lower conf-space dimension.
      // Primarily meant for loading:
      //   - 1x2v for a 2x2v sim.
      //   - 2x2v for a 3x2v sim.
      assert(pdim_do == pdim-1);
      for (int d=0; d<cdim_do-1; d++) {
        assert(grid_do.lower[d] == grid.lower[d]);
        assert(grid_do.upper[d] == grid.upper[d]);
        assert(grid_do.cells[d] == grid.cells[d]);
        assert(grid_do.dx[d] == grid.dx[d]);
      }
      assert(grid_do.lower[cdim_do-1] == grid.lower[cdim-1]);
      assert(grid_do.upper[cdim_do-1] == grid.upper[cdim-1]);
      assert(grid_do.cells[cdim_do-1] == grid.cells[cdim-1]);
      assert(grid_do.dx[cdim_do-1] == grid.dx[cdim-1]);
      for (int d=0; d<vdim; d++) {
        assert(grid_do.lower[cdim_do+d] == grid.lower[cdim+d]);
        assert(grid_do.upper[cdim_do+d] == grid.upper[cdim+d]);
        assert(grid_do.cells[cdim_do+d] == grid.cells[cdim+d]);
        assert(grid_do.dx[cdim_do+d] == grid.dx[cdim+d]);
      }
    }

    struct gyrokinetic_output_meta meta =
      gk_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }
      );
    assert(strcmp(app->basis.id, meta.basis_type_nm) == 0);
    assert(poly_order == meta.poly_order);
    gkyl_grid_sub_array_header_release(&hdr);
  }

  // Donor basis.
  struct gkyl_basis basis_do;
  switch (app->basis.b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (poly_order > 1)
        gkyl_cart_modal_serendip(&basis_do, pdim_do, poly_order);
      else if (poly_order == 1)
        gkyl_cart_modal_gkhybrid(&basis_do, cdim_do, vdim_do); // p=2 in vparallel
      break;
    default:
      assert(false);
      break;
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
    for (int d=0; d<cdim_do-1; d++) cuts_do[d] = cuts_tar[d];
    cuts_do[cdim_do-1] = cuts_tar[cdim-1];
  }
  else {
    for (int d=0; d<cdim; d++) cuts_do[d] = cuts_tar[d];
  }
  // Set velocity space cuts to 1 as we do not use MPI in vel-space.
  for (int d=0; d<vdim; d++) cuts_do[cdim_do+d] = 1;

  struct gkyl_rect_decomp *decomp_do = gkyl_rect_decomp_new_from_cuts(pdim_do, cuts_do, &global_do);
  struct gkyl_comm* comm_do = gkyl_comm_split_comm(gks->comm, 0, decomp_do);

  // Donor local range.
  int my_rank = 0;
  gkyl_comm_get_rank(comm_do, &my_rank);

  struct gkyl_range local_ext_do, local_do;
  gkyl_create_ranges(&decomp_do->ranges[my_rank], ghost_do, &local_ext_do, &local_do);

  // Donor array.
  struct gkyl_array *fdo = mkarr(app->use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_host = app->use_gpu? mkarr(false, basis_do.num_basis, local_ext_do.volume)
                                            : gkyl_array_acquire(fdo);

  // Read donor field.
  struct gkyl_app_restart_status rstat;
  rstat.io_status = gkyl_comm_array_read(comm_do, &grid_do, &local_do, fdo_host, inp.file_name);
  if (app->use_gpu)
    gkyl_array_copy(fdo, fdo_host);

  if (pdim_do == pdim-1) {
    struct gkyl_translate_dim_gyrokinetic* transdim = gkyl_translate_dim_gyrokinetic_new(cdim_do,
      basis_do, cdim, app->basis, app->use_gpu);
    gkyl_translate_dim_gyrokinetic_advance(transdim, &local_do, &gks->local, fdo, gks->f);
    gkyl_translate_dim_gyrokinetic_release(transdim);
  }
  else {
    if (same_res) {
      gkyl_array_copy(gks->f, fdo);
    }
    else {
      // Interpolate the donor distribution to the target grid.
      struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(app->cdim, &app->basis,
        &grid_do, &grid, &local_do, &gks->local, ghost_do, app->use_gpu);
      gkyl_dg_interpolate_advance(interp, fdo, gks->f);
      gkyl_dg_interpolate_release(interp);
    }
  }

  if (inp.type == GKYL_IC_IMPORT_AF || inp.type == GKYL_IC_IMPORT_AF_B) {
    // Scale f by a conf-space factor.
    gkyl_proj_on_basis *proj_conf_scale = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 1, inp.conf_scale, inp.conf_scale_ctx);
    struct gkyl_array *xfac = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_array *xfac_ho = app->use_gpu? mkarr(false, app->confBasis.num_basis, app->local_ext.volume)
                                             : gkyl_array_acquire(xfac_ho);
    gkyl_proj_on_basis_advance(proj_conf_scale, 0.0, &app->local, xfac_ho);
    gkyl_array_copy(xfac, xfac_ho);
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, gks->f, xfac, gks->f, &app->local, &gks->local);
    gkyl_proj_on_basis_release(proj_conf_scale);
    gkyl_array_release(xfac_ho);
    gkyl_array_release(xfac);
  }
  if (inp.type == GKYL_IC_IMPORT_F_B || inp.type == GKYL_IC_IMPORT_AF_B) {
    // Add a phase factor to f.
    struct gk_proj proj_phase_add;
    gk_species_projection_init(app, gks, inp.phase_add, &proj_phase_add);
    gk_species_projection_calc(app, gks, &proj_phase_add, gks->fnew, 0.0);
    gkyl_array_accumulate_range(gks->f, 1.0, gks->fnew, &gks->local);
    gk_species_projection_release(app, &proj_phase_add);
  }

  gkyl_rect_decomp_release(decomp_do);
  gkyl_comm_release(comm_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_host);
}

void
gk_species_init(struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = gk_app_inp->cells[d];
    lower[d] = gk_app_inp->lower[d];
    upper[d] = gk_app_inp->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    cells[cdim+d] = gks->info.cells[d];
    lower[cdim+d] = gks->info.lower[d];
    upper[cdim+d] = gks->info.upper[d];
    ghost[cdim+d] = 0; // No ghost-cells in velocity space.

    // Only velocity space.
    cells_vel[d] = gks->info.cells[d];
    lower_vel[d] = gks->info.lower[d];
    upper_vel[d] = gks->info.upper[d];
    ghost_vel[d] = 0; // No ghost-cells in velocity space.
  }
  // Full phase space grid.
  gkyl_rect_grid_init(&gks->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&gks->grid, ghost, &gks->global_ext, &gks->global);
  
  // Velocity space grid.
  gkyl_rect_grid_init(&gks->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&gks->grid_vel, ghost_vel, &gks->local_ext_vel, &gks->local_vel);

  // Phase-space communicator.
  gks->comm = gkyl_comm_extend_comm(app->comm, &gks->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // Local = conf-local X local_vel.
  gkyl_range_ten_prod(&local, &app->local, &gks->local_vel);
  gkyl_create_ranges(&local, ghost, &gks->local_ext, &gks->local);

  // Velocity space mapping.
  gks->vel_map = gkyl_velocity_map_new(gks->info.mapc2p, gks->grid, gks->grid_vel,
    gks->local, gks->local_ext, gks->local_vel, gks->local_ext_vel, app->use_gpu);

  // Write out the velocity space mapping and its Jacobian.
  gkyl_velocity_map_write(gks->vel_map, gks->comm, app->name, gks->info.name);

  // Determine field-type.
  gks->gkfield_id = app->field->gkfield_id;
  if (gks->info.no_by) {
    gks->gkmodel_id = GKYL_GK_MODEL_NO_BY;
  }
  else {
    gks->gkmodel_id = GKYL_GK_MODEL_GEN_GEO;
  }

  // Allocate distribution function arrays.
  gks->f = mkarr(app->use_gpu, app->basis.num_basis, gks->local_ext.volume);

  gks->f_host = gks->f;
  if (app->use_gpu)
    gks->f_host = mkarr(false, app->basis.num_basis, gks->local_ext.volume);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (app->poly_order > 1) {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, app->poly_order);
  }
  else {
    if (vdim>1) {
      gkyl_cart_modal_gkhybrid(&surf_basis, cdim-1, vdim); // p=2 in vparallel
      gkyl_cart_modal_gkhybrid(&surf_quad_basis, cdim-1, vdim); 
    }
    else {
      gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2); // p=2 in vparallel
      gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2); 
    }
  }
  int alpha_surf_sz = (cdim+1)*surf_basis.num_basis;
  int sgn_alpha_surf_sz = (cdim+1)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // Allocate arrays to store fields:
  // 1. alpha_surf (surface phase space flux)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  gks->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, gks->local_ext.volume);
  gks->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, gks->local_ext.volume);
  gks->const_sgn_alpha = mk_int_arr(app->use_gpu, (cdim+1), gks->local_ext.volume);
  // 4. EM fields: phi and (if EM GK) Apar and d/dt Apar  
  gks->phi = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gks->apar = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gks->apardot = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);    

  gks->calc_gk_vars = gkyl_dg_calc_gyrokinetic_vars_new(&gks->grid, &app->confBasis, &app->basis, 
    gks->info.charge, gks->info.mass, gks->gkmodel_id, app->gk_geom, gks->vel_map, app->use_gpu);

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // modify it in GK_IWL BCs without modifying the app's.
  gks->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<gks->num_periodic_dir; ++d)
    gks->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) gks->bc_is_np[d] = true;
  for (int d=0; d<gks->num_periodic_dir; ++d)
    gks->bc_is_np[gks->periodic_dirs[d]] = false;

  bool is_zero_flux[2*GKYL_MAX_DIM] = {false};
  if (!gks->info.is_static) {
    // Store the BCs from the input file.
    for (int dir=0; dir<app->cdim; ++dir) {
      gks->lower_bc[dir].type = gks->upper_bc[dir].type = GKYL_SPECIES_COPY;
      if (gks->bc_is_np[dir]) {
        const struct gkyl_gyrokinetic_bcs *bc;
        if (dir == 0)
          bc = &gks->info.bcx;
        else if (dir == 1)
          bc = &gks->info.bcy;
        else
          bc = &gks->info.bcz;
  
        gks->lower_bc[dir] = bc->lower;
        gks->upper_bc[dir] = bc->upper;
      }
    }
  
    // Determine which directions are zero-flux.  By default
    // we do not have zero-flux boundary conditions in any direction.
    for (int dir=0; dir<app->cdim; ++dir) {
      if (gks->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir] = true;
      }
      if (gks->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir+pdim] = true;
      }
    }
  }

  struct gkyl_dg_gyrokinetic_auxfields aux_inp = { .alpha_surf = gks->alpha_surf, 
    .sgn_alpha_surf = gks->sgn_alpha_surf, .const_sgn_alpha = gks->const_sgn_alpha, 
    .phi = gks->phi, .apar = gks->apar, .apardot = gks->apardot };
  // Create collisionless solver.
  gks->slvr = gkyl_dg_updater_gyrokinetic_new(&gks->grid, &app->confBasis, &app->basis, 
    &app->local, &gks->local, is_zero_flux, gks->info.charge, gks->info.mass,
    gks->gkmodel_id, app->gk_geom, gks->vel_map, &aux_inp, app->use_gpu);

  // Acquire equation object.
  gks->eqn_gyrokinetic = gkyl_dg_updater_gyrokinetic_acquire_eqn(gks->slvr);

  // allocate data for density (for use in charge density accumulation and weak division for upar)
  gk_species_moment_init(app, gks, &gks->m0, "M0");

  // allocate data for diagnostic moments
  int ndm = gks->info.num_diag_moments;
  gks->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_species_moment_init(app, gks, &gks->moms[m], gks->info.diag_moments[m]);

  // allocate data for integrated moments
  gk_species_moment_init(app, gks, &gks->integ_moms, "Integrated");
  if (app->use_gpu) {
    gks->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
    gks->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[vdim+2]));
  } else {
    gks->red_integ_diag = gkyl_malloc(sizeof(double[vdim+2]));
    gks->red_integ_diag_global = gkyl_malloc(sizeof(double[vdim+2]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  gks->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  gks->is_first_integ_write_call = true;

  // Objects for L2 norm diagnostic.
  gks->integ_wfsq_op = gkyl_array_integrate_new(&gks->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED, app->use_gpu);
  if (app->use_gpu) {
    gks->L2norm_local = gkyl_cu_malloc(sizeof(double));
    gks->L2norm_global = gkyl_cu_malloc(sizeof(double));
  } else {
    gks->L2norm_local = gkyl_malloc(sizeof(double));
    gks->L2norm_global = gkyl_malloc(sizeof(double));
  }
  gks->L2norm = gkyl_dynvec_new(GKYL_DOUBLE, 1); // L2 norm.
  gks->is_first_L2norm_write_call = true;

  // initialize projection routine for initial conditions
  if (gks->info.init_from_file.type == 0)
    gk_species_projection_init(app, gks, gks->info.projection, &gks->proj_init);
  else
    gk_species_file_import_init(app, gks, gks->info.init_from_file);

  // set species source id
  gks->src = (struct gk_source) { };

  // Create skin/ghost ranges for applying BCs.
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&gks->lower_skin[dir], &gks->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &gks->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&gks->upper_skin[dir], &gks->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &gks->local_ext, ghost);
  }

  // Global skin and ghost ranges, only valid (i.e. volume>0) in ranges
  // abutting boundaries.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&gks->global_lower_skin[dir], &gks->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE, &gks->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&gks->global_upper_skin[dir], &gks->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE, &gks->global_ext, ghost);

    gkyl_sub_range_intersect(&gks->global_lower_skin[dir], &gks->local_ext, &gks->global_lower_skin[dir]);
    gkyl_sub_range_intersect(&gks->global_upper_skin[dir], &gks->local_ext, &gks->global_upper_skin[dir]);

    gkyl_sub_range_intersect(&gks->global_lower_ghost[dir], &gks->local_ext, &gks->global_lower_ghost[dir]);
    gkyl_sub_range_intersect(&gks->global_upper_ghost[dir], &gks->local_ext, &gks->global_upper_ghost[dir]);
  }

  // Initialize boundary fluxes for diagnostics and/or neut recycle bcs 
  gks->bflux = (struct gk_boundary_fluxes) { };
  gk_species_bflux_init(app, gks, &gks->bflux);

  if (!gks->info.is_static) {
    gk_species_new_dynamic(gk_app_inp, app, gks);
  }
  else {
    gk_species_new_static(gk_app_inp, app, gks);
  }

}

void
gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double t0)
{
  if (gks->info.init_from_file.type == 0)
    gk_species_projection_calc(app, gks, &gks->proj_init, gks->f, t0);

  // We are pre-computing source for now as it is time-independent.
  gk_species_source_calc(app, gks, &gks->src, t0);

  gkyl_dg_calc_gyrokinetic_vars_alpha_surf(gks->calc_gk_vars, 
    &app->local, &gks->local, &gks->local_ext, gks->phi,
    gks->alpha_surf, gks->sgn_alpha_surf, gks->const_sgn_alpha);
  gk_species_bflux_rhs(app, gks, &gks->bflux, gks->f, gks->f);
}

void
gk_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_species *gks_self, double t0)
{
  // IC setup step that depends on the IC of other species.

  if (app->field->init_phi_pol && gks_self->info.charge > 0.0) {
    // Scale the distribution function so its guiding center density is computed
    // given the polarization density and the guiding center density of other species.

    // Compute the polarization density.
    struct gkyl_array *npol = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(app->confBasis, app->grid, app->use_gpu);
    gkyl_gyrokinetic_pol_density_advance(npol_op, &app->local, app->field->epsilon, app->field->phi_pol, npol);
    gkyl_gyrokinetic_pol_density_release(npol_op);

    // Calculate the guiding center density of this species: (npol - q_other*n^G_other)/q_self.
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *gks = &app->species[i];
      gk_species_moment_calc(&gks->m0, gks->local, app->local, gks->f);
      if (strcmp(gks->info.name, gks_self->info.name)) {
        gkyl_array_accumulate(npol, -gks->info.charge, gks->m0.marr);
      }
    }
    gkyl_array_scale(npol, 1./gks_self->info.charge);

    // Scale the distribution function so it has this guiding center density.
    struct gkyl_dg_bin_op_mem *div_mem = app->use_gpu? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis)
                                                     : gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
    struct gkyl_array *den_mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    gkyl_dg_div_op_range(div_mem, app->confBasis,
      0, den_mod, 0, npol, 0, gks_self->m0.marr, &app->local);
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, gks_self->f,
      den_mod, gks_self->f, &app->local_ext, &gks_self->local_ext);

    gkyl_array_release(den_mod);
    gkyl_dg_bin_op_mem_release(div_mem);
    gkyl_array_release(npol);
  }

}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  return species->rhs_func(app, species, fin, rhs);
}
  
// Accummulate function for forward euler method.
void
gk_species_step_f(struct gk_species *species, struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  species->step_f_func(out, a, inp);
}

// Combine function for rk3 updates.
void
gk_species_combine(struct gk_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng)
{
  species->combine_func(out, c1, arr1, c2, arr2, rng);
}

// Copy function for rk3 updates.
void
gk_species_copy_range(struct gk_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  species->copy_func(out, inp, range);
}

// Apply boundary conditions to the distribution function.
void
gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  species->bc_func(app, species, f);
}

void
gk_species_coll_tm(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].lbo.collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_gyrokinetic_tm tm =
        gkyl_dg_updater_lbo_gyrokinetic_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
gk_species_tm(gkyl_gyrokinetic_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_dg_updater_gyrokinetic_tm tm =
      gkyl_dg_updater_gyrokinetic_get_tm(app->species[i].slvr);
    app->stat.species_rhs_tm += tm.gyrokinetic_tm;
  }
}

// write functions
void
gk_species_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (frame == 0) {
    gk_species_write_dynamic(app, gks, tm, frame);
  }
  else
    gks->write_func(app, gks, tm, frame);
}

void
gk_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (frame == 0) {
    gk_species_write_mom_dynamic(app, gks, tm, frame);
  }
  else
    gks->write_mom_func(app, gks, tm, frame);
}

void
gk_species_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  gks->calc_integrated_mom_func(app, gks, tm);
}

void
gk_species_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  gks->write_integrated_mom_func(app, gks);
}

void
gk_species_calc_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  gks->calc_L2norm_func(app, gks, tm);
}

void
gk_species_write_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  gks->write_L2norm_func(app, gks);
}

// release resources for species
void
gk_species_release(const gkyl_gyrokinetic_app* app, const struct gk_species *s)
{
  // release various arrays and species objects
  gkyl_array_release(s->f);

  if (s->info.init_from_file.type == 0)
    gk_species_projection_release(app, &s->proj_init);

  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  gkyl_velocity_map_release(s->vel_map);

  gkyl_array_release(s->alpha_surf);
  gkyl_array_release(s->sgn_alpha_surf);
  gkyl_array_release(s->const_sgn_alpha);
  gkyl_array_release(s->phi);
  gkyl_array_release(s->apar);
  gkyl_array_release(s->apardot);

  gkyl_dg_calc_gyrokinetic_vars_release(s->calc_gk_vars);

  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_gyrokinetic);
  gkyl_dg_updater_gyrokinetic_release(s->slvr);

  // release moment data
  gk_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    gk_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);

  // Release integrated moment memory.
  gk_species_moment_release(app, &s->integ_moms); 
  gkyl_dynvec_release(s->integ_diag);
  if (app->use_gpu) {
    gkyl_cu_free(s->red_integ_diag);
    gkyl_cu_free(s->red_integ_diag_global);
  }
  else {
    gkyl_free(s->red_integ_diag);
    gkyl_free(s->red_integ_diag_global);
  }

  gk_species_source_release(app, &s->src);

  gk_species_bflux_release(app, &s->bflux);

  s->release_func(app, s);
}
