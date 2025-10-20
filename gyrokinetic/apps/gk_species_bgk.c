#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

void 
gk_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_bgk_collisions *bgk)
{
  bgk->collision_id = s->info.collisions.collision_id;
  bgk->write_diagnostics = s->info.collisions.write_diagnostics;

  int cdim = app->cdim, vdim = app->vdim;
  // allocate nu and initialize it
  bgk->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  bgk->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->basis.num_basis, app->local_ext.volume);

  bgk->num_cross_collisions = s->info.collisions.num_cross_collisions;
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(bgk->self_nu, self_nu);
  gkyl_array_copy(bgk->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  bgk->spitzer_calc = 0;
  bgk->normNu = false;
  if (s->info.collisions.normNu) {
    bgk->normNu = true;
    double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
    double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
    double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
    double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;

    double bmag_mid = s->info.collisions.bmag_mid ? s->info.collisions.bmag_mid : app->bmag_ref;

    // Compute a minimum representable temperature based on the smallest dv in the grid.
    double dv_min[vdim];
    gkyl_velocity_map_reduce_dv_range(s->vel_map, GKYL_MIN, dv_min, s->vel_map->local_vel);

    double tpar_min = (s->info.mass/6.0)*pow(dv_min[0],2);
    double tperp_min = vdim>1 ? (bmag_mid/3.0)*dv_min[1] : tpar_min;
    bgk->vtsq_min = (tpar_min + 2.0*tperp_min)/(3.0*s->info.mass);

    bgk->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->basis, app->poly_order+1,
      nuFrac, 1.0, 1.0, app->use_gpu);
    bgk->self_nu_fac = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, s->info.collisions.n_ref, 
      s->info.mass, s->info.mass, s->info.charge, s->info.charge, 
      s->info.collisions.T_ref, s->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);

    // Create arrays for scaling collisionality by normalization factor
    // norm_nu is computed from Spitzer calc and is the normalization factor for the local
    // density and thermal velocity, norm_nu_sr = n/(vth_s^2 + vth_r^2)^(3/2)
    // nu_init is the inital collisionality profile, which must be stored so that at every time
    // time step the collisionality profile is properly scaled and the effects are not cumulative
    bgk->norm_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    bgk->nu_init = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_array_copy(bgk->nu_init, bgk->self_nu);
  }

  // Host-side copy for I/O
  bgk->nu_sum_host = bgk->nu_sum;
  if (app->use_gpu) {
    bgk->nu_sum_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }
  // Density and T/m for Spitzer nu
  bgk->m0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  bgk->vtsq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  bgk->nu_fmax = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
  // BGK updater (also computes stable timestep)
  bgk->up_bgk = gkyl_bgk_collisions_new(&app->basis, &s->basis, app->use_gpu);
}

void 
gk_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_bgk_collisions *bgk)
{  
  // set pointers to species we cross-collide with
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    bgk->collide_with[i] = gk_find_species(app, s->info.collisions.collide_with[i]);
    if (s->info.collisions.normNu) {
      double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
      double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
      double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
      double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
      double bmag_mid = s->info.collisions.bmag_mid ? s->info.collisions.bmag_mid : app->bmag_ref;
      bgk->cross_nu_fac[i] = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, bgk->collide_with[i]->info.collisions.n_ref, 
        s->info.mass, bgk->collide_with[i]->info.mass, s->info.charge, bgk->collide_with[i]->info.charge, 
        s->info.collisions.T_ref, bgk->collide_with[i]->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);
    }    
    bgk->other_m[i] = bgk->collide_with[i]->info.mass;
    bgk->other_moms[i] = bgk->collide_with[i]->lte.moms.marr; // other species LTE moment array
    bgk->other_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    bgk->cross_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    
    if (bgk->other_m[i] > s->info.mass) {
      gkyl_array_set(bgk->cross_nu[i], sqrt(2.), bgk->self_nu);
      gkyl_array_set(bgk->other_nu[i], sqrt(2.)*(s->info.mass)/(bgk->other_m[i]), bgk->self_nu);
    } else {
      gkyl_array_set(bgk->cross_nu[i], sqrt(2.)*(bgk->other_m[i])/(s->info.mass), bgk->collide_with[i]->bgk.self_nu);
      gkyl_array_set(bgk->other_nu[i], sqrt(2.), bgk->collide_with[i]->bgk.self_nu);
    }
    gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[i]);

    bgk->cross_moms[i] = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
    // Host-side copy for I/O of cross moments
    bgk->cross_moms_host[i] = bgk->cross_moms[i];
    if (app->use_gpu) {
      bgk->cross_moms_host[i] = mkarr(false, 3*app->basis.num_basis, app->local_ext.volume);
    }    
  }

  bgk->betaGreenep1 = 1.0;

  bgk->cross_bgk = gkyl_gyrokinetic_cross_prim_moms_bgk_new(&s->basis, &app->basis, app->use_gpu);
}

// computes moments, boundary corrections, and primitive moments
void
gk_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // compute needed Maxwellian moments (J*n, u_par, T/m) 
  gk_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  // divide out the Jacobian from the density
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, app->basis, 
    0, species->lte.moms.marr, 0, species->lte.moms.marr, 0, 
    app->gk_geom->geo_int.jacobgeo, &app->local);  

  // Calculate self_nu if using spitzer nu
  if (bgk->normNu) {
    gkyl_array_clear(bgk->nu_sum, 0.0);
    // Fetch n and T/m from Maxwellian moments computed for BGK update
    gkyl_array_set_offset(bgk->m0, 1.0, species->lte.moms.marr, 0*app->basis.num_basis);
    gkyl_array_set_offset(bgk->vtsq, 1.0, species->lte.moms.marr, 2*app->basis.num_basis);
    gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, 
      bgk->vtsq, bgk->vtsq_min, bgk->m0, bgk->vtsq, bgk->vtsq_min, bgk->self_nu_fac, bgk->self_nu);
    gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->self_nu);
  }
  
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

void
gk_species_bgk_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();
  
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    // Calculate cross_nu if using spitzer nu
    if (bgk->normNu) {
      gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, 
        bgk->vtsq, bgk->vtsq_min, bgk->collide_with[i]->bgk.m0, 
        bgk->collide_with[i]->bgk.vtsq, bgk->collide_with[i]->bgk.vtsq_min, 
        bgk->cross_nu_fac[i], bgk->cross_nu[i]);
      gkyl_array_set(bgk->other_nu[i], (species->info.mass)/(bgk->other_m[i]), bgk->cross_nu[i]);
      gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[i]);
    }

    gkyl_gyrokinetic_cross_prim_moms_bgk_advance(bgk->cross_bgk, &app->local, bgk->betaGreenep1, 
      species->info.mass, species->lte.moms.marr, bgk->other_m[i], bgk->other_moms[i], 
      bgk->cross_nu[i], bgk->other_nu[i], bgk->cross_moms[i]);

  }

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

void
gk_species_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute the self-collisions Maxwellian.
  gk_species_lte_from_moms(app, species, &species->lte, species->lte.moms.marr);

  // Multiply the Maxwellian by the configuration-space Jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, species->lte.f_lte, 
    app->gk_geom->geo_int.jacobgeo, species->lte.f_lte, &app->local, &species->local);

  // Obtain and accumulate the self-collisions nu*fmax
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, bgk->nu_fmax, 
    bgk->self_nu, species->lte.f_lte, &app->local, &species->local);

  // Cross-collisions nu*fmax.
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    // Compute the cross-collisions Maxwellian.
    gk_species_lte_from_moms(app, species, &species->lte, bgk->cross_moms[i]);

    // Multiply the Maxwellian by the configuration-space Jacobian.
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, species->lte.f_lte, 
      app->gk_geom->geo_int.jacobgeo, species->lte.f_lte, &app->local, &species->local);

    // Compute and accumulate nu*fmax.
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, species->lte.f_lte, 
      bgk->cross_nu[i], species->lte.f_lte, &app->local, &species->local);
    gkyl_array_accumulate(bgk->nu_fmax, 1.0, species->lte.f_lte);
  }

  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &species->local, 
    bgk->nu_sum, bgk->nu_fmax, fin, bgk->implicit_step, bgk->dt_implicit, rhs, species->cflrate);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_bgk_write_cross_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }, GKYL_GK_META_NONE, 0
  );

  if (gks->bgk.num_cross_collisions && gks->bgk.write_diagnostics) {
    // Compute self and cross BGK moments
    struct timespec wst = gkyl_wall_clock();
    gk_species_bgk_moms(app, gks, &gks->bgk, gks->f);
    gk_species_bgk_cross_moms(app, gks, &gks->bgk, gks->f);
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

    // Loop over number of cross collisions and write out cross moments for each cross collision
    for (int i=0; i<gks->bgk.num_cross_collisions; ++i) {
      struct timespec wtm = gkyl_wall_clock();
      // Construct the file handles for cross moments
      const char *fmt_cross = "%s-%s_cross_moms_%d_%d.gkyl";
      int sz_cross = gkyl_calc_strlen(fmt_cross, app->name, gks->info.name, frame);
      char fileNm_cross[sz_cross+1]; // ensures no buffer overflow
      snprintf(fileNm_cross, sizeof fileNm_cross, fmt_cross, app->name, gks->info.name, i, frame);

      if (app->use_gpu)
        gkyl_array_copy(gks->bgk.cross_moms_host[i], gks->bgk.cross_moms[i]);

      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
        gks->bgk.cross_moms_host[i], fileNm_cross);
      app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.n_diag_io += 1;    
    }
  }

  gk_array_meta_release(mt); 
}

void 
gk_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk)
{
  gkyl_array_release(bgk->self_nu);
  gkyl_array_release(bgk->nu_sum);
  gkyl_array_release(bgk->m0);
  gkyl_array_release(bgk->vtsq);

  if (app->use_gpu) {
    gkyl_array_release(bgk->nu_sum_host);
  }

  if (bgk->normNu) {
    gkyl_array_release(bgk->norm_nu);
    gkyl_array_release(bgk->nu_init);
    gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);
  }

  if (bgk->num_cross_collisions) {
    for (int i=0; i<bgk->num_cross_collisions; ++i) {
      gkyl_array_release(bgk->cross_nu[i]);
      gkyl_array_release(bgk->other_nu[i]);
      gkyl_array_release(bgk->cross_moms[i]);
      if (app->use_gpu) {
        gkyl_array_release(bgk->cross_moms_host[i]);
      }
    }
    gkyl_gyrokinetic_cross_prim_moms_bgk_release(bgk->cross_bgk);
  }

  gkyl_array_release(bgk->nu_fmax);
  gkyl_bgk_collisions_release(bgk->up_bgk);
 }
