#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

void 
gk_species_lbo_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_lbo_collisions *lbo)
{
  lbo->collision_id = s->info.collisions.collision_id;
  lbo->write_diagnostics = s->info.collisions.write_diagnostics;

  int cdim = app->cdim, vdim = app->vdim;

  // Allocate nu and initialize it.
  lbo->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  lbo->num_cross_collisions = s->info.collisions.num_cross_collisions;
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(lbo->self_nu, self_nu);
  gkyl_array_copy(lbo->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  lbo->spitzer_calc = 0;
  lbo->normNu = false;
  if (s->info.collisions.normNu) {
    lbo->normNu = true;
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
    lbo->vtsq_min = (tpar_min + 2.0*tperp_min)/(3.0*s->info.mass);

    lbo->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->confBasis, app->poly_order+1,
      nuFrac, 1.0, 1.0, app->use_gpu);
    lbo->self_nu_fac = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, s->info.collisions.n_ref,
      s->info.mass, s->info.mass, s->info.charge, s->info.charge, s->info.collisions.T_ref,
      s->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);

    // Create arrays for scaling collisionality by normalization factor
    // norm_nu is computed from Spitzer calc and is the normalization factor for the local
    // density and thermal velocity, norm_nu_sr = n/(vtSq_s + vtSq_r)^(3/2)
    // nu_init is the inital collisionality profile, which must be stored so that at every time
    // time step the collisionality profile is properly scaled and the effects are not cumulative
    lbo->norm_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->nu_init = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_copy(lbo->nu_init, lbo->self_nu);
  }

  // Allocate needed arrays (boundary corrections, primitive moments, and nu*primitive moments)
  lbo->boundary_corrections = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);

  // Primitive moments in GK are (u_par, vtSq)
  lbo->prim_moms = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_prim_moms = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
  lbo->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->vtsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->m2self = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // Host-side copy for I/O
  lbo->nu_sum_host = lbo->nu_sum;
  lbo->prim_moms_host = lbo->prim_moms;
  lbo->nu_prim_moms_host = lbo->nu_prim_moms;
  if (app->use_gpu) {
    lbo->nu_sum_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    lbo->prim_moms_host = mkarr(false, 2*app->confBasis.num_basis, app->local_ext.volume);
    lbo->nu_prim_moms_host = mkarr(false, 2*app->confBasis.num_basis, app->local_ext.volume);    
  }

  // Allocate moments needed for LBO update.
  gk_species_moment_init(app, s, &lbo->moms, "ThreeMoments");

  lbo->dg_div_mem = 0; // Memory for weak division.
  if (app->use_gpu)
    lbo->dg_div_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
  else
    lbo->dg_div_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);

  // Edge of velocity space corrections to momentum and energy 
  lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(&s->grid, 
    &app->confBasis, &app->basis, s->info.mass, s->vel_map, app->use_gpu);
  
  // primitive moment calculator
  lbo->coll_pcalc = gkyl_prim_lbo_gyrokinetic_calc_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, app->use_gpu);

  // LBO updater
  struct gkyl_dg_lbo_gyrokinetic_drag_auxfields drag_inp = { .nuSum = lbo->nu_sum, 
    .nuPrimMomsSum = lbo->nu_prim_moms, .m2self = lbo->m2self };
  struct gkyl_dg_lbo_gyrokinetic_diff_auxfields diff_inp = { .nuSum = lbo->nu_sum, 
    .nuPrimMomsSum = lbo->nu_prim_moms, .m2self = lbo->m2self };
  lbo->coll_slvr = gkyl_dg_updater_lbo_gyrokinetic_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, &drag_inp, &diff_inp, s->info.mass,
    app->gk_geom, s->vel_map,  app->use_gpu);
}

void 
gk_species_lbo_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_lbo_collisions *lbo)
{
  lbo->cross_calc = gkyl_prim_lbo_gyrokinetic_cross_calc_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, app->use_gpu);
  
  lbo->cross_nu_prim_moms = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);

  // set pointers to species we cross-collide with
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    lbo->collide_with[i] = gk_find_species(app, s->info.collisions.collide_with[i]);
    if (s->info.collisions.normNu) {
      double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
      double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
      double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
      double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
      double bmag_mid = s->info.collisions.bmag_mid ? s->info.collisions.bmag_mid : app->bmag_ref;
      lbo->cross_nu_fac[i] = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, lbo->collide_with[i]->info.collisions.n_ref,
        s->info.mass, lbo->collide_with[i]->info.mass, s->info.charge, lbo->collide_with[i]->info.charge,
       	s->info.collisions.T_ref, lbo->collide_with[i]->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);
    }

    lbo->other_m[i] = lbo->collide_with[i]->info.mass;
    lbo->other_prim_moms[i] = lbo->collide_with[i]->lbo.prim_moms;
    lbo->other_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_prim_moms[i] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_num[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_den[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_factor[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    
    if (lbo->other_m[i] > s->info.mass) {
      gkyl_array_set(lbo->cross_nu[i], sqrt(2.0), lbo->self_nu);
      gkyl_array_set(lbo->other_nu[i], sqrt(2.0)*(s->info.mass)/(lbo->other_m[i]), lbo->self_nu);
    } else {
      gkyl_array_set(lbo->cross_nu[i], sqrt(2.0)*(lbo->other_m[i])/(s->info.mass), lbo->collide_with[i]->lbo.self_nu);
      gkyl_array_set(lbo->other_nu[i], sqrt(2.0), lbo->collide_with[i]->lbo.self_nu);
    }
    
    gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);

    lbo->other_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->other_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    lbo->self_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->self_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    gkyl_array_set(lbo->self_mnu[i], s->info.mass, lbo->cross_nu[i]);
    gkyl_array_set(lbo->other_mnu[i], lbo->other_m[i], lbo->other_nu[i]);
  }
    
  lbo->betaGreenep1 = 1.0;
}

// computes moments, boundary corrections, and primitive moments
void
gk_species_lbo_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute needed moments.
  gk_species_moment_calc(&lbo->moms, species->local, app->local, fin);

  gkyl_dg_div_op_range(lbo->dg_div_mem, app->confBasis, 0, lbo->m0,
    0, lbo->moms.marr, 0, app->gk_geom->jacobgeo, &app->local);  

  gkyl_array_set_offset_range(lbo->m2self, 1.0, lbo->moms.marr, 2*app->confBasis.num_basis, &app->local);
  
  if (app->use_gpu) {
    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance_cu(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // construct primitive moments  
    gkyl_prim_lbo_calc_advance_cu(lbo->coll_pcalc, &app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->prim_moms);
  } 
  else {
    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // construct primitive moments  
    gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, &app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->prim_moms);
  }

  // Calculate self_nu if using spitzer-nu
  if (lbo->normNu) {
    gkyl_array_set_offset(lbo->vtsq, 1.0, lbo->prim_moms, 1*app->confBasis.num_basis);
    gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, lbo->vtsq, lbo->vtsq_min,
      lbo->m0, lbo->vtsq, lbo->vtsq_min, lbo->self_nu_fac, lbo->self_nu);
    gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);
  }

  // Scale upar and vtSq by nu.
  for (int d=0; d<2; d++)
    gkyl_dg_mul_op(app->confBasis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// computes moments from cross-species collisions
void
gk_species_lbo_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();
  
  wst = gkyl_wall_clock();  
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->self_mnu_m0[i], 0,
      lbo->self_mnu[i], 0, lbo->m0, &app->local);
    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->other_mnu_m0[i], 0,
      lbo->other_mnu[i], 0, lbo->collide_with[i]->lbo.m0, &app->local);

    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->greene_num[i], 0,
      lbo->other_mnu_m0[i], 0, lbo->m0, &app->local);

    gkyl_array_set(lbo->greene_den[i], 1.0, lbo->self_mnu_m0[i]);
    gkyl_array_accumulate(lbo->greene_den[i], 1.0, lbo->other_mnu_m0[i]);

    gkyl_dg_div_op_range(lbo->dg_div_mem, app->confBasis, 0, lbo->greene_factor[i], 0,
      lbo->greene_num[i], 0, lbo->greene_den[i], &app->local);
    gkyl_array_scale(lbo->greene_factor[i], 2*lbo->betaGreenep1);

    if (app->use_gpu)
      gkyl_prim_lbo_cross_calc_advance_cu(lbo->cross_calc,
        &app->local, 
        lbo->greene_factor[i], 
        species->info.mass, lbo->moms.marr, lbo->prim_moms, 
        lbo->other_m[i], lbo->collide_with[i]->lbo.moms.marr, lbo->other_prim_moms[i],
        lbo->boundary_corrections, 
        lbo->cross_prim_moms[i]);
    else 
      gkyl_prim_lbo_cross_calc_advance(lbo->cross_calc,
        &app->local, 
        lbo->greene_factor[i], 
        species->info.mass, lbo->moms.marr, lbo->prim_moms, 
        lbo->other_m[i], lbo->collide_with[i]->lbo.moms.marr, lbo->other_prim_moms[i],
        lbo->boundary_corrections, 
        lbo->cross_prim_moms[i]);

    // Calculate cross nu if using spitzer nu
    if (lbo->normNu) {
      gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, lbo->vtsq, lbo->vtsq_min,
        lbo->collide_with[i]->lbo.m0, lbo->collide_with[i]->lbo.vtsq, lbo->collide_with[i]->lbo.vtsq_min,
       	lbo->cross_nu_fac[i], lbo->cross_nu[i]);
      gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);
    }
    
    // Scale upar_{rs} and vtSq_{rs} by nu_{rs}
    for (int d=0; d<2; d++)
      gkyl_dg_mul_op(app->confBasis, d, lbo->cross_nu_prim_moms, d, lbo->cross_prim_moms[i], 0, lbo->cross_nu[i]);

    gkyl_array_accumulate(lbo->nu_prim_moms, 1.0, lbo->cross_nu_prim_moms);

  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// updates the collision terms in the rhs
void
gk_species_lbo_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // accumulate update due to collisions onto rhs
  gkyl_dg_updater_lbo_gyrokinetic_advance(lbo->coll_slvr, &species->local,
    fin, species->cflrate, rhs);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_lbo_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (gks->lbo.collision_id == GKYL_LBO_COLLISIONS && gks->lbo.write_diagnostics) {
    struct timespec wst = gkyl_wall_clock();

    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    // Construct the file handles for collision frequency and primitive moments
    const char *fmt_prim = "%s-%s_prim_moms_%d.gkyl";
    int sz_prim = gkyl_calc_strlen(fmt_prim, app->name, gks->info.name, frame);
    char fileNm_prim[sz_prim+1]; // ensures no buffer overflow
    snprintf(fileNm_prim, sizeof fileNm_prim, fmt_prim, app->name, gks->info.name, frame);

    // Compute primitive moments
    const struct gkyl_array *fin[app->num_species];
    gk_species_lbo_moms(app, gks, &gks->lbo, gks->f);

    // copy data from device to host before writing it out
    if (app->use_gpu) {  
      gkyl_array_copy(gks->lbo.prim_moms_host, gks->lbo.prim_moms);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->lbo.prim_moms_host, fileNm_prim);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;

    // Uncomment the following to write out nu_sum and nu_prim_moms
    /*const char *fmt = "%s-%s_nu_sum_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
    
    const char *fmt_nu_prim = "%s-%s_nu_prim_moms_%d.gkyl";
    int sz_nu_prim = gkyl_calc_strlen(fmt_nu_prim, app->name, gks->info.name, frame);
    char fileNm_nu_prim[sz_nu_prim+1]; // ensures no buffer overflow
    snprintf(fileNm_nu_prim, sizeof fileNm_nu_prim, fmt_nu_prim, app->name, gks->info.name, frame);
    
    if (gks->lbo.num_cross_collisions)
      gk_species_lbo_cross_moms(app, gk_s, &gks->lbo, gks->f);
    
    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gks->lbo.nu_sum_host, gks->lbo.nu_sum);
      gkyl_array_copy(gks->lbo.nu_prim_moms_host, gks->lbo.nu_prim_moms);
    }
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->lbo.nu_sum_host, fileNm);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->lbo.nu_prim_moms_host, fileNm_nu_prim);*/

    gk_array_meta_release(mt); 
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void 
gk_species_lbo_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lbo_collisions *lbo)
{
  gkyl_array_release(lbo->boundary_corrections);
  gkyl_array_release(lbo->prim_moms);
  gkyl_array_release(lbo->self_nu);
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_prim_moms);
  gkyl_array_release(lbo->m0);
  gkyl_array_release(lbo->vtsq);
  gkyl_array_release(lbo->m2self);

  if (app->use_gpu) {
    gkyl_array_release(lbo->nu_sum_host);
    gkyl_array_release(lbo->prim_moms_host);
    gkyl_array_release(lbo->nu_prim_moms_host);    
  }

  gk_species_moment_release(app, &lbo->moms);
  gkyl_dg_bin_op_mem_release(lbo->dg_div_mem);

  gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
  gkyl_prim_lbo_calc_release(lbo->coll_pcalc);

  if (lbo->normNu) {
    gkyl_array_release(lbo->norm_nu);
    gkyl_array_release(lbo->nu_init);
    gkyl_spitzer_coll_freq_release(lbo->spitzer_calc);
  }

  if (lbo->num_cross_collisions) {
    gkyl_array_release(lbo->cross_nu_prim_moms);
    for (int i=0; i<lbo->num_cross_collisions; ++i) {
      gkyl_array_release(lbo->cross_prim_moms[i]);
      gkyl_array_release(lbo->cross_nu[i]);
      gkyl_array_release(lbo->other_nu[i]);
      gkyl_array_release(lbo->self_mnu[i]);
      gkyl_array_release(lbo->self_mnu_m0[i]);
      gkyl_array_release(lbo->other_mnu[i]);
      gkyl_array_release(lbo->other_mnu_m0[i]);
      gkyl_array_release(lbo->greene_num[i]);
      gkyl_array_release(lbo->greene_den[i]);
      gkyl_array_release(lbo->greene_factor[i]);
    }
    gkyl_prim_lbo_cross_calc_release(lbo->cross_calc);
  }
  gkyl_dg_updater_lbo_gyrokinetic_release(lbo->coll_slvr);
}
