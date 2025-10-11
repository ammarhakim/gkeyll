#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

void
gklbo_self_nu_calc_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
}

void
gklbo_self_nu_calc_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  // Calculate nu_ss(x,t).
  gk_species_moment_calc(&lbo->maxwellian_moms, species->local, app->local, fin);
  gkyl_array_set_offset(lbo->vtsq, 1.0, lbo->maxwellian_moms.marr, 2*app->basis.num_basis);
  gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, lbo->vtsq, lbo->vtsq_min,
    lbo->m0, lbo->vtsq, lbo->vtsq_min, lbo->self_norm_nu_fac, lbo->self_nu);

  gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);

  // Multiply moments and boundary corrections by nu.
  for (int d=0; d<3; d++)
    gkyl_dg_mul_op(app->basis, d, lbo->moms_buff, d, lbo->moms.marr, 0, lbo->self_nu);
  for (int d=0; d<2; d++)
    gkyl_dg_mul_op(app->basis, d, lbo->boundary_corrections_buff, d, lbo->boundary_corrections, 0, lbo->self_nu);
}

void
gklbo_cross_nu_calc_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_lbo_collisions *lbo)
{
}

void
gklbo_cross_nu_calc_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_lbo_collisions *lbo)
{
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    // Calculate nu_sr(x,t).
    gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, lbo->vtsq, lbo->vtsq_min,
      lbo->collide_with[i]->lbo.m0, lbo->collide_with[i]->lbo.vtsq, lbo->collide_with[i]->lbo.vtsq_min,
     	lbo->cross_norm_nu_fac[i], lbo->cross_nu[i]);

    gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);

    gkyl_array_set(lbo->self_mnu[i], s->info.mass, lbo->cross_nu[i]);

    // Multiply moments and boundary corrections by nu.
    for (int d=0; d<3; d++)
      gkyl_dg_mul_op(app->basis, d, lbo->moms_buff, d, lbo->moms.marr, 0, lbo->cross_nu[i]);
    for (int d=0; d<2; d++)
      gkyl_dg_mul_op(app->basis, d, lbo->boundary_corrections_buff, d, lbo->boundary_corrections, 0, lbo->cross_nu[i]);
  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

void
gklbo_cross_greene_num_constNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, int cross_coll_idx)
{
  // Calculate greene_num = n_s * m_r * nu_rs * n_r.
  gkyl_dg_mul_op_range(app->basis, 0, lbo->other_mnu_m0[cross_coll_idx], 0,
    lbo->other_mnu[cross_coll_idx], 0, lbo->collide_with[cross_coll_idx]->lbo.m0, &app->local);

  gkyl_dg_mul_op_range(app->basis, 0, lbo->greene_num, 0,
    lbo->other_mnu_m0[cross_coll_idx], 0, lbo->m0, &app->local);
}

void
gklbo_cross_greene_num_normNu(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, int cross_coll_idx)
{
  // Calculate greene_num = n_s * nu_sr * m_r * n_r * nu_rs.
  gkyl_array_set(lbo->other_mnu[cross_coll_idx], 1.0,
    lbo->collide_with[cross_coll_idx]->lbo.self_mnu[lbo->my_idx_in_other[cross_coll_idx]]);

  gkyl_dg_mul_op_range(app->basis, 0, lbo->other_mnu_m0[cross_coll_idx], 0,
    lbo->other_mnu[cross_coll_idx], 0, lbo->collide_with[cross_coll_idx]->lbo.m0, &app->local);

  gkyl_dg_mul_op_range(app->basis, 0, lbo->greene_num, 0,
    lbo->other_mnu_m0[cross_coll_idx], 0, lbo->m0, &app->local);

  gkyl_dg_mul_op_range(app->basis, 0, lbo->greene_num, 0,
    lbo->cross_nu[cross_coll_idx], 0, lbo->greene_num, &app->local);
}

void 
gk_species_lbo_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_lbo_collisions *lbo)
{
  lbo->collision_id = s->info.collisions.collision_id;
  lbo->write_diagnostics = s->info.collisions.write_diagnostics;
  lbo->num_cross_collisions = s->info.collisions.num_cross_collisions;
  
  int cdim = app->cdim, vdim = app->vdim;

  // Allocate nu and initialize it.
  lbo->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  lbo->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Allocate moments needed for LBO update.
  gk_species_moment_init(app, s, &lbo->moms, GKYL_F_MOMENT_M0M1M2, false);

  // Allocate needed arrays (boundary corrections, primitive moments, and nu*primitive moments)
  lbo->boundary_corrections = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);

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

    lbo->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->basis, app->poly_order+1,
      nuFrac, 1.0, 1.0, app->use_gpu);
    lbo->self_norm_nu_fac = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, s->info.collisions.n_ref,
      s->info.mass, s->info.mass, s->info.charge, s->info.charge, s->info.collisions.T_ref,
      s->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);

    // Allocate moments app used to compute vtsq.
    gk_species_moment_init(app, s, &lbo->maxwellian_moms, GKYL_F_MOMENT_MAXWELLIAN, false);

    lbo->vtsq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    // Set pointers to functions chosen at runtime.
    lbo->self_nu_calc = gklbo_self_nu_calc_normNu;
  }
  else {
    // Project the self-collisions collision frequency.
    struct gkyl_array *self_nu_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
    gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu_ho);
    gkyl_proj_on_basis_release(proj);
    gkyl_array_copy(lbo->self_nu, self_nu_ho);
    gkyl_array_release(self_nu_ho);

    gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);

    // Set pointers to functions chosen at runtime.
    lbo->self_nu_calc = gklbo_self_nu_calc_constNu;
  }

  if (!lbo->num_cross_collisions) {
    lbo->moms_buff = gkyl_array_acquire(lbo->moms.marr);
    lbo->boundary_corrections_buff = gkyl_array_acquire(lbo->boundary_corrections);
  }

  // Primitive moments in GK are (u_par, vtsq).
  lbo->prim_moms = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
  lbo->nu_prim_moms = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
  lbo->m0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  lbo->m2self = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Host-side copy for I/O.
  if (lbo->write_diagnostics) {
    lbo->nu_sum_host = lbo->nu_sum;
    lbo->prim_moms_host = lbo->prim_moms;
    lbo->nu_prim_moms_host = lbo->nu_prim_moms;
    if (app->use_gpu) {
      lbo->nu_sum_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      lbo->prim_moms_host = mkarr(false, 2*app->basis.num_basis, app->local_ext.volume);
      lbo->nu_prim_moms_host = mkarr(false, 2*app->basis.num_basis, app->local_ext.volume);    
    }
  }

  lbo->dg_div_mem = 0; // Memory for weak division.
  if (app->use_gpu)
    lbo->dg_div_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis);
  else
    lbo->dg_div_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);

  // Edge of velocity space corrections to momentum and energy 
  lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(&s->grid, 
    &app->basis, &s->basis, s->info.mass, s->vel_map, app->use_gpu);
  
  // Primitive moment calculator.
  lbo->coll_pcalc = gkyl_prim_lbo_gyrokinetic_calc_new(&s->grid, 
    &app->basis, &s->basis, &app->local, app->use_gpu);

  // LBO updater.
  struct gkyl_dg_lbo_gyrokinetic_drag_auxfields drag_inp = { .nuSum = lbo->nu_sum, 
    .nuPrimMomsSum = lbo->nu_prim_moms, .m2self = lbo->m2self };
  struct gkyl_dg_lbo_gyrokinetic_diff_auxfields diff_inp = { .nuSum = lbo->nu_sum, 
    .nuPrimMomsSum = lbo->nu_prim_moms, .m2self = lbo->m2self };
  lbo->coll_slvr = gkyl_dg_updater_lbo_gyrokinetic_new(&s->grid, 
    &app->basis, &s->basis, &app->local, &drag_inp, &diff_inp, s->info.mass, 
    s->info.skip_cell_threshold, app->gk_geom, s->vel_map,  app->use_gpu);
}

void 
gk_species_lbo_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_lbo_collisions *lbo)
{
  lbo->cross_nu_calc = gklbo_cross_nu_calc_constNu; // This method is empty.

  if (s->lbo.num_cross_collisions) {
    lbo->betaGreenep1 = 1.0; // Greene's beta factor + 1.
      
    lbo->cross_calc = gkyl_prim_lbo_gyrokinetic_cross_calc_new(&s->grid, 
      &app->basis, &s->basis, &app->local, app->use_gpu);

    lbo->prim_cross_m0deltas_op = gkyl_prim_cross_m0deltas_new(s->info.collisions.normNu,
      &app->basis, &app->local, lbo->betaGreenep1, app->use_gpu);
    
    lbo->cross_nu_prim_moms = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);

    lbo->my_idx_in_other = gkyl_malloc(lbo->num_cross_collisions*sizeof(int));
    for (int i=0; i<lbo->num_cross_collisions; ++i) {
      // Set pointers to species we cross-collide with.
      lbo->collide_with[i] = gk_find_species(app, s->info.collisions.collide_with[i]);
      lbo->my_idx_in_other[i] = -1;
      for (int j=0; j<lbo->collide_with[i]->lbo.num_cross_collisions; ++j) {
        if (0 == strcmp(s->info.name, lbo->collide_with[i]->info.collisions.collide_with[j])) {
          lbo->my_idx_in_other[i] = j;
          break;
        }
      }

      if (s->info.collisions.normNu) {
        double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
        double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
        double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
        double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
        double bmag_mid = s->info.collisions.bmag_mid ? s->info.collisions.bmag_mid : app->bmag_ref;
        lbo->cross_norm_nu_fac[i] = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, lbo->collide_with[i]->info.collisions.n_ref,
          s->info.mass, lbo->collide_with[i]->info.mass, s->info.charge, lbo->collide_with[i]->info.charge,
         	s->info.collisions.T_ref, lbo->collide_with[i]->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);
      }

      lbo->other_m[i] = lbo->collide_with[i]->info.mass;
      lbo->other_prim_moms[i] = lbo->collide_with[i]->lbo.prim_moms;
      lbo->other_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      lbo->cross_prim_moms[i] = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
      lbo->cross_nu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      
      if (lbo->other_m[i] > s->info.mass) {
        gkyl_array_set(lbo->cross_nu[i], sqrt(2.0), lbo->self_nu);
        gkyl_array_set(lbo->other_nu[i], sqrt(2.0)*(s->info.mass)/(lbo->other_m[i]), lbo->self_nu);
      } else {
        gkyl_array_set(lbo->cross_nu[i], sqrt(2.0)*(lbo->other_m[i])/(s->info.mass), lbo->collide_with[i]->lbo.self_nu);
        gkyl_array_set(lbo->other_nu[i], sqrt(2.0), lbo->collide_with[i]->lbo.self_nu);
      }
      
      gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);

      lbo->other_mnu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      lbo->other_mnu_m0[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      lbo->self_mnu[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      gkyl_array_set(lbo->self_mnu[i], s->info.mass, lbo->cross_nu[i]);
      gkyl_array_set(lbo->other_mnu[i], lbo->other_m[i], lbo->other_nu[i]);
    }

    lbo->greene_num = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    lbo->greene_den = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    lbo->greene_factor = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    if (s->info.collisions.normNu) {
      lbo->moms_buff = mkarr(app->use_gpu, lbo->moms.marr->ncomp, lbo->moms.marr->size);
      lbo->boundary_corrections_buff = mkarr(app->use_gpu, lbo->boundary_corrections->ncomp, lbo->boundary_corrections->size);

      // Set pointers to functions chosen at runtime.
      lbo->cross_nu_calc = gklbo_cross_nu_calc_normNu;
      lbo->cross_greene_num = gklbo_cross_greene_num_normNu;
    }
    else {
      lbo->moms_buff = gkyl_array_acquire(lbo->moms.marr);
      lbo->boundary_corrections_buff = gkyl_array_acquire(lbo->boundary_corrections);

      // Set pointers to functions chosen at runtime.
      lbo->cross_nu_calc = gklbo_cross_nu_calc_constNu;
      lbo->cross_greene_num = gklbo_cross_greene_num_constNu;
    }
  }
}

void
gk_species_lbo_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute J*M0, J*M1, J*M2 moments and separate our M0 and M2.
  gk_species_moment_calc(&lbo->moms, species->local, app->local, fin);
  gkyl_dg_div_op_range(lbo->dg_div_mem, app->basis, 0, lbo->m0,
    0, lbo->moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local);  
  gkyl_array_set_offset_range(lbo->m2self, 1.0, lbo->moms.marr, 2*app->basis.num_basis, &app->local);
  
  // Construct boundary corrections.
  gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
    &species->local, &app->local, fin, lbo->boundary_corrections);

  // Calculate nu_ss and multibly moms and corrections by it.
  lbo->self_nu_calc(app, species, lbo, fin);

  // Construct primitive moments.
  gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, &app->local, 
    lbo->moms_buff, lbo->boundary_corrections_buff, lbo->self_nu, lbo->prim_moms);

  // Scale upar and vtSq by nu.
  for (int d=0; d<2; d++)
    gkyl_dg_mul_op(app->basis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

void
gk_species_lbo_cross_nu(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_lbo_collisions *lbo)
{
  lbo->cross_nu_calc(app, s, lbo);
}

void
gk_species_lbo_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  // Compute primitive moments for cross-species collisions.
  struct timespec wst = gkyl_wall_clock();
  
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    // Calculate greene_num = n_s * nu_sr * m_r * n_r * nu_rs (normNu).
    //                      = n_s * m_r * n_r * nu_rs (constNu).
    lbo->cross_greene_num(app, species, lbo, i);

    // greene_den = m_s * nu_sr * n_s + m_r * nu_rs * n_r.
    gkyl_dg_mul_op_range(app->basis, 0, lbo->greene_den, 0,
      lbo->self_mnu[i], 0, lbo->m0, &app->local);
    gkyl_array_accumulate(lbo->greene_den, 1.0, lbo->other_mnu_m0[i]);

    // greene_fac = 2*(beta+1) * greene_num/greene_den
    //            = 2*(beta+1) * n_s * nu_sr * m_r * n_r * nu_rs / (m_s * nu_sr * n_s + m_r * nu_rs * n_r) (normNu).
    //            = 2*(beta+1) * n_s * m_r * n_r * nu_rs / (m_s * nu_sr * n_s + m_r * nu_rs * n_r) (constNu).
    gkyl_prim_cross_m0deltas_advance(lbo->prim_cross_m0deltas_op, 
      species->info.mass, lbo->m0, lbo->cross_nu[i],
      lbo->other_m[i], lbo->collide_with[i]->lbo.m0, lbo->collide_with[i]->lbo.cross_nu[lbo->my_idx_in_other[i]],
      lbo->greene_factor);

    // Compute cross primitive moments.
    gkyl_prim_lbo_cross_calc_advance(lbo->cross_calc, &app->local, lbo->greene_factor, species->info.mass,
      lbo->moms_buff, lbo->prim_moms, lbo->other_m[i], lbo->collide_with[i]->lbo.moms_buff,
      lbo->other_prim_moms[i], lbo->boundary_corrections_buff, lbo->cross_nu[i], lbo->cross_prim_moms[i]);

    // Scale upar_{sr} and vtSq_{sr} by nu_{sr}
    for (int d=0; d<2; d++)
      gkyl_dg_mul_op(app->basis, d, lbo->cross_nu_prim_moms, d, lbo->cross_prim_moms[i], 0, lbo->cross_nu[i]);

    gkyl_array_accumulate(lbo->nu_prim_moms, 1.0, lbo->cross_nu_prim_moms);

  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

void
gk_species_lbo_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // Accumulate update due to collisions onto rhs.
  gkyl_dg_updater_lbo_gyrokinetic_advance(lbo->coll_slvr, &species->local,
    fin, species->cflrate, rhs);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_lbo_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (gks->lbo.collision_id == GKYL_LBO_COLLISIONS && gks->lbo.write_diagnostics) {
    struct timespec wst = gkyl_wall_clock();
    // Compute primitive moments.
    const struct gkyl_array *fin[app->num_species];
    gk_species_lbo_moms(app, gks, &gks->lbo, gks->f);

    // Compute cross primitive moments.
    if (gks->lbo.num_cross_collisions)
      gk_species_lbo_cross_moms(app, gks, &gks->lbo, gks->f);
    
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

    struct timespec wtm = gkyl_wall_clock();
    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id
      }, GKYL_GK_META_NONE, 0
    );

    // Construct the file handles for collision frequency and primitive moments.
    const char *fmt_prim = "%s-%s_prim_moms_%d.gkyl";
    int sz_prim = gkyl_calc_strlen(fmt_prim, app->name, gks->info.name, frame);
    char fileNm_prim[sz_prim+1]; // Ensures no buffer overflow.
    snprintf(fileNm_prim, sizeof fileNm_prim, fmt_prim, app->name, gks->info.name, frame);

    // Copy data from device to host before writing it out.
    if (app->use_gpu)  
      gkyl_array_copy(gks->lbo.prim_moms_host, gks->lbo.prim_moms);

    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->lbo.prim_moms_host, fileNm_prim);
    app->stat.n_diag_io += 1;   

    // Write out nu_sum and nu_prim_moms.
    const char *fmt = "%s-%s_nu_sum_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
    
    const char *fmt_nu_prim = "%s-%s_nu_prim_moms_%d.gkyl";
    int sz_nu_prim = gkyl_calc_strlen(fmt_nu_prim, app->name, gks->info.name, frame);
    char fileNm_nu_prim[sz_nu_prim+1]; // ensures no buffer overflow
    snprintf(fileNm_nu_prim, sizeof fileNm_nu_prim, fmt_nu_prim, app->name, gks->info.name, frame);
    
    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gks->lbo.nu_sum_host, gks->lbo.nu_sum);
      gkyl_array_copy(gks->lbo.nu_prim_moms_host, gks->lbo.nu_prim_moms);
    }
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->lbo.nu_sum_host, fileNm);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->lbo.nu_prim_moms_host, fileNm_nu_prim);
    app->stat.n_diag_io += 2;

    gk_array_meta_release(mt); 
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
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
  gkyl_array_release(lbo->m2self);

  if (lbo->write_diagnostics) {
    if (app->use_gpu) {
      gkyl_array_release(lbo->nu_sum_host);
      gkyl_array_release(lbo->prim_moms_host);
      gkyl_array_release(lbo->nu_prim_moms_host);    
    }
  }

  gk_species_moment_release(app, &lbo->moms);
  gkyl_dg_bin_op_mem_release(lbo->dg_div_mem);

  gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
  gkyl_prim_lbo_calc_release(lbo->coll_pcalc);

  if (lbo->normNu) {
    gkyl_spitzer_coll_freq_release(lbo->spitzer_calc);
    gk_species_moment_release(app, &lbo->maxwellian_moms);
    gkyl_array_release(lbo->vtsq);
  }
  gkyl_array_release(lbo->boundary_corrections_buff);
  gkyl_array_release(lbo->moms_buff);

  if (lbo->num_cross_collisions) {
    gkyl_free(lbo->my_idx_in_other);
    gkyl_array_release(lbo->cross_nu_prim_moms);
    for (int i=0; i<lbo->num_cross_collisions; ++i) {
      gkyl_array_release(lbo->cross_prim_moms[i]);
      gkyl_array_release(lbo->cross_nu[i]);
      gkyl_array_release(lbo->other_nu[i]);
      gkyl_array_release(lbo->self_mnu[i]);
      gkyl_array_release(lbo->other_mnu[i]);
      gkyl_array_release(lbo->other_mnu_m0[i]);
    }
    gkyl_array_release(lbo->greene_num);
    gkyl_array_release(lbo->greene_den);
    gkyl_array_release(lbo->greene_factor);
    gkyl_prim_cross_m0deltas_release(lbo->prim_cross_m0deltas_op);
    gkyl_prim_lbo_cross_calc_release(lbo->cross_calc);
  }
  gkyl_dg_updater_lbo_gyrokinetic_release(lbo->coll_slvr);
}
