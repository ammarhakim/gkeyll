#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_neut_species_rrs_cross_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_recycle_react_scale *rrs, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  struct timespec wst = gkyl_wall_clock();    

  struct gk_species *gks_elc = &app->species[rrs->elc_idx]; 

  // Compute electron Maxwellian moments (J*n, u_par, T/m).
  gk_species_moment_calc(&gks_elc->lte.moms, 
    gks_elc->local, app->local, fin[rrs->elc_idx]);

  // Divide the electron density by the Jacobian.
  gkyl_dg_div_op_range(gks_elc->lte.moms.mem_geo, app->basis, 0, gks_elc->lte.moms.marr,
    0, gks_elc->lte.moms.marr, 0, app->gk_geom->jacobgeo, &app->local); 

  // Compute ionization reactivity <sigma v>_iz.
  gkyl_dg_iz_coll(rrs->iz_react_calc, gks_elc->lte.moms.marr, 
    rrs->dndt_react, rrs->dndt_react, rrs->reactivity, 0);
  
  app->stat.neut_species_react_mom_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_rrs_cross_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_recycle_react_scale *rrs, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  // Do nothing.
}

static void
gk_neut_species_rrs_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_react_scale *rrs, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();  

  struct gk_species *gks_elc = &app->species[rrs->elc_idx];

  // Compute (J*n_neut)*n_elc*<sigma v>_iz.
  gkyl_dg_mul_op_range(app->basis, 0, rrs->dndt_react,
    0, rrs->Jm0_init, 0, rrs->reactivity, &app->local);  
  gkyl_dg_mul_op_range(app->basis, 0, rrs->dndt_react,
    0, gks_elc->lte.moms.marr, 0, rrs->dndt_react, &app->local);  

  // Volume integrate the reaction contribution.
  gkyl_array_integrate_advance(rrs->integrate_op, rrs->dndt_react, 1.0, 0, &app->local, 0, rrs->react_vol_integ_local);
  // Reduce over MPI processes.
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, 
    rrs->react_vol_integ_local, rrs->react_vol_integ_global);
  if (app->use_gpu)
    gkyl_cu_memcpy(&rrs->react_vol_integ, rrs->react_vol_integ_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&rrs->react_vol_integ, rrs->react_vol_integ_global, sizeof(double));

  app->stat.neut_species_react_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_rrs_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_react_scale *rrs, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // Do nothing.
}

static void
gk_neut_species_rrs_apply_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_recycle_react_scale *rrs, struct gkyl_array *fin)
{
  struct gk_species *gks_ion = &app->species[rrs->ion_idx];

  if (rrs->react_vol_integ > 0.0) {
    gkyl_array_clear(rrs->dndt_react, 0.0);

    double bflux_intm0_local_ho = 0.0;
    for (int j=0; j<rrs->num_boundaries; ++j) {
      // Add integrated M0 moments of boundary fluxes.
      gk_species_bflux_get_flux_mom(&gks_ion->bflux, rrs->boundaries_dir[j], rrs->boundaries_edge[j],
        GKYL_F_MOMENT_M0, rrs->dndt_react, &rrs->boundaries_conf_ghost[j]);
      gkyl_array_integrate_advance(rrs->integrate_op, rrs->dndt_react, 1.0, 0,
        &rrs->boundaries_conf_ghost[j], 0, rrs->bflux_m0_vol_integ_local);

      double bflux_m0_vol_integ_local_ho;
      if (app->use_gpu)
        gkyl_cu_memcpy(&bflux_m0_vol_integ_local_ho, rrs->bflux_m0_vol_integ_local, sizeof(double), GKYL_CU_MEMCPY_D2H);
      else
        memcpy(&bflux_m0_vol_integ_local_ho, rrs->bflux_m0_vol_integ_local, sizeof(double));

      bflux_intm0_local_ho += bflux_m0_vol_integ_local_ho;
    }

    double bflux_intm0_global_ho;
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &bflux_intm0_local_ho, &bflux_intm0_global_ho);

    double neut_scaling_fac = rrs->recycling_coeff * bflux_intm0_global_ho / rrs->react_vol_integ;

    // Divide by the present J*rho, and multiply by neut_scaling_fac*mass*Jm0_init.
    gkyl_array_set_offset_range(rrs->dndt_react, 1.0, fin, 0, &app->local);
    for (int i=0; i<ns->num_moments;++i)
      gkyl_dg_div_op_range(gks_ion->lte.moms.mem_geo, app->basis, i, fin,
        i, fin, 0, rrs->dndt_react, &app->local); 
  
    for (int i=0; i<ns->num_moments;++i)
      gkyl_dg_mul_op_range(app->basis, i, fin,
        i, fin, 0, rrs->Jm0_init, &app->local);  

    gkyl_array_scale(fin, neut_scaling_fac*ns->info.mass);
  }
}

static void
gk_neut_species_rrs_apply_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_recycle_react_scale *rrs, struct gkyl_array *fin)
{
  // Do nothing.
}

static void
gk_neut_species_rrs_write_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_recycle_react_scale *rrs, int ridx, double tm, int frame)
{
  // React diagnostics usually written from gk_species.
  // In the case of static gk_species, write_diagnostics flag
  // can be used to check reaction rates from gk_neut_species.
  struct timespec wst = gkyl_wall_clock();
//  // Compute reaction rate
//  const struct gkyl_array *fin[app->num_species];
//  const struct gkyl_array *fin_neut[app->num_neut_species];
//  for (int i=0; i<app->num_species; ++i) {
//    fin[i] = app->species[i].f;
//  }
//  for (int i=0; i<app->num_neut_species; ++i) {
//    fin_neut[i] = app->neut_species[i].f;
//  }
//  gk_neut_species_recycle_react_scale_cross_moms(app, gkns, gkr, fin, fin_neut);
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  
  struct timespec wtm = gkyl_wall_clock();
  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gk_neut_species_rrs_write_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_recycle_react_scale *rrs, int ridx, double tm, int frame)
{
  // Do nothing
}

void 
gk_neut_species_recycle_react_scale_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns, 
  struct gk_recycle_react_scale *rrs)
{
  struct gkyl_gyrokinetic_recycling_reaction_scaling_inp *rrs_inp = &ns->info.recycling_reaction_scaling;
  rrs->num_boundaries = rrs_inp->num_boundaries;
  rrs->recycling_coeff = rrs_inp->recycling_coeff;

  if (rrs->num_boundaries > 0) {
    assert(ns->is_fluid);
    rrs->write_diagnostics = rrs_inp->write_diagnostics;

    // Initial number density times conf-space Jacobian.
    rrs->Jm0_init = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    // Create an updater that integrates an array.
    rrs->integrate_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    int num_mom = 1;
    if (app->use_gpu){
      rrs->react_vol_integ_local  = gkyl_cu_malloc(sizeof(double[num_mom]));
      rrs->react_vol_integ_global = gkyl_cu_malloc(sizeof(double[num_mom]));
      rrs->bflux_m0_vol_integ_local = gkyl_cu_malloc(sizeof(double));
    }
    else {
      rrs->react_vol_integ_local  = gkyl_malloc(sizeof(double[num_mom]));
      rrs->react_vol_integ_global = gkyl_malloc(sizeof(double[num_mom]));
      rrs->bflux_m0_vol_integ_local = gkyl_malloc(sizeof(double));
    }

    rrs->cross_moms_func = gk_neut_species_rrs_cross_moms_enabled;
    rrs->rhs_func = gk_neut_species_rrs_rhs_enabled;
    rrs->apply_func = gk_neut_species_rrs_apply_enabled;
    if (rrs->write_diagnostics)
      rrs->write_func = gk_neut_species_rrs_write_enabled;
    else
      rrs->write_func = gk_neut_species_rrs_write_disabled;
  }
  else {
    rrs->cross_moms_func = gk_neut_species_rrs_cross_moms_disabled;
    rrs->rhs_func = gk_neut_species_rrs_rhs_disabled;
    rrs->apply_func = gk_neut_species_rrs_apply_disabled;
    rrs->write_func = gk_neut_species_rrs_write_disabled;
  }
}

void 
gk_neut_species_recycle_react_scale_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_recycle_react_scale *rrs)
{
  if (rrs->num_boundaries > 0) {
    struct gkyl_gyrokinetic_recycling_reaction_scaling_inp *rrs_inp = &ns->info.recycling_reaction_scaling;

    // Fetch index of species for indexing arrays.
    rrs->elc_idx = gk_find_species_idx(app, rrs_inp->electron_name);
    rrs->ion_idx = gk_find_species_idx(app, rrs_inp->impacting_ion_name);

    struct gk_species *gks_ion = &app->species[rrs->ion_idx];

    for (int j=0; j<rrs->num_boundaries; ++j) {
      int dir  = rrs_inp->boundaries_dir[j];
      int edge = rrs_inp->boundaries_edge[j];
  
      // Source adaptation on periodic, zero flux, or reflect boundary is not allowed.
      assert(gks_ion->bc_is_np[dir]);
      if (edge == GKYL_LOWER_EDGE) {
        assert(gks_ion->lower_bc[dir].type != GKYL_SPECIES_ZERO_FLUX);
        assert(gks_ion->lower_bc[dir].type != GKYL_SPECIES_REFLECT);
      } else {
        assert(gks_ion->upper_bc[dir].type != GKYL_SPECIES_ZERO_FLUX);
        assert(gks_ion->upper_bc[dir].type != GKYL_SPECIES_REFLECT);
      }
  
      // Default scenario: we set the ranges to the full range of the ghost cells.
      rrs->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE ? app->lower_ghost[dir] : app->upper_ghost[dir];
      rrs->boundaries_dir[j]  = dir;
      rrs->boundaries_edge[j] = edge;
  
      // Specific scenario if we are in a inner wall limited case. We select only SOL range in parallel direction.
      if (edge == GKYL_LOWER_EDGE? gks_ion->lower_bc[dir].type == GKYL_SPECIES_GK_IWL
                                 : gks_ion->upper_bc[dir].type == GKYL_SPECIES_GK_IWL)
      {
        rrs->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE ? app->lower_ghost_par_sol : app->upper_ghost_par_sol;
      }
    }

    struct gkyl_dg_iz_inp iz_inp = {
      .cbasis = &app->basis, 
      .conf_rng = &app->local, 
      .type_ion = rrs_inp->impacting_ion_id, 
      .charge_state = 0,
      .type_self = GKYL_SELF_ION, // Could be GKYL_SELF_DONOR. It just can't be
                                  // GKYL_SELF_ELC because we don't need to
                                  // compute the ionization temperatures.
    };
    rrs->iz_react_calc = gkyl_dg_iz_new(&iz_inp, app->use_gpu);

    // Reaction contribution.
    rrs->reactivity = mkarr(app->use_gpu, ns->basis.num_basis, ns->local_ext.volume);
    rrs->dndt_react = mkarr(app->use_gpu, ns->basis.num_basis, ns->local_ext.volume);
  }
}

void 
gk_neut_species_recycle_react_scale_apply_ic_cross(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_recycle_react_scale *rrs)
{
  if (rrs->num_boundaries > 0) {
    // Store the initial particle number density.
    gkyl_array_set_offset(rrs->Jm0_init, 1.0/ns->info.mass, ns->f, 0);
  }
}

void
gk_neut_species_recycle_react_scale_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *ns,
  struct gk_recycle_react_scale *rrs, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  rrs->cross_moms_func(app, ns, rrs, fin, fin_neut);
}

void
gk_neut_species_recycle_react_scale_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_recycle_react_scale *rrs, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  rrs->rhs_func(app, ns, rrs, fin, rhs);
}

void
gk_neut_species_recycle_react_scale_apply(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_recycle_react_scale *rrs, struct gkyl_array *fin)
{
  rrs->apply_func(app, ns, rrs, fin);
}

void
gk_neut_species_recycle_react_scale_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_recycle_react_scale *rrs, int ridx, double tm, int frame)
{
  rrs->write_func(app, gkns, rrs, ridx, tm, frame);
}

void 
gk_neut_species_recycle_react_scale_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_recycle_react_scale *rrs)
{
  if (rrs->num_boundaries > 0) {
    gkyl_array_release(rrs->Jm0_init);
    gkyl_array_integrate_release(rrs->integrate_op);

    if (app->use_gpu){
      gkyl_cu_free(rrs->react_vol_integ_local );
      gkyl_cu_free(rrs->react_vol_integ_global);
      gkyl_cu_free(rrs->bflux_m0_vol_integ_local);
    }
    else {
      gkyl_free(rrs->react_vol_integ_local);
      gkyl_free(rrs->react_vol_integ_global);
      gkyl_free(rrs->bflux_m0_vol_integ_local);
    }
    gkyl_array_release(rrs->reactivity);
    gkyl_array_release(rrs->dndt_react);

    gkyl_dg_iz_release(rrs->iz_react_calc);
  } 
}
