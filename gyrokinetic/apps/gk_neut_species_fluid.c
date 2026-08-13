#include <gkyl_gk_neut_species_priv.h>

static void
gk_neut_species_fluid_diffusion_dirichlet_ghost(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  int dir, enum gkyl_edge_loc edge, double number_density,
  const struct gkyl_array *number_density_field,
  struct gkyl_array *density, struct gkyl_array *tensor)
{
  // Communication has already populated internal-rank ghosts. Only a rank
  // owning the physical edge supplies the exterior Dirichlet state.
  if ((edge == GKYL_LOWER_EDGE
      && species->local.lower[dir] != species->global.lower[dir])
    || (edge == GKYL_UPPER_EDGE
      && species->local.upper[dir] != species->global.upper[dir]))
    return;

  const struct gkyl_range *ghost = edge == GKYL_LOWER_EDGE
    ? &species->local_lower_ghost[dir] : &species->local_upper_ghost[dir];
  const int skin_idx_dir = edge == GKYL_LOWER_EDGE
    ? species->local.lower[dir] : species->local.upper[dir];
  const int nb = app->basis.num_basis;
  const int nK = app->cdim*app->cdim;
  const double rho_b_modal = pow(sqrt(2.0), app->cdim)
    *species->info.mass*number_density;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost);
  while (gkyl_range_iter_next(&iter)) {
    int skin_idx[GKYL_MAX_CDIM];
    for (int d=0; d<app->cdim; ++d)
      skin_idx[d] = iter.idx[d];
    skin_idx[dir] = skin_idx_dir;

    const double *rho_skin = gkyl_array_cfetch(density,
      gkyl_range_idx(&species->local_ext, skin_idx));
    double *rho_ghost = gkyl_array_fetch(density,
      gkyl_range_idx(&species->local_ext, iter.idx));

    // Mirror the skin polynomial and choose its exterior value so the
    // arithmetic trace at the physical face is rho_b:
    //   rho_g = 2 rho_b - R(rho_s).
    app->basis.flip_odd_sign(dir, rho_skin, rho_ghost);
    for (int k=0; k<nb; ++k)
      rho_ghost[k] = -rho_ghost[k];
    if (number_density_field) {
      const double *n_b = gkyl_array_cfetch(number_density_field,
        gkyl_range_idx(&species->local_ext, iter.idx));
      for (int k=0; k<nb; ++k)
        rho_ghost[k] += 2.0*species->info.mass*n_b[k];
    }
    else {
      rho_ghost[0] += 2.0*rho_b_modal;
    }

    const double *K_skin = gkyl_array_cfetch(tensor,
      gkyl_range_idx(&species->local_ext, skin_idx));
    double *K_ghost = gkyl_array_fetch(tensor,
      gkyl_range_idx(&species->local_ext, iter.idx));
    for (int k=0; k<nK; ++k)
      app->basis.flip_odd_sign(dir, &K_skin[k*nb], &K_ghost[k*nb]);
  }
}

static void
gk_neut_species_fluid_recycling_density(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  int dir, enum gkyl_edge_loc edge)
{
  const int e = edge == GKYL_LOWER_EDGE ? 0 : 1;
  const struct gkyl_gyrokinetic_bc *bc = edge == GKYL_LOWER_EDGE
    ? &species->lower_bc[dir] : &species->upper_bc[dir];
  const struct gkyl_range *ghost = edge == GKYL_LOWER_EDGE
    ? &species->local_lower_ghost[dir] : &species->local_upper_ghost[dir];
  gkyl_array_clear_range(species->diffusion_bc_density, 0.0, ghost);

  // A unit-density isotropic Maxwellian emits the one-sided particle flux
  // Gamma_M,1=sqrt(T_n/(2*pi*m_n)). Match the kinetic recycling convention:
  // n_b*Gamma_M,1=R*Gamma_i, locally along the target.
  const double unit_flux = sqrt(bc->emission.emission_temp
    /(2.0*M_PI*species->info.mass));
  for (int i=0; i<species->diffusion_recycling_num_species[e]; ++i) {
    struct gk_species *ion =
      &app->species[species->diffusion_recycling_ion_idx[e][i]];
    assert(ion->bflux.current_moms);
    gkyl_array_clear_range(species->diffusion_bc_density_tmp, 0.0, ghost);
    gk_species_bflux_get_flux_mom(&ion->bflux, dir, edge,
      GKYL_F_MOMENT_M0, ion->bflux.current_moms,
      species->diffusion_bc_density_tmp, ghost);
    gkyl_array_accumulate_range(species->diffusion_bc_density,
      bc->emission.recycling_frac/unit_flux,
      species->diffusion_bc_density_tmp, ghost);
  }
}

static void
gk_neut_species_fluid_apply_bc_dynamic(gkyl_gyrokinetic_app *app,
  const struct gk_neut_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    app->num_periodic_dir, app->periodic_dirs, f);
  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static double
gk_neut_species_fluid_rhs_dynamic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  double omega_cfl = 1/DBL_MAX;   
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);
  
  // Collisionless terms.
  struct timespec wst = gkyl_wall_clock();
  if (species->has_diffusion) {
    if (species->use_reaction_rate_diffusion) {
      struct gk_react *react = &species->react_neut;
      int iz_idx = species->ionization_react_idx;
      int cx_idx = species->diffusion_cx_react_idx;
      struct gk_species *elc = &app->species[react->elc_idx[iz_idx]];

      // denominator = ne*<sigma v>_cx.
      gkyl_dg_mul_op_range(&app->basis, 0, species->diffusion_moment_ratio,
        0, elc->lte.moms.marr, 0, react->coeff_react[cx_idx], &species->local);

      // Isotropic scalar DG coefficient: D = vti^2/(ne*<sigma v>_cx).
      gkyl_dg_div_op_range(species->diffusion_div_mem, &app->basis,
        0, species->diffusion_coeff, 0, react->vt_sq_ion[cx_idx],
        0, species->diffusion_moment_ratio, &species->local);

      // The variable-coefficient surface terms read neighboring cells.
      // Populate coefficient ghosts independently of the state BC update.
      gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
        app->num_periodic_dir, app->periodic_dirs, species->diffusion_coeff);
      gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext,
        species->diffusion_coeff);
    }

    // The evolved density is N=J*rho. Form the physical density rho=N/J and
    // the mapped tensor K^{ij}=J*D*g_neut^{ij} used by
    //   d_t N = d_i(K^{ij} d_j rho).
    // Use weak division by J instead of multiplication by its separately
    // projected inverse. This preserves a spatially constant physical density
    // exactly when the stored state is N=J*rho.
    gkyl_dg_div_op_range(species->diffusion_div_mem, &app->basis,
      0, species->diffusion_density, 0, fin, 0,
      app->gk_geom->geo_int.jacobgeo, &species->local);
    gkyl_dg_mul_op_range(&app->basis, 0, species->diffusion_geom_factor,
      0, species->diffusion_coeff, 0,
      app->gk_geom->geo_int.jacobgeo, &species->local);

    const int metric_idx_1x[1][1] = { { 5 } };
    const int metric_idx_2x[2][2] = { { 0, 2 }, { 2, 5 } };
    const int metric_idx_3x[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
    for (int i=0; i<app->cdim; ++i) {
      for (int j=0; j<app->cdim; ++j) {
        int gij_idx = app->cdim == 1 ? metric_idx_1x[i][j]
          : app->cdim == 2 ? metric_idx_2x[i][j] : metric_idx_3x[i][j];
        gkyl_dg_mul_op_range(&app->basis, i*app->cdim+j,
          species->diffusion_tensor, 0, species->diffusion_geom_factor,
          gij_idx, app->gk_geom->geo_int.gij_neut, &species->local);
      }
    }
    gkyl_comm_array_per_sync(species->comm, &species->local,
      &species->local_ext, app->num_periodic_dir, app->periodic_dirs,
      species->diffusion_density);
    gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext,
      species->diffusion_density);
    gkyl_comm_array_per_sync(species->comm, &species->local,
      &species->local_ext, app->num_periodic_dir, app->periodic_dirs,
      species->diffusion_tensor);
    gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext,
      species->diffusion_tensor);

    // Fixed-density or locally recycled open-field-line boundaries.
    const int par_dir = app->cdim-1;
    if (species->lower_bc[par_dir].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_fluid_recycling_density(app, species, par_dir,
        GKYL_LOWER_EDGE);
      gk_neut_species_fluid_diffusion_dirichlet_ghost(app, species, par_dir,
        GKYL_LOWER_EDGE, 0.0, species->diffusion_bc_density,
        species->diffusion_density, species->diffusion_tensor);
    }
    else if (species->info.diffusion.lower_bc_type
        == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET)
      gk_neut_species_fluid_diffusion_dirichlet_ghost(app, species, par_dir,
        GKYL_LOWER_EDGE, species->info.diffusion.lower_bc_density,
        0,
        species->diffusion_density, species->diffusion_tensor);
    if (species->upper_bc[par_dir].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_fluid_recycling_density(app, species, par_dir,
        GKYL_UPPER_EDGE);
      gk_neut_species_fluid_diffusion_dirichlet_ghost(app, species, par_dir,
        GKYL_UPPER_EDGE, 0.0, species->diffusion_bc_density,
        species->diffusion_density, species->diffusion_tensor);
    }
    else if (species->info.diffusion.upper_bc_type
        == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET)
      gk_neut_species_fluid_diffusion_dirichlet_ghost(app, species, par_dir,
        GKYL_UPPER_EDGE, species->info.diffusion.upper_bc_density,
        0,
        species->diffusion_density, species->diffusion_tensor);

    gkyl_dg_updater_conf_diffusion_advance(species->diffusion_slvr,
      &species->local, species->diffusion_density, species->cflrate, rhs);

    // Density is the only independently diffused quantity. Propagate its
    // change to the dependent conserved moments while holding their specific
    // values (moment/rho) fixed.
    for (int m=1; m<species->num_moments; ++m) {
      gkyl_dg_div_op_range(species->diffusion_div_mem, &app->basis,
        m, species->diffusion_moment_ratio, m, fin, 0, fin, &species->local);
      gkyl_dg_mul_op_range(&app->basis, m, rhs,
        m, species->diffusion_moment_ratio, 0, rhs, &species->local);
    }
  }
  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);

  // Ionization loss and other enabled reactions with charged species.
  gk_neut_species_react_rhs(app, species, &species->react_neut, fin, rhs);

  // Add the explicit ionization frequency ne*<sigma v>_iz to the CFL rate.
  if (species->ionization_react_idx >= 0) {
    struct gk_react *react = &species->react_neut;
    for (int i=0; i<react->num_react; ++i) {
      if (react->react_id[i] != GKYL_REACT_IZ)
        continue;
      struct gk_species *elc = &app->species[react->elc_idx[i]];
      gkyl_dg_mul_op_range(&app->basis, 0, species->ionization_rate,
        0, elc->lte.moms.marr, 0, react->coeff_react[i], &species->local);
      gkyl_dg_calc_average_range(&app->basis, 0, species->ionization_cflrate,
        0, species->ionization_rate, species->local);
      gkyl_array_accumulate_range(species->cflrate, 1.0,
        species->ionization_cflrate, &species->local);
    }
  }

  // Compute volume-integrated reactions in sca.
  gk_neut_species_scaling_rhs(app, species, &species->sca, fin, rhs);

  app->stat.n_neut_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);
  
  double omega_cfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omega_cfl_ho[0] = species->omega_cfl[0];

  omega_cfl = omega_cfl_ho[0];
  
  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl/omega_cfl;
}

static double
gk_neut_species_fluid_rhs_implicit_dynamic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt)
{ 
  double omega_cfl = 1/DBL_MAX;
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // No implicit terms yet.

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

static void
gk_neut_species_fluid_release_dynamic(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *ns)
{
  // Release memory allocated for dynamic neutrals.
  gkyl_array_release(ns->cflrate);

  if (ns->has_diffusion) {
    gkyl_array_release(ns->diffusion_coeff);
    gkyl_array_release(ns->diffusion_tensor);
    gkyl_array_release(ns->diffusion_geom_factor);
    gkyl_array_release(ns->diffusion_density);
    if (ns->diffusion_bc_density)
      gkyl_array_release(ns->diffusion_bc_density);
    if (ns->diffusion_bc_density_tmp)
      gkyl_array_release(ns->diffusion_bc_density_tmp);
    gkyl_array_release(ns->diffusion_moment_ratio);
    gkyl_dg_bin_op_mem_release(ns->diffusion_div_mem);
    gkyl_dg_updater_conf_diffusion_release(ns->diffusion_slvr);
  }
  if (ns->ionization_react_idx >= 0) {
    gkyl_array_release(ns->ionization_rate);
    gkyl_array_release(ns->ionization_cflrate);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(ns->omega_cfl);
  }
  else {
    gkyl_free(ns->omega_cfl);
  }

  // Release integrated mom data.
  gk_neut_species_moment_release(app, &ns->integ_moms); 

  // Release integrated mom diag data.
  gkyl_dynvec_release(ns->integ_diag);
  
  if (app->use_gpu) {
    gkyl_cu_free(ns->red_integ_diag);
    gkyl_cu_free(ns->red_integ_diag_global);
  }
  else {
    gkyl_free(ns->red_integ_diag);
    gkyl_free(ns->red_integ_diag_global);
  }
}

static void
gk_neut_species_fluid_release(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *ns)
{
  // Release resources for fluid neutral species.
  gkyl_msgpack_map_elem_release(ns->io_meta_basic_len, ns->io_meta_basic);
  gkyl_msgpack_map_elem_release(ns->io_meta_conf_len , ns->io_meta_conf );

  gkyl_array_release(ns->f);
  gkyl_array_release(ns->f1);
  gkyl_array_release(ns->fnew);
  gkyl_array_release(ns->f_host);

  if (ns->info.init_from_file.type == 0) {
    gk_neut_species_projection_release(app, &ns->proj_init);
  }
  gkyl_comm_release(ns->comm);

  for (int i=0; i<ns->info.num_diag_moments; ++i)
    gk_neut_species_moment_release(app, &ns->moms[i]);
  gkyl_free(ns->moms);

  gk_neut_species_bgk_release(app, &ns->bgk);

  gk_neut_species_positivity_release(app, &ns->positivity);

  gk_neut_species_react_release(app, &ns->react_neut);

  // Free boundary flux memory.
  gk_neut_species_bflux_release(app, ns, &ns->bflux);

  gk_neut_species_lte_release(app, &ns->lte);

  // Free memory for the object that scales the species according to a balance
  // between recycling and reactions.
  gk_neut_species_scaling_release(app, &ns->sca);

  ns->release_is_static_func(app, ns);
}

static void
gk_neut_species_fluid_init_dynamic(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns)
{
  int cdim = app->cdim;
  
  // Allocate additional moment arrays for time stepping.
  ns->f1 = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);
  ns->fnew = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);
  
  // Allocate cflrate (scalar array).
  ns->cflrate = mkarr(app->use_gpu, 1, ns->local_ext.volume);

  ns->omega_cfl = app->use_gpu? gkyl_cu_malloc(sizeof(double))
                              : gkyl_malloc(sizeof(double));

  // Initialize diffusion of neutral density (component 0) only. The same
  // scalar coefficient is used in every configuration-space direction.
  ns->use_reaction_rate_diffusion = ns->info.diffusion.use_reaction_rates;
  ns->has_diffusion = ns->info.diffusion.D > 0.0 || ns->use_reaction_rate_diffusion;
  ns->ionization_react_idx = -1;
  ns->diffusion_cx_react_idx = -1;
  ns->diffusion_coeff = 0;
  ns->diffusion_tensor = 0;
  ns->diffusion_geom_factor = 0;
  ns->diffusion_density = 0;
  ns->diffusion_bc_density = 0;
  ns->diffusion_bc_density_tmp = 0;
  ns->diffusion_moment_ratio = 0;
  ns->diffusion_div_mem = 0;
  ns->diffusion_slvr = 0;
  ns->diffusion_recycling_num_species[0] = 0;
  ns->diffusion_recycling_num_species[1] = 0;
  ns->ionization_rate = 0;
  ns->ionization_cflrate = 0;

  for (int i=0; i<ns->info.react_neut.num_react; ++i) {
    enum gkyl_react_id react_id = ns->info.react_neut.react_type[i].react_id;
    if (react_id == GKYL_REACT_IZ && ns->ionization_react_idx < 0)
      ns->ionization_react_idx = i;
    else if (react_id == GKYL_REACT_CX && ns->diffusion_cx_react_idx < 0)
      ns->diffusion_cx_react_idx = i;
  }
  if (ns->ionization_react_idx >= 0) {
    ns->ionization_rate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    ns->ionization_cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  }
  if (ns->has_diffusion) {
    // Configuration-space restriction of the single-integration-by-parts plus
    // recovery-DG anomalous-diffusion scheme. Generated CPU kernels currently
    // cover 1x-3x serendipity p1.
    assert(cdim > 0 && cdim <= 3);
    assert(app->basis.poly_order == 1);
    assert(!app->use_gpu);

    ns->diffusion_coeff = mkarr(false, app->basis.num_basis,
      app->local_ext.volume);
    ns->diffusion_tensor = mkarr(false,
      cdim*cdim*app->basis.num_basis, app->local_ext.volume);
    ns->diffusion_geom_factor = mkarr(false, app->basis.num_basis,
      app->local_ext.volume);
    ns->diffusion_density = mkarr(false, app->basis.num_basis,
      app->local_ext.volume);

    const int par_dir = cdim-1;
    const bool recycle_lo = ns->lower_bc[par_dir].type
      == GKYL_BC_GK_SPECIES_RECYCLE;
    const bool recycle_up = ns->upper_bc[par_dir].type
      == GKYL_BC_GK_SPECIES_RECYCLE;
    if (recycle_lo || recycle_up) {
      ns->diffusion_bc_density = mkarr(false, app->basis.num_basis,
        app->local_ext.volume);
      ns->diffusion_bc_density_tmp = mkarr(false, app->basis.num_basis,
        app->local_ext.volume);
      gkyl_array_clear(ns->diffusion_bc_density, 0.0);
      gkyl_array_clear(ns->diffusion_bc_density_tmp, 0.0);
    }
    for (int e=0; e<2; ++e) {
      const struct gkyl_gyrokinetic_bc *bc = e == 0
        ? &ns->lower_bc[par_dir] : &ns->upper_bc[par_dir];
      if (bc->type != GKYL_BC_GK_SPECIES_RECYCLE)
        continue;
      assert(bc->emission.emission_temp > 0.0);
      assert(bc->emission.num_species > 0);
      ns->diffusion_recycling_num_species[e] = bc->emission.num_species;
      for (int i=0; i<bc->emission.num_species; ++i) {
        ns->diffusion_recycling_ion_idx[e][i] = gk_find_species_idx(app,
          bc->emission.in_species[i]);
        assert(ns->diffusion_recycling_ion_idx[e][i] >= 0);
      }
    }
    gkyl_array_clear(ns->diffusion_coeff, 0.0);
    gkyl_array_clear(ns->diffusion_tensor, 0.0);
    gkyl_array_clear(ns->diffusion_geom_factor, 0.0);
    gkyl_array_clear(ns->diffusion_density, 0.0);
    if (!ns->use_reaction_rate_diffusion) {
      // Modal constant coefficient in cdim dimensions.
      gkyl_array_shiftc(ns->diffusion_coeff,
        pow(sqrt(2.0),cdim)*ns->info.diffusion.D, 0);
    }

    if (ns->use_reaction_rate_diffusion) {
      assert(ns->ionization_react_idx >= 0);
      assert(ns->diffusion_cx_react_idx >= 0);
    }

    ns->diffusion_moment_ratio = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);
    ns->diffusion_div_mem = app->use_gpu
      ? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis)
      : gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);

    int diffusion_zero_flux_flags[2*GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<cdim; ++d)
      diffusion_zero_flux_flags[d] = diffusion_zero_flux_flags[d+cdim]
        = ns->bc_is_np[d];
    if (recycle_lo || ns->info.diffusion.lower_bc_type
        == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET)
      diffusion_zero_flux_flags[par_dir] = 0;
    if (recycle_up || ns->info.diffusion.upper_bc_type
        == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET)
      diffusion_zero_flux_flags[par_dir+cdim] = 0;
    ns->diffusion_slvr = gkyl_dg_updater_conf_diffusion_new(&app->grid,
      &app->basis, &app->local_ext, ns->diffusion_tensor,
      app->gk_geom->geo_int.jacobgeo_inv, diffusion_zero_flux_flags,
      app->use_gpu);
  }

  // Allocate data for integrated moments.
  gk_neut_species_moment_init(app, ns, &ns->integ_moms, GKYL_F_MOMENT_M0M1M2, true);

  // Allocate data for integrated diagnostics.
  if (app->use_gpu) {
    ns->red_integ_diag = gkyl_cu_malloc(sizeof(double[ns->integ_moms.num_mom]));
    ns->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[ns->integ_moms.num_mom]));
  } else {
    ns->red_integ_diag = gkyl_malloc(sizeof(double[ns->integ_moms.num_mom]));
    ns->red_integ_diag_global = gkyl_malloc(sizeof(double[ns->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  ns->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, ns->integ_moms.num_mom);
  ns->is_first_integ_write_call = true;

  // Set function pointers
  ns->rhs_func = gk_neut_species_fluid_rhs_dynamic;
  ns->rhs_implicit_func = gk_neut_species_fluid_rhs_implicit_dynamic;
  ns->bc_func = gk_neut_species_fluid_apply_bc_dynamic;
  ns->release_func = gk_neut_species_fluid_release;
  ns->release_is_static_func = gk_neut_species_fluid_release_dynamic;
  ns->step_f_func = gk_neut_species_step_f_dynamic;
  ns->combine_func = gk_neut_species_combine_dynamic;
  ns->copy_func = gk_neut_species_copy_range_dynamic;
  ns->write_func = gk_neut_species_write_dynamic;
  ns->write_mom_func = gk_neut_species_write_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->write_integrated_mom_func = gk_neut_species_write_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

static void
gk_neut_species_fluid_init_static(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  // Set pointers for RK methods.
  s->f1 = gkyl_array_acquire(s->f);
  s->fnew = gkyl_array_acquire(s->f);

  // Set function pointers
  s->rhs_func = gk_neut_species_rhs_static;
  s->rhs_implicit_func = gk_neut_species_rhs_implicit_static;
  s->bc_func = gk_neut_species_apply_bc_static;
  s->release_func = gk_neut_species_fluid_release;
  s->release_is_static_func = gk_neut_species_release_static;
  s->step_f_func = gk_neut_species_step_f_static;
  s->combine_func = gk_neut_species_combine_static;
  s->copy_func = gk_neut_species_copy_range_static;
  s->write_func = gk_neut_species_write_init_only;
  s->write_mom_func = gk_neut_species_write_mom_init_only;
  s->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_static;
  s->write_integrated_mom_func = gk_neut_species_write_integrated_mom_static;
  s->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

void
gk_neut_species_fluid_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns)
{
  ns->is_fluid = true; // Fluid neutrals.
  assert(ns->info.vdim == 0); // Ensure user provided vdim=0 in input file, or didn't provide it at all.

  ns->model_id = GKYL_MODEL_DEFAULT;
  ns->field_id = GKYL_FIELD_NULL;

  int cdim = app->cdim;

  // Number of moments depends on eqn_type.
  ns->num_moments = 5; // rho, rho*ux, rho*uy, rho*uy, totE

  // Use the same basis as conf-space.
  ns->basis = app->basis;
  ns->basis_on_dev = app->basis_on_dev;
  
  // Use the same grid as conf-space.
  ns->grid = app->grid;
  ns->global_ext = app->global_ext;
  ns->global = app->global;
  
  // Use the same communicator as conf-space.
  ns->comm = gkyl_comm_acquire(app->comm);

  // Use the same local range as conf-space.
  ns->local_ext = app->local_ext;
  ns->local = app->local;

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // add the parallel direction in case TS BCs are needed.
  ns->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<ns->num_periodic_dir; ++d)
    ns->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) ns->bc_is_np[d] = true;
  for (int d=0; d<ns->num_periodic_dir; ++d)
    ns->bc_is_np[ns->periodic_dirs[d]] = false;

  // Store the BCs from the input file.
  for (int d=0; d<app->cdim; ++d) {
    struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(ns->info.bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
    if (bc_lo != 0)
      ns->lower_bc[d] = *bc_lo;
    else
      ns->lower_bc[d].type = GKYL_BC_GK_SKIP;

    struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(ns->info.bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
    if (bc_up != 0)
      ns->upper_bc[d] = *bc_up;
    else
      ns->upper_bc[d].type = GKYL_BC_GK_SKIP;
  }

  // Species properties metadata.
  struct gkyl_msgpack_map_elem io_meta_sprop[] = {
    { .key = "mass", .elem_type = GKYL_MP_DOUBLE, .dval = ns->info.mass },
    { .key = "charge", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0 },
    { .key = "gas_gamma", .elem_type = GKYL_MP_DOUBLE, .dval = ns->info.gas_gamma },
    { .key = "vdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = ns->info.vdim },
  };

  // Metadata for integrated quantities.
  const struct gkyl_msgpack_map_elem *io_meta_basic_union[] = {app->io_meta_basic, io_meta_sprop};
  int io_meta_basic_union_len[] = {app->io_meta_basic_len, sizeof(io_meta_sprop)/sizeof(io_meta_sprop[0])};
  ns->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_basic_union)/sizeof(io_meta_basic_union[0]),
    io_meta_basic_union_len, io_meta_basic_union, &ns->io_meta_basic_len);

  // Metadata for conf-space quantities.
  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = app->basis.id },
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0 },
    { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
  };
  const struct gkyl_msgpack_map_elem *io_meta_conf_union[] = {app->io_meta_basic, io_meta_sprop, io_meta_conf};
  int io_meta_conf_union_len[] = {app->io_meta_basic_len, sizeof(io_meta_sprop)/sizeof(io_meta_sprop[0]),
    sizeof(io_meta_conf)/sizeof(io_meta_conf[0])};
  ns->io_meta_conf = gkyl_msgpack_map_elem_union(sizeof(io_meta_conf_union)/sizeof(io_meta_conf_union[0]),
    io_meta_conf_union_len, io_meta_conf_union, &ns->io_meta_conf_len);

  // Metadata for phase-space quantities.
  ns->io_meta_phase = ns->io_meta_conf;
  ns->io_meta_phase_len = ns->io_meta_conf_len;

  // Allocate distribution function array for initialization and I/O.
  ns->f = mkarr(app->use_gpu, ns->num_moments*ns->basis.num_basis, ns->local_ext.volume);
  ns->f_host = app->use_gpu? mkarr(false, ns->f->ncomp, ns->f->size)
                           : gkyl_array_acquire(ns->f);

  // Create skin/ghost ranges fir applying BCs. Only used for dynamic neutrals but included here to avoid
  // code duplication since the "ghost" array is needed.
  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d)
    ghost[d] = 1;

  // Create skin/ghost ranges.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&ns->local_lower_skin[dir], &ns->local_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &ns->local_ext, ghost);
    gkyl_skin_ghost_ranges(&ns->local_upper_skin[dir], &ns->local_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &ns->local_ext, ghost);
    gkyl_skin_ghost_ranges(&ns->global_lower_skin[dir], &ns->global_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &ns->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&ns->global_upper_skin[dir], &ns->global_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &ns->local_ext, ghost);
  }

  // Initialize projection routine for initial conditions.
  if (ns->info.init_from_file.type == 0) {
    gk_neut_species_projection_init(app, ns, ns->info.projection, &ns->proj_init);
  }

  // Allocate objects for computing diagnostic moments.
  int ndm = ns->info.num_diag_moments;
  ns->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_neut_species_moment_init(app, ns, &ns->moms[m], ns->info.diag_moments[m], false);

  // Initialize boundary fluxes.
  ns->bflux = (struct gk_boundary_fluxes) { };
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp) { };
  // Set the operation type for the bflux app.
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  gk_neut_species_bflux_init(app, ns, &ns->bflux, bflux_type, add_bflux_moms_inp);

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  ns->lte = (struct gk_lte) { };
  struct correct_all_moms_inp corr_inp = { .correct_all_moms = false,
    .max_iter = 0, .iter_eps = 10, .use_last_converged = false };
  gk_neut_species_lte_init(app, ns, &ns->lte, corr_inp);

  // Initialize the object that scales the species according to a balance
  // between recycling and reactions.
  ns->sca = (struct gk_scaling) { };
  gk_neut_species_scaling_init(app, ns, &ns->sca);

  // Initialize BGK collisions with null type (not applicable to fluids).
  ns->bgk = (struct gk_bgk_collisions) { };
  ns->info.collisions.collision_id = 0;
  gk_neut_species_bgk_init(app, ns, &ns->bgk);

  // Initialize positivity enforcing operator with null type (NYI for fluids).
  ns->positivity = (struct gk_positivity) { };
  gk_neut_species_positivity_init(app, ns, &ns->positivity);

  // Initialize reactions with charged species. Fluid neutrals use the
  // ionization reaction as a density loss; CX supplies the diffusion rate.
  ns->react_neut = (struct gk_react) { };
  gk_neut_species_react_init(app, ns, ns->info.react_neut, &ns->react_neut);

  ns->src = (struct gk_source) { };
  if (!ns->info.is_static) {
    gk_neut_species_fluid_init_dynamic(gk, app, ns);
  }
  else {
    gk_neut_species_fluid_init_static(gk, app, ns);
  }
}
