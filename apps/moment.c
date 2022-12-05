#include <gkyl_moment_priv.h>

static inline int
int_max(int a, int b)
{
  return a>b ? a : b;
}

gkyl_moment_app*
gkyl_moment_app_new(struct gkyl_moment *mom)
{
  disable_denorm_float();
  
  struct gkyl_moment_app *app = gkyl_malloc(sizeof(gkyl_moment_app));

  int ndim = app->ndim = mom->ndim;
  strcpy(app->name, mom->name);
  app->tcurr = 0.0; // reset on init

  app->scheme_type = mom->scheme_type;
  app->mp_recon = mom->mp_recon;
  app->use_hybrid_flux_kep = mom->use_hybrid_flux_kep;
  
  if (app->scheme_type == GKYL_MOMENT_WAVE_PROP)
    app->update_func = moment_update_one_step;
  else if (app->scheme_type == GKYL_MOMENT_MP) 
    app->update_func = moment_update_ssp_rk3;
  else if (app->scheme_type == GKYL_MOMENT_KEP)
    app->update_func = moment_update_ssp_rk3;

  int ghost[3] = { 2, 2, 2 }; // 2 ghost-cells for wave
  if (mom->scheme_type != GKYL_MOMENT_WAVE_PROP)
    for (int d=0; d<3; ++d) ghost[d] = 3; // 3 for MP scheme and KEP
  
  gkyl_rect_grid_init(&app->grid, ndim, mom->lower, mom->upper, mom->cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);
  
  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  app->c2p_ctx = app->mapc2p = 0;  
  app->has_mapc2p = mom->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = mom->c2p_ctx;
    app->mapc2p = mom->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, ndim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(false, ndim*basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p = gkyl_eval_on_nodes_new(&app->grid, &basis, ndim, mom->mapc2p, mom->c2p_ctx);
    gkyl_eval_on_nodes_advance(ev_c2p, 0.0, &app->local_ext, c2p);

    // write DG projection of mapc2p to file
    cstr fileNm = cstr_from_fmt("%s-mapc2p.gkyl", app->name);
    gkyl_grid_sub_array_write(&app->grid, &app->local, c2p, fileNm.str);
    cstr_drop(&fileNm);

    gkyl_array_release(c2p);
    gkyl_eval_on_nodes_release(ev_c2p);
  }

  // create geometry object
  app->geom = gkyl_wave_geom_new(&app->grid, &app->local_ext,
    app->mapc2p, app->c2p_ctx);

  double cfl_frac = mom->cfl_frac == 0 ? 0.95 : mom->cfl_frac;
  app->cfl = 1.0*cfl_frac;
  if (app->scheme_type == GKYL_MOMENT_MP)
    app->cfl = 0.4*cfl_frac; // this should be 1/(1+alpha) = 0.2 but is set to a larger value

  app->num_periodic_dir = mom->num_periodic_dir;
  for (int d=0; d<ndim; ++d)
    app->periodic_dirs[d] = mom->periodic_dirs[d];

  // construct list of directions to skip
  for (int d=0; d<3; ++d)
    app->is_dir_skipped[d] = 0;
  for (int i=0; i<mom->num_skip_dirs; ++i)
    app->is_dir_skipped[mom->skip_dirs[i]] = 1;

  app->has_field = 0;
  // initialize field if we have one
  if (mom->field.init) {
    app->has_field = 1;
    moment_field_init(mom, &mom->field, app, &app->field);
  }

  int ns = app->num_species = mom->num_species;
  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct moment_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i)
    moment_species_init(mom, &mom->species[i], app, &app->species[i]);

  // specify collision parameters FIXME move to a better place
  app->has_collision = mom->has_collision;
  app->gas_gamma = mom->gas_gamma;
  int num_entries = app->num_species * (app->num_species-1) / 2;
  for (int s=0; s<app->num_species; ++s)
    for (int r=0; r<app->num_species; ++r)
      app->nu_base[s][r] = mom->nu_base[s][r];

  app->has_user_sources = mom->has_user_sources;

  // check if we should update sources
  app->update_sources = 0;
  if (app->has_field && ns>0) {
    app->update_sources = 1; // only update if field and species are present
    moment_coupling_init(app, &app->sources);
  }

  app->update_mhd_source = false;
  if (ns==1 && mom->species[0].equation->type==GKYL_EQN_MHD) {
    app->update_mhd_source = true;
    mhd_src_init(app, &mom->species[0], &app->mhd_source);
  }

  // allocate work array for use in MP scheme
  if (app->scheme_type == GKYL_MOMENT_MP || app->scheme_type == GKYL_MOMENT_KEP) {
    int max_eqn = 0;
    for (int i=0; i<ns; ++i)
      max_eqn = int_max(max_eqn, app->species[i].num_equations);
    if (app->has_field)
      max_eqn = int_max(max_eqn, 8); // maxwell equations have 8 components
    app->ql = mkarr(false, max_eqn, app->local_ext.volume);
    app->qr = mkarr(false, max_eqn, app->local_ext.volume);

    app->amdq = mkarr(false, max_eqn, app->local_ext.volume);
    app->apdq = mkarr(false, max_eqn, app->local_ext.volume);
  }

  // initialize stat object to all zeros
  app->stat = (struct gkyl_moment_stat) {
  };

  return app;
}

double
gkyl_moment_app_max_dt(gkyl_moment_app* app)
{
  double max_dt = DBL_MAX;
  for (int i=0;  i<app->num_species; ++i) 
    max_dt = fmin(max_dt, moment_species_max_dt(app, &app->species[i]));

  if (app->has_field && app->field.evolve)
    max_dt = fmin(max_dt, moment_field_max_dt(app, &app->field));

  return max_dt;
}

void
gkyl_moment_app_apply_ic(gkyl_moment_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_moment_app_apply_ic_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_moment_app_apply_ic_species(app, i, t0);
}

void
gkyl_moment_app_apply_ic_field(gkyl_moment_app* app, double t0)
{
  if (app->has_field != 1) return;
  
  app->tcurr = t0;
  int num_quad = app->scheme_type == GKYL_MOMENT_MP ? 4 : 2;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, num_quad, 8, app->field.init, app->field.ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->field.fcurr);
  gkyl_fv_proj_release(proj);

  moment_field_apply_bc(app, t0, &app->field, app->field.fcurr);
}

void
gkyl_moment_app_apply_ic_species(gkyl_moment_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  int num_quad = app->scheme_type == GKYL_MOMENT_MP ? 4 : 2;  
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, num_quad, app->species[sidx].num_equations,
    app->species[sidx].init, app->species[sidx].ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->species[sidx].fcurr);
  gkyl_fv_proj_release(proj);

  moment_species_apply_bc(app, t0, &app->species[sidx], app->species[sidx].fcurr);
}

void
gkyl_moment_app_write(const gkyl_moment_app* app, double tm, int frame)
{
  gkyl_moment_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_moment_app_write_species(app, i, tm, frame);
}

void
gkyl_moment_app_write_field(const gkyl_moment_app* app, double tm, int frame)
{
  if (app->has_field != 1) return;

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  gkyl_grid_sub_array_write(&app->grid, &app->local, app->field.fcurr, fileNm.str);
  cstr_drop(&fileNm);

  // write external EM field if it is present
  if (app->field.ext_em) {
    cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "ext_em_field", frame);
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->field.ext_em, fileNm.str);
    cstr_drop(&fileNm);
  }
}

void
gkyl_moment_app_write_field_energy(gkyl_moment_app *app)
{
  if (app->has_field) {
    // write out field energy
    cstr fileNm = cstr_from_fmt("%s-field-energy.gkyl", app->name);

    if (app->field.is_first_energy_write_call) {
      // write to a new file (this ensure previous output is removed)
      gkyl_dynvec_write(app->field.integ_energy, fileNm.str);
      app->field.is_first_energy_write_call = false;
    }
    else {
      // append to existing file
      gkyl_dynvec_awrite(app->field.integ_energy, fileNm.str);
    }
    gkyl_dynvec_clear(app->field.integ_energy);
    
    cstr_drop(&fileNm);
  }
}

void
gkyl_moment_app_write_integrated_mom(gkyl_moment_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    // write out diagnostic moments
    cstr fileNm = cstr_from_fmt("%s-%s-%s.gkyl", app->name, app->species[i].name,
      "imom");
    
    if (app->species[i].is_first_q_write_call) {
      gkyl_dynvec_write(app->species[i].integ_q, fileNm.str);
      app->species[i].is_first_q_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(app->species[i].integ_q, fileNm.str);
    }
    gkyl_dynvec_clear(app->species[i].integ_q);

    cstr_drop(&fileNm);
  }
}

void
gkyl_moment_app_write_species(const gkyl_moment_app* app, int sidx, double tm, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->species[sidx].name, frame);
  gkyl_grid_sub_array_write(&app->grid, &app->local, app->species[sidx].fcurr, fileNm.str);
  cstr_drop(&fileNm);

  if (app->scheme_type == GKYL_MOMENT_KEP) {
    cstr fileNm = cstr_from_fmt("%s-%s-alpha_%d.gkyl", app->name, app->species[sidx].name, frame);
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->species[sidx].alpha, fileNm.str);
    cstr_drop(&fileNm);
  }
}

struct gkyl_update_status
gkyl_moment_update(gkyl_moment_app* app, double dt)
{
  app->stat.nup += 1;
  
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status = app->update_func(app, dt);
  app->tcurr += status.dt_actual;
  
  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  
  return status;
}

void
gkyl_moment_app_calc_field_energy(gkyl_moment_app* app, double tm)
{
  if (app->has_field) {
    double energy[6] = { 0.0 };
    calc_integ_quant(app->field.maxwell, app->grid.cellVolume, app->field.fcurr, app->geom,
      app->local, energy);
    gkyl_dynvec_append(app->field.integ_energy, tm, energy);
  }
}

void
gkyl_moment_app_calc_integrated_mom(gkyl_moment_app *app, double tm)
{
  for (int sidx=0; sidx<app->num_species; ++sidx) {

    int num_diag = app->species[sidx].equation->num_diag;
    double q_integ[num_diag];

    calc_integ_quant(app->species[sidx].equation, app->grid.cellVolume, app->species[sidx].fcurr, app->geom,
      app->local, q_integ);
    
    gkyl_dynvec_append(app->species[sidx].integ_q, tm, q_integ);
  }
}

struct gkyl_moment_stat
gkyl_moment_app_stat(gkyl_moment_app* app)
{
  return app->stat;
}

void
gkyl_moment_app_stat_write(const gkyl_moment_app* app)
{
  cstr fileNm = cstr_from_fmt("%s-%s", app->name, "stat.json");

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  // compute total number of cells updated in simulation
  long tot_cells_up = app->local.volume*app->num_species*app->ndim*app->stat.nup;

  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  with_file (fp, fileNm.str, "a") {
    fprintf(fp, "{\n");

    if (strftime(buff, sizeof buff, "%c", &curr_tm))
      fprintf(fp, " date : %s\n", buff);

    fprintf(fp, " nup : %ld,\n", app->stat.nup);
    fprintf(fp, " nfail : %ld,\n", app->stat.nfail);
    fprintf(fp, " total_tm : %lg,\n", app->stat.total_tm);

    if (app->scheme_type == GKYL_MOMENT_WAVE_PROP) {    
      fprintf(fp, " species_tm : %lg,\n", app->stat.species_tm);
      fprintf(fp, " field_tm : %lg,\n", app->stat.field_tm);
      fprintf(fp, " sources_tm : %lg\n", app->stat.sources_tm);
    }
    else if (app->scheme_type == GKYL_MOMENT_MP || app->scheme_type == GKYL_MOMENT_KEP) {
      fprintf(fp, " nfeuler : %ld,\n", app->stat.nfeuler);
      fprintf(fp, " nstage_2_fail : %ld,\n", app->stat.nstage_2_fail);
      fprintf(fp, " nstage_3_fail : %ld,\n", app->stat.nstage_3_fail);

      fprintf(fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
        app->stat.stage_2_dt_diff[0], app->stat.stage_2_dt_diff[1]);
      fprintf(fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
        app->stat.stage_3_dt_diff[0], app->stat.stage_3_dt_diff[1]);
      
      fprintf(fp, " total_tm : %lg,\n", app->stat.total_tm);
      fprintf(fp, " init_species_tm : %lg,\n", app->stat.init_species_tm);
      if (app->has_field)
        fprintf(fp, " init_field_tm : %lg,\n", app->stat.init_field_tm);
      
      fprintf(fp, " species_rhs_tm : %lg,\n", app->stat.species_rhs_tm);

      if (app->has_field)
        fprintf(fp, " field_rhs_tm : %lg,\n", app->stat.field_rhs_tm);
    }
    
    for (int i=0; i<app->num_species; ++i) {
      long tot_bad_cells = 0L;
      if (app->scheme_type == GKYL_MOMENT_WAVE_PROP) {
        for (int d=0; d<app->ndim; ++d) {
          struct gkyl_wave_prop_stats wvs = gkyl_wave_prop_stats(app->species[i].slvr[d]);
          fprintf(fp, " %s_n_bad_1D_sweeps[%d] = %ld\n", app->species[i].name, d, wvs.n_bad_advance_calls);
          fprintf(fp, " %s_n_bad_cells[%d] = %ld\n", app->species[i].name, d, wvs.n_bad_cells);
          fprintf(fp, " %s_n_max_bad_cells[%d] = %ld\n", app->species[i].name, d, wvs.n_max_bad_cells);

          tot_bad_cells += wvs.n_bad_cells;
        }
      }
      fprintf(fp, " %s_bad_cell_frac = %lg\n", app->species[i].name, (double) tot_bad_cells/tot_cells_up );
    }
  
    fprintf(fp, "}\n");
  }

  cstr_drop(&fileNm);
}

void
gkyl_moment_app_release(gkyl_moment_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    moment_species_release(&app->species[i]);
  gkyl_free(app->species);

  if (app->has_field)
    moment_field_release(&app->field);

  if (app->update_sources)
    moment_coupling_release(&app->sources);

  if (app->update_mhd_source)
    mhd_src_release(&app->mhd_source);

  gkyl_wave_geom_release(app->geom);

  if (app->scheme_type == GKYL_MOMENT_MP || app->scheme_type == GKYL_MOMENT_KEP) {
    gkyl_array_release(app->ql);
    gkyl_array_release(app->qr);
    gkyl_array_release(app->amdq);
    gkyl_array_release(app->apdq);
  }

  gkyl_free(app);
}

