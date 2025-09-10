#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_neut_species_lte_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, struct gk_lte *lte, 
  struct correct_all_moms_inp corr_inp)
{
  int cdim = app->cdim, vdim = 3;

  // allocate moments needed for lte update
  gk_neut_species_moment_init(app, s, &lte->moms, GKYL_F_MOMENT_LTE, false);

  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = &s->grid,
    .vel_grid = &s->grid_vel, 
    .conf_basis = &app->basis,
    .phase_basis = &s->basis,
    .conf_range =  &app->local,
    .conf_range_ext = &app->local_ext,
    .vel_range = &s->local_vel,
    .phase_range = &s->local,
    .h_ij = s->g_ij,
    .h_ij_inv = s->gij,
    .det_h = app->gk_geom->jacobgeo,
    .hamil = s->hamil,
    .model_id = s->model_id,
    .use_gpu = app->use_gpu,
  };
  lte->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

  lte->correct_all_moms = corr_inp.correct_all_moms;
  int max_iter = corr_inp.max_iter > 0 ? corr_inp.max_iter : 50;
  double iter_eps = corr_inp.iter_eps > 0 ? corr_inp.iter_eps  : 1e-10;
  bool use_last_converged = corr_inp.use_last_converged;
  
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_inp inp_corr = {
      .phase_grid = &s->grid,
      .vel_grid = &s->grid_vel, 
      .conf_basis = &app->basis,
      .phase_basis = &s->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel,
      .phase_range = &s->local,
      .h_ij = s->g_ij,
      .h_ij_inv = s->gij,
      .det_h = app->gk_geom->jacobgeo,
      .hamil = s->hamil,
      .model_id = s->model_id,
      .use_gpu = app->use_gpu,
      .max_iter = max_iter,
      .eps = iter_eps,
      .use_last_converged = use_last_converged, 
    };
    lte->n_iter = 0;
    lte->corr_lte = gkyl_vlasov_lte_correct_inew( &inp_corr );

    lte->corr_stat = gkyl_dynvec_new(GKYL_DOUBLE, 7);
    lte->is_first_corr_status_write_call = true;
  }

  lte->f_lte = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
}

// Compute f_lte from input LTE moments
void
gk_neut_species_lte_from_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_lte *lte, const struct gkyl_array *moms_lte)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(lte->f_lte, 0.0);

  // Project the LTE distribution function to obtain f_lte.
  // Projection routine also corrects the density of the projected distribution function.
  gkyl_vlasov_lte_proj_on_basis_advance(lte->proj_lte, &species->local, &app->local, 
    moms_lte, lte->f_lte);

  // Correct all the moments of the projected LTE distribution function.
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_status status_corr;
    status_corr = gkyl_vlasov_lte_correct_all_moments(lte->corr_lte, lte->f_lte, moms_lte,
      &species->local, &app->local);
    double corr_vec[7] = { 0.0 };
    corr_vec[0] = status_corr.num_iter;
    corr_vec[1] = status_corr.iter_converged;
    for (int i=0; i<5; ++i) {
      corr_vec[2+i] = status_corr.error[i];
    }
    double corr_vec_global[7] = { 0.0 };
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 7, corr_vec, corr_vec_global);    
    gkyl_dynvec_append(lte->corr_stat, app->tcurr, corr_vec_global);

    lte->n_iter += status_corr.num_iter;
  } 

  app->stat.neut_species_lte_tm += gkyl_time_diff_now_sec(wst);   
}

// Compute equivalent f_lte from fin
void
gk_neut_species_lte(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_lte *lte, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();
  gk_neut_species_moment_calc(&lte->moms, species->local, app->local, fin);

  // divide out the Jacobian from the density
  gkyl_dg_div_op_range(lte->moms.mem_geo, app->basis, 
    0, lte->moms.marr, 0, lte->moms.marr, 0, 
    app->gk_geom->jacobgeo, &app->local);  
  app->stat.neut_species_lte_tm += gkyl_time_diff_now_sec(wst);   

  gk_neut_species_lte_from_moms(app, species, lte, lte->moms.marr);
}

void
gk_neut_species_lte_write_max_corr_status(gkyl_gyrokinetic_app* app, struct gk_neut_species *gk_ns)
{
  if (gk_ns->lte.correct_all_moms) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out correction status 
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gk_ns->info.name, "corr-max-stat");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_ns->info.name, "corr-max-stat");

      if (gk_ns->lte.is_first_corr_status_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(gk_ns->lte.corr_stat, fileNm);
        gk_ns->lte.is_first_corr_status_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(gk_ns->lte.corr_stat, fileNm);
      }
    }
    gkyl_dynvec_clear(gk_ns->lte.corr_stat);

    app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_neut_diag_io += 1;
  }
}

void 
gk_neut_species_lte_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lte *lte)
{
  gkyl_array_release(lte->f_lte);

  gk_neut_species_moment_release(app, &lte->moms);

  gkyl_vlasov_lte_proj_on_basis_release(lte->proj_lte);
  if (lte->correct_all_moms) {
    gkyl_vlasov_lte_correct_release(lte->corr_lte);
    gkyl_dynvec_release(lte->corr_stat);
  }
}
