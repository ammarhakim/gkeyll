#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_neut_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, struct gk_bgk_collisions *bgk)
{
  int cdim = app->cdim, vdim = 3;
  // allocate nu and initialize it
  bgk->nu_sum = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  bgk->self_nu = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->basis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(bgk->self_nu, self_nu);
  gkyl_array_copy(bgk->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  // Host-side copy for I/O
  bgk->nu_sum_host = bgk->nu_sum;
  if (app->use_gpu) {
    bgk->nu_sum_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  bgk->nu_fmax = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
  // BGK updater (also computes stable timestep)
  bgk->up_bgk = gkyl_bgk_collisions_new(&app->basis, &s->basis, app->use_gpu);
}

// computes moments
void
gk_neut_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  gk_neut_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  // divide out the Jacobian from the density
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, app->basis, 
    0, species->lte.moms.marr, 0, species->lte.moms.marr, 0, 
    app->gk_geom->jacobgeo, &app->local);  

  app->stat.neut_species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// updates the collision terms in the rhs
void
gk_neut_species_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(bgk->nu_fmax, 0.0);

  // Project the LTE distribution function from the computed LTE moments
  gk_neut_species_lte_from_moms(app, species, &species->lte, species->lte.moms.marr);

  // Multiply the Maxwellian by the configuration-space Jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, species->lte.f_lte, 
    app->gk_geom->jacobgeo, species->lte.f_lte, &app->local, &species->local);

  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, bgk->nu_fmax, 
    bgk->self_nu, species->lte.f_lte, &app->local, &species->local);

  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &species->local, 
    bgk->nu_sum, bgk->nu_fmax, fin, bgk->implicit_step, bgk->dt_implicit, rhs, species->cflrate);

  app->stat.neut_species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_neut_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk)
{
  gkyl_array_release(bgk->self_nu);
  gkyl_array_release(bgk->nu_sum);

  if (app->use_gpu) {
    gkyl_array_release(bgk->nu_sum_host);
  }

  gkyl_array_release(bgk->nu_fmax);
  gkyl_bgk_collisions_release(bgk->up_bgk);
}
