#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_neut_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->proj_id = inp.proj_id;
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid,
        .basis = &s->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = s->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.func,
        .ctx = inp.ctx_func,
      }
    );
    if (app->use_gpu) {
      proj->proj_host = mkarr(false, s->basis.num_basis, s->local_ext.volume);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    int vdim = app->vdim+1; // neutral species are 3v otherwise
    proj->dens = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->udrift = mkarr(false, vdim*app->basis.num_basis, app->local_ext.volume);
    proj->vtsq = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->prim_moms_host = mkarr(false, (2+vdim)*app->basis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(app->use_gpu, (2+vdim)*app->basis.num_basis, app->local_ext.volume);

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      s->basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_udrift = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      s->basis.poly_order+1, vdim, inp.udrift, inp.ctx_udrift);
    proj->proj_temp = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      s->basis.poly_order+1, 1, inp.temp, inp.ctx_temp);


    struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
      .phase_grid = &s->grid,
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
    proj->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

    proj->correct_all_moms = false; 
    if (inp.correct_all_moms) {
      proj->correct_all_moms = true;

      struct gkyl_vlasov_lte_correct_inp inp_corr = {
        .phase_grid = &s->grid,
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
        .max_iter = 100,
        .eps = 1e-12,
      };
      proj->corr_lte = gkyl_vlasov_lte_correct_inew( &inp_corr );
    }
  }
}

void
gk_neut_species_projection_calc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    if (app->use_gpu) {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local_ext, proj->proj_host);
      gkyl_array_copy(f, proj->proj_host);
    }
    else {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local_ext, f);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    int vdim = app->vdim+1;
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_udrift, tm, &app->local, proj->udrift);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local, proj->vtsq);
    gkyl_array_scale(proj->vtsq, 1/s->info.mass);

    // Projection routines expect the LTE moments as a single array.
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->udrift, 1*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsq, (vdim+1)*app->basis.num_basis);

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);

    // Multiply density by the conf-space jacobian.
    gkyl_dg_mul_op_range(app->basis, 0, proj->prim_moms, 
      0, app->gk_geom->jacobgeo, 0, proj->prim_moms, &app->local);

    // Project the Maxwellian distribution function.
    // Projection routine also corrects the density of the projected distribution function.
    gkyl_vlasov_lte_proj_on_basis_advance(proj->proj_lte, &s->local, &app->local, 
      proj->prim_moms, f);

    // Correct all the moments of the projected Maxwellian distribution function.
    if (proj->correct_all_moms) {
      struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(proj->corr_lte, 
        f, proj->prim_moms, &s->local, &app->local);
    } 
  }
}

void
gk_neut_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    gkyl_proj_on_basis_release(proj->proj_func);
    if (app->use_gpu) {
      gkyl_array_release(proj->proj_host);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) { 
    gkyl_array_release(proj->dens);
    gkyl_array_release(proj->udrift); 
    gkyl_array_release(proj->vtsq);
    gkyl_array_release(proj->prim_moms_host);
    gkyl_array_release(proj->prim_moms);

    gkyl_proj_on_basis_release(proj->proj_dens);
    gkyl_proj_on_basis_release(proj->proj_udrift);
    gkyl_proj_on_basis_release(proj->proj_temp);
    
    gkyl_vlasov_lte_proj_on_basis_release(proj->proj_lte);
    if (proj->correct_all_moms) {
      gkyl_vlasov_lte_correct_release(proj->corr_lte);
    }
  } 
}
