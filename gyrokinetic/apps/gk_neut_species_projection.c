#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void read_neut_ics(struct gkyl_gyrokinetic_app *app, struct gkyl_array *out, cstr fileNm)
{

  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &app->local, app->poly_order);
  struct gkyl_array* nnodal = mkarr(false, 1, nrange.volume);

  // See if we are in upper half
  bool reverse = false;
  char *p; 
  p = strstr(fileNm.str, "b3");
  if (p!= NULL) {
    p[1] = '1';
    reverse=true;
  }
  p = strstr(fileNm.str, "b4");
  if (p!= NULL) {
    p[1] = '0';
    reverse=true;
  }
  p = strstr(fileNm.str, "b5");
  if (p!= NULL) {
    p[1] = '9';
    reverse=true;
  }
  p = strstr(fileNm.str, "b8");
  if (p!= NULL) {
    p[1] = '6';
    reverse=true;
  }

  FILE *ptr = fopen(fileNm.str,"r");
  size_t status;

  int nr = gkyl_range_shape(&nrange, 0);
  int nz = gkyl_range_shape(&nrange, 1);
  int idx[2];
  if (reverse) {
    for(int ir = 0; ir < nr; ir++){
      idx[0] = ir;
      for(int iz = nz-1; iz >= 0; iz--){
        idx[1] = iz;
        // set psi
        double *nnodal_n = gkyl_array_fetch(nnodal, gkyl_range_idx(&nrange, idx));
        status = fscanf(ptr,"%lf", nnodal_n);
      }
    }
  }
  else {
    for(int ir = 0; ir < nr; ir++){
      idx[0] = ir;
      for(int iz = 0; iz < nz; iz++){
        idx[1] = iz;
        // set psi
        double *nnodal_n = gkyl_array_fetch(nnodal, gkyl_range_idx(&nrange, idx));
        status = fscanf(ptr,"%lf", nnodal_n);
      }
    }
  }

  fclose(ptr);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
  gkyl_nodal_ops_n2m(n2m, &app->basis, &app->grid, &nrange, &app->local, 1, nnodal, out, false);
  gkyl_array_release(nnodal);


  // Reflect for symmetric blocks
  // Reflect DG coeffs rather than nodal data to avoid symmetry errors in n2m conversion
  bool reflect = false;
  if (strstr(fileNm.str, "b7") != NULL)
    reflect = true;
  if (strstr(fileNm.str, "b2") != NULL)
    reflect = true;
  if (strstr(fileNm.str, "b10") != NULL)
    reflect = true;
  if (strstr(fileNm.str, "b11") != NULL)
    reflect = true;
  if (reflect) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[1] < gkyl_range_shape(&app->local,1)/2 +1 ) {
        int idx_change[2] = {iter.idx[0], gkyl_range_shape(&app->local, 1) - iter.idx[1]+1};
        const double *coeffs_ref = gkyl_array_cfetch(out, gkyl_range_idx(&app->local, iter.idx));
        double *coeffs  = gkyl_array_fetch(out, gkyl_range_idx(&app->local, idx_change));
        app->basis.flip_odd_sign( 1, coeffs_ref, coeffs);
      }
    }
  }
 
}

static void half_read_neut_ics(struct gkyl_gyrokinetic_app *app, struct gkyl_array *out, cstr fileNm)
{

  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &app->local, app->poly_order);
  struct gkyl_array* nnodal = mkarr(false, 1, nrange.volume);

  FILE *ptr = fopen(fileNm.str,"r");
  size_t status;

  int nr = gkyl_range_shape(&nrange, 0);
  int nz = gkyl_range_shape(&nrange, 1);
  int idx[2];

  for(int ir = 0; ir < nr; ir++){
    idx[0] = ir;
    for(int iz = 0; iz < nz; iz++){
      idx[1] = iz;
      // set psi
      double *nnodal_n = gkyl_array_fetch(nnodal, gkyl_range_idx(&nrange, idx));
      status = fscanf(ptr,"%lf", nnodal_n);
    }
  }

  fclose(ptr);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
  gkyl_nodal_ops_n2m(n2m, &app->basis, &app->grid, &nrange, &app->local, 1, nnodal, out, false);
  gkyl_array_release(nnodal);

}

static void
gk_neut_species_projection_kinetic_calc(gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
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
    int vdim = s->info.vdim;
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
      0, app->gk_geom->geo_int.jacobgeo, 0, proj->prim_moms, &app->local);

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

static void
gk_neut_species_projection_kinetic_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
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

static void 
gk_neut_species_projection_kinetic_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
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
    int vdim = s->info.vdim;
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
      .vel_map = s->vel_map,
      .phase_range = &s->local,
      .h_ij = s->g_ij,
      .h_ij_inv = s->gij,
      .det_h = app->gk_geom->geo_int.jacobgeo,
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
        .vel_map = s->vel_map,
        .phase_range = &s->local,
        .h_ij = s->g_ij,
        .h_ij_inv = s->gij,
        .det_h = app->gk_geom->geo_int.jacobgeo,
        .hamil = s->hamil,	
        .model_id = s->model_id,
        .use_gpu = app->use_gpu,
        .max_iter = 100,
        .eps = 1e-12,
      };
      proj->corr_lte = gkyl_vlasov_lte_correct_inew( &inp_corr );
    }
  }

  proj->neut_calc_func = gk_neut_species_projection_kinetic_calc;
  proj->release_func = gk_neut_species_projection_kinetic_release;
}

static void
gk_neut_species_projection_fluid_calc(gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
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
    cstr fileNm = cstr_from_fmt("gkeyll_text_input/%s-D0_M0.txt", app->name);
    half_read_neut_ics(app, proj->dens, fileNm);
    
    fileNm = cstr_from_fmt("gkeyll_text_input/%s-D0_ux.txt", app->name);
    half_read_neut_ics(app, proj->vtsq, fileNm);
    gkyl_array_set_offset(proj->udrift, 1.0, proj->vtsq, 0*app->basis.num_basis);
    fileNm = cstr_from_fmt("gkeyll_text_input/%s-D0_uy.txt", app->name);
    half_read_neut_ics(app, proj->vtsq, fileNm);
    gkyl_array_set_offset(proj->udrift, 1.0, proj->vtsq, 1*app->basis.num_basis);
    fileNm = cstr_from_fmt("gkeyll_text_input/%s-D0_uz.txt", app->name);
    half_read_neut_ics(app, proj->vtsq, fileNm);
    gkyl_array_set_offset(proj->udrift, 1.0, proj->vtsq, 2*app->basis.num_basis);

    fileNm = cstr_from_fmt("gkeyll_text_input/%s-D0_Temp.txt", app->name);
    half_read_neut_ics(app, proj->vtsq, fileNm);


    // f[0] = mass*dens
    gkyl_array_set_offset_range(s->f_host, s->info.mass, proj->dens, 0, &app->local);

    // f[1] = f[0]*udrift[0]
    // f[2] = f[0]*udrift[1]
    // f[3] = f[0]*udrift[2]
    for (int d=0; d<3; d++) {
      gkyl_dg_mul_op_range(app->basis, d+1, s->f_host, 0, s->f_host, d, proj->udrift, &app->local);
    }

    // f[4] = 0.5*rho*u^2 + p/(gas_gamma-1)
    //      = 0.5*(rho*ux^2+rho*uy^2+rho*uz^2) + dens*temp/(gas_gamma-1)
    //      = 0.5*(f[1].udrift[0]+f[2].udrift[1]+f[3].udrift[2]) + dens*temp/(gas_gamma-1)
    gkyl_dg_mul_op_range(app->basis, 0, proj->vtsq, 0, proj->dens, 0, proj->vtsq, &app->local);
    gkyl_array_set_offset_range(s->f_host, 1.0/(s->info.gas_gamma-1.0), proj->vtsq, 4*app->basis.num_basis, &app->local);
    for (int d=0; d<3; d++) {
      gkyl_dg_mul_op_range(app->basis, 0, proj->dens, d+1, s->f_host, d, proj->udrift, &app->local);
      gkyl_array_accumulate_offset_range(s->f_host, 0.5, proj->dens, 4*app->basis.num_basis, &app->local);
    }

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(f, s->f_host);

    // Multiply moments by the conf-space Jacobian.
    for (int d=0; d<s->num_moments; d++) {
      gkyl_dg_mul_op_range(app->basis, d, f, 0, app->gk_geom->geo_int.jacobgeo, d, f, &app->local);
    }

  }
}

static void
gk_neut_species_projection_fluid_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
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

    gkyl_proj_on_basis_release(proj->proj_dens);
    gkyl_proj_on_basis_release(proj->proj_udrift);
    gkyl_proj_on_basis_release(proj->proj_temp);
  } 
}

static void 
gk_neut_species_projection_fluid_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->proj_id = inp.proj_id;
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid,
        .basis = &s->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = s->basis.poly_order+1,
        .num_ret_vals = s->num_moments,
        .eval = inp.func,
        .ctx = inp.ctx_func,
      }
    );
    if (app->use_gpu) {
      proj->proj_host = mkarr(false, s->num_moments*s->basis.num_basis, s->local_ext.volume);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    int udim = s->num_moments-2;
    proj->dens = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->udrift = mkarr(false, udim*app->basis.num_basis, app->local_ext.volume);
    proj->vtsq = mkarr(false, app->basis.num_basis, app->local_ext.volume);

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      s->basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_udrift = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      s->basis.poly_order+1, udim, inp.udrift, inp.ctx_udrift);
    proj->proj_temp = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      s->basis.poly_order+1, 1, inp.temp, inp.ctx_temp);
  }

  proj->neut_calc_func = gk_neut_species_projection_fluid_calc;
  proj->release_func = gk_neut_species_projection_fluid_release;
}

void 
gk_neut_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  if (s->is_fluid)
    gk_neut_species_projection_fluid_init(app, s, inp, proj);
  else
    gk_neut_species_projection_kinetic_init(app, s, inp, proj);
}

void
gk_neut_species_projection_calc(gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  proj->neut_calc_func(app, s, proj, f, tm);
}

void
gk_neut_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
{
  proj->release_func(app, proj);
}
