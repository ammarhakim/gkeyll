#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_anomalous_diff_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_species_anomalous_diff_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_dg_updater_gk_anomalous_diffusion_advance(gkad->slvr, &species->local, 
    fin, species->cflrate, rhs);

  app->stat.species_diffusion_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_anomalous_diff_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, double tm, int frame)
{
}

static void
gk_species_anomalous_diff_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_anomalous_diff_write_conf_array(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, int frame, double stime, char* file_suffix,
  struct gkyl_array *arrout, struct gkyl_array *arrout_host)
{
  // Write out a conf-space array.
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = stime,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );
  // Construct the file handles for collision frequency and primitive moments.
  const char *fmt = "%s-%s_%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, file_suffix, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, file_suffix, frame);

  struct gkyl_array *arr_ho;
  if (app->use_gpu) {  
    if (arrout_host)
      arr_ho = gkyl_array_acquire(arrout_host);
    else {
      arr_ho = mkarr(false, arrout->ncomp, arrout->size);
    }
    // Copy data from device to host before writing it out.
    gkyl_array_copy(arr_ho, arrout);
  }
  else {
    if (arrout_host)
      arr_ho = gkyl_array_acquire(arrout_host);
    else
      arr_ho = gkyl_array_acquire(arrout);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, arr_ho, fileNm);
  gk_array_meta_release(mt); 
  gkyl_array_release(arr_ho);
}

static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}

void 
gk_species_anomalous_diff_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_anomalous_diff *gkad)
{
  gkad->anom_diff_id = gks->info.anomalous_diffusion.anomalous_diff_id;
  gkad->write_diagnostics = gks->info.anomalous_diffusion.write_diagnostics;

  gkad->write_diags_func = gk_species_anomalous_diff_write_diags_disabled;
  gkad->rhs_func = gk_species_anomalous_diff_rhs_disabled;

  if (gkad->anom_diff_id) {

    // Project the diffusivity.
    gkad->diffD = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *diffD_ho = app->use_gpu? mkarr(false, gkad->diffD->ncomp, gkad->diffD->size)
                                              : gkyl_array_acquire(gkad->diffD);
    struct gkyl_eval_on_nodes *diffDproj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .num_ret_vals = 1,
      .eval = gks->info.anomalous_diffusion.D_profile,
      .ctx = gks->info.anomalous_diffusion.D_profile_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
    });
    gkyl_eval_on_nodes_advance(diffDproj, 0.0, &app->local, diffD_ho);
    gkyl_array_copy(gkad->diffD, diffD_ho);

    gkyl_eval_on_nodes_release(diffDproj);
    gkyl_array_release(diffD_ho);

    // Multiply diffD by g^xx*jacobgeo.
    gkyl_dg_mul_op(app->basis, 0, gkad->diffD, 0, app->gk_geom->gxxj, 0, gkad->diffD);

    // Use zero-flux BCs for the anomalous diffusion updater. NOTE: this
    // doesn't actually mean there is no flux due to this term, we simply use
    // it as a way to trigger use of boundary flux kernels.
    int pdim = app->cdim + app->vdim;
    bool is_zero_flux[2*GKYL_MAX_DIM] = {false};
    is_zero_flux[0] = true;
    is_zero_flux[0+pdim] = true;

    gkad->slvr = gkyl_dg_updater_gk_anomalous_diffusion_new(&gks->grid, &gks->basis, &app->basis,
      &app->local, is_zero_flux, gks->info.skip_cell_threshold,
      gkad->diffD, app->gk_geom->jacobgeo_inv, app->use_gpu);

    if (gkad->write_diagnostics) {
      // Write out the diffusivity.
      gk_anomalous_diff_write_conf_array(app, gks, gkad, 0, 0.0, "anomalous_diffusivity", gkad->diffD, 0);
    }

    // Methods chosen at runtime.
    gkad->rhs_func = gk_species_anomalous_diff_rhs_enabled;
    if (gkad->write_diagnostics) {
      gkad->write_diags_func = gk_species_anomalous_diff_write_diags_enabled;
    }
  }
}

void
gk_species_anomalous_diff_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkad->rhs_func(app, species, gkad, fin, rhs);
}

void
gk_species_anomalous_diff_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, double tm, int frame)
{
  gkad->write_diags_func(app, gks, gkad, tm, frame);
}

void
gk_species_anomalous_diff_release(const struct gkyl_gyrokinetic_app *app, const struct gk_anomalous_diff *gkad)
{
  if (gkad->anom_diff_id) {

    gkyl_array_release(gkad->diffD);
    gkyl_dg_updater_gk_anomalous_diffusion_release(gkad->slvr);

//    if (gkad->write_diagnostics) {
//    }
  }
}
