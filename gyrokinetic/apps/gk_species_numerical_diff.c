#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_numerical_diff_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_numerical_diff *gknd, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_species_numerical_diff_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_numerical_diff *gknd, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_dg_updater_gk_numerical_diffusion_advance(gknd->slvr, &species->local, 
    fin, species->cflrate, rhs);

  app->stat.species_nume_diffusion_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_numerical_diff_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_numerical_diff *gknd, double tm, int frame)
{
}

static void
gk_species_numerical_diff_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_numerical_diff *gknd, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_numerical_diff_write_conf_array(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_numerical_diff *gknd, int frame, double stime, char* file_suffix, char* description,
  struct gkyl_array *arrout, struct gkyl_array *arrout_host)
{
  // Write out a conf-space array.
  
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(app->io_meta_grid_len, app->io_meta_grid, "time", stime);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_grid_len, app->io_meta_grid, "frame", frame);
  struct gkyl_msgpack_map_elem desc[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = description }
  };
  int io_meta_len[] = {app->io_meta_grid_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_grid, app->gk_geom->io_meta_basic, desc};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

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
  gkyl_msgpack_data_release(mt); 
  gkyl_array_release(arr_ho);
}

static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}

void 
gk_species_numerical_diff_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_numerical_diff *gknd)
{
  gknd->num_diff_id = gks->info.numerical_diffusion.type;
  gknd->write_diagnostics = gks->info.numerical_diffusion.write_diagnostics;

  gknd->write_diags_func = gk_species_numerical_diff_write_diags_disabled;
  gknd->rhs_func = gk_species_numerical_diff_rhs_disabled;

  if (gknd->num_diff_id) {

    bool is_diff_dir[GKYL_MAX_CDIM] = {false};

    int num_diff_dir = gks->info.numerical_diffusion.num_dirs;
    assert(num_diff_dir > 0 && num_diff_dir <= app->cdim);

    int szD = num_diff_dir*app->basis.num_basis;
    gknd->diffD = mkarr(app->use_gpu, szD, app->local_ext.volume);

    struct gkyl_array *diffD1 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *diffD1_ho = app->use_gpu? mkarr(false, diffD1->ncomp, diffD1->size)
                                              : gkyl_array_acquire(diffD1);

    // Project the diffusivity.
    for (int d=0; d<num_diff_dir; ++d) {
      int dir = gks->info.numerical_diffusion.dirs[d];
      is_diff_dir[dir] = true;

      struct gkyl_eval_on_nodes *diffDproj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp){
        .grid = &app->grid,
        .basis = &app->basis,
        .num_ret_vals = 1,
        .eval = gks->info.numerical_diffusion.D_profile[d],
        .ctx = gks->info.numerical_diffusion.D_profile_ctx,
        .c2p_func = eval_on_nodes_c2p_position_func,
        .c2p_func_ctx = app->position_map,
      });
      gkyl_eval_on_nodes_advance(diffDproj, 0.0, &app->local, diffD1_ho);
      gkyl_eval_on_nodes_release(diffDproj);

      gkyl_array_copy(diffD1, diffD1_ho);
      gkyl_array_accumulate_offset(gknd->diffD, 1.0, diffD1, d*app->basis.num_basis);

      // Multiply diffD by g^{ii}*jacobgeo.
      struct gkyl_array *giiJ;
      if (dir == 0)
        giiJ = app->gk_geom->geo_int.gxxj;
      else if (dir == 1)
        giiJ = app->gk_geom->geo_int.gyyj;
      else if (dir == 2) {
        fprintf(stderr, "Numerical diffusion in z not available. Need g^{zz}.\n");
        assert(false);
//        giiJ = app->gk_geom->geo_int.gzzj;
      }
      gkyl_dg_mul_op(app->basis, d, gknd->diffD, 0, giiJ, d, gknd->diffD);

    }

    gkyl_array_release(diffD1_ho);
    gkyl_array_release(diffD1);

    // Sync diffusivity.
    int num_periodic_dir = app->num_periodic_dir;
    gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
      num_periodic_dir, app->periodic_dirs, gknd->diffD); 

    // ABSORB and FIXED_FUNC BCs need to fill the ghost cell. MF 2025/11/03: The
    // only sensible way to fill the ghost cell is to use the value at the
    // boundary.
    for (int d=0; d<num_diff_dir; ++d) {
      int dir = gks->info.numerical_diffusion.dirs[d];

      for (int b=0; b<2; ++b) {
        if ((b == 0 && ((gks->lower_bc[0].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) || (gks->lower_bc[0].type == GKYL_BC_GK_SPECIES_ABSORB))) || 
            (b == 1 && ((gks->upper_bc[0].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) || (gks->upper_bc[0].type == GKYL_BC_GK_SPECIES_ABSORB))) ) {
          enum gkyl_edge_loc edge = b==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
          struct gkyl_range *skin_r = b==0? &app->local_lower_skin[dir] : &app->local_upper_skin[dir];
          struct gkyl_range *ghost_r = b==0? &app->local_lower_ghost[dir] : &app->local_upper_ghost[dir];

          long vol = skin_r->volume;
          long buff_sz = 1;
          buff_sz = buff_sz > vol ? buff_sz : vol;
          struct gkyl_array *bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

          struct gkyl_bc_basic_gyrokinetic *gfss_bc_op = gkyl_bc_basic_gyrokinetic_new(dir, edge,
            GKYL_BC_GK_FIELD_BOUNDARY_VALUE, app->basis_on_dev, skin_r, ghost_r, szD, app->cdim, app->use_gpu);
          gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, bc_buffer, gknd->diffD);
          gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, bc_buffer, app->gk_geom->geo_int.jacobgeo_inv);

          gkyl_bc_basic_gyrokinetic_release(gfss_bc_op);
          gkyl_array_release(bc_buffer);
        }
      }
    }
    gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, gknd->diffD);

    int diffusion_order = gks->info.numerical_diffusion.order ? gks->info.numerical_diffusion.order : 4;
    assert(diffusion_order == 2 || diffusion_order == 4);

    enum gkyl_gyrokinetic_bc_type bc_type_lower[GKYL_MAX_CDIM] = {0}, bc_type_upper[GKYL_MAX_CDIM] = {0};
    for (int d=0; d<num_diff_dir; ++d) {
      int dir = gks->info.numerical_diffusion.dirs[d];
      bc_type_lower[dir] = gks->lower_bc[dir].type;
      bc_type_upper[dir] = gks->upper_bc[dir].type;
    }

    // Create solver.
    gknd->slvr = gkyl_dg_updater_gk_numerical_diffusion_new(&gks->grid, &gks->basis, &app->basis,
      &app->local, bc_type_lower, bc_type_upper, is_diff_dir, diffusion_order, 
      gknd->diffD, app->gk_geom->geo_int.jacobgeo_inv, app->use_gpu);

    if (gknd->write_diagnostics) {
      // Write out the diffusivity.
      gk_numerical_diff_write_conf_array(app, gks, gknd, 0, 0.0, "nume_diff", "Numerical diffusivity", gknd->diffD, 0);
    }

    // Methods chosen at runtime.
    gknd->rhs_func = gk_species_numerical_diff_rhs_enabled;
    if (gknd->write_diagnostics) {
      gknd->write_diags_func = gk_species_numerical_diff_write_diags_enabled;
    }
  }
}

void
gk_species_numerical_diff_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_numerical_diff *gknd, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gknd->rhs_func(app, species, gknd, fin, rhs);
}

void
gk_species_numerical_diff_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_numerical_diff *gknd, double tm, int frame)
{
  gknd->write_diags_func(app, gks, gknd, tm, frame);
}

void
gk_species_numerical_diff_release(const struct gkyl_gyrokinetic_app *app, const struct gk_numerical_diff *gknd)
{
  if (gknd->num_diff_id) {

    gkyl_array_release(gknd->diffD);
    gkyl_dg_updater_gk_numerical_diffusion_release(gknd->slvr);

//    if (gknd->write_diagnostics) {
//    }
  }
}
