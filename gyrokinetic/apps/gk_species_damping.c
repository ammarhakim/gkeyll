#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_array_dg_find_peaks.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>

static void
proj_on_basis_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
}

void
gk_species_damping_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
}

void
gk_species_damping_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // DG metadata for damping rate.
  struct gkyl_msgpack_map_elem mpe_drate[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
  };
  int mpe_drate_len = sizeof(mpe_drate)/sizeof(mpe_drate[0]);
  // Update app basic metada with time/frame.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  // Package metadata.
  int io_meta_len[] = {app->io_meta_basic_len, mpe_drate_len, app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, mpe_drate, app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Write out the damping rate.
  const char *fmt = "%s-%s_damping_rate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->damping.rate_host, gks->damping.rate);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->damping.rate_host, fileNm);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_damping_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_damping_write_enabled(app, gks, tm, frame);
  gks->damping.write_func = gk_species_damping_write_disabled;
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

void
gk_species_damping_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_damping *damp)
{
  damp->type = gks->info.damping.type;
  damp->evolve = false; // Whether the rate is time dependent.
  damp->is_tandem = false; // Default to single mirror.

  int num_quad = gks->info.damping.num_quad? gks->info.damping.num_quad : 1; // Default is a p=0 mask.
  assert(num_quad == 1); // MF 2025/06/11: Limited to this for now.

  // Default function pointers.
  damp->write_func = gk_species_damping_write_disabled;


  damp->proj_on_basis_c2p_ctx.cdim = app->cdim;
  damp->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
  damp->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
  damp->proj_on_basis_c2p_ctx.pos_map = app->position_map;


  if (damp->type) {
    // Allocate rate array.
    damp->rate = mkarr(app->use_gpu, num_quad==1? 1 : gks->basis.num_basis, gks->local_ext.volume);
    damp->rate_host = damp->rate;
    if (app->use_gpu)
      damp->rate_host = mkarr(false, damp->rate->ncomp, damp->rate->size); 

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
      proj_on_basis_c2p_ctx.cdim = app->cdim;
      proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
      proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &gks->basis,
          .num_quad = num_quad,
          .num_ret_vals = 1,
          .eval = gks->info.damping.rate_profile,
          .ctx = gks->info.damping.rate_profile_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &proj_on_basis_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, damp->rate_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(damp->rate, damp->rate_host);

      if (num_quad == 1)
        gkyl_array_scale_range(damp->rate, 1.0/pow(sqrt(2.0),gks->grid.ndim), &gks->local);
    }
    else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
      damp->evolve = true; // Since the loss cone boundary is proportional to phi(t).

      // Create peak finder for bmag to find the mirror throat.
      // Search along the parallel (z) direction, which is the last configuration space dimension.
      int search_dir = app->cdim - 1;
      struct gkyl_array_dg_find_peaks_inp peak_inp = {
        .basis = &app->basis,
        .grid = &app->grid,
        .range = &app->local,
        .range_ext = &app->local_ext,
        .search_dir = search_dir,
        .use_gpu = app->use_gpu,
      };
      // Pass a global bmag_int into the peak finder
      struct gkyl_array *bmag_int_global = mkarr(false, 
        app->gk_geom->geo_int.bmag->ncomp, app->gk_geom->geo_int.bmag->size);
      gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->geo_int.bmag, bmag_int_global);
      damp->bmag_peak_finder = gkyl_array_dg_find_peaks_new(&peak_inp, bmag_int_global);
      gkyl_array_dg_find_peaks_advance(damp->bmag_peak_finder, app->gk_geom->geo_int.bmag);
      gkyl_array_release(bmag_int_global);
      
      // Get the LOCAL_MAX peak (bmag maximum along z direction).
      int num_peaks = gkyl_array_dg_find_peaks_num_peaks(damp->bmag_peak_finder);
      damp->bmag_max_peak_idx = num_peaks-2; // Edge is num_peaks-1, so maximum is one less
      damp->bmag_max = gkyl_array_dg_find_peaks_acquire_vals(damp->bmag_peak_finder, damp->bmag_max_peak_idx);
      damp->bmag_max_z_coord = gkyl_array_dg_find_peaks_acquire_coords(damp->bmag_peak_finder, damp->bmag_max_peak_idx);
      damp->bmag_wall = gkyl_array_dg_find_peaks_acquire_vals(damp->bmag_peak_finder, num_peaks-1);
      damp->bmag_wall_z_coord = gkyl_array_dg_find_peaks_acquire_coords(damp->bmag_peak_finder, num_peaks-1);
      damp->bmag_max_basis = gkyl_array_dg_find_peaks_get_basis(damp->bmag_peak_finder);
      damp->bmag_max_range = gkyl_array_dg_find_peaks_get_range(damp->bmag_peak_finder);
      damp->bmag_max_range_ext = gkyl_array_dg_find_peaks_get_range_ext(damp->bmag_peak_finder);

      damp->phi_at_bmag_max = mkarr(app->use_gpu, damp->bmag_max_basis->num_basis, 
        damp->bmag_max_range_ext->volume);
      damp->phi_at_bmag_tandem = mkarr(app->use_gpu, damp->bmag_max_basis->num_basis, 
        damp->bmag_max_range_ext->volume);
      // phi is defined as 0 at the wall

      bool is_symmetric;
      int cdim = app->cdim;
      if (gkyl_compare_double(-app->grid.lower[cdim-1], app->grid.upper[cdim-1], 1e-12)) {
        is_symmetric = true;
      }
      else if (gkyl_compare_double(app->grid.lower[cdim-1], 0.0, 1e-12)) {
        is_symmetric = false;
      }
      else {
        assert(false); // Needs either the lower bound at 0 or symmetric grid
      }

      if ( (is_symmetric && num_peaks == 5) || (!is_symmetric && num_peaks == 3) ) {
        damp->is_tandem = false;
      }
      else if ((is_symmetric && num_peaks == 9) || (!is_symmetric && num_peaks == 5)) {
        damp->is_tandem = true;
      }
      else {
        assert(false); // Unsupported number of extrema for loss-cone damping
      }

      if (damp->is_tandem) {
        damp->bmag_tandem_peak_idx = num_peaks-4;
      } else {
        damp->bmag_tandem_peak_idx = num_peaks-2;
      }
      damp->bmag_tandem = gkyl_array_dg_find_peaks_acquire_vals(damp->bmag_peak_finder, damp->bmag_tandem_peak_idx);
      damp->bmag_tandem_z_coord = gkyl_array_dg_find_peaks_acquire_coords(damp->bmag_peak_finder, damp->bmag_tandem_peak_idx);


      // Operator that projects the loss cone mask.
      struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
        .phase_grid = &gks->grid,
        .conf_basis = &app->basis,
        .phase_basis = &gks->basis,
        .conf_range =  &app->local,
        .conf_range_ext = &app->local_ext,
        .vel_range = &gks->local_vel, 
        .vel_map = gks->vel_map,
        .bmag = app->gk_geom->geo_int.bmag,
        .bmag_max = damp->bmag_max,
        .bmag_max_z_coord = damp->bmag_max_z_coord,
        .bmag_wall = damp->bmag_wall,
        .bmag_wall_z_coord = damp->bmag_wall_z_coord,
        .bmag_tandem = damp->bmag_tandem,
        .bmag_tandem_z_coord = damp->bmag_tandem_z_coord,
        .is_tandem = damp->is_tandem,
        .bmag_max_basis = damp->bmag_max_basis,
        .bmag_max_range = damp->bmag_max_range,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
        .num_quad = num_quad,
        .c2p_pos_func = proj_on_basis_c2p_position_func,
        .c2p_pos_func_ctx = &damp->proj_on_basis_c2p_ctx,
        .use_gpu = app->use_gpu,
      };
      damp->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew( &inp_proj );

      // Project the conf-space rate profile provided.
      struct gkyl_array *scale_prof_high_order = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
      struct gkyl_array *scale_prof_high_order_ho = app->use_gpu? mkarr(false, scale_prof_high_order->ncomp, scale_prof_high_order->size)
                                                     : gkyl_array_acquire(scale_prof_high_order);
      
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_new(&gks->grid, &gks->basis, num_quad, 1, 
        gks->info.damping.rate_profile, gks->info.damping.rate_profile_ctx);
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, scale_prof_high_order_ho);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(scale_prof_high_order, scale_prof_high_order_ho);

      damp->scale_prof = mkarr(app->use_gpu, num_quad == 1? 1 : gks->basis.num_basis, gks->local_ext.volume);
      gkyl_array_set_offset(damp->scale_prof, pow(sqrt(2.0),gks->grid.ndim), scale_prof_high_order, 0);

      gkyl_array_release(scale_prof_high_order_ho);
      gkyl_array_release(scale_prof_high_order);

      // Compute the initial damping rate (assuming phi=0 because phi hasn't been computed).
      // Find the potential at the mirror throat.
      gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, app->field->phi_smooth,
        damp->bmag_max_peak_idx, damp->phi_at_bmag_max);
      
      if (damp->is_tandem) {
        gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, app->field->phi_smooth,
          damp->bmag_tandem_peak_idx, damp->phi_at_bmag_tandem);
        // Project the loss cone mask.
        gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
          app->field->phi_smooth, damp->phi_at_bmag_max, damp->phi_at_bmag_tandem, damp->rate);
      } else {
        // Project the loss cone mask using the phi_m array.
        gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
          app->field->phi_smooth, damp->phi_at_bmag_max, damp->phi_at_bmag_max, damp->rate);
      }
      // Multiply by the user's scaling profile.
      gkyl_array_scale_by_cell(damp->rate, damp->scale_prof);
    }

    // Set function pointers chosen at runtime.
    if (damp->evolve) {
      damp->write_func = gk_species_damping_write_enabled;
    }
    else {
      damp->write_func = gk_species_damping_write_init_only;
    }
  }
}

void
gk_species_damping_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_damping *damp, 
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  if (damp->type) {
    struct timespec wst = gkyl_wall_clock();
    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      gkyl_array_set(f_buffer, 1.0, fin);
      gkyl_array_scale_by_cell(f_buffer, damp->rate);
      gkyl_array_accumulate(rhs, -1.0, f_buffer);
    }
    else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
      // Find the potential at all peak locations (including the mirror throat).
      gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, phi,
        damp->bmag_max_peak_idx, damp->phi_at_bmag_max);
      // Allgather on phi_at_bmag_max. It's not an allgather.
      // One process has the correct one, but the others do not. Is it a bcast or a sync?

      if (damp->is_tandem) {
        gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, phi,
          damp->bmag_tandem_peak_idx, damp->phi_at_bmag_tandem);
        // Allgather on phi_at_bmag_tandem. It's not an allgather.
        // One process has the correct one, but the others do not. Is it a bcast or a sync?
        gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
          phi, damp->phi_at_bmag_max, damp->phi_at_bmag_tandem, damp->rate);
      } else {
        // Project the loss cone mask using the phi_m array.
        gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
          phi, damp->phi_at_bmag_max, damp->phi_at_bmag_max, damp->rate);
      }

      // Assemble the damping term -scale_prof * mask * f.
      gkyl_array_set(f_buffer, 1.0, fin);
      gkyl_array_scale_by_cell(damp->rate, damp->scale_prof);
      gkyl_array_scale_by_cell(f_buffer, damp->rate);
      gkyl_array_accumulate(rhs, -1.0, f_buffer);

    }

    // Add the frequency to the CFL frequency.
    gkyl_array_accumulate(cflrate, 1.0, damp->rate);

    app->stat.species_damp_tm += gkyl_time_diff_now_sec(wst);
  }
}

void
gk_species_damping_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->damping.write_func(app, gks, tm, frame);
}

void
gk_species_damping_release(const struct gkyl_gyrokinetic_app *app, const struct gk_damping *damp)
{
  if (damp->type) {
    gkyl_array_release(damp->rate);
    if (app->use_gpu) {
      gkyl_array_release(damp->rate_host);
    }

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      // Nothing to release.
    }
    else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
      gkyl_array_release(damp->bmag_max);
      gkyl_array_release(damp->bmag_max_z_coord);
      gkyl_array_release(damp->bmag_wall);
      gkyl_array_release(damp->bmag_wall_z_coord);
      gkyl_array_release(damp->bmag_tandem);
      gkyl_array_release(damp->bmag_tandem_z_coord);

      gkyl_array_release(damp->phi_at_bmag_max);
      gkyl_array_release(damp->phi_at_bmag_tandem);

      gkyl_array_dg_find_peaks_release(damp->bmag_peak_finder);
      gkyl_loss_cone_mask_gyrokinetic_release(damp->lcm_proj_op);
      gkyl_array_release(damp->scale_prof);
    }
  }
}
