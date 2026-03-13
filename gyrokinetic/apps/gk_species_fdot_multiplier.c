#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_array_dg_find_peaks.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>

void
gk_species_fdot_multiplier_write_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
}

void
gk_species_fdot_multiplier_write_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // DG metadata for multiplier.
  struct gkyl_msgpack_map_elem mpe_mult[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
  };
  int mpe_mult_len = sizeof(mpe_mult) / sizeof(mpe_mult[0]);
  // Update app basic metada with time/frame.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  // Package metadata.
  int io_meta_len[] = { app->io_meta_basic_len, mpe_mult_len, app->gk_geom->io_meta_len };
  const struct gkyl_msgpack_map_elem *io_meta[] = { app->io_meta_basic, mpe_mult,
                                                    app->gk_geom->io_meta };
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int),
    io_meta_len, io_meta);

  // Write out the multiplicative function.
  const char *fmt = "%s-%s_fdot_multiplier_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->fdot_mult.multiplier_host, gks->fdot_mult.multiplier);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->fdot_mult.multiplier_host,
    fileNm);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_fdot_multiplier_write_init_only(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
  gk_species_fdot_multiplier_write_enabled(app, gks, tm, frame);
  gks->fdot_mult.write_func = gk_species_fdot_multiplier_write_disabled;
}

void
gk_species_fdot_multiplier_advance_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  // Multiply out by the multplier.
  gkyl_array_scale_by_cell(out, fdmul->multiplier);
}

void
gk_species_fdot_multiplier_advance_omegaH_mult(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double *out)
{
  // Multiply out by the multplier.
  out[0] = out[0] / gks->collisionless.scale_fac;
}

void
gk_species_fdot_multiplier_advance_omegaH_disabled(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double *out)
{
}

void
gk_species_fdot_multiplier_advance_loss_cone_mult(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, phi, fdmul->phi_smooth_global);
  // Find the potential at bmag_max
  gkyl_array_dg_find_peaks_project_on_peak_idx(fdmul->bmag_peak_finder, fdmul->phi_smooth_global,
    fdmul->bmag_max_peak_idx, fdmul->phi_at_bmag_max);

  if (fdmul->is_tandem) {
    gkyl_array_dg_find_peaks_project_on_peak_idx(fdmul->bmag_peak_finder, fdmul->phi_smooth_global,
      fdmul->bmag_tandem_peak_idx, fdmul->phi_at_bmag_tandem);
    gkyl_loss_cone_mask_gyrokinetic_advance(fdmul->lcm_proj_op, &gks->local, &app->local,
      phi, fdmul->phi_at_bmag_max, fdmul->phi_at_bmag_tandem,
      fdmul->multiplier);
  }
  else {
    gkyl_loss_cone_mask_gyrokinetic_advance(fdmul->lcm_proj_op, &gks->local, &app->local,
      phi, fdmul->phi_at_bmag_max, fdmul->phi_at_bmag_max, fdmul->multiplier);
  }

  // Multiply out by the multiplier.
  gkyl_array_scale_by_cell(out, fdmul->multiplier);
}

void
gk_species_fdot_multiplier_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

static void
proj_on_basis_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
}

void
gk_species_fdot_multiplier_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul)
{
  fdmul->type = gks->info.time_rate_multiplier.type;
  fdmul->write_diagnostics = gks->info.time_rate_multiplier.write_diagnostics;

  // Default function pointers.
  fdmul->write_func = gk_species_fdot_multiplier_write_disabled;
  fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_disabled;
  fdmul->advance_times_omegaH_func = gk_species_fdot_multiplier_advance_omegaH_disabled;
  fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_disabled;

  if (fdmul->type) {
    bool cellwise_const = gks->info.time_rate_multiplier.cellwise_const;
    assert(cellwise_const); // MF 2025/06/11: Limited to this for now.

    // Create a basis for the multiplier.
    struct gkyl_basis basis_mult;
    if (cellwise_const)
      gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);
    else
      basis_mult = gks->basis;

    // Allocate multiplier array.
    fdmul->multiplier = mkarr(app->use_gpu, basis_mult.num_basis, gks->local_ext.volume);
    fdmul->multiplier_host = app->use_gpu? mkarr(false, fdmul->multiplier->ncomp,
      fdmul->multiplier->size)
                                         : gkyl_array_acquire(fdmul->multiplier);

    // Context for c2p function passed to proj_on_basis.
    fdmul->proj_on_basis_c2p_ctx.cdim = app->cdim;
    fdmul->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
    fdmul->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
    fdmul->proj_on_basis_c2p_ctx.pos_map = app->position_map;

    if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_USER_INPUT) {

      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
        .grid = &gks->grid,
        .basis = &basis_mult,
        .num_quad = basis_mult.poly_order + 1,
        .num_ret_vals = 1,
        .eval = gks->info.time_rate_multiplier.profile,
        .ctx = gks->info.time_rate_multiplier.profile_ctx,
        .c2p_func = proj_on_basis_c2p_phase_func,
        .c2p_func_ctx = &fdmul->proj_on_basis_c2p_ctx,
      });
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, fdmul->multiplier_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(fdmul->multiplier, fdmul->multiplier_host);

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_mult;
      fdmul->advance_times_omegaH_func = gk_species_fdot_multiplier_advance_omegaH_mult;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      if (fdmul->write_diagnostics)
        fdmul->write_func = gk_species_fdot_multiplier_write_init_only;
      else
        gkyl_array_release(fdmul->multiplier_host);
    }
    else if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE) {
      // Available options:
      // A) num_quad=1, qtype=GKYL_GAUSS_QUAD. Output: ncomp=1 array.
      // B) num_quad>1, qtype=GKYL_GAUSS_QUAD or GKYL_GAUSS_LOBATTO_QUAD, cellwise_const=true. Output: ncomp=1 array.
      enum gkyl_quad_type qtype = GKYL_GAUSS_LOBATTO_QUAD;
      int num_quad = gks->basis.poly_order + 1; // This can be p+1 or 1. Must be
                                                // at least p+1 for Gauss-Lobatto.

      // Create peak finder for bmag to find the mirror throat.
      // Search along the parallel (z) direction, which is the last configuration space dimension.
      int search_dir = app->cdim - 1;
      struct gkyl_array_dg_find_peaks_inp peak_inp = {
        .basis = &app->basis,
        .grid = &app->grid,
        .range = &app->global,
        .range_ext = &app->global_ext,
        .search_dir = search_dir,
        .use_gpu = app->use_gpu,
      };
      // Pass a global bmag_int into the peak finder
      struct gkyl_array *bmag_int_global = mkarr(app->use_gpu,
        app->gk_geom->geo_int.bmag->ncomp, app->global_ext.volume);
      fdmul->phi_smooth_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

      gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->geo_int.bmag,
        bmag_int_global);

      fdmul->bmag_peak_finder = gkyl_array_dg_find_peaks_new(&peak_inp, bmag_int_global);
      gkyl_array_dg_find_peaks_advance(fdmul->bmag_peak_finder, bmag_int_global);
      gkyl_array_release(bmag_int_global);

      // Get the LOCAL_MAX peak (bmag maximum along z direction).
      int num_peaks = gkyl_array_dg_find_peaks_num_peaks(fdmul->bmag_peak_finder);
      fdmul->bmag_max_peak_idx = num_peaks - 2; // Edge is num_peaks-1, so maximum is one less
      fdmul->bmag_max = gkyl_array_dg_find_peaks_acquire_vals(fdmul->bmag_peak_finder,
        fdmul->bmag_max_peak_idx);
      fdmul->bmag_max_z_coord = gkyl_array_dg_find_peaks_acquire_coords(fdmul->bmag_peak_finder,
        fdmul->bmag_max_peak_idx);
      fdmul->bmag_wall = gkyl_array_dg_find_peaks_acquire_vals(fdmul->bmag_peak_finder,
        num_peaks - 1);
      fdmul->bmag_wall_z_coord = gkyl_array_dg_find_peaks_acquire_coords(fdmul->bmag_peak_finder,
        num_peaks - 1);
      fdmul->bmag_max_basis = gkyl_array_dg_find_peaks_get_basis(fdmul->bmag_peak_finder);
      fdmul->bmag_max_range = gkyl_array_dg_find_peaks_get_range(fdmul->bmag_peak_finder);
      fdmul->bmag_max_range_ext = gkyl_array_dg_find_peaks_get_range_ext(fdmul->bmag_peak_finder);

      fdmul->phi_at_bmag_max = mkarr(app->use_gpu, fdmul->bmag_max_basis->num_basis,
        fdmul->bmag_max_range_ext->volume);
      fdmul->phi_at_bmag_tandem = mkarr(app->use_gpu, fdmul->bmag_max_basis->num_basis,
        fdmul->bmag_max_range_ext->volume);
      // phi is defined as 0 at the wall

      bool is_symmetric, is_tandem;
      int cdim = app->cdim;
      if (gkyl_compare_double(-app->grid.lower[cdim - 1], app->grid.upper[cdim - 1], 1e-12)) {
        is_symmetric = true;
      }
      else if (gkyl_compare_double(app->grid.lower[cdim - 1], 0.0, 1e-12)) {
        is_symmetric = false;
      }
      else {
        assert(false); // Needs either the lower bound at 0 or symmetric grid
      }

      if ( (is_symmetric && num_peaks == 5) || (!is_symmetric && num_peaks == 3) ) {
        is_tandem = false;
      }
      else if ((is_symmetric && num_peaks == 9) || (!is_symmetric && num_peaks == 5)) {
        is_tandem = true;
      }
      else {
        assert(false); // Unsupported number of extrema for loss-cone multiplier
      }

      if (is_tandem) {
        fdmul->bmag_tandem_peak_idx = num_peaks - 4;
      }
      else {
        fdmul->bmag_tandem_peak_idx = num_peaks - 2;
      }
      fdmul->bmag_tandem = gkyl_array_dg_find_peaks_acquire_vals(fdmul->bmag_peak_finder,
        fdmul->bmag_tandem_peak_idx);
      fdmul->bmag_tandem_z_coord = gkyl_array_dg_find_peaks_acquire_coords(fdmul->bmag_peak_finder,
        fdmul->bmag_tandem_peak_idx);

      // Operator that projects the loss cone mask.
      struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
        .phase_grid = &gks->grid,
        .conf_basis = &app->basis,
        .phase_basis = &gks->basis,
        .conf_range = &app->local,
        .conf_range_ext = &app->local_ext,
        .vel_range = &gks->local_vel,
        .vel_map = gks->vel_map,
        .bmag = app->gk_geom->geo_int.bmag,
        .bmag_max = fdmul->bmag_max,
        .bmag_max_z_coord = fdmul->bmag_max_z_coord,
        .bmag_wall = fdmul->bmag_wall,
        .bmag_wall_z_coord = fdmul->bmag_wall_z_coord,
        .bmag_tandem = fdmul->bmag_tandem,
        .bmag_tandem_z_coord = fdmul->bmag_tandem_z_coord,
        .is_tandem = is_tandem,
        .bmag_max_basis = fdmul->bmag_max_basis,
        .bmag_max_range = fdmul->bmag_max_range,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
        .qtype = qtype,
        .num_quad = num_quad,
        .cellwise_trap_loss = cellwise_const,
        .use_gpu = app->use_gpu,
      };
      fdmul->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_loss_cone_mult;
      fdmul->advance_times_omegaH_func = gk_species_fdot_multiplier_advance_omegaH_mult;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      if (fdmul->write_diagnostics) {
        fdmul->write_func = gk_species_fdot_multiplier_write_enabled;
      }
      else {
        gkyl_array_release(fdmul->multiplier_host);
      }
    }
  }
}

void
gk_species_fdot_multiplier_advance_times_cfl(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_cfl_func(app, gks, fdmul, phi, out);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_fdot_multiplier_advance_times_omegaH(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double *out)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_omegaH_func(app, gks, fdmul, out);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_fdot_multiplier_advance_times_rate(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_rate_func(app, gks, fdmul, phi, out);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_fdot_multiplier_write(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm,
  int frame)
{
  gks->fdot_mult.write_func(app, gks, tm, frame);
}

void
gk_species_fdot_multiplier_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_fdot_multiplier *fdmul)
{
  if (fdmul->type) {
    gkyl_array_release(fdmul->multiplier);
    if (fdmul->write_diagnostics) {
      gkyl_array_release(fdmul->multiplier_host);
    }

    if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_USER_INPUT) {
      // Nothing to release.
    }
    else if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE) {
      gkyl_array_release(fdmul->bmag_max);
      gkyl_array_release(fdmul->bmag_max_z_coord);
      gkyl_array_release(fdmul->bmag_wall);
      gkyl_array_release(fdmul->bmag_wall_z_coord);
      gkyl_array_release(fdmul->bmag_tandem);
      gkyl_array_release(fdmul->bmag_tandem_z_coord);

      gkyl_array_release(fdmul->phi_at_bmag_max);
      gkyl_array_release(fdmul->phi_at_bmag_tandem);

      gkyl_array_release(fdmul->phi_smooth_global);
      gkyl_array_dg_find_peaks_release(fdmul->bmag_peak_finder);
      gkyl_loss_cone_mask_gyrokinetic_release(fdmul->lcm_proj_op);
    }
  }
}

void
gk_species_fdot_multiplier_reset(gkyl_gyrokinetic_app *app, double tm, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp)
{
  gk_species_fdot_multiplier_release(app, fdmul);

  gks->info.time_rate_multiplier = fdot_mult_inp;
  gk_species_fdot_multiplier_init(app, gks, fdmul);
}
