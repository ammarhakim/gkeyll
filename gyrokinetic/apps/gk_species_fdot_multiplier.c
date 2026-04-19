#include <assert.h>
#include <gkyl_alloc.h>
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
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, phi, fdmul->phi_global);

  gkyl_loss_cone_mask_gyrokinetic_advance(fdmul->lcm_gk, &gks->local, &app->global,
    fdmul->bmag_global, fdmul->phi_global, fdmul->multiplier);

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
      // Operator that projects the loss cone mask.
      struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
        .conf_basis = &app->basis,
        .vel_map = gks->vel_map,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
      };
      fdmul->lcm_gk = gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);
      fdmul->bmag_global = mkarr(app->use_gpu, app->gk_geom->geo_corn.bmag->ncomp,
        app->global_ext.volume);
      fdmul->phi_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

      gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->geo_corn.bmag,
        fdmul->bmag_global);

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
      gkyl_array_release(fdmul->bmag_global);
      gkyl_array_release(fdmul->phi_global);
      gkyl_loss_cone_mask_gyrokinetic_release(fdmul->lcm_gk);
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
