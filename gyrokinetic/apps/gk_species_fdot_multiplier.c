#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>

void
gk_species_fdot_multiplier_write_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double tm, int frame)
{
}

void
gk_species_fdot_multiplier_write_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double tm, int frame)
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
  const char *fmt = gks->num_fdot_mult > 1 ?
    "%s-%s_fdot_multiplier_c%d_%d.gkyl" : "%s-%s_fdot_multiplier_%d.gkyl";
  int sz = gks->num_fdot_mult > 1 ?
    gkyl_calc_strlen(fmt, app->name, gks->info.name, fdmul->component_id, frame) :
    gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  if (gks->num_fdot_mult > 1) {
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, fdmul->component_id,
      frame);
  }
  else {
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
  }

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(fdmul->multiplier_host, fdmul->multiplier);
  }

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, fdmul->multiplier_host,
    fileNm);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_fdot_multiplier_write_init_only(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double tm, int frame)
{
  gk_species_fdot_multiplier_write_enabled(app, gks, fdmul, tm, frame);
  fdmul->write_func = gk_species_fdot_multiplier_write_disabled;
}

void
gk_species_fdot_multiplier_advance_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f,
  struct gkyl_array *cflrate)
{
  // Multiply cflrate by the multplier.
  gkyl_array_scale_by_cell(cflrate, fdmul->multiplier);
}

void
gk_species_fdot_multiplier_advance_omegaH_mult(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double *out)
{
  // Multiply out by the multplier.
  out[0] = out[0] / gks->collisionless.scale_fac / fdmul->time_dilation_scale_const;
}

void
gk_species_fdot_multiplier_advance_omegaH_disabled(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, double *out)
{
}

void
gk_species_fdot_multiplier_advance_loss_cone_mult(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, phi, fdmul->phi_global);

  gkyl_loss_cone_mask_gyrokinetic_advance(fdmul->lcm_gk, &gks->local, &app->global,
    fdmul->bmag_global, fdmul->phi_global, fdmul->multiplier);

  // Offset by a constant factor
  gkyl_array_scale(fdmul->multiplier, fdmul->time_dilation_scale_const);

  // Multiply cflrate by the multplier.
  gkyl_array_scale_by_cell(cflrate, fdmul->multiplier);
}

// Set the multiplier array based on cflrate and time dilation parameters.
// See https://arxiv.org/html/2510.09756
static double
compute_global_array_max(const gkyl_gyrokinetic_app *app, const struct gk_fdot_multiplier *fdmul,
  const struct gkyl_array *array)
{
  // Compute the global maximum of the array on the host
  double local_max;
  if (app->use_gpu) {
#ifdef GKYL_HAVE_CUDA
    gkyl_array_reduce(fdmul->omega_max_local_cu, array, GKYL_MAX);
    gkyl_cu_memcpy(&local_max, fdmul->omega_max_local_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
#endif
  }
  else {
    gkyl_array_reduce(&local_max, array, GKYL_MAX);
  }
  double global_max = DBL_MAX;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &local_max, &global_max);
  return global_max;
}

static void
clamp_cflrate_by_omega_max(const struct gk_fdot_multiplier *fdmul, const double omega_max,
  struct gkyl_array *multiplier, struct gkyl_array *cflrate)
{
  // multiplier = max(1.0, omega_max/(cflrate)) * time_dilation_scale_const.
  gkyl_array_copy(multiplier, cflrate);
  gkyl_array_invert_by_cell(multiplier);
  gkyl_array_scale(multiplier, omega_max);
  gkyl_array_min_by_cell(multiplier, 1.0);
  gkyl_array_scale(multiplier, fdmul->time_dilation_scale_const);
  gkyl_array_scale_by_cell(cflrate, multiplier);
}

void
gk_species_fdot_multiplier_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f,
  struct gkyl_array *cflrate)
{
}

void
gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_omegaH(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  double omega_max = 1.0 / gks->dt_omegaH;
  clamp_cflrate_by_omega_max(fdmul, omega_max, fdmul->multiplier, cflrate);
}

void
gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_user_specified(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  double omega_max = 1.0 / fdmul->cfl_dt_min_value;
  clamp_cflrate_by_omega_max(fdmul, omega_max, fdmul->multiplier, cflrate);
}

void
gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_set_by_species(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  double omega_max = 1.0 / fdmul->species_dt_is_set_from->dt_cfl_global_ho[0];
  clamp_cflrate_by_omega_max(fdmul, omega_max, fdmul->multiplier, cflrate);
}

void
gk_species_fdot_multiplier_advance_time_dilation_cfl_factor_user_specified(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  double omega_max = compute_global_array_max(app, fdmul, cflrate);
  omega_max = fdmul->cfl_factor_times_omega_max * omega_max;
  clamp_cflrate_by_omega_max(fdmul, omega_max, fdmul->multiplier, cflrate);
}

void
gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_global(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  gkyl_array_reduce(fdmul->local_max_f, f, GKYL_MAX);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, app->basis.num_basis,
    fdmul->local_max_f, fdmul->global_max_f);
  gkyl_dg_array_mask_advance_threshold(fdmul->cfl_mask, fdmul->global_max_f[0]);

  gkyl_dg_array_mask_advance(fdmul->cfl_mask, f);
  gkyl_dg_array_mask_scale_by_cell(fdmul->cfl_mask, cflrate);
  const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(fdmul->cfl_mask);
  double omega_max = compute_global_array_max(app, fdmul, mask_array);
  clamp_cflrate_by_omega_max(fdmul, omega_max, fdmul->multiplier, cflrate);
}

void
gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_local(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  gkyl_dg_array_mask_advance(fdmul->cfl_mask, f);
  gkyl_dg_array_mask_scale_by_cell(fdmul->cfl_mask, cflrate);
  const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(fdmul->cfl_mask);
  double omega_max = compute_global_array_max(app, fdmul, mask_array);
  clamp_cflrate_by_omega_max(fdmul, omega_max, fdmul->multiplier, cflrate);
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
gk_species_fdot_multiplier_init(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul,
  const struct gkyl_gyrokinetic_fdot_multiplier *fdot_mult_inp, int component_id)
{
  fdmul->type = fdot_mult_inp->type;
  fdmul->component_id = component_id;
  fdmul->write_diagnostics = fdot_mult_inp->write_diagnostics;

  // Default function pointers.
  fdmul->write_func = gk_species_fdot_multiplier_write_disabled;
  fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_disabled;
  fdmul->advance_times_omegaH_func = gk_species_fdot_multiplier_advance_omegaH_disabled;
  fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_disabled;

  if (fdmul->type) {
    bool cellwise_const = fdot_mult_inp->cellwise_const;
    assert(cellwise_const); // MF 2025/06/11: Limited to this for now.

    // Create a basis for the multiplier.
    struct gkyl_basis basis_mult;
    if (cellwise_const) {
      gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);
    }
    else {
      basis_mult = gks->basis;
    }

    // Allocate multiplier array.
    fdmul->multiplier = mkarr(app->use_gpu, basis_mult.num_basis, gks->local_ext.volume);
    fdmul->multiplier_host = app->use_gpu? mkarr(false, fdmul->multiplier->ncomp,
      fdmul->multiplier->size) : gkyl_array_acquire(fdmul->multiplier);

    // Context for c2p function passed to proj_on_basis.
    fdmul->proj_on_basis_c2p_ctx.cdim = app->cdim;
    fdmul->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
    fdmul->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
    fdmul->proj_on_basis_c2p_ctx.pos_map = app->position_map;

    fdmul->time_dilation_scale_const = fdot_mult_inp->time_dilation_scale_const;
    if (fdmul->time_dilation_scale_const <= 0.0) {
      fdmul->time_dilation_scale_const = 1.0;
    }

    if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_USER_INPUT) {

      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
        .grid = &gks->grid,
        .basis = &basis_mult,
        .num_quad = basis_mult.poly_order + 1,
        .num_ret_vals = 1,
        .eval = fdot_mult_inp->profile,
        .ctx = fdot_mult_inp->profile_ctx,
        .c2p_func = proj_on_basis_c2p_phase_func,
        .c2p_func_ctx = &fdmul->proj_on_basis_c2p_ctx,
      });
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, fdmul->multiplier_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(fdmul->multiplier, fdmul->multiplier_host);

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_mult;
      fdmul->advance_times_omegaH_func = gk_species_fdot_multiplier_advance_omegaH_mult;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      if (fdmul->write_diagnostics) {
        fdmul->write_func = gk_species_fdot_multiplier_write_init_only;
      }
      else {
        gkyl_array_release(fdmul->multiplier_host);
      }
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
    else if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_DT_SET_BY_SPECIES) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_CONSTANT)) {

      // Copy input parameters to struct.
      fdmul->cfl_dt_min_value = fdot_mult_inp->cfl_dt_min_value;
      fdmul->f_threshold = fdot_mult_inp->f_threshold;
      fdmul->cfl_factor_times_omega_max = fdot_mult_inp->cfl_factor_times_omega_max;

      if (app->use_gpu) {
#ifdef GKYL_HAVE_CUDA
        fdmul->omega_max_local_cu = (double *)gkyl_cu_malloc(sizeof(double));
        fdmul->local_max_f = (double *)gkyl_cu_malloc(sizeof(double) * gks->basis.num_basis);
        fdmul->global_max_f = (double *)gkyl_cu_malloc(sizeof(double) * gks->basis.num_basis);
#endif
      }
      else {
        fdmul->local_max_f = gkyl_malloc(sizeof(double) * gks->basis.num_basis);
        fdmul->global_max_f = gkyl_malloc(sizeof(double) * gks->basis.num_basis);
      }

      enum gkyl_dg_array_mask_types mask_type = GKYL_DG_ARRAY_MASK_NONE;
      switch (fdmul->type) {
        case GKYL_GK_FDOT_MULTIPLIER_CONSTANT:
          mask_type = GKYL_DG_ARRAY_MASK_NONE;
          fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_mult;
          break;
        case GKYL_GK_FDOT_MULTIPLIER_FIXED_DT:
          mask_type = GKYL_DG_ARRAY_MASK_NONE;
          fdmul->advance_times_cfl_func =
            gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_user_specified;
          break;
        case GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX:
          mask_type = GKYL_DG_ARRAY_MASK_NONE;
          fdmul->advance_times_cfl_func =
            gk_species_fdot_multiplier_advance_time_dilation_cfl_factor_user_specified;
          break;
        case GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH:
          mask_type = GKYL_DG_ARRAY_MASK_NONE;
          fdmul->advance_times_cfl_func =
            gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_omegaH;
          break;
        case GKYL_GK_FDOT_MULTIPLIER_DT_SET_BY_SPECIES:
          mask_type = GKYL_DG_ARRAY_MASK_NONE;
          fdmul->species_dt_is_set_from = gk_find_species(app, fdot_mult_inp->dt_set_by_species);
          fdmul->advance_times_cfl_func =
            gk_species_fdot_multiplier_advance_time_dilation_cfl_dt_set_by_species;
          break;
        case GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD:
          mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER;
          fdmul->advance_times_cfl_func =
            gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_local;
          break;
        case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL:
          mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF;
          fdmul->advance_times_cfl_func =
            gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_local;
          break;
        case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL:
          mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC;
          fdmul->advance_times_cfl_func =
            gk_species_fdot_multiplier_advance_time_dilation_cfl_f_frac_global;
          break;
        default:
          mask_type = GKYL_DG_ARRAY_MASK_NONE;
          fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_disabled;
          break;
      }

      struct gkyl_dg_array_mask_inp cfl_mask_inp = {
        .type = mask_type,
        .threshold = fdmul->f_threshold,
        .phase_rng = &gks->local,
        .phase_rng_ext = &gks->local_ext,
        .conf_rng = &app->local,
        .conf_rng_ext = &app->local_ext,
        .vel_rng = &gks->local_vel,
        .use_gpu = app->use_gpu,
      };
      fdmul->cfl_mask = gkyl_dg_array_mask_new(cfl_mask_inp);

      // Initialize multiplier to 1.0.
      gkyl_array_clear(fdmul->multiplier, 1.0);

      if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_CONSTANT) {
        gkyl_array_clear(fdmul->multiplier, fdmul->time_dilation_scale_const);
      }
      
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
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_cfl_func(app, gks, fdmul, phi, f, cflrate);

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
  const struct gk_species *gks, struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi,
  const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_rate_func(app, gks, fdmul, phi, f, cflrate);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_fdot_multiplier_write(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm,
  int frame)
{
  for (int i = 0; i < gks->num_fdot_mult; ++i) {
    struct gk_fdot_multiplier *fdmul = &gks->fdot_mult[i];
    fdmul->write_func(app, gks, fdmul, tm, frame);
  }
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
    else if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_DT_SET_BY_SPECIES) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
      (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL)) {
      gkyl_dg_array_mask_release(fdmul->cfl_mask);
      if (app->use_gpu) {
#ifdef GKYL_HAVE_CUDA
        gkyl_cu_free(fdmul->omega_max_local_cu);
#endif
      }
      else {
        gkyl_free(fdmul->local_max_f);
        gkyl_free(fdmul->global_max_f);
      }
    }
  }
}

void
gk_species_fdot_multiplier_reset(gkyl_gyrokinetic_app *app, double tm, struct gk_species *gks,
  struct gkyl_gyrokinetic_fdot_multipliers fdot_mult_inp)
{
  for (int i = 0; i < gks->num_fdot_mult; ++i) {
    gk_species_fdot_multiplier_release(app, &gks->fdot_mult[i]);
  }

  gks->info.time_rate_multipliers = fdot_mult_inp;
  gks->num_fdot_mult = fdot_mult_inp.num_multipliers;
  
  for (int i = 0; i < gks->num_fdot_mult; ++i) {
    gk_species_fdot_multiplier_init(app, gks, &gks->fdot_mult[i],
      &fdot_mult_inp.multiplier[i], i);
  }
}
