#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_dg_array_mask.h>

// Disabled write function - does nothing.
static void
gk_species_time_dilation_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
}

// Enabled write function - writes the multiplier array to file.
static void
gk_species_time_dilation_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = 0,
      .basis_type = "tensor",
    }, GKYL_GK_META_NONE, 0
  );

  // Write out the multiplicative function.
  const char *fmt = "%s-%s_time_dilation_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(gks->time_dilation.multiplier_host, gks->time_dilation.multiplier);
  }

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->time_dilation.multiplier_host, fileNm);
  app->stat.n_io += 1;

  gk_array_meta_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

// Write only on init - writes once then disables further writes.
static void
gk_species_time_dilation_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_time_dilation_write_enabled(app, gks, tm, frame);
  gks->time_dilation.write_func = gk_species_time_dilation_write_disabled;
}

// Set the multiplier array based on cflrate and time dilation parameters.
// See https://arxiv.org/html/2510.09756
void
gk_species_time_dilation_advance_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, const struct gkyl_array *f, const struct gkyl_array *cflrate)
{
  // Compute omega_max - a ceiling on omega is a floor on dt.
  // WARNING: dt_omegaH is DBL_MAX for boltzmann and adiabatic fields!
  double omega_max = DBL_MAX;

  if (time_dilation->enable_cfl_dt_floor) {
    // Use omega_H based CFL dt flooring.
    if (time_dilation->cfl_dt_min_omegaH) {
      omega_max = (gks->dt_omegaH > 1e-30) ? 1.0 / gks->dt_omegaH : DBL_MAX;
    }
    // Also apply user-specified minimum dt if provided.
    if (time_dilation->cfl_dt_min_value > 0.0) {
      double omega_from_user = 1.0 / time_dilation->cfl_dt_min_value;
      omega_max = fmin(omega_max, omega_from_user); // Take the largest timestep.
    }
  } else if (time_dilation->enable_mask_based_omega) {
    // Use mask-based approach to find omega_max from masked cells.
    gkyl_dg_array_mask_advance(time_dilation->cfl_mask, f);
    gkyl_dg_array_mask_scale_by_cell(time_dilation->cfl_mask, cflrate);
    const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(time_dilation->cfl_mask);

    double omega_max_local;
    if (app->use_gpu) {
    #ifdef GKYL_HAVE_CUDA
      gkyl_array_reduce(time_dilation->omega_max_local_cu, mask_array, GKYL_MAX);
      gkyl_cu_memcpy(&omega_max_local, time_dilation->omega_max_local_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
    #endif
    } else {
      gkyl_array_reduce(&omega_max_local, mask_array, GKYL_MAX);
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &omega_max_local, &omega_max);
  }

  // Compute multiplier = min(1.0, omega_max / omega_cfl).
  gkyl_array_copy(time_dilation->multiplier, cflrate);
  gkyl_array_invert_by_cell(time_dilation->multiplier); // 1/omega_cfl
  gkyl_array_scale(time_dilation->multiplier, omega_max); // omega_max / omega_cfl
  gkyl_array_min_by_cell(time_dilation->multiplier, 1.0); // min(1.0, omega_max / omega_cfl)
}

// Multiply f by the multiplier (f = f * multiplier).
static void
gk_species_time_dilation_mul_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, struct gkyl_array *f)
{
  gkyl_array_scale_by_cell(f, time_dilation->multiplier);
}

// Divide f by the multiplier (f = f / multiplier).
static void
gk_species_time_dilation_div_enabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, struct gkyl_array *f)
{
  gkyl_array_divide_by_cell(f, time_dilation->multiplier);
}

// Disabled advance functions - do nothing.
static void
gk_species_time_dilation_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, const struct gkyl_array *f, const struct gkyl_array *cflrate){}
static void
gk_species_time_dilation_mul_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, struct gkyl_array *f){}
static void
gk_species_time_dilation_div_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, struct gkyl_array *f){}


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
gk_species_time_dilation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_time_dilation *time_dilation)
{
  time_dilation->type = gks->info.time_dilation.type;
  time_dilation->write_diagnostics = gks->info.time_dilation.write_diagnostics;

  // Default function pointers - disabled.
  time_dilation->write_func = gk_species_time_dilation_write_disabled;
  time_dilation->advance_func = gk_species_time_dilation_advance_disabled;
  time_dilation->mul_func = gk_species_time_dilation_mul_disabled;
  time_dilation->div_func = gk_species_time_dilation_div_disabled;

  if (time_dilation->type) {

    // Create a basis for the multiplier.
    struct gkyl_basis basis_mult;
    gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);

    // Allocate multiplier array.
    time_dilation->multiplier = mkarr(app->use_gpu, 1, gks->local_ext.volume);
    time_dilation->multiplier_host = app->use_gpu ? mkarr(false, 1, gks->local_ext.volume)
                                         : gkyl_array_acquire(time_dilation->multiplier);

    // Context for c2p function passed to proj_on_basis.
    time_dilation->proj_on_basis_c2p_ctx.cdim = app->cdim;
    time_dilation->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
    time_dilation->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
    time_dilation->proj_on_basis_c2p_ctx.pos_map = app->position_map;

    if (time_dilation->type == GKYL_GK_F_MULTIPLIER_USER_INPUT) {
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &basis_mult,
          .num_quad = basis_mult.poly_order+1,
          .num_ret_vals = 1,
          .eval = gks->info.time_dilation.profile,
          .ctx = gks->info.time_dilation.profile_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &time_dilation->proj_on_basis_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, time_dilation->multiplier_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(time_dilation->multiplier, time_dilation->multiplier_host);

      time_dilation->mul_func = gk_species_time_dilation_mul_enabled;
      time_dilation->div_func = gk_species_time_dilation_div_enabled;
      if (time_dilation->write_diagnostics) {
        time_dilation->write_func = gk_species_time_dilation_write_init_only;
      } else {
        gkyl_array_release(time_dilation->multiplier_host);
      }
    }
    else if (time_dilation->type == GKYL_GK_CFL_DRIVEN_TIME_DILATION) {
      // Copy input parameters to struct.
      time_dilation->cfl_dt_min_omegaH = gks->info.time_dilation.cfl_dt_min_omegaH;
      time_dilation->cfl_dt_min_value = gks->info.time_dilation.cfl_dt_min_value;
      time_dilation->time_dilation_f_threshold = gks->info.time_dilation.time_dilation_f_threshold;
      time_dilation->time_dilation_f_frac = gks->info.time_dilation.time_dilation_f_frac;
      time_dilation->time_dilation_spatial_frac = gks->info.time_dilation.time_dilation_spatial_frac;

      // Determine which time dilation features are enabled.
      bool has_cfl_dt_floor = time_dilation->cfl_dt_min_omegaH || 
                              (time_dilation->cfl_dt_min_value > 0.0);
      bool has_mask_threshold = (time_dilation->time_dilation_f_threshold > 0.0) ||
                                (time_dilation->time_dilation_f_frac > 0.0);

      time_dilation->enable_cfl_dt_floor = has_cfl_dt_floor;
      time_dilation->enable_mask_based_omega = has_mask_threshold && !has_cfl_dt_floor;

      // Allocate GPU scratch space for reduce operation if using GPU.
      if (app->use_gpu) {
      #ifdef GKYL_HAVE_CUDA
        time_dilation->omega_max_local_cu = (double*) gkyl_cu_malloc(sizeof(double));
      #endif
      }

      // Create mask object if using mask-based time dilation.
      if (time_dilation->enable_mask_based_omega) {
        enum gkyl_dg_array_mask_types mask_type = GKYL_DG_ARRAY_MASK_NONE;
        if (time_dilation->time_dilation_f_frac > 0.0) {
          mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD;
        } else if (time_dilation->time_dilation_f_threshold > 0.0) {
          mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD;
        }
        if (time_dilation->time_dilation_spatial_frac) {
          mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL;
        }

        struct gkyl_dg_array_mask_inp cfl_mask_inp = {
          .type = mask_type,
          .val_threshold = time_dilation->time_dilation_f_threshold,
          .frac_threshold = time_dilation->time_dilation_f_frac,
          .phase_rng = gks->local_ext,
          .config_rng = app->local_ext,
          .vel_rng = gks->local_ext_vel,
          .use_gpu = app->use_gpu,
        };
        time_dilation->cfl_mask = gkyl_dg_array_mask_new(cfl_mask_inp);
      }

      // Initialize multiplier to 1.0.
      gkyl_array_clear(time_dilation->multiplier, 1.0);

      time_dilation->advance_func = gk_species_time_dilation_advance_enabled;
      time_dilation->mul_func = gk_species_time_dilation_mul_enabled;
      time_dilation->div_func = gk_species_time_dilation_div_enabled;
      if (time_dilation->write_diagnostics) {
        time_dilation->write_func = gk_species_time_dilation_write_enabled;
      } else {
        gkyl_array_release(time_dilation->multiplier_host);
      }
    }
  }
}

void
gk_species_time_dilation_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, const struct gkyl_array *f, const struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  time_dilation->advance_func(app, gks, time_dilation, f, cflrate);

  app->stat.species_f_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_time_dilation_div(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  time_dilation->div_func(app, gks, time_dilation, f);

  app->stat.species_f_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_time_dilation_mul(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *time_dilation, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  time_dilation->mul_func(app, gks, time_dilation, f);

  app->stat.species_f_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_time_dilation_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->time_dilation.write_func(app, gks, tm, frame);
}

void
gk_species_time_dilation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_time_dilation *time_dilation)
{
  if (time_dilation->type) {
    gkyl_array_release(time_dilation->multiplier);
    if (time_dilation->write_diagnostics) {
      gkyl_array_release(time_dilation->multiplier_host);
    }

    if (time_dilation->type == GKYL_GK_CFL_DRIVEN_TIME_DILATION) {

      // Free GPU scratch space if it was allocated.
      if (app->use_gpu) {
      #ifdef GKYL_HAVE_CUDA
        gkyl_cu_free(time_dilation->omega_max_local_cu);
      #endif
      }

      if (time_dilation->enable_mask_based_omega) {
        gkyl_dg_array_mask_release(time_dilation->cfl_mask);
      }
    }
  }
}