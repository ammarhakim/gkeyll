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
    gkyl_array_copy(gks->f_mult.multiplier_host, gks->f_mult.multiplier);
  }

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->f_mult.multiplier_host, fileNm);
  app->stat.n_io += 1;

  gk_array_meta_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

// Write only on init - writes once then disables further writes.
static void
gk_species_time_dilation_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_time_dilation_write_enabled(app, gks, tm, frame);
  gks->f_mult.write_func = gk_species_time_dilation_write_disabled;
}

// Multiply f by the multiplier (f = f * multiplier).
static void
gk_species_time_dilation_advance_mul_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *fmul, struct gkyl_array *f)
{
  gkyl_array_scale_by_cell(f, fmul->multiplier);
}

// Divide f by the multiplier (f = f / multiplier).
static void
gk_species_time_dilation_advance_div_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *fmul, struct gkyl_array *f)
{
  gkyl_array_divide_by_cell(f, fmul->multiplier);
}

// Disabled advance functions - do nothing.
static void
gk_species_time_dilation_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *fmul, struct gkyl_array *f)
{
}

// Multiply f by the time dilation mask (f = f * mask).
// The mask is computed based on f using gkyl_dg_array_mask.
static void
gk_species_time_dilation_advance_mul_mask(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *fmul, struct gkyl_array *f)
{
  // Compute the mask based on current f.
  gkyl_dg_array_mask_advance(fmul->f_mask, f);
  
  // Multiply f by the mask.
  const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(fmul->f_mask);
  gkyl_array_scale_by_cell(f, mask_array);
}

// Divide f by the time dilation mask (f = f / mask).
// The mask is computed based on f using gkyl_dg_array_mask.
static void
gk_species_time_dilation_advance_div_mask(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *fmul, struct gkyl_array *f)
{
  // Compute the mask based on current f.
  gkyl_dg_array_mask_advance(fmul->f_mask, f);
  
  // Divide f by the mask.
  const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(fmul->f_mask);
  gkyl_array_divide_by_cell(f, mask_array);
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
gk_species_time_dilation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_time_dilation *fmul)
{
  fmul->type = gks->info.time_dilation.type;
  fmul->write_diagnostics = gks->info.time_dilation.write_diagnostics;

  // Default function pointers - disabled.
  fmul->write_func = gk_species_time_dilation_write_disabled;
  fmul->advance_mul_func = gk_species_time_dilation_advance_disabled;
  fmul->advance_div_func = gk_species_time_dilation_advance_disabled;

  if (fmul->type) {

    // Create a basis for the multiplier.
    struct gkyl_basis basis_mult;
    gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);

    // Allocate multiplier array.
    fmul->multiplier = mkarr(app->use_gpu, basis_mult.num_basis, gks->local_ext.volume);
    fmul->multiplier_host = app->use_gpu ? mkarr(false, basis_mult.num_basis, gks->local_ext.volume)
                                         : gkyl_array_acquire(fmul->multiplier);

    // Context for c2p function passed to proj_on_basis.
    fmul->proj_on_basis_c2p_ctx.cdim = app->cdim;
    fmul->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
    fmul->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
    fmul->proj_on_basis_c2p_ctx.pos_map = app->position_map;

    if (fmul->type == GKYL_GK_F_MULTIPLIER_USER_INPUT) {
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &basis_mult,
          .num_quad = basis_mult.poly_order+1,
          .num_ret_vals = 1,
          .eval = gks->info.time_dilation.profile,
          .ctx = gks->info.time_dilation.profile_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &fmul->proj_on_basis_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, fmul->multiplier_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(fmul->multiplier, fmul->multiplier_host);

      fmul->advance_mul_func = gk_species_time_dilation_advance_mul_mult;
      fmul->advance_div_func = gk_species_time_dilation_advance_div_mult;
      if (fmul->write_diagnostics) {
        fmul->write_func = gk_species_time_dilation_write_init_only;
      } else {
        gkyl_array_release(fmul->multiplier_host);
      }
    }
    else if (fmul->type == GKYL_GK_F_MULTIPLIER_TIME_DILATION) {
      // Copy input parameters to struct.
      fmul->time_dilation_f_threshold = gks->info.time_dilation.time_dilation_f_threshold;
      fmul->time_dilation_f_frac = gks->info.time_dilation.time_dilation_f_frac;
      fmul->time_dilation_spatial_frac = gks->info.time_dilation.time_dilation_spatial_frac;

      // Determine mask type based on input parameters.
      enum gkyl_dg_array_mask_types mask_type = GKYL_DG_ARRAY_MASK_NONE;
      if (fmul->time_dilation_f_frac > 0.0) {
        mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD;
      } else if (fmul->time_dilation_f_threshold > 0.0) {
        mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD;
      }
      if (fmul->time_dilation_spatial_frac) {
        mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL;
      }

      // Create mask object for time dilation.
      struct gkyl_dg_array_mask_inp f_mask_inp = {
        .type = mask_type,
        .val_threshold = fmul->time_dilation_f_threshold,
        .frac_threshold = fmul->time_dilation_f_frac,
        .phase_rng = gks->local_ext,
        .config_rng = app->local_ext,
        .vel_rng = gks->local_ext_vel,
        .use_gpu = app->use_gpu,
      };
      fmul->f_mask = gkyl_dg_array_mask_new(f_mask_inp);

      // Initialize multiplier to 1.0.
      gkyl_array_clear(fmul->multiplier, 1.0);

      fmul->advance_mul_func = gk_species_time_dilation_advance_mul_mask;
      fmul->advance_div_func = gk_species_time_dilation_advance_div_mask;
      if (fmul->write_diagnostics) {
        fmul->write_func = gk_species_time_dilation_write_enabled;
      } else {
        gkyl_array_release(fmul->multiplier_host);
      }
    }
  }
}

void
gk_species_time_dilation_advance_div(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *fmul, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  fmul->advance_div_func(app, gks, fmul, f);

  app->stat.species_f_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_time_dilation_advance_mul(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_time_dilation *fmul, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  fmul->advance_mul_func(app, gks, fmul, f);

  app->stat.species_f_mult_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_time_dilation_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->f_mult.write_func(app, gks, tm, frame);
}

void
gk_species_time_dilation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_time_dilation *fmul)
{
  if (fmul->type) {
    gkyl_array_release(fmul->multiplier);
    if (fmul->write_diagnostics) {
      gkyl_array_release(fmul->multiplier_host);
    }

    if (fmul->type == GKYL_GK_F_MULTIPLIER_USER_INPUT) {
      // Nothing to release.
    }
    else if (fmul->type == GKYL_GK_F_MULTIPLIER_TIME_DILATION) {
      gkyl_dg_array_mask_release(fmul->f_mask);
    }
  }
}