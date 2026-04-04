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

// Damping state synchronization helpers.

void

gk_species_damping_set_fbar_to_f_enabled(const struct gk_species *gks, struct gk_damping *damp,
  const struct gkyl_array *f)
{
  gkyl_array_set_offset(damp->fbar, 1.0, f, 0);
  gkyl_array_set_offset(damp->fbar1, 1.0, f, 0);
  gkyl_array_set_offset(damp->fbarnew, 1.0, f, 0);
}

void
gk_species_damping_set_fbar_to_f_disabled(const struct gk_species *gks, struct gk_damping *damp,
  const struct gkyl_array *f)
{
}

void
gk_species_damping_set_fbar_to_f(const struct gk_species *gks, struct gk_damping *damp,
  const struct gkyl_array *f)
{
  damp->set_fbar_to_f_func(gks, damp, f);
}

// Damping diagnostics write helpers.

void
gk_species_damping_write_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm,
  int frame)
{
}

static void
gk_species_damping_write_rate_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
}

static void
gk_species_damping_write_rate_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
  // DG metadata for damping rate.
  struct gkyl_msgpack_map_elem mpe_drate[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
  };
  int mpe_drate_len = sizeof(mpe_drate) / sizeof(mpe_drate[0]);

  // Update app basic metadata with time/frame.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  // Package metadata.
  int io_meta_len[] = { app->io_meta_basic_len, mpe_drate_len, app->gk_geom->io_meta_len };
  const struct gkyl_msgpack_map_elem *io_meta[] = { app->io_meta_basic, mpe_drate,
                                                    app->gk_geom->io_meta };
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int),
    io_meta_len, io_meta);

  // Write out the damping rate.
  const char *fmt = "%s-%s_damping_rate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->damping.rate_host, gks->damping.rate);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->damping.rate_host, fileNm);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt);
}

static void
gk_species_damping_write_fbar_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
}

static void
gk_species_damping_write_fbar_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
  // Metadata from app/species/geometry. fbar is always stored as a p0 field.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);

  struct gkyl_msgpack_map_elem mpe_fbar_p0[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
  };
  int io_meta_fbar_len[] = {
    app->io_meta_basic_len,
    (int)(sizeof(mpe_fbar_p0) / sizeof(mpe_fbar_p0[0])),
    app->gk_geom->io_meta_len
  };
  const struct gkyl_msgpack_map_elem *io_meta_fbar[] = {
    app->io_meta_basic,
    mpe_fbar_p0,
    app->gk_geom->io_meta
  };
  struct gkyl_msgpack_data *mt_fbar = gkyl_msgpack_create_union(
    sizeof(io_meta_fbar_len) / sizeof(int), io_meta_fbar_len, io_meta_fbar);

  const char *fmt_fbar = "%s-%s_damping_fbar_%d.gkyl";
  int sz_fbar = gkyl_calc_strlen(fmt_fbar, app->name, gks->info.name, frame);
  char fileNm_fbar[sz_fbar + 1]; // ensures no buffer overflow
  snprintf(fileNm_fbar, sizeof fileNm_fbar, fmt_fbar, app->name, gks->info.name, frame);

  // Copy fbar from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->damping.fbar_host, gks->damping.fbar);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt_fbar, gks->damping.fbar_host,
    fileNm_fbar);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt_fbar);
}

void
gk_species_damping_write_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm,
  int frame)
{
  struct timespec wst = gkyl_wall_clock();
  gks->damping.write_rate_func(app, gks, tm, frame);

  gks->damping.write_fbar_func(app, gks, tm, frame);

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_damping_write_init_only(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm,
  int frame)
{
  gk_species_damping_write_enabled(app, gks, tm, frame);
  gks->damping.write_func = gk_species_damping_write_disabled;
}

// Damping rate projection helpers.

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

static void
gk_species_damping_project_phase_rate(const struct gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, int num_quad,
  void (*rate_profile)(double t, const double *xn, double *fout, void *ctx),
  void *rate_profile_ctx, struct gkyl_array *rate_host, struct gkyl_array *rate)
{
  struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
  proj_on_basis_c2p_ctx.cdim = app->cdim;
  proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
  proj_on_basis_c2p_ctx.vel_map = gks->vel_map;

  struct gkyl_basis basis_mult;
  gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);

  gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
    .grid = &gks->grid,
    .basis = &basis_mult,
    .num_quad = num_quad,
    .num_ret_vals = 1,
    .eval = rate_profile,
    .ctx = rate_profile_ctx,
    .c2p_func = proj_on_basis_c2p_phase_func,
    .c2p_func_ctx = &proj_on_basis_c2p_ctx,
  });
  gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, rate_host);
  gkyl_proj_on_basis_release(projup);
  gkyl_array_copy(rate, rate_host);

  if (num_quad == 1)
    gkyl_array_scale_range(rate, 1.0 / pow(sqrt(2.0), gks->grid.ndim), &gks->local);
}

// Damping RHS assembly dispatch helpers.

void
gk_species_damping_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
}

void
gk_species_damping_advance_user_input(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  gkyl_array_set(f_buffer, 1.0, fin);
  gkyl_array_scale_by_cell(f_buffer, damp->rate);

  // Add damping to f and the CFL frequency.
  gkyl_array_accumulate(rhs, -1.0, f_buffer);
  gkyl_array_accumulate(cflrate, 1.0, damp->rate);
}

void
gk_species_damping_advance_loss_cone(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, phi, damp->phi_smooth_global);

  // Find phi at bmag_max and (if present) at tandem mirror peak.
  gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, damp->phi_smooth_global,
    damp->bmag_max_peak_idx, damp->phi_at_bmag_max);

  if (damp->is_tandem) {
    gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, damp->phi_smooth_global,
      damp->bmag_tandem_peak_idx, damp->phi_at_bmag_tandem);
    gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
      phi,  damp->phi_at_bmag_max, damp->phi_at_bmag_tandem, damp->rate);
  }
  else {
    gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
      phi, damp->phi_at_bmag_max, damp->phi_at_bmag_max, damp->rate);
  }

  // Assemble the damping term -scale_prof * mask * f.
  gkyl_array_set(f_buffer, 1.0, fin);
  gkyl_array_scale_by_cell(damp->rate, damp->scale_prof);
  gkyl_array_scale_by_cell(f_buffer, damp->rate);

  // Add damping to f and the CFL frequency.
  gkyl_array_accumulate(rhs, -1.0, f_buffer);
  gkyl_array_accumulate(cflrate, 1.0, damp->rate);
}

void
gk_species_damping_advance_low_pass_filter(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  // Compute f - fbar, where fbar is stored as a cellwise-constant p0 field.
  gkyl_array_set(f_buffer, 1.0, fin); // f_buffer = fin
  gkyl_array_accumulate_offset(f_buffer, -1.0, fbar_in, 0); // f_buffer = fin - fbar
  gkyl_array_scale_by_cell(f_buffer, damp->rate); // f_buffer = rate * (fin - fbar)

  // Add damping term to RHS: df/dt -= rate * (f - fbar)
  // Add to the CFL frequency.
  gkyl_array_accumulate(rhs, -1.0, f_buffer);
  gkyl_array_accumulate(cflrate, 1.0, damp->rate);
}

void
gk_species_damping_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  damp->advance_func(app, gks, damp, phi, fin, fbar_in, f_buffer, rhs, cflrate);

  app->stat.species_damp_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_damping_write(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  gks->damping.write_func(app, gks, tm, frame);
}

// Low-pass filter fbar RHS helpers.

void
gk_species_damping_calc_fbar_rhs_disabled(const struct gk_damping *damp,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in, struct gkyl_array *rhs_fbar)
{
}

void
gk_species_damping_calc_fbar_rhs_enabled(const struct gk_damping *damp,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in, struct gkyl_array *rhs_fbar)
{
  // rhs_fbar = rate * (f - fbar)
  gkyl_array_set_offset(rhs_fbar, 1.0, fin, 0);
  gkyl_array_accumulate(rhs_fbar, -1.0, fbar_in);
  gkyl_array_scale_by_cell(rhs_fbar, damp->rate);
}

void
gk_species_damping_calc_fbar_rhs(const struct gk_damping *damp,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in, struct gkyl_array *rhs_fbar)
{
  damp->calc_fbar_rhs_func(damp, fin, fbar_in, rhs_fbar);
}

// Low-pass filter stage update helpers.

static void
gk_species_damping_forward_euler_disabled(struct gk_species *gks,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *fbar_out, double dt)
{
}

static void
gk_species_damping_forward_euler_enabled(struct gk_species *gks,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *fbar_out, double dt)
{
  gk_species_damping_calc_fbar_rhs(&gks->damping, fin, fbar_in, fbar_out);
  gk_species_step_f(gks, fbar_out, dt, fbar_in);
}

void
gk_species_damping_forward_euler(struct gk_species *gks,
  const struct gkyl_array *fin, const struct gkyl_array *fbar_in,
  struct gkyl_array *fbar_out, double dt)
{
  gks->damping.forward_euler_func(gks, fin, fbar_in, fbar_out, dt);
}

// Low-pass filter RK stage combine/copy helpers.

static void
gk_species_damping_combine_disabled(struct gk_species *gks, struct gkyl_array *fout,
  double c1, const struct gkyl_array *f1, double c2, const struct gkyl_array *f2,
  const struct gkyl_range *rng)
{
}

static void
gk_species_damping_combine_enabled(struct gk_species *gks, struct gkyl_array *fout,
  double c1, const struct gkyl_array *f1, double c2, const struct gkyl_array *f2,
  const struct gkyl_range *rng)
{
  gk_species_combine(gks, fout, c1, f1, c2, f2, rng);
}

void
gk_species_damping_combine(struct gk_species *gks, struct gkyl_array *fout, double c1,
  const struct gkyl_array *f1, double c2, const struct gkyl_array *f2,
  const struct gkyl_range *rng)
{
  gks->damping.combine_func(gks, fout, c1, f1, c2, f2, rng);
}

static void
gk_species_damping_copy_range_disabled(struct gk_species *gks, struct gkyl_array *fout,
  const struct gkyl_array *fin, const struct gkyl_range *range)
{
}

static void
gk_species_damping_copy_range_enabled(struct gk_species *gks, struct gkyl_array *fout,
  const struct gkyl_array *fin, const struct gkyl_range *range)
{
  gk_species_copy_range(gks, fout, fin, range);
}

void
gk_species_damping_copy_range(struct gk_species *gks, struct gkyl_array *fout,
  const struct gkyl_array *fin, const struct gkyl_range *range)
{
  gks->damping.copy_func(gks, fout, fin, range);
}

// Damping object lifecycle.

void
gk_species_damping_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_damping *damp)
{
  damp->type = gks->info.damping.type;
  damp->evolve = false; // Whether the rate is time dependent.
  damp->is_tandem = false; // Default to single mirror.
  damp->write_rate = gks->info.damping.write_rate;
  damp->write_fbar = gks->info.damping.write_fbar;
  const double rate_const = gks->info.damping.rate_const;

  int num_quad = gks->info.damping.num_quad? gks->info.damping.num_quad : 1; // Default is a p=0 mask.
  assert(num_quad == 1); // MF 2025/06/11: Limited to this for now.

  // Default function pointers.
  damp->write_func = gk_species_damping_write_disabled;
  damp->write_rate_func = gk_species_damping_write_rate_disabled;
  damp->write_fbar_func = gk_species_damping_write_fbar_disabled;
  damp->advance_func = gk_species_damping_advance_disabled;
  damp->set_fbar_to_f_func = gk_species_damping_set_fbar_to_f_disabled;
  damp->calc_fbar_rhs_func = gk_species_damping_calc_fbar_rhs_disabled;
  damp->forward_euler_func = gk_species_damping_forward_euler_disabled;
  damp->combine_func = gk_species_damping_combine_disabled;
  damp->copy_func = gk_species_damping_copy_range_disabled;

  damp->proj_on_basis_c2p_ctx.cdim = app->cdim;
  damp->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
  damp->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
  damp->proj_on_basis_c2p_ctx.pos_map = app->position_map;

  if (!damp->type)
    return;

  if (damp->write_rate)
    damp->write_rate_func = gk_species_damping_write_rate_enabled;

  // Damping rate is hardwired cellwise-constant (p0).
  damp->rate = mkarr(app->use_gpu, 1, gks->local_ext.volume);
  damp->rate_host = damp->rate;
  if (app->use_gpu)
    damp->rate_host = mkarr(false, damp->rate->ncomp, damp->rate->size);

  if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
    // USER_INPUT supports either a projected profile, or a uniform constant if no profile is provided.
    if (gks->info.damping.rate_profile) {
      gk_species_damping_project_phase_rate(app, gks, num_quad, gks->info.damping.rate_profile,
        gks->info.damping.rate_profile_ctx, damp->rate_host, damp->rate);
      gkyl_array_scale(damp->rate, rate_const == 0.0 ? 1.0 : rate_const);
    }
    else {
      gkyl_array_clear(damp->rate, rate_const);
    }

    damp->advance_func = gk_species_damping_advance_user_input;
  }
  else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
    damp->evolve = true; // Since the loss cone boundary is proportional to phi(t).

    enum gkyl_quad_type qtype = GKYL_GAUSS_LOBATTO_QUAD;

    // Create peak finder for bmag to find the mirror throat.
    // Search along the parallel (z) direction, which is the last configuration-space dimension.
    int search_dir = app->cdim - 1;
    struct gkyl_array_dg_find_peaks_inp peak_inp = {
      .basis = &app->basis,
      .grid = &app->grid,
      .range = &app->global,
      .range_ext = &app->global_ext,
      .search_dir = search_dir,
      .use_gpu = app->use_gpu,
    };

    struct gkyl_array *bmag_int_global = mkarr(app->use_gpu,
      app->gk_geom->geo_int.bmag->ncomp, app->global_ext.volume);
    damp->phi_smooth_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->geo_int.bmag,
      bmag_int_global);
    damp->bmag_peak_finder = gkyl_array_dg_find_peaks_new(&peak_inp, bmag_int_global);
    gkyl_array_dg_find_peaks_advance(damp->bmag_peak_finder, bmag_int_global);
    gkyl_array_release(bmag_int_global);

    int num_peaks = gkyl_array_dg_find_peaks_num_peaks(damp->bmag_peak_finder);
    damp->bmag_max_peak_idx = num_peaks - 2; // Edge is num_peaks-1, so maximum is one less.
    damp->bmag_max = gkyl_array_dg_find_peaks_acquire_vals(damp->bmag_peak_finder,
      damp->bmag_max_peak_idx);
    damp->bmag_max_z_coord = gkyl_array_dg_find_peaks_acquire_coords(damp->bmag_peak_finder,
      damp->bmag_max_peak_idx);
    damp->bmag_wall = gkyl_array_dg_find_peaks_acquire_vals(damp->bmag_peak_finder,
      num_peaks - 1);
    damp->bmag_wall_z_coord = gkyl_array_dg_find_peaks_acquire_coords(damp->bmag_peak_finder,
      num_peaks - 1);
    damp->bmag_max_basis = gkyl_array_dg_find_peaks_get_basis(damp->bmag_peak_finder);
    damp->bmag_max_range = gkyl_array_dg_find_peaks_get_range(damp->bmag_peak_finder);
    damp->bmag_max_range_ext = gkyl_array_dg_find_peaks_get_range_ext(damp->bmag_peak_finder);

    damp->phi_at_bmag_max = mkarr(app->use_gpu, damp->bmag_max_basis->num_basis,
      damp->bmag_max_range_ext->volume);
    damp->phi_at_bmag_tandem = mkarr(app->use_gpu, damp->bmag_max_basis->num_basis,
      damp->bmag_max_range_ext->volume);

    bool is_symmetric;
    int cdim = app->cdim;
    if (gkyl_compare_double(-app->grid.lower[cdim - 1], app->grid.upper[cdim - 1], 1e-12)) {
      is_symmetric = true;
    }
    else if (gkyl_compare_double(app->grid.lower[cdim - 1], 0.0, 1e-12)) {
      is_symmetric = false;
    }
    else {
      assert(false); // Needs either the lower bound at 0 or symmetric grid.
    }

    if ((is_symmetric && num_peaks == 5) || (!is_symmetric && num_peaks == 3)) {
      damp->is_tandem = false;
    }
    else if ((is_symmetric && num_peaks == 9) || (!is_symmetric && num_peaks == 5)) {
      damp->is_tandem = true;
    }
    else {
      assert(false); // Unsupported number of extrema for loss-cone damping.
    }

    if (damp->is_tandem)
      damp->bmag_tandem_peak_idx = num_peaks - 4;
    else
      damp->bmag_tandem_peak_idx = num_peaks - 2;

    damp->bmag_tandem = gkyl_array_dg_find_peaks_acquire_vals(damp->bmag_peak_finder,
      damp->bmag_tandem_peak_idx);
    damp->bmag_tandem_z_coord = gkyl_array_dg_find_peaks_acquire_coords(damp->bmag_peak_finder,
      damp->bmag_tandem_peak_idx);

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
      .qtype = qtype,
      .num_quad = num_quad,
      .cellwise_trap_loss = true,
      .use_gpu = app->use_gpu,
    };
    damp->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

    // Project conf-space scaling profile (if provided), otherwise use 1.
    struct gkyl_array *scale_prof_high_order = mkarr(app->use_gpu, 1, gks->local_ext.volume);
    struct gkyl_array *scale_prof_high_order_ho = app->use_gpu ?
      mkarr(false, scale_prof_high_order->ncomp, scale_prof_high_order->size) :
      gkyl_array_acquire(scale_prof_high_order);

    if (gks->info.damping.rate_profile) {
      struct gkyl_basis basis_mult;
      gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
        .grid = &gks->grid,
        .basis = &basis_mult,
        .num_quad = num_quad,
        .num_ret_vals = 1,
        .eval = gks->info.damping.rate_profile,
        .ctx = gks->info.damping.rate_profile_ctx,
        .c2p_func = proj_on_basis_c2p_position_func,
        .c2p_func_ctx = &damp->proj_on_basis_c2p_ctx,
      });
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, scale_prof_high_order_ho);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(scale_prof_high_order, scale_prof_high_order_ho);
    }
    else {
      gkyl_array_clear(scale_prof_high_order, 1.0);
    }

    gkyl_array_scale(scale_prof_high_order, rate_const == 0.0 ? 1.0 : rate_const);

    damp->scale_prof = mkarr(app->use_gpu, 1, gks->local_ext.volume);
    gkyl_array_set_offset(damp->scale_prof, pow(sqrt(2.0), gks->grid.ndim), scale_prof_high_order,
      0);

    gkyl_array_release(scale_prof_high_order_ho);
    gkyl_array_release(scale_prof_high_order);

    // Compute initial damping rate (phi at init time).
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->field->phi_smooth,
      damp->phi_smooth_global);
    gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, damp->phi_smooth_global,
      damp->bmag_max_peak_idx, damp->phi_at_bmag_max);

    if (damp->is_tandem) {
      gkyl_array_dg_find_peaks_project_on_peak_idx(damp->bmag_peak_finder, damp->phi_smooth_global,
        damp->bmag_tandem_peak_idx, damp->phi_at_bmag_tandem);
      gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
        damp->phi_smooth_global, damp->phi_at_bmag_max, damp->phi_at_bmag_tandem, damp->rate);
    }
    else {
      gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
        damp->phi_smooth_global, damp->phi_at_bmag_max, damp->phi_at_bmag_max, damp->rate);
    }
    gkyl_array_scale_by_cell(damp->rate, damp->scale_prof);

    damp->advance_func = gk_species_damping_advance_loss_cone;
  }
  else if (damp->type == GKYL_GK_DAMPING_LOW_PASS_FILTER) {
    // LOW_PASS_FILTER supports either a projected phase-space profile or a uniform constant.
    damp->evolve = true; // Since fbar must evolve in time.
    if (gks->info.damping.rate_profile) {
      gk_species_damping_project_phase_rate(app, gks, num_quad, gks->info.damping.rate_profile,
        gks->info.damping.rate_profile_ctx, damp->rate_host, damp->rate);
      gkyl_array_scale(damp->rate, rate_const == 0.0 ? 1.0 : rate_const);
    }
    else {
      gkyl_array_clear(damp->rate, rate_const);
    }

    // Allocate filtered distribution function array as p0.
    damp->fbar = mkarr(app->use_gpu, 1, gks->local_ext.volume);
    damp->fbar1 = mkarr(app->use_gpu, 1, gks->local_ext.volume);
    damp->fbarnew = mkarr(app->use_gpu, 1, gks->local_ext.volume);
    damp->fbar_host = 0;
    if (damp->write_fbar) {
      damp->write_fbar_func = gk_species_damping_write_fbar_enabled;
      damp->fbar_host = damp->fbar;
      if (app->use_gpu)
        damp->fbar_host = mkarr(false, damp->fbar->ncomp, damp->fbar->size);
    }

    // Initialize fbar from the projection of the initial distribution
    // (will be set from the initial f in the main app loop).
    gkyl_array_clear(damp->fbar, 0.0);

    damp->advance_func = gk_species_damping_advance_low_pass_filter;
    damp->set_fbar_to_f_func = gk_species_damping_set_fbar_to_f_enabled;
    damp->calc_fbar_rhs_func = gk_species_damping_calc_fbar_rhs_enabled;
    damp->forward_euler_func = gk_species_damping_forward_euler_enabled;
    damp->combine_func = gk_species_damping_combine_enabled;
    damp->copy_func = gk_species_damping_copy_range_enabled;
  }

  // Set function pointers chosen at runtime.
  if (damp->evolve) {
    damp->write_func = gk_species_damping_write_enabled;
  }
  else {
    damp->write_func = gk_species_damping_write_init_only;
  }
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

      gkyl_array_release(damp->phi_smooth_global);
      gkyl_array_dg_find_peaks_release(damp->bmag_peak_finder);
      gkyl_loss_cone_mask_gyrokinetic_release(damp->lcm_proj_op);
      gkyl_array_release(damp->scale_prof);
    }
    else if (damp->type == GKYL_GK_DAMPING_LOW_PASS_FILTER) {
      gkyl_array_release(damp->fbar);
      gkyl_array_release(damp->fbar1);
      gkyl_array_release(damp->fbarnew);
      if (app->use_gpu && damp->write_fbar)
        gkyl_array_release(damp->fbar_host);
    }
  }
}

// Damping runtime reset.

void
gk_species_damping_reset(gkyl_gyrokinetic_app *app, double tm, struct gk_species *gks,
  struct gk_damping *damp, struct gkyl_gyrokinetic_damping damp_inp)
{
  gk_species_damping_release(app, damp);

  gks->info.damping = damp_inp;
  gk_species_damping_init(app, gks, damp);
  gk_species_damping_set_fbar_to_f(gks, damp, gks->f);
}
