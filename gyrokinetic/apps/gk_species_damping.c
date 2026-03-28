#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_alloc.h>
#include <gkyl_dg_basis_ops.h>

// Damping state synchronization helpers.

void
gk_species_damping_set_fbar_to_f_enabled(const struct gk_species *gks, struct gk_damping *damp,
  const struct gkyl_array *f)
{
  gkyl_array_set(damp->fbar, 1.0, f);
  gkyl_array_set(damp->fbar1, 1.0, f);
  gkyl_array_set(damp->fbarnew, 1.0, f);
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
gk_species_damping_write_fbar_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
}

static void
gk_species_damping_write_fbar_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  double tm, int frame)
{
  // Metadata from app, species, and geometry (same pattern as species f write).
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  int io_meta_fbar_len[] = {app->io_meta_basic_len, gks->io_meta_len, app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem *io_meta_fbar[] = {app->io_meta_basic, gks->io_meta,
    app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt_fbar = gkyl_msgpack_create_union(
    sizeof(io_meta_fbar_len) / sizeof(int), io_meta_fbar_len, io_meta_fbar);

  const char *fmt_fbar = "%s-%s_fbar_%d.gkyl";
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
  // DG metadata for damping rate.
  struct gkyl_msgpack_map_elem mpe_drate[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
  };
  int mpe_drate_len = sizeof(mpe_drate) / sizeof(mpe_drate[0]);
  // Update app basic metada with time/frame.
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

  gks->damping.write_fbar_func(app, gks, tm, frame);

  gkyl_msgpack_data_release(mt);
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

  gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
    .grid = &gks->grid,
    .basis = &gks->basis,
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
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
}

void
gk_species_damping_advance_user_input(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
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
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  // Find the potential at the mirror throat.
  gkyl_dg_basis_ops_eval_array_at_coord_comp(phi, damp->bmag_max_coord,
    app->basis_on_dev, &app->grid, &app->local, damp->phi_m);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, damp->phi_m, damp->phi_m_global);

  // Project the loss cone mask.
  gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
    phi, damp->phi_m_global, damp->rate);

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
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  // Match fbar state to the RK stage of fin.
  const struct gkyl_array *fbar_in = damp->fbar;
  if (fin == gks->f1)
    fbar_in = damp->fbar1;
  else if (fin == gks->fnew)
    fbar_in = damp->fbarnew;

  // Compute f - fbar and scale by the damping rate: rate * (f - fbar)
  gkyl_array_set(f_buffer, 1.0, fin);      // f_buffer = f
  gkyl_array_accumulate(f_buffer, -1.0, fbar_in); // f_buffer = f - fbar
  gkyl_array_scale_by_cell(f_buffer, damp->rate);    // f_buffer = rate * (f - fbar)

  // Add damping term to RHS: df/dt -= rate * (f - fbar)
  // Add to the CFL frequency.
  gkyl_array_accumulate(rhs, -1.0, f_buffer);
  gkyl_array_accumulate(cflrate, 1.0, damp->rate);
}

void
gk_species_damping_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_damping *damp,
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  damp->advance_func(app, gks, damp, phi, fin, f_buffer, rhs, cflrate);

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
  gkyl_array_set(rhs_fbar, 1.0, fin);
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
  // Reuse species scratch storage to assemble fbar RHS for this FE substep.
  gk_species_damping_calc_fbar_rhs(&gks->damping, fin, fbar_in, gks->lte.f_lte);
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
  damp->write_fbar = gks->info.damping.write_fbar;
  const double rate_const = gks->info.damping.rate_const;

  int num_quad = gks->info.damping.num_quad? gks->info.damping.num_quad : 1; // Default is a p=0 mask.
  assert(num_quad == 1); // MF 2025/06/11: Limited to this for now.

  // Default function pointers.
  damp->write_func = gk_species_damping_write_disabled;
  damp->write_fbar_func = gk_species_damping_write_fbar_disabled;
  damp->advance_func = gk_species_damping_advance_disabled;
  damp->set_fbar_to_f_func = gk_species_damping_set_fbar_to_f_disabled;
  damp->calc_fbar_rhs_func = gk_species_damping_calc_fbar_rhs_disabled;
  damp->forward_euler_func = gk_species_damping_forward_euler_disabled;
  damp->combine_func = gk_species_damping_combine_disabled;
  damp->copy_func = gk_species_damping_copy_range_disabled;

  if (damp->type) {
    // Allocate rate array.
    damp->rate = mkarr(app->use_gpu, num_quad == 1? 1 : gks->basis.num_basis,
      gks->local_ext.volume);
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

      // Maximum bmag and its location.
      // NOTE: if the same max bmag occurs at multiple locations,
      // bmag_max_coord may have different values on different MPI processes.
      double bmag_max_coord_ho[GKYL_MAX_CDIM];
      double bmag_max_ho = gkyl_gk_geometry_reduce_arg_bmag(app->gk_geom, GKYL_MAX,
        bmag_max_coord_ho);
      double bmag_max_local = bmag_max_ho;
      double bmag_max_global;
      gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &bmag_max_local,
        &bmag_max_global);
      double bmag_max_coord_local[app->cdim], bmag_max_coord_global[app->cdim];
      if (fabs(bmag_max_ho - bmag_max_global) < 1e-16) {
        for (int d = 0; d < app->cdim; d++) {
          bmag_max_coord_local[d] = bmag_max_coord_ho[d];
        }
      }
      else {
        for (int d = 0; d < app->cdim; d++) {
          bmag_max_coord_local[d] = -DBL_MAX;
        }
      }
      gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, app->cdim, bmag_max_coord_local,
        bmag_max_coord_global);

      if (app->use_gpu) {
        damp->bmag_max = gkyl_cu_malloc(sizeof(double));
        damp->bmag_max_coord = gkyl_cu_malloc(app->cdim * sizeof(double));
        gkyl_cu_memcpy(damp->bmag_max, &bmag_max_global, sizeof(double), GKYL_CU_MEMCPY_H2D);
        gkyl_cu_memcpy(damp->bmag_max_coord, bmag_max_coord_ho, app->cdim * sizeof(double),
          GKYL_CU_MEMCPY_H2D);
      }
      else {
        damp->bmag_max = gkyl_malloc(sizeof(double));
        damp->bmag_max_coord = gkyl_malloc(app->cdim * sizeof(double));
        memcpy(damp->bmag_max, &bmag_max_global, sizeof(double));
        memcpy(damp->bmag_max_coord, bmag_max_coord_ho, app->cdim * sizeof(double));
      }

      // Electrostatic potential at bmag_max_coord.
      if (app->use_gpu) {
        damp->phi_m = gkyl_cu_malloc(sizeof(double));
        damp->phi_m_global = gkyl_cu_malloc(sizeof(double));
      }
      else {
        damp->phi_m = gkyl_malloc(sizeof(double));
        damp->phi_m_global = gkyl_malloc(sizeof(double));
      }

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
        .bmag_max_loc = damp->bmag_max_coord,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
        .num_quad = num_quad,
        .use_gpu = app->use_gpu,
      };
      damp->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

      // Project the conf-space rate profile provided.
      struct gkyl_array *scale_prof_high_order = mkarr(app->use_gpu, gks->basis.num_basis,
        gks->local_ext.volume);
      struct gkyl_array *scale_prof_high_order_ho = app->use_gpu? mkarr(false,
        scale_prof_high_order->ncomp, scale_prof_high_order->size)
                                                     : gkyl_array_acquire(scale_prof_high_order);

      if (gks->info.damping.rate_profile) {
        gkyl_proj_on_basis *projup = gkyl_proj_on_basis_new(&gks->grid, &gks->basis, num_quad, 1,
          gks->info.damping.rate_profile, gks->info.damping.rate_profile_ctx);
        gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, scale_prof_high_order_ho);
        gkyl_proj_on_basis_release(projup);
        gkyl_array_copy(scale_prof_high_order, scale_prof_high_order_ho);
      }
      else {
        gkyl_array_clear(scale_prof_high_order, 1.0);
      }

      gkyl_array_scale(scale_prof_high_order, rate_const == 0.0 ? 1.0 : rate_const);

      damp->scale_prof = mkarr(app->use_gpu, num_quad == 1? 1 : gks->basis.num_basis,
        gks->local_ext.volume);
      gkyl_array_set_offset(damp->scale_prof, pow(sqrt(2.0), gks->grid.ndim), scale_prof_high_order,
        0);

      gkyl_array_release(scale_prof_high_order_ho);
      gkyl_array_release(scale_prof_high_order);

      // Compute the initial damping rate (assuming phi=0 because phi hasn't been computed).
      // Find the potential at the mirror throat.
      gkyl_dg_basis_ops_eval_array_at_coord_comp(app->field->phi_smooth, damp->bmag_max_coord,
        app->basis_on_dev, &app->grid, &app->local, damp->phi_m);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, damp->phi_m, damp->phi_m_global);
      // Project the loss cone mask.
      gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
        app->field->phi_smooth, damp->phi_m_global, damp->rate);
      // Multiply by the user's scaling profile.
      gkyl_array_scale_by_cell(damp->rate, damp->scale_prof);

      damp->advance_func = gk_species_damping_advance_loss_cone;
    }
    else if (damp->type == GKYL_GK_DAMPING_LOW_PASS_FILTER) {
      // LOW_PASS_FILTER supports either a projected phase-space profile or a uniform constant.
      if (gks->info.damping.rate_profile) {
        gk_species_damping_project_phase_rate(app, gks, num_quad, gks->info.damping.rate_profile,
          gks->info.damping.rate_profile_ctx, damp->rate_host, damp->rate);
        gkyl_array_scale(damp->rate, rate_const == 0.0 ? 1.0 : rate_const);
      }
      else {
        gkyl_array_clear(damp->rate, rate_const);
      }

      // Allocate filtered distribution function array
      damp->fbar = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
      damp->fbar1 = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
      damp->fbarnew = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
      damp->fbar_host = 0;
      if (damp->write_fbar) {
        damp->write_fbar_func = gk_species_damping_write_fbar_enabled;
        damp->fbar_host = damp->fbar;
        if (app->use_gpu)
          damp->fbar_host = mkarr(false, damp->fbar->ncomp, damp->fbar->size);
      }

      // Initialize fbar from the projection of the initial distribution
      // (will be set from the initial f in the main app loop)
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
}

void
gk_species_damping_release(const struct gkyl_gyrokinetic_app *app, const struct gk_damping *damp)
{
  if (damp->type) {
    gkyl_array_release(damp->rate);
    if (app->use_gpu)
      gkyl_array_release(damp->rate_host);

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      // Nothing to release.
    }
    else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
      if (app->use_gpu) {
        gkyl_cu_free(damp->bmag_max);
        gkyl_cu_free(damp->bmag_max_coord);
        gkyl_cu_free(damp->phi_m);
        gkyl_cu_free(damp->phi_m_global);
      }
      else {
        gkyl_free(damp->bmag_max);
        gkyl_free(damp->bmag_max_coord);
        gkyl_free(damp->phi_m);
        gkyl_free(damp->phi_m_global);
      }
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
