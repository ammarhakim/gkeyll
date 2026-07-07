#include <gkyl_moment_priv.h>
#include <gkyl_wv_embed_geo.h>

// ---- Backend behavior variants (installed by moment_spacetime_init) ----
// Three sets of {update, max_dt, copy, calc_products, apply_bc}: a no-op set
// for apps with no spacetime at all (nothing allocated), a no-op/analytic
// set for a static background, and a dynamic-Einstein set.
// moment_spacetime_init picks one; the public functions at the bottom just
// dispatch through the stored pointers.

// Analytic (static) backend.
static struct gkyl_update_status spacetime_update_static(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr, double dt);
static double spacetime_max_dt_static(const gkyl_moment_app *app,
  const struct moment_spacetime *sp);
static void spacetime_copy_static(const struct moment_spacetime *sp,
  struct gkyl_array *dst, const struct gkyl_array *src);
static void spacetime_calc_products_static(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr);
static void spacetime_apply_bc_static(gkyl_moment_app *app, double tcurr,
  const struct moment_spacetime *sp, struct gkyl_array *f);

// Dynamic Einstein backend.
static struct gkyl_update_status spacetime_update_dynamic(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr, double dt);
static double spacetime_max_dt_dynamic(const gkyl_moment_app *app,
  const struct moment_spacetime *sp);
static void spacetime_copy_dynamic(const struct moment_spacetime *sp,
  struct gkyl_array *dst, const struct gkyl_array *src);
static void spacetime_calc_products_dynamic(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr);
static void spacetime_apply_bc_dynamic(gkyl_moment_app *app, double tcurr,
  const struct moment_spacetime *sp, struct gkyl_array *f);

// ICs, IO, and restart. Only the dynamic backend carries an evolving
// Einstein state to project/write/read; the static backend's IC work is the
// initial products fill, and its write/restart are no-ops (the analytic
// background is reconstructed from its callbacks).
static void spacetime_apply_ic_static(gkyl_moment_app *app,
  struct moment_spacetime *sp, double t0);
static void spacetime_apply_ic_dynamic(gkyl_moment_app *app,
  struct moment_spacetime *sp, double t0);
static void spacetime_write_dynamic(const gkyl_moment_app *app,
  const struct moment_spacetime *sp, double tm, int frame);
static struct gkyl_app_restart_status spacetime_read_from_frame_dynamic(
  gkyl_moment_app *app, struct moment_spacetime *sp, int frame);

// No-spacetime set: every method is a no-op; nothing is allocated, so an app
// without a spacetime pays nothing.
static void no_spacetime_calc_products(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr);
static void no_spacetime_apply_ic(gkyl_moment_app *app,
  struct moment_spacetime *sp, double t0);
static void spacetime_write_static(const gkyl_moment_app *app,
  const struct moment_spacetime *sp, double tm, int frame);
static struct gkyl_app_restart_status spacetime_read_from_frame_static(
  gkyl_moment_app *app, struct moment_spacetime *sp, int frame);
static void no_spacetime_release(const struct moment_spacetime *sp);

// Per-backend init/release: each installs its own function set (incl. release)
// and allocates/frees only the memory that backend uses, so neither the
// initializer nor the releaser has to re-test which backend is active.
static void spacetime_init_static(const struct gkyl_moment *mom,
  const struct gkyl_moment_spacetime *mom_st, struct gkyl_moment_app *app,
  struct moment_spacetime *sp);
static void spacetime_init_dynamic(const struct gkyl_moment *mom,
  const struct gkyl_moment_spacetime *mom_st, struct gkyl_moment_app *app,
  struct moment_spacetime *sp);
static void spacetime_release_static(const struct moment_spacetime *sp);
static void spacetime_release_dynamic(const struct moment_spacetime *sp);

// Initialize the no-spacetime object: zero the struct, install the no-op
// method set, allocate nothing (not even the products array).
static void
no_spacetime_init(const struct gkyl_moment *mom, struct gkyl_moment_app *app,
  struct moment_spacetime *sp)
{
  *sp = (struct moment_spacetime) {
    .ndim = mom->ndim,
    .is_static = true,
  };

  sp->update_func        = spacetime_update_static;   // no-op
  sp->max_dt_func        = spacetime_max_dt_static;   // DBL_MAX
  sp->copy_func          = spacetime_copy_static;     // no-op
  sp->calc_products_func = no_spacetime_calc_products;
  sp->apply_bc_func      = spacetime_apply_bc_static;
  sp->apply_ic_func      = no_spacetime_apply_ic;
  sp->write_func         = spacetime_write_static;
  sp->read_func    = spacetime_read_from_frame_static;
  sp->release_func       = no_spacetime_release;
}

// Initialize the spacetime component. Mirrors moment_field_init's shape but
// for the GR background. Allocates the shared spacetime-products array
// (consumed by the modular GR fluid equation objects via auxfields), then
// dispatches to the backend initializer: static-analytic, or dynamic
// Bona-Masso (einstein_eqn != NULL), which also allocates the Einstein-state
// arrays, wave-prop solvers, and BC handles.
void
moment_spacetime_init(const struct gkyl_moment *mom,
  const struct gkyl_moment_spacetime *mom_st,
  struct gkyl_moment_app *app, struct moment_spacetime *sp)
{
  // No spacetime configured: install the no-op set and allocate nothing.
  if (mom_st->analytic_spacetime == NULL && mom_st->einstein_eqn == NULL) {
    no_spacetime_init(mom, app, sp);
    return;
  }

  sp->ndim             = mom->ndim;
  sp->is_static        = mom_st->is_static;
  sp->has_tetrad       = mom_st->has_tetrad;
  sp->ctx              = mom_st->ctx;
  sp->init             = mom_st->init;
  sp->has_einstein_eqn = (mom_st->einstein_eqn != NULL);
  sp->spacetime_gauge  = mom_st->spacetime_gauge;
  sp->reinit_freq      = mom_st->reinit_freq;

  // Exactly one backend must be configured.
  assert((mom_st->analytic_spacetime != NULL) ^ (mom_st->einstein_eqn != NULL));

  // The shared products array, common to both backends. Mod fluid equations
  // point at this via gkyl_gr_euler_tetrad_set_auxfields. Zero-initialised here;
  // calc_products fills it at IC time and (dynamic) as needed each step.
  sp->prods_ncomp = mom_st->has_tetrad
    ? GKYL_GR_SP_NCOMP_TETRAD
    : GKYL_GR_SP_NCOMP_BASE;
  sp->prods = mkarr(false, sp->prods_ncomp, app->local_ext.volume);
  gkyl_array_clear(sp->prods, 0.0);

  // The per-interface tetrad cache is allocated lazily on the first
  // calc_products call. Zeroed here so equation objects fall back to per-call
  // averaging until the cache is wired in.
  sp->wave_spacetime = 0;

  // Dispatch to the backend initializer: each acquires/allocates only what its
  // backend needs and installs the full method set (incl. release_func), so the
  // initializer/releaser never re-test which backend is active.
  if (mom_st->einstein_eqn != NULL)
    spacetime_init_dynamic(mom, mom_st, app, sp);
  else
    spacetime_init_static(mom, mom_st, app, sp);
}

// Static-analytic backend: acquire the analytic-spacetime reference, leave
// the Einstein-state plumbing NULL, install the no-op/analytic method set.
static void
spacetime_init_static(const struct gkyl_moment *mom,
  const struct gkyl_moment_spacetime *mom_st, struct gkyl_moment_app *app,
  struct moment_spacetime *sp)
{
  // Acquire our own reference on the backend object: the creating wrapper
  // (Lua or a C driver) holds the only other reference and may release it
  // (or have it garbage-collected) while the app is still using the object.
  sp->analytic_spacetime = gkyl_gr_spacetime_acquire(mom_st->analytic_spacetime);
  sp->einstein_eqn = 0;

  // No evolving Einstein state for a fixed background.
  for (int d = 0; d < 3; d++) {
    sp->slvr[d]     = 0;
    sp->lower_bc[d] = 0;
    sp->upper_bc[d] = 0;
  }
  for (int d = 0; d < 4; d++) sp->f[d] = 0;
  sp->fdup       = 0;
  sp->fcurr      = 0;
  sp->bc_buffer  = 0;
  sp->embed_mask = 0;

  sp->update_func        = spacetime_update_static;
  sp->max_dt_func        = spacetime_max_dt_static;
  sp->copy_func          = spacetime_copy_static;
  sp->calc_products_func = spacetime_calc_products_static;
  sp->apply_bc_func      = spacetime_apply_bc_static;
  sp->apply_ic_func      = spacetime_apply_ic_static;
  sp->write_func         = spacetime_write_static;
  sp->read_func    = spacetime_read_from_frame_static;
  sp->release_func       = spacetime_release_static;
}

// Dynamic Bona-Masso backend: acquire the Einstein equation, allocate the
// Einstein-state arrays, wave-prop solvers, and BC handles, install the dynamic
// method set.
static void
spacetime_init_dynamic(const struct gkyl_moment *mom,
  const struct gkyl_moment_spacetime *mom_st, struct gkyl_moment_app *app,
  struct moment_spacetime *sp)
{
  sp->analytic_spacetime = 0;
  sp->einstein_eqn = gkyl_wv_eqn_acquire(mom_st->einstein_eqn);

  sp->update_func        = spacetime_update_dynamic;
  sp->max_dt_func        = spacetime_max_dt_dynamic;
  sp->copy_func          = spacetime_copy_dynamic;
  sp->calc_products_func = spacetime_calc_products_dynamic;
  sp->apply_bc_func      = spacetime_apply_bc_dynamic;
  sp->apply_ic_func      = spacetime_apply_ic_dynamic;
  sp->write_func         = spacetime_write_dynamic;
  sp->read_func    = spacetime_read_from_frame_dynamic;
  sp->release_func       = spacetime_release_dynamic;

  // BC handles are created only for non-periodic directions below; zero them
  // first so periodic directions are reliably null (the one case the releaser
  // must guard).
  for (int d = 0; d < 3; d++) { sp->lower_bc[d] = 0; sp->upper_bc[d] = 0; }

  int ndim    = mom->ndim;
  int ncomp   = sp->einstein_eqn->num_equations;
  enum gkyl_wave_limiter limiter =
    mom_st->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_st->limiter;

  for (int d = 0; d < ndim; d++) {
    sp->slvr[d] = gkyl_wave_prop_new(&(struct gkyl_wave_prop_inp) {
        .grid        = &app->grid,
        .equation    = sp->einstein_eqn,
        .split_type  = GKYL_WAVE_QWAVE,
        .limiter     = limiter,
        .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
        .update_dirs = { d },
        // Mirror mom_species: test each updated cell's admissibility via the
        // equation's check_inv_func and fall back to low-order flux where the
        // high-order update would leave the invariant domain. Without this
        // fallback, stiff collapses (e.g. conformal Gowdy) diverge.
        .check_inv_domain = true,
        .force_low_order_flux = mom_st->force_low_order_flux,
        .cfl  = app->cfl,
        .geom = app->geom,
        .comm = app->comm,
      });
  }
  sp->fdup = mkarr(false, ncomp, app->local_ext.volume);
  for (int d = 0; d < ndim + 1; d++)
    sp->f[d] = mkarr(false, ncomp, app->local_ext.volume);
  sp->fcurr = sp->f[0];

  // Embedding/excision mask the wave-prop reads (mirrors mom_species): 1.0
  // everywhere unless the Einstein equation carries an embed_geo.
  sp->embed_mask = mkarr(false, 1, app->local_ext.volume);
  gkyl_array_clear(sp->embed_mask, 1.0);
  if (sp->einstein_eqn->embed_geo)
    gkyl_wv_embed_geo_new_mask(sp->einstein_eqn->embed_geo, &app->grid,
      &app->local, sp->embed_mask);

  // Boundary conditions on the Einstein state array — mirrors the BC setup
  // in moment_field_init. Non-periodic directions only.
  int num_periodic_dir = app->num_periodic_dir;
  int is_np[3] = { 1, 1, 1 };
  for (int d = 0; d < num_periodic_dir; d++)
    is_np[app->periodic_dirs[d]] = 0;

  int nghost[3] = { 2, 2, 2 };
  for (int dir = 0; dir < ndim; dir++) {
    if (!is_np[dir]) continue;

    const enum gkyl_field_bc_type *bc;
    wv_bc_func_t bc_lower_func, bc_upper_func;
    if (dir == 0)      { bc = mom_st->bcx; bc_lower_func = mom_st->bcx_func[0]; bc_upper_func = mom_st->bcx_func[1]; }
    else if (dir == 1) { bc = mom_st->bcy; bc_lower_func = mom_st->bcy_func[0]; bc_upper_func = mom_st->bcy_func[1]; }
    else               { bc = mom_st->bcz; bc_lower_func = mom_st->bcz_func[0]; bc_upper_func = mom_st->bcz_func[1]; }

    sp->lower_bct[dir] = bc[0];
    sp->upper_bct[dir] = bc[1];

    switch (bc[0]) {
      case GKYL_FIELD_FUNC:
        sp->lower_bc[dir] = gkyl_wv_apply_bc_new(&app->grid, sp->einstein_eqn,
          app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_lower_func, mom_st->ctx);
        break;
      case GKYL_FIELD_COPY:
      case GKYL_FIELD_WEDGE:
        sp->lower_bc[dir] = gkyl_wv_apply_bc_new(&app->grid, sp->einstein_eqn,
          app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
        break;
      case GKYL_FIELD_SKIP:
        sp->lower_bc[dir] = gkyl_wv_apply_bc_new(&app->grid, sp->einstein_eqn,
          app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_skip, 0);
        break;
      default:
        assert(false);
        break;
    }
    switch (bc[1]) {
      case GKYL_FIELD_FUNC:
        sp->upper_bc[dir] = gkyl_wv_apply_bc_new(&app->grid, sp->einstein_eqn,
          app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_upper_func, mom_st->ctx);
        break;
      case GKYL_FIELD_COPY:
      case GKYL_FIELD_WEDGE:
        sp->upper_bc[dir] = gkyl_wv_apply_bc_new(&app->grid, sp->einstein_eqn,
          app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
        break;
      case GKYL_FIELD_SKIP:
        sp->upper_bc[dir] = gkyl_wv_apply_bc_new(&app->grid, sp->einstein_eqn,
          app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_skip, 0);
        break;
      default:
        assert(false);
        break;
    }
  }

  // BC sync buffer.
  long buff_sz = 0;
  for (int d = 0; d < ndim; d++) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  sp->bc_buffer = mkarr(false, ncomp, buff_sz);
}

// ---- apply_bc ----------------------------------------------------------------

// No Einstein state to evolve (static-analytic backend or no spacetime).
static void
spacetime_apply_bc_static(gkyl_moment_app *app, double tcurr,
  const struct moment_spacetime *sp, struct gkyl_array *f) { }

static void
spacetime_apply_bc_dynamic(gkyl_moment_app *app, double tcurr,
  const struct moment_spacetime *sp, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir;
  int ndim = app->ndim;
  int is_non_periodic[3] = { 1, 1, 1 };
  for (int d = 0; d < num_periodic_dir; d++)
    is_non_periodic[app->periodic_dirs[d]] = 0;

  for (int d = 0; d < ndim; d++) {
    if (!is_non_periodic[d]) continue;
    if (sp->lower_bct[d] != GKYL_FIELD_WEDGE)
      gkyl_wv_apply_bc_advance(sp->lower_bc[d], tcurr, &app->local, f);
    if (sp->upper_bct[d] != GKYL_FIELD_WEDGE)
      gkyl_wv_apply_bc_advance(sp->upper_bc[d], tcurr, &app->local, f);

    if (sp->lower_bct[d] == GKYL_FIELD_WEDGE)
      moment_apply_wedge_bc(app, tcurr, &app->local, sp->bc_buffer,
        d, sp->lower_bc[d], sp->upper_bc[d], f);
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, app->periodic_dirs, f);
}

void
moment_spacetime_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_spacetime *sp, struct gkyl_array *f)
{
  sp->apply_bc_func(app, tcurr, sp, f);
}

// ---- max_dt ----------------------------------------------------------------

// Analytic backend imposes no CFL constraint of its own.
static double
spacetime_max_dt_static(const gkyl_moment_app *app, const struct moment_spacetime *sp)
{
  return DBL_MAX;
}

// Einstein backend: CFL of the hyperbolic Einstein wave-prop.
static double
spacetime_max_dt_dynamic(const gkyl_moment_app *app, const struct moment_spacetime *sp)
{
  double max_dt = DBL_MAX;
  for (int d = 0; d < app->ndim; d++)
    max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(sp->slvr[d], &app->local, sp->f[0]));
  return max_dt;
}

double
moment_spacetime_max_dt(const gkyl_moment_app *app, const struct moment_spacetime *sp)
{
  return sp->max_dt_func(app, sp);
}

// ---- update ----------------------------------------------------------------

// Analytic backend: nothing to evolve.
static struct gkyl_update_status
spacetime_update_static(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr, double dt)
{
  // Nothing to evolve for a fixed background.
  return (struct gkyl_update_status) { .success = true, .dt_suggested = DBL_MAX };
}

// Einstein backend: hyperbolic wave-prop of the Einstein state + BC each
// stage, then refresh the derived geometry the fluid step consumes.
static struct gkyl_update_status
spacetime_update_dynamic(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr, double dt)
{
  struct gkyl_wave_prop_status stat = { true, DBL_MAX };
  int ndim = sp->ndim;
  for (int d = 0; d < ndim; d++) {
    stat = gkyl_wave_prop_advance(sp->slvr[d], tcurr, dt, &app->local,
      sp->embed_mask, sp->f[d], sp->f[d + 1]);
    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false, .dt_suggested = stat.dt_suggested,
      };
    moment_spacetime_apply_bc(app, tcurr, sp, sp->f[d + 1]);
  }

  // Refresh products + tetrad cache from the advanced state so the upcoming
  // fluid step reads the new geometry.
  spacetime_calc_products_dynamic(app, sp, tcurr);

  return (struct gkyl_update_status) {
    .success = true, .dt_suggested = stat.dt_suggested,
  };
}

struct gkyl_update_status
moment_spacetime_update(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr, double dt)
{
  return sp->update_func(app, sp, tcurr, dt);
}

// ---- copy (time-stepper backup/commit/restore) -----------------------------

// Analytic backend has no evolving state; the passed-in arrays are NULL.
static void
spacetime_copy_static(const struct moment_spacetime *sp,
  struct gkyl_array *dst, const struct gkyl_array *src)
{
  // No evolving state; the passed-in arrays are NULL.
}

static void
spacetime_copy_dynamic(const struct moment_spacetime *sp,
  struct gkyl_array *dst, const struct gkyl_array *src)
{
  gkyl_array_copy(dst, src);
}

void
moment_spacetime_copy(const struct moment_spacetime *sp,
  struct gkyl_array *dst, const struct gkyl_array *src)
{
  sp->copy_func(sp, dst, src);
}

// ---- stepper protocol (backup/restore/commit) --------------------------------
// The spacetime always uses the wave-prop state layout (f[0] -> f[ndim]);
// the static/no-spacetime copy method is a no-op, so these are safe to call
// unconditionally.

void
moment_spacetime_step_backup(const struct moment_spacetime *sp)
{
  sp->copy_func(sp, sp->fdup, sp->f[0]);
}

void
moment_spacetime_step_restore(const struct moment_spacetime *sp)
{
  sp->copy_func(sp, sp->f[0], sp->fdup);
}

void
moment_spacetime_step_commit(const struct moment_spacetime *sp)
{
  sp->copy_func(sp, sp->f[0], sp->f[sp->ndim]);
}

// ---- calc_products ---------------------------------------------------------
// Recompute the derived spacetime quantities the fluid solver consumes:
// cell-center products (over interior + ghost cells) and the per-interface
// tetrad cache. The two are inseparable — refreshing one without the other
// would leave the cache stale — so they live in one call.

// Fill the products array from the analytic background — closed-form, per cell,
// over the interior + ghost cells (wave_prop reads spacetime at boundary-ghost
// interfaces too). This is the static-backend derivation; it owns the geometry
// here rather than going through the spacetime-coupling object.
static void
spacetime_fill_products_analytic(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&app->local_ext, iter.idx);
    double *prods_row = gkyl_array_fetch(sp->prods, cidx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    double x = app->ndim >= 1 ? xc[0] : 0.0;
    double y = app->ndim >= 2 ? xc[1] : 0.0;
    double z = app->ndim >= 3 ? xc[2] : 0.0;

    gkyl_moment_spacetime_coupling_fill_products_analytic(sp->analytic_spacetime,
      tcurr, x, y, z, prods_row);
  }
}

// Post-process filled products (shared by both backends): mirror-copy skin →
// ghost on each non-periodic boundary so the ghost spacetime equals the
// interior spacetime (matching the BC semantics of the self-contained
// wv_gr_euler formulation), sync periodic directions, then refresh the
// per-interface tetrad cache (created once at app construction by
// moment_spacetime_create_tetrad_cache).
static void
spacetime_finish_products(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr)
{
  int num_periodic_dir = app->num_periodic_dir;
  int is_non_periodic[3] = { 1, 1, 1 };
  for (int d = 0; d < num_periodic_dir; d++)
    is_non_periodic[app->periodic_dirs[d]] = 0;
  long ncomp = sp->prods_ncomp;
  for (int d = 0; d < app->ndim; d++) {
    if (!is_non_periodic[d]) continue;
    int lo = app->local.lower[d], up = app->local.upper[d];
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->skin_ghost.lower_skin[d]);
    while (gkyl_range_iter_next(&iter)) {
      int skin_idx[GKYL_MAX_DIM], ghost_idx[GKYL_MAX_DIM];
      gkyl_copy_int_arr(app->ndim, iter.idx, skin_idx);
      gkyl_copy_int_arr(app->ndim, iter.idx, ghost_idx);
      ghost_idx[d] = 2 * lo - skin_idx[d] - 1;  // mirror across lo - 0.5
      long sloc = gkyl_range_idx(&app->local_ext, skin_idx);
      long gloc = gkyl_range_idx(&app->local_ext, ghost_idx);
      memcpy(gkyl_array_fetch(sp->prods, gloc),
             gkyl_array_cfetch(sp->prods, sloc), ncomp * sizeof(double));
    }
    gkyl_range_iter_init(&iter, &app->skin_ghost.upper_skin[d]);
    while (gkyl_range_iter_next(&iter)) {
      int skin_idx[GKYL_MAX_DIM], ghost_idx[GKYL_MAX_DIM];
      gkyl_copy_int_arr(app->ndim, iter.idx, skin_idx);
      gkyl_copy_int_arr(app->ndim, iter.idx, ghost_idx);
      ghost_idx[d] = 2 * up - skin_idx[d] + 1;  // mirror across up + 0.5
      long sloc = gkyl_range_idx(&app->local_ext, skin_idx);
      long gloc = gkyl_range_idx(&app->local_ext, ghost_idx);
      memcpy(gkyl_array_fetch(sp->prods, gloc),
             gkyl_array_cfetch(sp->prods, sloc), ncomp * sizeof(double));
    }
  }
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, app->periodic_dirs, sp->prods);

  if (sp->wave_spacetime)
    gkyl_wave_spacetime_refresh(sp->wave_spacetime, &app->grid, app->geom,
      sp->analytic_spacetime, sp->prods, tcurr);
}

// Create the per-interface tetrad cache if any species consumes it. Called
// once from gkyl_moment_app_new, after the products array has been filled
// for the first time (the cache is built from it). The tetrad equation
// objects receive the cache pointer through moment_coupling_init's auxfield
// wiring -- all cross-object pointers are handed out at construction --
// and calc_products refreshes the cache contents whenever the products
// change.
void
moment_spacetime_create_tetrad_cache(gkyl_moment_app *app,
  struct moment_spacetime *sp, double tcurr)
{
  if (sp->wave_spacetime) return;

  bool needed = false;
  for (int i = 0; i < app->num_species; i++)
    if (app->species[i].eqn_type == GKYL_EQN_GR_EULER_TETRAD) {
      needed = true; break;
    }
  if (!needed) return;

  sp->wave_spacetime = gkyl_wave_spacetime_new(&app->grid, &app->local_ext,
    app->geom, sp->analytic_spacetime, sp->prods, tcurr, /*use_gpu=*/false);
}

static void
no_spacetime_calc_products(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr) { }

static void
spacetime_calc_products_static(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr)
{
  spacetime_fill_products_analytic(app, sp, tcurr);
  spacetime_finish_products(app, sp, tcurr);
}

static void
spacetime_calc_products_dynamic(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr)
{
  // TODO: derive the products from the evolved Einstein state (sp->fcurr);
  // for now only the shared post-processing (BC mirror + tetrad cache) runs.
  spacetime_finish_products(app, sp, tcurr);
}

void
moment_spacetime_calc_products(gkyl_moment_app *app, struct moment_spacetime *sp,
  double tcurr)
{
  sp->calc_products_func(app, sp, tcurr);
}

// ---- initial conditions ------------------------------------------------------

static void
no_spacetime_apply_ic(gkyl_moment_app *app, struct moment_spacetime *sp,
  double t0) { }

// Static-analytic backend: no evolving state to project; the IC work is
// filling the derived spacetime quantities (products + tetrad cache) the
// fluid solver consumes.
static void
spacetime_apply_ic_static(gkyl_moment_app *app, struct moment_spacetime *sp,
  double t0)
{
  moment_spacetime_calc_products(app, sp, t0);
}

// Dynamic backend: project the user-supplied IC into the Einstein state
// array (analogous to the field IC projection), apply BCs, then fill the
// derived spacetime quantities. The products are also refreshed inside
// moment_spacetime_update each step.
static void
spacetime_apply_ic_dynamic(gkyl_moment_app *app, struct moment_spacetime *sp,
  double t0)
{
  if (sp->init) {
    int num_quad = app->scheme_type == GKYL_MOMENT_MP ? 4 : 2;
    int ncomp = sp->einstein_eqn->num_equations;
    gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, num_quad, ncomp,
      sp->init, sp->ctx);
    gkyl_fv_proj_advance(proj, t0, &app->local, sp->fcurr);
    gkyl_fv_proj_release(proj);
    moment_spacetime_apply_bc(app, t0, sp, sp->fcurr);
  }

  moment_spacetime_calc_products(app, sp, t0);
}

void
moment_spacetime_apply_ic(gkyl_moment_app *app, struct moment_spacetime *sp,
  double t0)
{
  sp->apply_ic_func(app, sp, t0);
}

// ---- write (frame IO) ----------------------------------------------------------

static void
spacetime_write_static(const gkyl_moment_app *app, const struct moment_spacetime *sp,
  double tm, int frame) { }

static void
spacetime_write_dynamic(const gkyl_moment_app *app, const struct moment_spacetime *sp,
  double tm, int frame)
{
  struct gkyl_msgpack_data *mt = moment_array_meta_new( (struct moment_output_meta) {
      .frame = frame,
      .stime = tm
    }
  );

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "spacetime", frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, sp->fcurr, fileNm.str);
  cstr_drop(&fileNm);

  moment_array_meta_release(mt);
}

void
moment_spacetime_write(const gkyl_moment_app *app, const struct moment_spacetime *sp,
  double tm, int frame)
{
  sp->write_func(app, sp, tm, frame);
}

// ---- restart read ----------------------------------------------------------------

// No evolving Einstein state to read back. For the static-analytic backend
// the products are filled at construction and by apply_ic; the analytic
// background needs nothing from the restart files.
static struct gkyl_app_restart_status
spacetime_read_from_frame_static(gkyl_moment_app *app, struct moment_spacetime *sp,
  int frame)
{
  return (struct gkyl_app_restart_status) {
    .io_status = GKYL_ARRAY_RIO_SUCCESS,
    .frame = 0,
    .stime = 0.0
  };
}

// Dynamic backend: read the Einstein state written by spacetime_write_dynamic,
// apply BCs, and rebuild the derived geometry (products + tetrad cache) at
// the restart time so the first fluid step reads a consistent spacetime.
static struct gkyl_app_restart_status
spacetime_read_from_frame_dynamic(gkyl_moment_app *app, struct moment_spacetime *sp,
  int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "spacetime", frame);
  struct gkyl_app_restart_status rstat = moment_app_header_from_file(app, fileNm.str);

  if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, sp->fcurr, fileNm.str);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      moment_spacetime_apply_bc(app, rstat.stime, sp, sp->fcurr);
      moment_spacetime_calc_products(app, sp, rstat.stime);
    }
  }
  cstr_drop(&fileNm);

  return rstat;
}

struct gkyl_app_restart_status
moment_spacetime_read_from_frame(gkyl_moment_app *app, struct moment_spacetime *sp,
  int frame)
{
  return sp->read_func(app, sp, frame);
}

// ---- release ---------------------------------------------------------------

// No-spacetime set: nothing was allocated.
static void
no_spacetime_release(const struct moment_spacetime *sp) { }

// Static backend: release the acquired analytic-spacetime reference.
static void
spacetime_release_static(const struct moment_spacetime *sp)
{
  gkyl_gr_spacetime_release(sp->analytic_spacetime);
}

// Dynamic backend: release the Einstein equation reference and the wave-prop /
// state-array / BC plumbing this backend allocates. Everything here is
// unconditionally created by spacetime_init_dynamic, except the per-direction
// BC handles, which exist only for non-periodic directions.
static void
spacetime_release_dynamic(const struct moment_spacetime *sp)
{
  gkyl_wv_eqn_release(sp->einstein_eqn);

  for (int d = 0; d < sp->ndim; d++) {
    if (sp->lower_bc[d]) gkyl_wv_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d]) gkyl_wv_apply_bc_release(sp->upper_bc[d]);
    gkyl_wave_prop_release(sp->slvr[d]);
    gkyl_array_release(sp->f[d]);
  }
  gkyl_array_release(sp->f[sp->ndim]);
  gkyl_array_release(sp->fdup);
  gkyl_array_release(sp->bc_buffer);
  gkyl_array_release(sp->embed_mask);
}

// Release everything owned by the spacetime component: the shared products /
// tetrad cache (common to both backends), then the backend-specific memory.
void
moment_spacetime_release(const struct moment_spacetime *sp)
{
  if (sp->prods) gkyl_array_release(sp->prods);
  if (sp->wave_spacetime) gkyl_wave_spacetime_release(sp->wave_spacetime);

  sp->release_func(sp);
}
