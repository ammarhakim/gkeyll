#include <gkyl_moment_priv.h>

// Initialize the spacetime component. Mirrors moment_field_init's shape but
// for the GR-mod background. Allocates the shared spacetime-products array
// (consumed by the modular GR fluid equation objects via auxfields). For
// the dynamic Bona-Masso backend (einstein_eqn != NULL), also allocates the
// Einstein-state arrays, wave-prop solvers, and BC handles — Phase A only
// exercises the static-analytic path, but the dynamic plumbing is wired in
// the same place so the Phase B transition is a local change.
void
moment_spacetime_init(const struct gkyl_moment *mom,
  const struct gkyl_moment_spacetime *mom_st,
  struct gkyl_moment_app *app, struct moment_spacetime *sp)
{
  sp->ndim       = mom->ndim;
  sp->is_static  = mom_st->is_static;
  sp->has_tetrad = mom_st->has_tetrad;
  sp->ctx        = mom_st->ctx;
  sp->init       = mom_st->init;

  sp->analytic_spacetime = mom_st->analytic_spacetime;
  sp->einstein_eqn       = mom_st->einstein_eqn;
  sp->has_einstein_eqn   = (mom_st->einstein_eqn != NULL);
  sp->spacetime_gauge    = mom_st->spacetime_gauge;
  sp->reinit_freq        = mom_st->reinit_freq;

  // Exactly one backend must be configured.
  assert((mom_st->analytic_spacetime != NULL) ^ (mom_st->einstein_eqn != NULL));

  sp->prods_ncomp = mom_st->has_tetrad
    ? GKYL_GR_SP_NCOMP_TETRAD
    : GKYL_GR_SP_NCOMP_BASE;

  // The shared products array. Mod fluid equations point at this via
  // gkyl_gr_euler_mod_set_auxfields after construction. Zero-initialised
  // here; the coupling object fills it at IC time and as needed each step.
  sp->prods = mkarr(false, sp->prods_ncomp, app->local_ext.volume);
  gkyl_array_clear(sp->prods, 0.0);

  // Dynamic-case plumbing: Einstein-state arrays + wave-prop solvers + BCs.
  // Skipped entirely when running with the static-analytic backend.
  for (int d = 0; d < 3; d++) sp->slvr[d] = NULL;
  for (int d = 0; d < 4; d++) sp->f[d] = NULL;
  sp->fdup  = NULL;
  sp->fcurr = NULL;
  sp->bc_buffer = NULL;
  for (int d = 0; d < 3; d++) {
    sp->lower_bc[d] = NULL;
    sp->upper_bc[d] = NULL;
  }

  if (!sp->has_einstein_eqn)
    return;

  // ---- Dynamic Bona-Masso backend (Phase B) ----
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
        .check_inv_domain = false,
        .cfl  = app->cfl,
        .geom = app->geom,
        .comm = app->comm,
      });
  }
  sp->fdup = mkarr(false, ncomp, app->local_ext.volume);
  for (int d = 0; d < ndim + 1; d++)
    sp->f[d] = mkarr(false, ncomp, app->local_ext.volume);
  sp->fcurr = sp->f[0];

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

// Apply BCs to the Einstein state array. No-op in the static-analytic case
// (there is no Einstein state to evolve).
void
moment_spacetime_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_spacetime *sp, struct gkyl_array *f)
{
  if (!sp->has_einstein_eqn) return;

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

// Maximum stable dt for the Einstein hyperbolic step. Returns DBL_MAX in the
// static case (the spacetime imposes no CFL constraint of its own).
double
moment_spacetime_max_dt(const gkyl_moment_app *app,
  const struct moment_spacetime *sp)
{
  if (!sp->has_einstein_eqn || sp->is_static) return DBL_MAX;

  double max_dt = DBL_MAX;
  for (int d = 0; d < app->ndim; d++)
    max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(sp->slvr[d], &app->local, sp->f[0]));
  return max_dt;
}

// Advance the Einstein state by dt. No-op when is_static or when running
// with the analytic backend (no Einstein state to evolve).
struct gkyl_update_status
moment_spacetime_update(gkyl_moment_app *app,
  const struct moment_spacetime *sp, double tcurr, double dt)
{
  struct gkyl_wave_prop_status stat = { true, DBL_MAX };

  if (sp->has_einstein_eqn && !sp->is_static) {
    int ndim = sp->ndim;
    for (int d = 0; d < ndim; d++) {
      stat = gkyl_wave_prop_advance(sp->slvr[d], tcurr, dt, &app->local,
        NULL, sp->f[d], sp->f[d + 1]);
      if (!stat.success)
        return (struct gkyl_update_status) {
          .success = false, .dt_suggested = stat.dt_suggested,
        };
      moment_spacetime_apply_bc(app, tcurr, sp, sp->f[d + 1]);
    }
  }

  return (struct gkyl_update_status) {
    .success = true, .dt_suggested = stat.dt_suggested,
  };
}

// Release everything owned by the spacetime component.
void
moment_spacetime_release(const struct moment_spacetime *sp)
{
  if (sp->prods) gkyl_array_release(sp->prods);

  if (!sp->has_einstein_eqn) return;

  for (int d = 0; d < sp->ndim; d++) {
    if (sp->lower_bc[d]) gkyl_wv_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d]) gkyl_wv_apply_bc_release(sp->upper_bc[d]);
  }
  for (int d = 0; d < sp->ndim; d++)
    if (sp->slvr[d]) gkyl_wave_prop_release(sp->slvr[d]);
  if (sp->fdup) gkyl_array_release(sp->fdup);
  for (int d = 0; d < sp->ndim + 1; d++)
    if (sp->f[d]) gkyl_array_release(sp->f[d]);
  if (sp->bc_buffer) gkyl_array_release(sp->bc_buffer);
}
