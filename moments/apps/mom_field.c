#include <gkyl_moment_priv.h>
#include <gkyl_wv_embed_geo.h>
#include <gkyl_wv_maxwell.h>

// ---- Scheme/behavior variants (installed by moment_field_init) ------------
// All scheme/static/presence dispatch is wired once, at init: the
// moment_field_* wrappers call through the fld->*_func pointers and never
// branch on scheme_type / is_static / field presence themselves. When no
// field is configured (mom_fld->init == 0), the no-op "none" set is
// installed and NOTHING is allocated: a fluid-only run pays no Maxwell
// state, solvers, or equation object.

// Wave-propagation scheme (stepped by the one-step Strang stepper).
static double field_max_dt_wave_prop(const gkyl_moment_app *app,
  const struct moment_field *fld);
static struct gkyl_update_status field_update_wave_prop(gkyl_moment_app *app,
  const struct moment_field *fld, double tcurr, double dt);
static void field_release_wave_prop(const struct moment_field *fld);

// MP scheme (RHS-based; also used for KEP runs -- there is no KEP Maxwell
// solver).
static double field_max_dt_mp(const gkyl_moment_app *app,
  const struct moment_field *fld);
static double field_rhs_mp(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs);
static void field_release_mp(const struct moment_field *fld);

// Static field: pass the state through the stepper's staging arrays
// unchanged / zero RHS (the sources still read the field; only its
// evolution is frozen).
static struct gkyl_update_status field_update_static(gkyl_moment_app *app,
  const struct moment_field *fld, double tcurr, double dt);
static double field_rhs_static(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

// MP full-step update: a complete SSP-RK3 step built on rhs_func.
static struct gkyl_update_status field_update_ssp_rk3(gkyl_moment_app *app,
  const struct moment_field *fld, double tcurr, double dt);

// Stepper protocol: Strang stage-state accessor and end-of-step commit, per
// state-array layout (backup/restore work through fcurr and need no
// per-scheme variants; without a field they degrade to no-ops through
// copy_func and a NULL fcurr).
static struct gkyl_array* field_stage_state_wave_prop(const struct moment_field *fld, int nstrang);
static void field_step_commit_wave_prop(const struct moment_field *fld);
static struct gkyl_array* field_stage_state_curr(const struct moment_field *fld, int nstrang);
static void field_step_commit_curr(const struct moment_field *fld);

// BCs and stepper state copy.
static void field_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_field *field, struct gkyl_array *f);
static void field_copy(const struct moment_field *fld,
  struct gkyl_array *dst, const struct gkyl_array *src);

// ICs, IO, energy diagnostics, and restart (one real variant each; the
// no-field set below no-ops them).
static void field_apply_ic(gkyl_moment_app *app, struct moment_field *fld,
  double t0);
static void field_write(const gkyl_moment_app *app,
  const struct moment_field *fld, double tm, int frame);
static void field_calc_energy(gkyl_moment_app *app,
  struct moment_field *fld, double tm);
static void field_write_energy(gkyl_moment_app *app,
  struct moment_field *fld);
static struct gkyl_app_restart_status field_from_file(gkyl_moment_app *app,
  struct moment_field *fld, const char *fname);

// No-field set: every method is a no-op; max_dt/update/rhs impose no
// constraint.
static double no_field_max_dt(const gkyl_moment_app *app,
  const struct moment_field *fld);
static struct gkyl_update_status no_field_update(gkyl_moment_app *app,
  const struct moment_field *fld, double tcurr, double dt);
static double no_field_rhs(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs);
static void no_field_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_field *field, struct gkyl_array *f);
static void no_field_copy(const struct moment_field *fld,
  struct gkyl_array *dst, const struct gkyl_array *src);
static void no_field_apply_ic(gkyl_moment_app *app, struct moment_field *fld,
  double t0);
static void no_field_write(const gkyl_moment_app *app,
  const struct moment_field *fld, double tm, int frame);
static void no_field_calc_energy(gkyl_moment_app *app,
  struct moment_field *fld, double tm);
static void no_field_write_energy(gkyl_moment_app *app,
  struct moment_field *fld);
static struct gkyl_app_restart_status no_field_read(gkyl_moment_app *app,
  struct moment_field *fld, const char *fname);
static void no_field_release(const struct moment_field *fld);

// Initialize the no-field object: zero the struct (so the coupling reads
// benign epsilon0/mu0/flags), install the no-op method set, allocate
// nothing.
static void
no_field_init(const struct gkyl_moment *mom, struct gkyl_moment_app *app,
  struct moment_field *fld)
{
  *fld = (struct moment_field) {
    .ndim = mom->ndim,
  };

  fld->max_dt_func = no_field_max_dt;
  fld->update_func = no_field_update;
  fld->rhs_func = no_field_rhs;
  fld->apply_bc_func = no_field_apply_bc;
  fld->copy_func = no_field_copy;
  fld->stage_state_func = field_stage_state_curr; // returns the NULL fcurr
  fld->step_commit_func = field_step_commit_curr; // no-op
  fld->apply_ic_func = no_field_apply_ic;
  fld->write_func = no_field_write;
  fld->calc_energy_func = no_field_calc_energy;
  fld->write_energy_func = no_field_write_energy;
  fld->read_func = no_field_read;
  fld->release_func = no_field_release;
}

// initialize field
void
moment_field_init(const struct gkyl_moment *mom, const struct gkyl_moment_field *mom_fld,
  struct gkyl_moment_app *app, struct moment_field *fld)
{
  if (!mom_fld->init) {
    no_field_init(mom, app, fld);
    return;
  }

  fld->ndim = mom->ndim;
  double epsilon0 = fld->epsilon0 = mom_fld->epsilon0;
  double mu0 = fld->mu0 = mom_fld->mu0;
  bool is_static = fld->is_static = mom_fld->is_static; 

  fld->ctx = mom_fld->ctx;
  fld->init = mom_fld->init;

  // The field follows the (resolved) app-level scheme. The CFL number is a
  // property of the scheme.
  fld->scheme_type = app->scheme_type;
  double cfl_frac = mom->cfl_frac == 0 ? 0.95 : mom->cfl_frac;
  fld->cfl = (fld->scheme_type == GKYL_MOMENT_MP ? 0.4 : 1.0) * cfl_frac;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_fld->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_fld->limiter;

  double c = 1/sqrt(epsilon0*mu0);
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_inew(&(struct gkyl_wv_maxwell_inp) {
      .c = c,
      .e_fact = mom_fld->elc_error_speed_fact,
      .b_fact = mom_fld->mag_error_speed_fact,
      .rp_type = WV_MAXWELL_RP_ROE,
      .embed_geo = mom_fld->embed_geo,
      .use_gpu = false,
    }
  );

  fld->maxwell = gkyl_wv_eqn_acquire(maxwell);
  
  int ndim = mom->ndim;

  if (fld->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    // create updaters for each directional update
    for (int d=0; d<ndim; ++d)
      fld->slvr[d] = gkyl_wave_prop_new( &(struct gkyl_wave_prop_inp) {
          .grid = &app->grid,
          .equation = maxwell,
          .split_type = GKYL_WAVE_QWAVE, // q-waves is fine for linear systems
          .limiter = limiter,
          .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
          .update_dirs = { d },
          .check_inv_domain = false,
          .cfl = fld->cfl,
          .geom = app->geom,
          .comm = app->comm
        }
      );

    // allocate arrays
    fld->fdup = mkarr(false, 8, app->local_ext.volume);
    for (int d=0; d<ndim+1; ++d)
      fld->f[d] = mkarr(false, 8, app->local_ext.volume);

    // set current solution so ICs and IO work properly
    fld->fcurr = fld->f[0];
  }
  else if (fld->scheme_type == GKYL_MOMENT_MP || fld->scheme_type == GKYL_MOMENT_KEP) {
    // NOTE: there is no KEP scheme for Maxwell, and we simply use MP scheme instead
    
    // determine directions to update
    int num_up_dirs = 0, update_dirs[GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<ndim; ++d)
      if (!app->is_dir_skipped[d]) {
        update_dirs[num_up_dirs] = d;
        num_up_dirs += 1;
      }

    // choose U3 as field reconstruction if using KEP
    enum gkyl_mp_recon mp_recon =
      fld->scheme_type == GKYL_MOMENT_KEP ? GKYL_MP_U3 : app->mp_recon;
    
    // single MP updater updates all directions
    fld->mp_slvr = gkyl_mp_scheme_new( &(struct gkyl_mp_scheme_inp) {
        .grid = &app->grid,
        .equation = maxwell,
        .mp_recon = mp_recon,
        .skip_mp_limiter = mom->skip_mp_limiter,
        .num_up_dirs = num_up_dirs,
        .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
        .cfl = fld->cfl,
        .geom = app->geom,
      }
    );

    // allocate arrays
    fld->fdup = mkarr(false, 8, app->local_ext.volume);
    fld->f0 = mkarr(false, 8, app->local_ext.volume);
    fld->f1 = mkarr(false, 8, app->local_ext.volume);
    fld->fnew = mkarr(false, 8, app->local_ext.volume);
    fld->cflrate = mkarr(false, 1, app->local_ext.volume);

    // set current solution so ICs and IO work properly
    fld->fcurr = fld->f0;
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int i=0; i<3; ++i) {
    fld->lower_bc[i] = 0;
    fld->upper_bc[i] = 0;
  }

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_field_bc_type *bc;
      if (dir == 0)
        bc = mom_fld->bcx;
      else if (dir == 1)
        bc = mom_fld->bcy;
      else
        bc = mom_fld->bcz;
      
      wv_bc_func_t bc_lower_func;
      if (dir == 0)
        bc_lower_func = mom_fld->bcx_func[0];
      else if (dir == 1)
        bc_lower_func = mom_fld->bcy_func[0];
      else
        bc_lower_func = mom_fld->bcz_func[0];

      wv_bc_func_t bc_upper_func;
      if (dir == 0)
        bc_upper_func = mom_fld->bcx_func[1];
      else if (dir == 1)
        bc_upper_func = mom_fld->bcy_func[1];
      else
        bc_upper_func = mom_fld->bcz_func[1];

      fld->lower_bct[dir] = bc[0];
      fld->upper_bct[dir] = bc[1];

      switch (bc[0]) {
        case GKYL_FIELD_PEC_WALL:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            maxwell->wall_bc_func, 0);
          break;

        case GKYL_FIELD_FUNC:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_lower_func, mom_fld->ctx);
          break;

        case GKYL_FIELD_COPY:
        case GKYL_FIELD_WEDGE:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_copy, 0);
          break;

        case GKYL_FIELD_SKIP:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_skip, 0);
          break;

        default:
          assert(false);
          break;
      }

      switch (bc[1]) {
        case GKYL_FIELD_PEC_WALL:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            maxwell->wall_bc_func, 0);
          break;

        case GKYL_FIELD_FUNC:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_upper_func, mom_fld->ctx);
          break;

        case GKYL_FIELD_COPY:
        case GKYL_FIELD_WEDGE:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE,
            nghost, bc_copy, 0);
          break;

        case GKYL_FIELD_SKIP:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE,
            nghost, bc_skip, 0);
          break;          
          
        default:
          assert(false);
          break;
      }
    }
  }

  fld->use_explicit_em_coupling = mom_fld->use_explicit_em_coupling;
  fld->use_gr_em_coupling = mom_fld->use_gr_em_coupling;

  fld->embed_mask = mkarr(false, 1, app->local_ext.volume);
  gkyl_array_clear(fld->embed_mask, 1.0);

  if (maxwell->embed_geo) {
    gkyl_wv_embed_geo_new_mask(maxwell->embed_geo, &app->grid,
      &app->local, fld->embed_mask);
  }

  fld->ext_em = mkarr(false, 6, app->local_ext.volume);
  gkyl_array_clear(fld->ext_em, 0.0);
  fld->has_ext_em = false;
  fld->ext_em_evolve = false;
  fld->t_ramp_E = mom_fld->t_ramp_E ? mom_fld->t_ramp_E : 0.0;
  // setup external electromagnetic field
  if (mom_fld->ext_em) {
    fld->has_ext_em = true;
    // Only set the external field to evolve if a user asks for 
    // dynamic external field and t_ramp_E = 0.0, otherwise
    // we project the external field once and any variation
    // in time is encoded linearly in t_ramp_E. 
    if (mom_fld->ext_em_evolve && mom_fld->t_ramp_E == 0.0) {
      fld->ext_em_evolve = mom_fld->ext_em_evolve;
    }
    fld->ext_em_proj = gkyl_fv_proj_new(&app->grid, 2, GKYL_MOM_APP_NUM_EXT_EM,
      mom_fld->ext_em, mom_fld->ext_em_ctx);    
  }  

  fld->app_current = mkarr(false, 3, app->local_ext.volume);
  gkyl_array_clear(fld->app_current, 0.0);
  if(mom_fld->use_explicit_em_coupling){
    fld->app_current1 = mkarr(false, 3, app->local_ext.volume);
    gkyl_array_clear(fld->app_current1, 0.0);
    fld->app_current2 = mkarr(false, 3, app->local_ext.volume);
    gkyl_array_clear(fld->app_current2, 0.0);
  }
  fld->has_app_current = false;
  fld->app_current_evolve = false;
  if (mom_fld->app_current) {
    fld->has_app_current = true;
    // Only set the applied current to evolve if a user asks for 
    // dynamic applied current and t_ramp_curr = 0.0, otherwise
    // we project the applied current once and any variation
    // in time is encoded linearly in t_ramp_curr. 
    if (mom_fld->app_current_evolve &&  mom_fld->t_ramp_curr == 0.0) {
      fld->app_current_evolve = mom_fld->app_current_evolve;
    }
    fld->app_current_proj = gkyl_fv_proj_new(&app->grid, 2, GKYL_MOM_APP_NUM_APPLIED_CURRENT,
      mom_fld->app_current, mom_fld->app_current_ctx);  
  }

  fld->volume_sources = mom_fld->volume_sources;

  // allocate buffer for applying BCs (used for periodic BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->ndim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  fld->bc_buffer = mkarr(false, 8, buff_sz);

  gkyl_wv_eqn_release(maxwell);

  fld->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  fld->is_first_energy_write_call = true;

  // ---- Method wiring -------------------------------------------------------
  // The single place scheme/static dispatch is decided; everything
  // downstream calls through these pointers.
  fld->apply_bc_func = field_apply_bc;
  fld->copy_func = field_copy;
  fld->apply_ic_func = field_apply_ic;
  fld->write_func = field_write;
  fld->calc_energy_func = field_calc_energy;
  fld->write_energy_func = field_write_energy;
  fld->read_func = field_from_file;
  if (fld->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    fld->max_dt_func = field_max_dt_wave_prop;
    fld->update_func = field_update_wave_prop;
    fld->rhs_func = NULL; // the wave-prop full step does not decompose into an RHS
    fld->stage_state_func = field_stage_state_wave_prop;
    fld->step_commit_func = field_step_commit_wave_prop;
    fld->release_func = field_release_wave_prop;
  }
  else { // MP and KEP runs both use the MP Maxwell solver
    fld->max_dt_func = field_max_dt_mp;
    fld->update_func = field_update_ssp_rk3;
    fld->rhs_func = field_rhs_mp;
    fld->stage_state_func = field_stage_state_curr;
    fld->step_commit_func = field_step_commit_curr;
    fld->release_func = field_release_mp;
  }
  // A static field freezes its hyperbolic update, whatever the scheme: the
  // solution stays in fcurr (the in-place stage-state/commit set) and the
  // RHS is zero for the RK stages.
  if (fld->is_static) {
    fld->update_func = field_update_static;
    fld->rhs_func = field_rhs_static;
    fld->stage_state_func = field_stage_state_curr;
    fld->step_commit_func = field_step_commit_curr;
  }
}

// ---- apply_bc ----------------------------------------------------------------

static void
no_field_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_field *field, struct gkyl_array *f) { }

static void
field_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_field *field, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};
  
  for (int d=0; d<num_periodic_dir; ++d)
    is_non_periodic[app->periodic_dirs[d]] = 0;

  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {
      // handle non-wedge BCs
      if (field->lower_bct[d] != GKYL_FIELD_WEDGE)
        gkyl_wv_apply_bc_advance(field->lower_bc[d], tcurr, &app->local, f);
      if (field->upper_bct[d] != GKYL_FIELD_WEDGE)      
        gkyl_wv_apply_bc_advance(field->upper_bc[d], tcurr, &app->local, f);

      // wedge BCs for upper/lower must be handled in one shot
      if (field->lower_bct[d] == GKYL_FIELD_WEDGE)
        moment_apply_wedge_bc(app, tcurr, &app->local,
          field->bc_buffer, d, field->lower_bc[d], field->upper_bc[d], f);
    }

  // sync interior ghost cells
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);
  // sync periodic ghost cells
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext, num_periodic_dir,
    app->periodic_dirs, f);

  app->stat.field_bc_tm += gkyl_time_diff_now_sec(wst);
}

// apply BCs to EM field
void
moment_field_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_field *field, struct gkyl_array *f)
{
  field->apply_bc_func(app, tcurr, field, f);
}

// ---- max_dt ------------------------------------------------------------------

static double
no_field_max_dt(const gkyl_moment_app *app, const struct moment_field *fld)
{
  return DBL_MAX; // no field, no constraint
}

static double
field_max_dt_wave_prop(const gkyl_moment_app *app, const struct moment_field *fld)
{
  double max_dt = DBL_MAX;
  for (int d=0; d<app->ndim; ++d)
    max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(fld->slvr[d], &app->local, fld->f[0]));
  return max_dt;
}

static double
field_max_dt_mp(const gkyl_moment_app *app, const struct moment_field *fld)
{
  return gkyl_mp_scheme_max_dt(fld->mp_slvr, &app->local, fld->f0);
}

double
moment_field_max_dt(const gkyl_moment_app *app, const struct moment_field *fld)
{
  return fld->max_dt_func(app, fld);
}

// ---- update (wave-prop path: input in fld->f[0], output in fld->f[ndim]) ------

static struct gkyl_update_status
no_field_update(gkyl_moment_app *app, const struct moment_field *fld,
  double tcurr, double dt)
{
  return (struct gkyl_update_status) { .success = true, .dt_suggested = DBL_MAX };
}

static struct gkyl_update_status
field_update_wave_prop(gkyl_moment_app *app, const struct moment_field *fld,
  double tcurr, double dt)
{
  int ndim = fld->ndim;
  struct gkyl_wave_prop_status stat = { true, DBL_MAX };

  for (int d=0; d<ndim; ++d) {
    // update solution
    stat = gkyl_wave_prop_advance(fld->slvr[d], tcurr, dt, &app->local, fld->embed_mask, fld->f[d], fld->f[d+1]);

    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = stat.dt_suggested
      };
    // apply BC
    moment_field_apply_bc(app, tcurr, fld, fld->f[d+1]);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = stat.dt_suggested
  };
}

// Static field: nothing to evolve. The solution stays in fcurr (the
// in-place stage-state/commit set), so the sources read and the stepper
// commits the live field state -- a static field survives the step (unlike
// the historical behavior, which committed never-written staging arrays and
// silently erased it).
static struct gkyl_update_status
field_update_static(gkyl_moment_app *app, const struct moment_field *fld,
  double tcurr, double dt)
{
  return (struct gkyl_update_status) { .success = true, .dt_suggested = DBL_MAX };
}

// One forward-Euler stage of the SSP-RK3 update (see the species analogue).
static double
field_rk3_stage(gkyl_moment_app *app, const struct moment_field *fld,
  double tcurr, double dt, const struct gkyl_array *fin, struct gkyl_array *fout)
{
  app->stat.nfeuler += 1;

  double dtmin = fld->rhs_func(app, (struct moment_field*) fld, fin, fout);
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < 0.01)
    dtmin = dt; // avoid retaking steps on very small dt differences

  if (dtmin < dt)
    return dtmin;

  gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dt, &(app->local)),
    1.0, fin, &(app->local));
  moment_field_apply_bc(app, tcurr, fld, fout);

  return dtmin;
}

// A full SSP-RK3 step for the field, built on rhs_func. f0 holds the
// pre-step state and is only overwritten by the final commit, so a stage
// failure reports back and the stepper's redo protocol retakes the Strang
// step at the suggested dt.
static struct gkyl_update_status
field_update_ssp_rk3(gkyl_moment_app *app, const struct moment_field *fld,
  double tcurr, double dt)
{
  // Stage 1: f1 = f0 + dt*L(f0).
  double dt1 = field_rk3_stage(app, fld, tcurr, dt, fld->f0, fld->f1);
  if (dt1 < dt)
    return (struct gkyl_update_status) { .success = false, .dt_suggested = dt1 };

  // Stage 2: fnew = f1 + dt*L(f1); f1 = 3/4 f0 + 1/4 fnew.
  double dt2 = field_rk3_stage(app, fld, tcurr+dt, dt, fld->f1, fld->fnew);
  if (dt2 < dt) {
    double dt_rel_diff = (dt-dt2)/dt2;
    app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0], dt_rel_diff);
    app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1], dt_rel_diff);
    app->stat.nstage_2_fail += 1;
    return (struct gkyl_update_status) { .success = false, .dt_suggested = dt2 };
  }
  array_combine(fld->f1, 3.0/4.0, fld->f0, 1.0/4.0, fld->fnew, &app->local_ext);

  // Stage 3: fnew = f1 + dt*L(f1); f1 = 1/3 f0 + 2/3 fnew; commit f0 = f1.
  double dt3 = field_rk3_stage(app, fld, tcurr+dt/2, dt, fld->f1, fld->fnew);
  if (dt3 < dt) {
    double dt_rel_diff = (dt-dt3)/dt3;
    app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0], dt_rel_diff);
    app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1], dt_rel_diff);
    app->stat.nstage_3_fail += 1;
    return (struct gkyl_update_status) { .success = false, .dt_suggested = dt3 };
  }
  array_combine(fld->f1, 1.0/3.0, fld->f0, 2.0/3.0, fld->fnew, &app->local_ext);
  gkyl_array_copy_range(fld->f0, fld->f1, &app->local_ext);

  return (struct gkyl_update_status) { .success = true, .dt_suggested = dt3 };
}


// update solution: initial solution is in fld->f[0] and updated
// solution in fld->f[ndim]
struct gkyl_update_status
moment_field_update(gkyl_moment_app *app,
  const struct moment_field *fld, double tcurr, double dt)
{
  return fld->update_func(app, fld, tcurr, dt);
}

// ---- rhs (MP path) -------------------------------------------------------------

static double
no_field_rhs(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  return DBL_MAX;
}

static double
field_rhs_mp(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec tm = gkyl_wall_clock();

  gkyl_array_clear(fld->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  gkyl_mp_scheme_advance(fld->mp_slvr, &app->local, fin,
    app->ql, app->qr, app->amdq, app->apdq,
    fld->cflrate, fld->embed_mask, rhs);

  double omegaCfl[1];
  gkyl_array_reduce_range(omegaCfl, fld->cflrate, GKYL_MAX, &(app->local));

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(tm);

  return fld->cfl/omegaCfl[0];
}

// Static field: zero RHS (forward Euler then gives emout = emin), no CFL
// constraint of its own.
static double
field_rhs_static(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_clear(rhs, 0.0);
  return DBL_MAX;
}


// Compute RHS of EM equations
double
moment_field_rhs(gkyl_moment_app *app, struct moment_field *fld,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  return fld->rhs_func(app, fld, fin, rhs);
}

// ---- copy (time-stepper backup/commit/restore) ---------------------------------

static void
no_field_copy(const struct moment_field *fld,
  struct gkyl_array *dst, const struct gkyl_array *src) { }

static void
field_copy(const struct moment_field *fld,
  struct gkyl_array *dst, const struct gkyl_array *src)
{
  gkyl_array_copy(dst, src);
}

void
moment_field_copy(const struct moment_field *fld,
  struct gkyl_array *dst, const struct gkyl_array *src)
{
  fld->copy_func(fld, dst, src);
}

// ---- stepper protocol (backup/restore/stage-state/commit) ---------------------
// The pre-step backup/restore work through fcurr and copy_func, so they
// need no per-scheme variants and degrade to no-ops without a field.

void
moment_field_step_backup(const struct moment_field *fld)
{
  fld->copy_func(fld, fld->fdup, fld->fcurr);
}

void
moment_field_step_restore(const struct moment_field *fld)
{
  fld->copy_func(fld, fld->fcurr, fld->fdup);
}

static struct gkyl_array*
field_stage_state_wave_prop(const struct moment_field *fld, int nstrang)
{
  return nstrang == 0 ? fld->f[0] : fld->f[fld->ndim];
}

static void
field_step_commit_wave_prop(const struct moment_field *fld)
{
  fld->copy_func(fld, fld->f[0], fld->f[fld->ndim]);
}

// In-place set (MP scheme, static fields, and the no-field object): the
// current solution lives in fcurr through the whole step (NULL without a
// field), so both stage states are fcurr and there is nothing to commit.
static struct gkyl_array*
field_stage_state_curr(const struct moment_field *fld, int nstrang)
{
  return fld->fcurr;
}

static void
field_step_commit_curr(const struct moment_field *fld)
{
}

struct gkyl_array*
moment_field_stage_state(const struct moment_field *fld, int nstrang)
{
  return fld->stage_state_func(fld, nstrang);
}

void
moment_field_step_commit(const struct moment_field *fld)
{
  fld->step_commit_func(fld);
}

// ---- initial conditions --------------------------------------------------------

static void
no_field_apply_ic(gkyl_moment_app *app, struct moment_field *fld, double t0) { }

static void
field_apply_ic(gkyl_moment_app *app, struct moment_field *fld, double t0)
{
  int num_quad = fld->scheme_type == GKYL_MOMENT_MP ? 4 : 2;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, num_quad, 8, fld->init, fld->ctx);

  gkyl_fv_proj_advance(proj, t0, &app->local, fld->fcurr);
  gkyl_fv_proj_release(proj);

  if (fld->has_ext_em) {
    gkyl_fv_proj_advance(fld->ext_em_proj, t0, &app->local, fld->ext_em);
  }
  if (fld->has_app_current) {
    gkyl_fv_proj_advance(fld->app_current_proj, t0, &app->local, fld->app_current);
  }

  moment_field_apply_bc(app, t0, fld, fld->fcurr);
}

void
moment_field_apply_ic(gkyl_moment_app *app, struct moment_field *fld, double t0)
{
  fld->apply_ic_func(app, fld, t0);
}

// ---- write (frame IO) ------------------------------------------------------------

static void
no_field_write(const gkyl_moment_app *app, const struct moment_field *fld,
  double tm, int frame) { }

static void
field_write(const gkyl_moment_app *app, const struct moment_field *fld,
  double tm, int frame)
{
  struct gkyl_msgpack_data *mt = moment_array_meta_new( (struct moment_output_meta) {
      .frame = frame,
      .stime = tm
    }
  );

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, fld->fcurr, fileNm.str);
  cstr_drop(&fileNm);

  // write external EM field if it is present
  if (fld->ext_em_proj) {
    if (fld->ext_em_evolve || frame == 0) {
      cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "ext_em", frame);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, fld->ext_em, fileNm.str);
      cstr_drop(&fileNm);
    }
  }

  // write applied current if it is present
  if (fld->app_current_proj) {
    if (fld->app_current_evolve || frame == 0) {
      cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "app_current", frame);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, fld->app_current, fileNm.str);
      cstr_drop(&fileNm);
    }
  }

  moment_array_meta_release(mt);
}

void
moment_field_write(const gkyl_moment_app *app, const struct moment_field *fld,
  double tm, int frame)
{
  fld->write_func(app, fld, tm, frame);
}

// ---- energy diagnostics ------------------------------------------------------------

static void
no_field_calc_energy(gkyl_moment_app *app, struct moment_field *fld, double tm) { }

static void
field_calc_energy(gkyl_moment_app *app, struct moment_field *fld, double tm)
{
  double energy[6] = { 0.0 };
  calc_integ_quant(fld->maxwell, app->grid.cellVolume, fld->fcurr, app->geom,
    app->local, energy);

  double energy_global[6] = { 0.0 };
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, energy, energy_global);

  gkyl_dynvec_append(fld->integ_energy, tm, energy_global);
}

static void
no_field_write_energy(gkyl_moment_app *app, struct moment_field *fld) { }

static void
field_write_energy(gkyl_moment_app *app, struct moment_field *fld)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // write out field energy
    cstr fileNm = cstr_from_fmt("%s-field-energy.gkyl", app->name);

    if (fld->is_first_energy_write_call) {
      // write to a new file (this ensure previous output is removed)
      gkyl_dynvec_write(fld->integ_energy, fileNm.str);
      fld->is_first_energy_write_call = false;
    }
    else {
      // append to existing file
      gkyl_dynvec_awrite(fld->integ_energy, fileNm.str);
    }
    cstr_drop(&fileNm);
  }
  gkyl_dynvec_clear(fld->integ_energy);
}

void
moment_field_calc_energy(gkyl_moment_app *app, struct moment_field *fld, double tm)
{
  fld->calc_energy_func(app, fld, tm);
}

void
moment_field_write_energy(gkyl_moment_app *app, struct moment_field *fld)
{
  fld->write_energy_func(app, fld);
}

// ---- restart read ------------------------------------------------------------------

static struct gkyl_app_restart_status
no_field_read(gkyl_moment_app *app, struct moment_field *fld, const char *fname)
{
  return (struct gkyl_app_restart_status) {
    .io_status = GKYL_ARRAY_RIO_SUCCESS,
    .frame = 0,
    .stime = 0.0
  };
}

static struct gkyl_app_restart_status
field_from_file(gkyl_moment_app *app, struct moment_field *fld, const char *fname)
{
  struct gkyl_app_restart_status rstat = moment_app_header_from_file(app, fname);

  if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, fld->fcurr, fname);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      moment_field_apply_bc(app, rstat.stime, fld, fld->fcurr);
    }
  }

  // Compute external EM field and applied current if present
  // Computation necessary in case external EM field or applied current
  // are time-independent and not computed in the time-stepping loop
  // since they are not read-in as part of restarts.
  if (fld->has_ext_em) {
    gkyl_fv_proj_advance(fld->ext_em_proj, rstat.stime, &app->local, fld->ext_em);
  }
  if (fld->has_app_current) {
    gkyl_fv_proj_advance(fld->app_current_proj, rstat.stime, &app->local, fld->app_current);
  }

  return rstat;
}

struct gkyl_app_restart_status
moment_field_from_file(gkyl_moment_app *app, struct moment_field *fld, const char *fname)
{
  return fld->read_func(app, fld, fname);
}

// ---- release --------------------------------------------------------------------

static void
no_field_release(const struct moment_field *fld) { }

// Release the members common to the wave-prop and MP variants.
static void
field_release_common(const struct moment_field *fld)
{
  gkyl_wv_eqn_release(fld->maxwell);

  for (int d=0; d<fld->ndim; ++d) {
    if (fld->lower_bc[d])
      gkyl_wv_apply_bc_release(fld->lower_bc[d]);
    if (fld->upper_bc[d])
      gkyl_wv_apply_bc_release(fld->upper_bc[d]);
  }

  gkyl_array_release(fld->ext_em);
  if (fld->has_ext_em) {
    gkyl_fv_proj_release(fld->ext_em_proj);
  }

  gkyl_array_release(fld->app_current);
  if(fld->use_explicit_em_coupling) {
    gkyl_array_release(fld->app_current1);
    gkyl_array_release(fld->app_current2);
  }
  if (fld->has_app_current) {
    gkyl_fv_proj_release(fld->app_current_proj);
  }

  gkyl_array_release(fld->embed_mask);

  gkyl_dynvec_release(fld->integ_energy);
  gkyl_array_release(fld->bc_buffer);
}

static void
field_release_wave_prop(const struct moment_field *fld)
{
  field_release_common(fld);

  for (int d=0; d<fld->ndim; ++d)
    gkyl_wave_prop_release(fld->slvr[d]);

  gkyl_array_release(fld->fdup);
  for (int d=0; d<fld->ndim+1; ++d)
    gkyl_array_release(fld->f[d]);
}

static void
field_release_mp(const struct moment_field *fld)
{
  field_release_common(fld);

  gkyl_mp_scheme_release(fld->mp_slvr);
  gkyl_array_release(fld->fdup);
  gkyl_array_release(fld->f0);
  gkyl_array_release(fld->f1);
  gkyl_array_release(fld->fnew);
  gkyl_array_release(fld->cflrate);
}

// free field
void
moment_field_release(const struct moment_field *fld)
{
  fld->release_func(fld);
}

