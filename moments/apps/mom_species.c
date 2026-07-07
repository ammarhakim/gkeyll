#include <gkyl_moment_priv.h>
#include <gkyl_util.h>
#include <gkyl_wv_embed_geo.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_ten_moment.h>

// ---- Scheme/behavior variants (installed by moment_species_init) ----------
// All scheme and species-type dispatch is wired once, at init: the
// moment_species_* wrappers at the bottom of this file (and the time
// steppers) call through the sp->*_func pointers and never branch on
// scheme_type / is_static / eqn_type themselves.

// Wave-propagation scheme (dimensionally-split second-order FV; stepped by
// the one-step Strang stepper through update_func).
static double species_max_dt_wave_prop(const gkyl_moment_app *app,
  const struct moment_species *sp);
static struct gkyl_update_status species_update_wave_prop(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt);
static void species_release_wave_prop(const struct moment_species *sp);

// MP and KEP schemes (RHS-based; stepped by the SSP-RK3 stepper through
// rhs_func).
static double species_max_dt_mp(const gkyl_moment_app *app,
  const struct moment_species *sp);
static double species_max_dt_kep(const gkyl_moment_app *app,
  const struct moment_species *sp);
static double species_rhs_mp(gkyl_moment_app *app, struct moment_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *rhs);
static double species_rhs_kep(gkyl_moment_app *app, struct moment_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *rhs);
static void species_release_mp(const struct moment_species *sp);
static void species_release_kep(const struct moment_species *sp);

// Static species: the hyperbolic update passes the state through unchanged
// and imposes no CFL constraint of its own (sources, if any, still apply in
// the operator-split source step).
static struct gkyl_update_status species_update_static(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt);
static double species_rhs_static(gkyl_moment_app *app, struct moment_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

// MP/KEP full-step update: a complete SSP-RK3 step built on rhs_func, run
// inside the same Strang-split stepper as the wave-prop one-step update.
static struct gkyl_update_status species_update_ssp_rk3(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt);

// Stepper protocol: Strang stage-state accessor and end-of-step commit, per
// state-array layout. Wave-prop steps f[0] -> f[ndim] and commits at the
// end of the Strang step; the RK schemes and static species hold the
// current solution in fcurr throughout (commit is a no-op). The pre-step
// backup/restore for the redo protocol works through fcurr and needs no
// per-scheme variants.
static struct gkyl_array* species_stage_state_wave_prop(const struct moment_species *sp, int nstrang);
static void species_step_commit_wave_prop(const struct moment_species *sp);
static struct gkyl_array* species_stage_state_curr(const struct moment_species *sp, int nstrang);
static void species_step_commit_curr(const struct moment_species *sp);

// Post-update hook, applied after each directional sweep. Composed from
// (species type x spacetime presence): the excision scrub for modular GR
// tetrad species under a spacetime, a no-op otherwise.
static void species_no_post_update(gkyl_moment_app *app,
  struct moment_species *sp, struct gkyl_array *f);
static void species_post_update_scrub_excised(gkyl_moment_app *app,
  struct moment_species *sp, struct gkyl_array *f);

// Stepper state copy (backup/commit/restore).
static void species_copy(const struct moment_species *sp,
  struct gkyl_array *dst, const struct gkyl_array *src);

// ICs, IO, integrated diagnostics, and restart (one variant each for now;
// wired through the vtable so species types can override them later).
static void species_apply_ic(gkyl_moment_app *app,
  struct moment_species *sp, double t0);
static void species_write(const gkyl_moment_app *app,
  const struct moment_species *sp, double tm, int frame);
static void species_calc_integrated_mom(gkyl_moment_app *app,
  struct moment_species *sp, double tm);
static void species_write_integrated_mom(gkyl_moment_app *app,
  struct moment_species *sp);
static struct gkyl_app_restart_status species_from_file(gkyl_moment_app *app,
  struct moment_species *sp, const char *fname);

// initialize species
void
moment_species_init(const struct gkyl_moment *mom, const struct gkyl_moment_species *mom_sp,
  struct gkyl_moment_app *app, struct moment_species *sp)
{
  sp->ndim = mom->ndim;
  strcpy(sp->name, mom_sp->name);
  sp->charge = mom_sp->charge;
  sp->mass = mom_sp->mass;
  sp->is_static = mom_sp->is_static; 
  sp->ctx = mom_sp->ctx;
  sp->init = mom_sp->init;

  sp->eqn_type = mom_sp->equation->type;
  sp->num_equations = mom_sp->equation->num_equations;
  sp->equation = gkyl_wv_eqn_acquire(mom_sp->equation);

  // Do we need to update source terms for this fluid?
  // Sources can be electromagnetic fields, closure-related, applied accelerations
  // volume expansion, friction, reactivity, etc.
  sp->update_sources = false;
  if (app->has_field) {
    sp->update_sources = true;
  }
  // Modular GR equation types always require source integration through the
  // spacetime-coupling object — turn sources on regardless of the rest of
  // the species configuration. Note: tetrad species deliberately do NOT set
  // has_gr_euler (the self-contained GR-Euler source flag), so the
  // EM-coupling object's GR-source path skips them; the spacetime-coupling
  // object owns them.
  if (sp->eqn_type == GKYL_EQN_GR_EULER_TETRAD) {
    sp->update_sources = true;
  }

  sp->k0 = 0.0;
  sp->has_grad_closure = false;
  sp->has_nn_closure = false;
  if (mom_sp->equation->type == GKYL_EQN_TEN_MOMENT) {
    sp->update_sources = true; 
    sp->k0 = gkyl_wv_ten_moment_k0(mom_sp->equation);
    sp->has_grad_closure = gkyl_wv_ten_moment_use_grad_closure(mom_sp->equation);
    sp->has_nn_closure = gkyl_wv_ten_moment_use_nn_closure(mom_sp->equation);
    sp->poly_order = gkyl_wv_ten_moment_poly_order(mom_sp->equation);
    sp->ann = gkyl_wv_ten_moment_ann(mom_sp->equation);
  }

  // check if we are running with Braginskii transport and fetch Braginskii type
  if (app->has_braginskii) {
    sp->update_sources = true;
    sp->type_brag = mom_sp->type_brag;
  }

  // Source-family parameter bundles travel by struct assignment; any enabled
  // family requires the operator-split source update.
  sp->friction = mom_sp->friction;
  sp->volume_sources = mom_sp->volume_sources;
  sp->reactivity = mom_sp->reactivity;
  sp->einstein_medium = mom_sp->einstein_medium;
  sp->gr_ultra_rel = mom_sp->gr_ultra_rel;
  sp->gr_euler = mom_sp->gr_euler;
  sp->gr_twofluid = mom_sp->gr_twofluid;
  sp->vacuum_einstein = mom_sp->vacuum_einstein;
  sp->vacuum_einstein_conformal = mom_sp->vacuum_einstein_conformal;
  sp->gr_mhd = mom_sp->gr_mhd;
  if (sp->friction.enabled || sp->volume_sources.enabled || sp->reactivity.enabled ||
    sp->einstein_medium.enabled || sp->gr_ultra_rel.enabled || sp->gr_euler.enabled ||
    sp->gr_twofluid.enabled || sp->vacuum_einstein.enabled ||
    sp->vacuum_einstein_conformal.enabled || sp->gr_mhd.enabled) {
    sp->update_sources = true;
  }

  // Resolve this species' scheme: the input may override the app-level
  // default (GKYL_MOMENT_DEFAULT inherits it). The CFL number is a property
  // of the scheme, so it is per-species too.
  sp->scheme_type = mom_sp->scheme_type == GKYL_MOMENT_DEFAULT
    ? app->scheme_type : mom_sp->scheme_type;
  double cfl_frac = mom->cfl_frac == 0 ? 0.95 : mom->cfl_frac;
  sp->cfl = (sp->scheme_type == GKYL_MOMENT_MP ? 0.4 : 1.0) * cfl_frac;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  enum gkyl_wave_split_type split_type = mom_sp->split_type;    

  int ndim = mom->ndim;
  int meqn = sp->num_equations;  

  if (sp->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    // create updaters for each directional update
    for (int d=0; d<ndim; ++d)
      sp->slvr[d] = gkyl_wave_prop_new( &(struct gkyl_wave_prop_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .split_type = split_type,
          .limiter = limiter,
          .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
          .force_low_order_flux = mom_sp->force_low_order_flux,
          .check_inv_domain = true,
          .update_dirs = { d },
          .cfl = sp->cfl,
          .geom = app->geom,
          .comm = app->comm
        }
      );
      
    sp->fdup = mkarr(false, meqn, app->local_ext.volume);
    // allocate arrays
    for (int d=0; d<ndim+1; ++d)
      sp->f[d] = mkarr(false, meqn, app->local_ext.volume);
  
    // set current solution so ICs and IO work properly
    sp->fcurr = sp->f[0];
  }
  else if ( sp->scheme_type == GKYL_MOMENT_MP || sp->scheme_type == GKYL_MOMENT_KEP ) {
    // determine directions to update
    int num_up_dirs = 0, update_dirs[GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<ndim; ++d)
      if (!app->is_dir_skipped[d]) {
        update_dirs[num_up_dirs] = d;
        num_up_dirs += 1;
      }

    if (sp->scheme_type == GKYL_MOMENT_MP)
      // single MP updater updates all directions
      sp->mp_slvr = gkyl_mp_scheme_new( &(struct gkyl_mp_scheme_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .mp_recon = app->mp_recon,
          .skip_mp_limiter = mom->skip_mp_limiter,
          .num_up_dirs = num_up_dirs,
          .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
          .cfl = sp->cfl,
          .geom = app->geom,
        }
      );
    else
      // single KEP updater updates all directions
      sp->kep_slvr = gkyl_kep_scheme_new( &(struct gkyl_kep_scheme_inp) {
          .grid = &app->grid,
          .equation = mom_sp->equation,
          .use_hybrid_flux = app->use_hybrid_flux_kep,
          .num_up_dirs = num_up_dirs,
          .update_dirs = { update_dirs[0], update_dirs[1], update_dirs[2] } ,
          .cfl = sp->cfl,
          .geom = app->geom,
        }
      );
    
    // allocate arrays
    sp->fdup = mkarr(false, meqn, app->local_ext.volume);
    sp->f0 = mkarr(false, meqn, app->local_ext.volume);
    sp->f1 = mkarr(false, meqn, app->local_ext.volume);
    sp->fnew = mkarr(false, meqn, app->local_ext.volume);
    sp->cflrate = mkarr(false, 1, app->local_ext.volume);
    sp->alpha = mkarr(false, 1, app->local_ext.volume);
    
    // set current solution so ICs and IO work properly
    sp->fcurr = sp->f0;
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int i=0; i<3; ++i) {
    sp->lower_bc[i] = 0;
    sp->upper_bc[i] = 0;
  }

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = mom_sp->bcx;
      else if (dir == 1)
        bc = mom_sp->bcy;
      else
        bc = mom_sp->bcz;

      wv_bc_func_t bc_lower_func;
      if (dir == 0)
        bc_lower_func = mom_sp->bcx_func[0];
      else if (dir == 1)
        bc_lower_func = mom_sp->bcy_func[0];
      else
        bc_lower_func = mom_sp->bcz_func[0];

      wv_bc_func_t bc_upper_func;
      if (dir == 0)
        bc_upper_func = mom_sp->bcx_func[1];
      else if (dir == 1)
        bc_upper_func = mom_sp->bcy_func[1];
      else
        bc_upper_func = mom_sp->bcz_func[1];

      sp->lower_bct[dir] = bc[0];
      sp->upper_bct[dir] = bc[1];

      // lower BCs
      switch (bc[0]) {
        case GKYL_SPECIES_REFLECT:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_lower_func, mom_sp->ctx);
          break;
        
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE: // wedge also uses bc_copy
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_copy, 0);
          break;

        case GKYL_SPECIES_SKIP:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            bc_skip, 0);
          break;

        default:
          assert(false);
          break;
      }
      
      // upper BCs
      switch (bc[1]) {
        case GKYL_SPECIES_REFLECT:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_upper_func, mom_sp->ctx);
          break;
          
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_copy, 0);
          break;

        case GKYL_SPECIES_SKIP:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            bc_skip, 0);
          break;

        default:
          assert(false);
          break;
      }
    }
  }

  // allocate array for applied acceleration/forces for each species
  sp->app_accel = mkarr(false, 3, app->local_ext.volume);
  gkyl_array_clear(sp->app_accel, 0.0);
  sp->has_app_accel = false;
  sp->app_accel_evolve = false;
  if (mom_sp->app_accel) {
    sp->update_sources = true; 
    sp->has_app_accel = true;
    if (mom_sp->app_accel_evolve) {
      sp->app_accel_evolve = mom_sp->app_accel_evolve;
    }
    sp->app_accel_proj = gkyl_fv_proj_new(&app->grid, 2, GKYL_MOM_APP_NUM_APPLIED_ACCELERATION,
      mom_sp->app_accel, mom_sp->app_accel_ctx);  
  }

  sp->embed_mask = mkarr(false, 1, app->local_ext.volume);
  gkyl_array_clear(sp->embed_mask, 1.0);

  if (mom_sp->equation->embed_geo) {
    gkyl_wv_embed_geo_new_mask(mom_sp->equation->embed_geo, &app->grid,
      &app->local, sp->embed_mask);
  }

  sp->nT_source = mkarr(false, 2, app->local_ext.volume);
  sp->nT_source_is_set = false;
  sp->proj_nT_source = 0;
  if (mom_sp->nT_source_func) {
    sp->update_sources = true; 
    void *ctx = sp->ctx;
    if (mom_sp->nT_source_ctx)
      ctx = mom_sp->nT_source_ctx;
    
    sp->proj_nT_source = gkyl_fv_proj_new(&app->grid, 2, GKYL_MOM_APP_NUM_NT_SOURCE,
        mom_sp->nT_source_func, ctx);
    sp->nT_source_set_only_once = mom_sp->nT_source_set_only_once;
  }

  // allocate buffer for applying BCs (used for periodic BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->ndim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  sp->bc_buffer = mkarr(false, meqn, buff_sz);

  // Size the integrated-diagnostics vector by the equation's diagnostic
  // count, which the equation declares (Euler: 6, KE and PE stored
  // independently; GR Euler: 5; vacuum Einstein: 1; ...). Sizing it by
  // num_equations instead reads past the num_diag-sized buffer that
  // calc_integ_quant fills whenever num_diag < num_equations.
  sp->integ_q = gkyl_dynvec_new(GKYL_DOUBLE, mom_sp->equation->num_diag);

  sp->is_first_q_write_call = true;

  // ---- Method wiring -------------------------------------------------------
  // The single place scheme/static/type dispatch is decided; everything
  // downstream calls through these pointers. The scheme is a per-species
  // property here (currently inherited from the app input), so species
  // updated by different schemes can coexist.
  sp->copy_func = species_copy;
  sp->apply_ic_func = species_apply_ic;
  sp->write_func = species_write;
  sp->calc_integrated_mom_func = species_calc_integrated_mom;
  sp->write_integrated_mom_func = species_write_integrated_mom;
  sp->read_func = species_from_file;
  switch (sp->scheme_type) {
    case GKYL_MOMENT_WAVE_PROP:
      sp->max_dt_func = species_max_dt_wave_prop;
      sp->update_func = species_update_wave_prop;
      sp->rhs_func = NULL; // the wave-prop full step does not decompose into an RHS
      sp->stage_state_func = species_stage_state_wave_prop;
      sp->step_commit_func = species_step_commit_wave_prop;
      sp->release_func = species_release_wave_prop;
      break;
    case GKYL_MOMENT_MP:
      sp->max_dt_func = species_max_dt_mp;
      sp->update_func = species_update_ssp_rk3;
      sp->rhs_func = species_rhs_mp;
      sp->stage_state_func = species_stage_state_curr;
      sp->step_commit_func = species_step_commit_curr;
      sp->release_func = species_release_mp;
      break;
    case GKYL_MOMENT_KEP:
      sp->max_dt_func = species_max_dt_kep;
      sp->update_func = species_update_ssp_rk3;
      sp->rhs_func = species_rhs_kep;
      sp->stage_state_func = species_stage_state_curr;
      sp->step_commit_func = species_step_commit_curr;
      sp->release_func = species_release_kep;
      break;
    default:
      assert(false); // species scheme must resolve to a concrete scheme
      break;
  }
  // A static species freezes its hyperbolic update, whatever the scheme:
  // the solution stays in fcurr (so the stage states and commit are the
  // in-place set) and the RHS is zero for the RK stages.
  if (sp->is_static) {
    sp->update_func = species_update_static;
    sp->rhs_func = species_rhs_static;
    sp->stage_state_func = species_stage_state_curr;
    sp->step_commit_func = species_step_commit_curr;
  }
  // Post-update hook: modular GR tetrad species under a spacetime must scrub
  // hydro in excised cells after each sweep; everything else is a no-op. The
  // spacetime is initialized before the species, so its presence is known.
  sp->post_update_func =
    (sp->eqn_type == GKYL_EQN_GR_EULER_TETRAD && app->has_spacetime)
    ? species_post_update_scrub_excised
    : species_no_post_update;
}

// apply BCs to species
void
moment_species_apply_bc(gkyl_moment_app *app, double tcurr,
  const struct moment_species *sp, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};

  for (int d=0; d<num_periodic_dir; ++d)
    is_non_periodic[app->periodic_dirs[d]] = 0;

  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {
      // handle non-wedge BCs
      if (sp->lower_bct[d] != GKYL_SPECIES_WEDGE)
        gkyl_wv_apply_bc_advance(sp->lower_bc[d], tcurr, &app->local, f);
      if (sp->upper_bct[d] != GKYL_SPECIES_WEDGE)
        gkyl_wv_apply_bc_advance(sp->upper_bc[d], tcurr, &app->local, f);

      // wedge BCs for upper/lower must be handled in one shot
      if (sp->lower_bct[d] == GKYL_SPECIES_WEDGE)
        moment_apply_wedge_bc(app, tcurr, &app->local,
          sp->bc_buffer, d, sp->lower_bc[d], sp->upper_bc[d], f);
    }

  // sync interior ghost cells
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);
  // sync periodic ghost cells
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext, num_periodic_dir,
    app->periodic_dirs, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// ---- max_dt ----------------------------------------------------------------

static double
species_max_dt_wave_prop(const gkyl_moment_app *app, const struct moment_species *sp)
{
  double max_dt = DBL_MAX;
  for (int d=0; d<app->ndim; ++d)
    max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(sp->slvr[d], &app->local, sp->f[0]));
  return max_dt;
}

static double
species_max_dt_mp(const gkyl_moment_app *app, const struct moment_species *sp)
{
  return gkyl_mp_scheme_max_dt(sp->mp_slvr, &app->local, sp->f0);
}

static double
species_max_dt_kep(const gkyl_moment_app *app, const struct moment_species *sp)
{
  return gkyl_kep_scheme_max_dt(sp->kep_slvr, &app->local, sp->f0);
}

// maximum stable time-step
double
moment_species_max_dt(const gkyl_moment_app *app, const struct moment_species *sp)
{
  return sp->max_dt_func(app, sp);
}

// ---- update (wave-prop path: input in sp->f[0], output in sp->f[ndim]) ------

static struct gkyl_update_status
species_update_wave_prop(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt)
{
  int ndim = sp->ndim;
  double dt_suggested = DBL_MAX;
  double max_speed = 0.0;
  struct gkyl_wave_prop_status stat;

  for (int d=0; d<ndim; ++d) {
    stat = gkyl_wave_prop_advance(sp->slvr[d], tcurr, dt, &app->local, sp->embed_mask, sp->f[d], sp->f[d+1]);

    double my_max_speed = stat.max_speed;
    max_speed = max_speed > my_max_speed ? max_speed : my_max_speed;

    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = stat.dt_suggested
      };

    dt_suggested = fmin(dt_suggested, stat.dt_suggested);
    moment_species_apply_bc(app, tcurr, sp, sp->f[d+1]);

    sp->post_update_func(app, sp, sp->f[d+1]);
  }

  for (int d=0; d<ndim; ++d) {
    struct gkyl_wv_eqn *eqn = sp->equation;
    if (eqn->type==GKYL_EQN_MHD) {
      gkyl_wv_mhd_set_glm_ch(eqn, max_speed);
    }
  }
  if (app->update_mhd_source) {
    // app->mhd_source.slvr->glm_ch = max_speed;
    // gkyl_mhd_src_set_glm_ch(app->mhd_source.slvr, max_speed);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

// Static species: nothing to evolve. The solution stays in fcurr (the
// in-place stage-state/commit set), so no state needs to pass through the
// staging arrays.
static struct gkyl_update_status
species_update_static(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt)
{
  return (struct gkyl_update_status) { .success = true, .dt_suggested = DBL_MAX };
}

// One forward-Euler stage of the SSP-RK3 update: rhs into fout, then
// fout = fin + dt*rhs, then BCs. Returns the stage's stable dt (with the
// same 1%-overrun tolerance the historical joint stepper applied); the
// caller fails the step if it is below the requested dt.
static double
species_rk3_stage(gkyl_moment_app *app, struct moment_species *sp,
  double tcurr, double dt, const struct gkyl_array *fin, struct gkyl_array *fout)
{
  app->stat.nfeuler += 1;

  double dtmin = sp->rhs_func(app, sp, fin, fout);
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < 0.01)
    dtmin = dt; // avoid retaking steps on very small dt differences

  if (dtmin < dt)
    return dtmin;

  gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dt, &(app->local)),
    1.0, fin, &(app->local));
  moment_species_apply_bc(app, tcurr, sp, fout);

  return dtmin;
}

// A full SSP-RK3 step for this species, built on rhs_func. f0 holds the
// pre-step state and is only overwritten by the final commit, so a stage
// whose stable dt is below the requested dt simply reports failure and the
// stepper's redo protocol retakes the whole Strang step at the suggested dt
// (the historical joint stepper instead restarted its stages internally --
// same trajectory, different bookkeeping).
static struct gkyl_update_status
species_update_ssp_rk3(gkyl_moment_app *app, struct moment_species *sp,
  double tcurr, double dt)
{
  // Stage 1: f1 = f0 + dt*L(f0).
  double dt1 = species_rk3_stage(app, sp, tcurr, dt, sp->f0, sp->f1);
  if (dt1 < dt)
    return (struct gkyl_update_status) { .success = false, .dt_suggested = dt1 };

  // Stage 2: fnew = f1 + dt*L(f1); f1 = 3/4 f0 + 1/4 fnew.
  double dt2 = species_rk3_stage(app, sp, tcurr+dt, dt, sp->f1, sp->fnew);
  if (dt2 < dt) {
    double dt_rel_diff = (dt-dt2)/dt2;
    app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0], dt_rel_diff);
    app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1], dt_rel_diff);
    app->stat.nstage_2_fail += 1;
    return (struct gkyl_update_status) { .success = false, .dt_suggested = dt2 };
  }
  array_combine(sp->f1, 3.0/4.0, sp->f0, 1.0/4.0, sp->fnew, &app->local_ext);

  // Stage 3: fnew = f1 + dt*L(f1); f1 = 1/3 f0 + 2/3 fnew; commit f0 = f1.
  double dt3 = species_rk3_stage(app, sp, tcurr+dt/2, dt, sp->f1, sp->fnew);
  if (dt3 < dt) {
    double dt_rel_diff = (dt-dt3)/dt3;
    app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0], dt_rel_diff);
    app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1], dt_rel_diff);
    app->stat.nstage_3_fail += 1;
    return (struct gkyl_update_status) { .success = false, .dt_suggested = dt3 };
  }
  array_combine(sp->f1, 1.0/3.0, sp->f0, 2.0/3.0, sp->fnew, &app->local_ext);
  gkyl_array_copy_range(sp->f0, sp->f1, &app->local_ext);

  return (struct gkyl_update_status) { .success = true, .dt_suggested = dt3 };
}


// update solution: initial solution is in sp->f[0] and updated
// solution in sp->f[ndim]
struct gkyl_update_status
moment_species_update(gkyl_moment_app *app,
  struct moment_species *sp, double tcurr, double dt)
{
  return sp->update_func(app, sp, tcurr, dt);
}

// ---- post-update hook --------------------------------------------------------

static void
species_no_post_update(gkyl_moment_app *app, struct moment_species *sp,
  struct gkyl_array *f) { }

// Scrub hydro to zero inside the excision region for modular GR species.
// The self-contained formulation does the same in gr_euler_impose_gauge,
// which zeroes the state in excised cells each step. Without this scrub,
// excised cells whose IC projection averaged in non-zero hydro from
// adjacent non-excised quadrature points would retain that hydro forever.
static void
species_post_update_scrub_excised(gkyl_moment_app *app, struct moment_species *sp,
  struct gkyl_array *f)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&app->local, iter.idx);
    const double *prods_row = gkyl_array_cfetch(app->spacetime.prods, cidx);
    if (prods_row[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
      double *q = gkyl_array_fetch(f, cidx);
      for (int k = 0; k < sp->num_equations; k++) q[k] = 0.0;
    }
  }
}

// ---- rhs (MP/KEP path) -------------------------------------------------------

static double
species_rhs_mp(gkyl_moment_app *app, struct moment_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec tm = gkyl_wall_clock();

  gkyl_array_clear(sp->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  gkyl_mp_scheme_advance(sp->mp_slvr, &app->local, fin,
    app->ql, app->qr, app->amdq, app->apdq,
    sp->cflrate, sp->embed_mask, rhs);

  double omegaCfl[1];
  gkyl_array_reduce_range(omegaCfl, sp->cflrate, GKYL_MAX, &(app->local));

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);

  return sp->cfl/omegaCfl[0];
}

static double
species_rhs_kep(gkyl_moment_app *app, struct moment_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec tm = gkyl_wall_clock();

  gkyl_array_clear(sp->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  gkyl_kep_scheme_advance(sp->kep_slvr, &app->local, fin, sp->alpha,
    sp->cflrate, rhs);

  double omegaCfl[1];
  gkyl_array_reduce_range(omegaCfl, sp->cflrate, GKYL_MAX, &(app->local));

  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);

  return sp->cfl/omegaCfl[0];
}

// Static species: zero RHS (forward Euler then gives fout = fin), no CFL
// constraint of its own.
static double
species_rhs_static(gkyl_moment_app *app, struct moment_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_clear(rhs, 0.0);
  return DBL_MAX;
}


// Compute RHS of moment equations
double
moment_species_rhs(gkyl_moment_app *app, struct moment_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  return species->rhs_func(app, species, fin, rhs);
}

// ---- copy (time-stepper backup/commit/restore) -------------------------------

static void
species_copy(const struct moment_species *sp,
  struct gkyl_array *dst, const struct gkyl_array *src)
{
  gkyl_array_copy(dst, src);
}

void
moment_species_copy(const struct moment_species *sp,
  struct gkyl_array *dst, const struct gkyl_array *src)
{
  sp->copy_func(sp, dst, src);
}

// ---- stepper protocol (backup/restore/stage-state/commit) ---------------------
// Hides each scheme's state-array layout from the stepper and the source
// coupling. The pre-step backup/restore work through fcurr (the current
// solution, whatever the scheme), so they need no per-scheme variants; the
// stage-state accessor and the commit do.

void
moment_species_step_backup(const struct moment_species *sp)
{
  sp->copy_func(sp, sp->fdup, sp->fcurr);
}

void
moment_species_step_restore(const struct moment_species *sp)
{
  sp->copy_func(sp, sp->fcurr, sp->fdup);
}

static struct gkyl_array*
species_stage_state_wave_prop(const struct moment_species *sp, int nstrang)
{
  return nstrang == 0 ? sp->f[0] : sp->f[sp->ndim];
}

static void
species_step_commit_wave_prop(const struct moment_species *sp)
{
  sp->copy_func(sp, sp->f[0], sp->f[sp->ndim]);
}

// In-place set (RK schemes and static species): the current solution lives
// in fcurr through the whole step, so both stage states are fcurr and there
// is nothing to commit (the SSP-RK3 update commits into fcurr itself).
static struct gkyl_array*
species_stage_state_curr(const struct moment_species *sp, int nstrang)
{
  return sp->fcurr;
}

static void
species_step_commit_curr(const struct moment_species *sp)
{
}

struct gkyl_array*
moment_species_stage_state(const struct moment_species *sp, int nstrang)
{
  return sp->stage_state_func(sp, nstrang);
}

void
moment_species_step_commit(const struct moment_species *sp)
{
  sp->step_commit_func(sp);
}

// ---- initial conditions --------------------------------------------------------

static void
species_apply_ic(gkyl_moment_app *app, struct moment_species *sp, double t0)
{
  int num_quad = app->scheme_type == GKYL_MOMENT_MP ? 4 : 2;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, num_quad, sp->num_equations,
    sp->init, sp->ctx);

  gkyl_fv_proj_advance(proj, t0, &app->local, sp->fcurr);
  gkyl_fv_proj_release(proj);

  if (sp->has_app_accel) {
    gkyl_fv_proj_advance(sp->app_accel_proj, t0, &app->local, sp->app_accel);
  }

  moment_species_apply_bc(app, t0, sp, sp->fcurr);
}

void
moment_species_apply_ic(gkyl_moment_app *app, struct moment_species *sp, double t0)
{
  sp->apply_ic_func(app, sp, t0);
}

// ---- write (frame IO) ------------------------------------------------------------

static void
species_write(const gkyl_moment_app *app, const struct moment_species *sp,
  double tm, int frame)
{
  struct gkyl_msgpack_data *mt = moment_array_meta_new( (struct moment_output_meta) {
      .frame = frame,
      .stime = tm
    }
  );

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, sp->name, frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, sp->fcurr, fileNm.str);
  cstr_drop(&fileNm);

  if (sp->scheme_type == GKYL_MOMENT_KEP) {
    cstr fileNm = cstr_from_fmt("%s-%s-alpha_%d.gkyl", app->name, sp->name, frame);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, sp->alpha, fileNm.str);
    cstr_drop(&fileNm);
  }

  if (sp->has_app_accel) {
    if (sp->app_accel_evolve || frame == 0) {
      cstr fileNm = cstr_from_fmt("%s-%s-app_accel_%d.gkyl", app->name, sp->name, frame);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, sp->app_accel, fileNm.str);
      cstr_drop(&fileNm);
    }
  }

  moment_array_meta_release(mt);
}

void
moment_species_write(const gkyl_moment_app *app, const struct moment_species *sp,
  double tm, int frame)
{
  sp->write_func(app, sp, tm, frame);
}

// ---- integrated diagnostics ------------------------------------------------------

static void
species_calc_integrated_mom(gkyl_moment_app *app, struct moment_species *sp,
  double tm)
{
  int num_diag = sp->equation->num_diag;
  double q_integ[num_diag];

  calc_integ_quant(sp->equation, app->grid.cellVolume, sp->fcurr, app->geom,
    app->local, q_integ);

  double q_integ_global[num_diag];
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_diag, q_integ, q_integ_global);
  gkyl_dynvec_append(sp->integ_q, tm, q_integ_global);
}

static void
species_write_integrated_mom(gkyl_moment_app *app, struct moment_species *sp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // write out diagnostic moments
    cstr fileNm = cstr_from_fmt("%s-%s-%s.gkyl", app->name, sp->name, "imom");

    if (sp->is_first_q_write_call) {
      gkyl_dynvec_write(sp->integ_q, fileNm.str);
      sp->is_first_q_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(sp->integ_q, fileNm.str);
    }
    cstr_drop(&fileNm);
  }
  gkyl_dynvec_clear(sp->integ_q);
}

void
moment_species_calc_integrated_mom(gkyl_moment_app *app, struct moment_species *sp,
  double tm)
{
  sp->calc_integrated_mom_func(app, sp, tm);
}

void
moment_species_write_integrated_mom(gkyl_moment_app *app, struct moment_species *sp)
{
  sp->write_integrated_mom_func(app, sp);
}

// ---- restart read ------------------------------------------------------------------

static struct gkyl_app_restart_status
species_from_file(gkyl_moment_app *app, struct moment_species *sp,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = moment_app_header_from_file(app, fname);

  if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, sp->fcurr, fname);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      moment_species_apply_bc(app, rstat.stime, sp, sp->fcurr);
    }
  }

  // Compute applied acceleration if present.
  // Computation necessary in case applied acceleration
  // is time-independent and not computed in the time-stepping loop
  // since it is not read-in as part of restarts.
  if (sp->has_app_accel) {
    gkyl_fv_proj_advance(sp->app_accel_proj, rstat.stime, &app->local, sp->app_accel);
  }

  return rstat;
}

struct gkyl_app_restart_status
moment_species_from_file(gkyl_moment_app *app, struct moment_species *sp,
  const char *fname)
{
  return sp->read_func(app, sp, fname);
}

// ---- release -----------------------------------------------------------------

// Release the state/solver memory specific to each scheme; the shared
// members are released by moment_species_release below.
static void
species_release_wave_prop(const struct moment_species *sp)
{
  for (int d=0; d<sp->ndim; ++d)
    gkyl_wave_prop_release(sp->slvr[d]);

  gkyl_array_release(sp->fdup);
  for (int d=0; d<sp->ndim+1; ++d)
    gkyl_array_release(sp->f[d]);
}

static void
species_release_mp(const struct moment_species *sp)
{
  gkyl_mp_scheme_release(sp->mp_slvr);
  gkyl_array_release(sp->fdup);
  gkyl_array_release(sp->f0);
  gkyl_array_release(sp->f1);
  gkyl_array_release(sp->fnew);
  gkyl_array_release(sp->cflrate);
  gkyl_array_release(sp->alpha);
}

static void
species_release_kep(const struct moment_species *sp)
{
  gkyl_kep_scheme_release(sp->kep_slvr);
  gkyl_array_release(sp->fdup);
  gkyl_array_release(sp->f0);
  gkyl_array_release(sp->f1);
  gkyl_array_release(sp->fnew);
  gkyl_array_release(sp->cflrate);
  gkyl_array_release(sp->alpha);
}

// free species
void
moment_species_release(const struct moment_species *sp)
{
  gkyl_wv_eqn_release(sp->equation);

  for (int d=0; d<sp->ndim; ++d) {
    if (sp->lower_bc[d])
      gkyl_wv_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d])
      gkyl_wv_apply_bc_release(sp->upper_bc[d]);
  }

  sp->release_func(sp);

  gkyl_array_release(sp->app_accel);
  if (sp->has_app_accel) {
    gkyl_fv_proj_release(sp->app_accel_proj);
  }

  gkyl_array_release(sp->nT_source);
  if (sp->proj_nT_source) {
    gkyl_fv_proj_release(sp->proj_nT_source);
  }

  gkyl_array_release(sp->embed_mask);

  gkyl_array_release(sp->bc_buffer);

  gkyl_dynvec_release(sp->integ_q);
}

/** mhd_src functions */

void
mhd_src_init(const struct gkyl_moment_app *app,
             const struct gkyl_moment_species *sp,
             struct mhd_src *src)
{
  double dxyz_min = DBL_MAX;
  for (int d=0; d<app->grid.ndim; ++d) {
    double dx = app->grid.dx[d];
    dxyz_min = dx < dxyz_min ? dx : dxyz_min;
  }

  struct gkyl_mhd_src_inp src_inp = {
    .grid = &app->grid,
    .divergence_constraint = gkyl_wv_mhd_divergence_constraint(sp->equation),
    .glm_ch = gkyl_wv_mhd_glm_ch(sp->equation),
    .glm_alpha = gkyl_wv_mhd_glm_ch(sp->equation),
    .dxyz_min = dxyz_min,
  };

  src->slvr = gkyl_mhd_src_new(src_inp, &app->local_ext);
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
void
mhd_src_update(gkyl_moment_app *app, struct mhd_src *src, int nstrang,
               double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  int i = 0; // mhd has only one 'species'
  struct gkyl_array *fluid = app->species[i].f[sidx[nstrang]];

  if (app->species[i].app_accel_proj)
    gkyl_fv_proj_advance(app->species[i].app_accel_proj, tcurr, &app->local,
                         app->species[i].app_accel);

  // FIXME presently needed for computing divB etc
  moment_species_apply_bc(app, tcurr, &app->species[i], fluid);

  gkyl_mhd_src_advance(src->slvr, dt, &app->local, fluid,
                       app->species[i].app_accel);

  moment_species_apply_bc(app, tcurr, &app->species[i], fluid);

}

void
mhd_src_release(const struct mhd_src *src)
{
  gkyl_mhd_src_release(src->slvr);
}
