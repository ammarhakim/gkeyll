#include <gkyl_moment_priv.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_iso_euler.h>

// initialize source solver: this should be called after all species
// and fields are initialized
void
moment_coupling_init(const struct gkyl_moment_app *app, struct moment_coupling *src)
{
  struct gkyl_moment_em_coupling_inp src_inp = {
    .grid = &app->grid,
    .nfluids = app->num_species,
    // if there is a field, need to update electric field too, otherwise just updating fluid
    .epsilon0 = app->field.epsilon0 ? app->field.epsilon0 : 0.0, 
    .mu0 = app->field.mu0 ? app->field.mu0 : 0.0, 
    // is the field static?
    .static_field = app->field.is_static, 
    // linear ramping function for slowing turning on applied accelerations, E fields, or currents
    .t_ramp_E = app->field.t_ramp_E ? app->field.t_ramp_E : 0.0,
    .t_ramp_curr = app->field.t_ramp_curr ? app->field.t_ramp_curr : 0.0,
  };

  for (int i=0; i<app->num_species; ++i)
    src_inp.param[i] = (struct gkyl_moment_em_coupling_data) {
      .type = app->species[i].eqn_type,
      .charge = app->species[i].charge,
      .mass = app->species[i].mass,
      // The gradient-based closure defines its heat flux through k0, so k0=0.0 in the source solve to avoid double-applying it.
      // The neural-network closure supplies the heat flux directly, so k0 is retained here as the integrating-factor relaxation rate.
      .k0 = (app->species[i].has_grad_closure) ? 0.0 : app->species[i].k0,
    };

  src_inp.has_collision = app->has_collision;
  for (int s=0; s<app->num_species; ++s)
    for (int r=0; r<app->num_species; ++r)
      src_inp.nu_base[s][r] = app->nu_base[s][r];

  src_inp.has_nT_sources = false;
  for (int i=0; i<app->num_species; ++i)
    if (app->species[i].proj_nT_source)
      src_inp.has_nT_sources = true;

  // check for relativistic-species
  bool use_rel = 0;
  for (int n=0;  n<app->num_species; ++n) 
    use_rel = use_rel || (app->species[n].eqn_type == GKYL_EQN_COLDFLUID_SR) || app->field.use_explicit_em_coupling;

  // Merge the per-species source-family bundles into the solver input. A
  // family is enabled if any species enables it; preserving the historical
  // merge, its parameters are taken from the last enabling species that set
  // any of them (the vacuum-Einstein families take the last enabling
  // species' parameters unconditionally).
  for (int i = 0; i < app->num_species; i++) {
    const struct moment_species *msp = &app->species[i];

    if (msp->friction.enabled) {
      bool use_explicit = src_inp.friction.use_explicit || msp->friction.use_explicit;
      if (msp->friction.Z != 0.0 || msp->friction.T_elc != 0.0 || msp->friction.Lambda_ee != 0.0)
        src_inp.friction = msp->friction;
      src_inp.friction.enabled = true;
      src_inp.friction.use_explicit = use_explicit;
    }

    if (msp->volume_sources.enabled) {
      if (msp->volume_sources.gas_gamma != 0.0 || msp->volume_sources.U0 != 0.0 || msp->volume_sources.R0 != 0.0)
        src_inp.volume_sources = msp->volume_sources;
      src_inp.volume_sources.enabled = true;
    }

    if (msp->reactivity.enabled) {
      if (msp->reactivity.gas_gamma != 0.0 || msp->reactivity.specific_heat_capacity != 0.0 ||
        msp->reactivity.energy_of_formation != 0.0 || msp->reactivity.ignition_temperature != 0.0 ||
        msp->reactivity.reaction_rate != 0.0)
        src_inp.reactivity = msp->reactivity;
      src_inp.reactivity.enabled = true;
    }

    if (msp->einstein_medium.enabled) {
      if (msp->einstein_medium.gas_gamma != 0.0 || msp->einstein_medium.kappa != 0.0)
        src_inp.einstein_medium = msp->einstein_medium;
      src_inp.einstein_medium.enabled = true;
    }

    if (msp->gr_ultra_rel.enabled) {
      if (msp->gr_ultra_rel.gas_gamma != 0.0)
        src_inp.gr_ultra_rel = msp->gr_ultra_rel;
      src_inp.gr_ultra_rel.enabled = true;
    }

    if (msp->gr_euler.enabled) {
      if (msp->gr_euler.gas_gamma != 0.0)
        src_inp.gr_euler = msp->gr_euler;
      src_inp.gr_euler.enabled = true;
    }

    if (msp->gr_twofluid.enabled) {
      if (msp->gr_twofluid.mass_elc != 0.0 || msp->gr_twofluid.mass_ion != 0.0 ||
        msp->gr_twofluid.charge_elc != 0.0 || msp->gr_twofluid.charge_ion != 0.0 ||
        msp->gr_twofluid.gas_gamma_elc != 0.0 || msp->gr_twofluid.gas_gamma_ion != 0.0 ||
        msp->gr_twofluid.e_fact != 0.0)
        src_inp.gr_twofluid = msp->gr_twofluid;
      src_inp.gr_twofluid.enabled = true;
    }

    if (msp->gr_mhd.enabled) {
      if (msp->gr_mhd.gas_gamma != 0.0)
        src_inp.gr_mhd = msp->gr_mhd;
      src_inp.gr_mhd.enabled = true;
    }

    if (msp->vacuum_einstein.enabled)
      src_inp.vacuum_einstein = msp->vacuum_einstein;
    if (msp->vacuum_einstein_conformal.enabled)
      src_inp.vacuum_einstein_conformal = msp->vacuum_einstein_conformal;
  }

  // The field can also carry volume-based geometrical sources.
  if (app->field.volume_sources.enabled) {
    if (app->field.volume_sources.gas_gamma != 0.0 || app->field.volume_sources.U0 != 0.0 ||
      app->field.volume_sources.R0 != 0.0)
      src_inp.volume_sources = app->field.volume_sources;
    src_inp.volume_sources.enabled = true;
  }

  // Explicit (special-)relativistic multi-fluid + Maxwell coupling on modular
  // GR-Euler fluids. The per-species adiabatic index comes from the EOS bundle
  // carried on the tetrad equation object (the same closure the wave step uses).
  src_inp.has_gr_em_coupling = app->field.use_gr_em_coupling;
  if (src_inp.has_gr_em_coupling) {
    for (int i = 0; i < app->num_species; i++) {
      if (app->species[i].eqn_type == GKYL_EQN_GR_EULER_TETRAD) {
        src_inp.gr_em_eos[i] = gkyl_wv_gr_euler_tetrad_eos(app->species[i].equation);
      }
      else {
        src_inp.gr_em_eos[i] = gkyl_gr_euler_eos_ideal(0.0);
      }
    }
  }

  // save the use-rel bool
  src_inp.use_rel = use_rel;

  // check for explicit em-coupling
  src_inp.use_explicit_em_coupling = 0;
  if (app->field.use_explicit_em_coupling)
    src_inp.use_explicit_em_coupling = 1;

  // create updater to solve for sources
  src->slvr = gkyl_moment_em_coupling_new(src_inp);

  // Several source paths (closures, Braginskii, volume sources, implicit
  // friction) read the EM/external-EM/applied-current arrays
  // unconditionally; when the app has no field, hand them zero scratch
  // arrays instead of the (nonexistent) field state. Allocated only when
  // the coupling exists without a field -- a source-free fluid run builds
  // no coupling and allocates nothing.
  src->nofield_em = NULL;
  src->nofield_ext_em = NULL;
  src->nofield_app_current = NULL;
  if (!app->has_field) {
    src->nofield_em = mkarr(false, 8, app->local_ext.volume);
    gkyl_array_clear(src->nofield_em, 0.0);
    src->nofield_ext_em = mkarr(false, 6, app->local_ext.volume);
    gkyl_array_clear(src->nofield_ext_em, 0.0);
    src->nofield_app_current = mkarr(false, 3, app->local_ext.volume);
    gkyl_array_clear(src->nofield_app_current, 0.0);
  }

  // Spacetime coupling for modular GR fluids. Only constructed when at
  // least one species is a modular GR variant; otherwise spacetime_slvr stays
  // NULL and the spacetime path costs nothing.
  src->spacetime_slvr = NULL;
  for (int i = 0; i < GKYL_MAX_SPECIES; i++) {
    src->gr_euler_prim_status_wave_prop[i]   = NULL;
    src->gr_euler_repair_status_wave_prop[i] = NULL;
    src->gr_euler_prim_status_source[i]      = NULL;
    src->gr_euler_repair_status_source[i]    = NULL;
  }
  if (app->has_spacetime) {
    bool any_mod = false;
    for (int i = 0; i < app->num_species; i++) {
      if (app->species[i].eqn_type == GKYL_EQN_GR_EULER_TETRAD) {
        any_mod = true;
        break;
      }
    }
    // Build the spacetime coupling when there are modular GR fluids feeling
    // gravity, OR a dynamic Einstein backend whose state needs the vacuum-
    // Einstein source — even with no fluid (the pure-spacetime case, analogous
    // to running moment_em_coupling for an applied current with no fluid).
    bool dynamic_einstein = (app->spacetime.einstein_eqn != NULL);
    if (any_mod || dynamic_einstein) {
      struct gkyl_moment_spacetime_coupling_inp st_inp = {
        .grid               = &app->grid,
        .nfluids            = app->num_species,
        .is_static          = app->spacetime.is_static,
        .has_tetrad         = app->spacetime.has_tetrad,
        .analytic_spacetime = app->spacetime.analytic_spacetime,
        .spacetime_gauge    = app->spacetime.spacetime_gauge,
        .reinit_freq        = app->spacetime.reinit_freq,
        .einstein_eqn       = app->spacetime.einstein_eqn,
      };
      for (int i = 0; i < app->num_species; i++) {
        enum gkyl_eqn_type t = app->species[i].eqn_type;
        if (t == GKYL_EQN_GR_EULER_TETRAD) {
          // Modular GR tetrad species: full EOS bundle (IDEAL or
          // APPROXIMATE_SYNGE)
          // lives on the equation object — read it back so the
          // source-step uses the same closure the wave-step does.
          st_inp.fluid_param[i] = (struct gkyl_moment_spacetime_coupling_data) {
            .type = t,
            .eos = gkyl_wv_gr_euler_tetrad_eos(app->species[i].equation),
          };
          st_inp.eqn[i] = app->species[i].equation;
        } else {
          // Non-modular species: leave the entry in a benign-default state;
          // the coupling's advance loop skips non-modular types.
          st_inp.fluid_param[i] = (struct gkyl_moment_spacetime_coupling_data) {
            .type = t,
            .eos = gkyl_gr_euler_eos_ideal(0.0),
          };
          st_inp.eqn[i] = NULL;
        }
      }
      src->spacetime_slvr = gkyl_moment_spacetime_coupling_new(st_inp);

      // Allocate per-species GR Euler tetrad instrumentation buckets:
      // one wave-prop side prim+repair status and one source side prim+
      // repair status. Attach wave-prop + source-repair pointers to the
      // equation's auxfields so the equation's wave-prop hooks and
      // repair_state callback can write into them.
      for (int i = 0; i < app->num_species; i++) {
        if (app->species[i].eqn_type != GKYL_EQN_GR_EULER_TETRAD)
          continue;
        src->gr_euler_prim_status_wave_prop[i] =
          gkyl_malloc(sizeof(struct gkyl_gr_euler_prim_status));
        src->gr_euler_repair_status_wave_prop[i] =
          gkyl_malloc(sizeof(struct gkyl_gr_euler_repair_status));
        src->gr_euler_prim_status_source[i] =
          gkyl_malloc(sizeof(struct gkyl_gr_euler_prim_status));
        src->gr_euler_repair_status_source[i] =
          gkyl_malloc(sizeof(struct gkyl_gr_euler_repair_status));
        memset(src->gr_euler_prim_status_wave_prop[i], 0,
          sizeof(struct gkyl_gr_euler_prim_status));
        memset(src->gr_euler_repair_status_wave_prop[i], 0,
          sizeof(struct gkyl_gr_euler_repair_status));
        memset(src->gr_euler_prim_status_source[i], 0,
          sizeof(struct gkyl_gr_euler_prim_status));
        memset(src->gr_euler_repair_status_source[i], 0,
          sizeof(struct gkyl_gr_euler_repair_status));

        // The single wiring point for the tetrad equation's cross-object
        // pointers: the products array and tetrad cache (both created by
        // gkyl_moment_app_new before this runs) plus the instrumentation
        // buckets above.
        gkyl_gr_euler_tetrad_set_auxfields(app->species[i].equation,
          (struct gkyl_wv_gr_euler_tetrad_auxfields){
            .prods                   = app->spacetime.prods,
            .wave_spacetime          = app->spacetime.wave_spacetime,
            .prim_status_wave_prop   = src->gr_euler_prim_status_wave_prop[i],
            .repair_status_wave_prop = src->gr_euler_repair_status_wave_prop[i],
            .prim_status_source      = src->gr_euler_prim_status_source[i],
            .repair_status_source    = src->gr_euler_repair_status_source[i],
          });
      }
    }
  }

  for (int n=0; n<app->num_species; ++n) {
    int meqn = app->species[n].num_equations;
    src->pr_rhs[n] = mkarr(false, meqn, app->local_ext.volume);
    src->non_ideal_cflrate[n] = mkarr(false, 1, app->local_ext.volume);
  }

  int ghost[3] = { 1, 1, 1 };
  // create non-ideal local extended range from local range
  // has one additional cell in each direction because non-ideal
  // variables are stored at cell vertices
  gkyl_create_vertex_ranges(&app->local, ghost, &src->non_ideal_local_ext,
    &src->non_ideal_local);

  // check if gradient-closure is present
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      int nadj[3] = { 1, 4, 8 }; // cells adjacent to a vertex
      src->non_ideal_vars[i] = mkarr(false, nadj[app->ndim - 1]*10,
        src->non_ideal_local_ext.volume);
      struct gkyl_ten_moment_grad_closure_inp grad_closure_inp = {
        .grid = &app->grid,
        .k0 = app->species[i].k0,
        .cfl = app->cfl,
        .comm = app->comm,
        .update_range = &app->local,
        .heat_flux_range = &src->non_ideal_local,
      };
      src->grad_closure_slvr[i] = gkyl_ten_moment_grad_closure_new(&grad_closure_inp);
    }
  }

  // Check whether neural network-based closure is present.
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_nn_closure) {
      struct gkyl_ten_moment_nn_closure_inp nn_closure_inp = {
        .grid = &app->grid,
        .k0 = app->species[i].k0,
        .poly_order = app->species[i].poly_order,
        .ann = app->species[i].ann,
      };
      src->nn_closure_slvr[i] = gkyl_ten_moment_nn_closure_new(nn_closure_inp);
    }
  }

  // check if braginskii terms are present
  if (app->has_braginskii) {
    struct gkyl_moment_braginskii_inp brag_inp = {
      .grid = &app->grid,
      .nfluids = app->num_species,
      .epsilon0 = app->field.epsilon0,
      // Check for multiplicative collisionality factor, default is 1.0
      .coll_fac = app->coll_fac == 0 ? 1.0 : app->coll_fac,
    };
    for (int i=0; i<app->num_species; ++i) {
      // Braginskii coefficients depend on pressure and coefficient to obtain
      // pressure is different for different equation systems (gasGamma, vt, Tr(P))
      double p_fac = 1.0;
      if (app->species[i].eqn_type == GKYL_EQN_EULER) {
        p_fac =  gkyl_wv_euler_gas_gamma(app->species[i].equation);
      }
      else if (app->species[i].eqn_type == GKYL_EQN_ISO_EULER) {
        p_fac =  gkyl_wv_iso_euler_vt(app->species[i].equation);
      }
      brag_inp.param[i] = (struct gkyl_moment_braginskii_data) {
        .type_eqn = app->species[i].eqn_type,
        .type_brag = app->species[i].type_brag,
        .charge = app->species[i].charge,
        .mass = app->species[i].mass,
        .p_fac = p_fac,
      };
    }
    src->brag_slvr = gkyl_moment_braginskii_new(brag_inp);
  }
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
struct gkyl_update_status
moment_coupling_update(gkyl_moment_app *app, struct moment_coupling *src,
  int nstrang, double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];
  const struct gkyl_array *pr_rhs_const[GKYL_MAX_SPECIES];
  const struct gkyl_array *nT_sources[GKYL_MAX_SPECIES];

  // Field state for this Strang half-step, through the field's stepper
  // protocol (the accessor hides the per-scheme state layout; it returns
  // NULL without a field, where the zero scratch arrays stand in because
  // several source paths read these arrays unconditionally).
  struct gkyl_array *em_state = app->has_field
    ? moment_field_stage_state(&app->field, nstrang) : src->nofield_em;
  const struct gkyl_array *em_closures = em_state;
  const struct gkyl_array *ext_em =
    app->has_field ? app->field.ext_em : src->nofield_ext_em;
  const struct gkyl_array *app_current =
    app->has_field ? app->field.app_current : src->nofield_app_current;

  double dt_suggested = DBL_MAX;
  struct gkyl_ten_moment_grad_closure_status stat;

  for (int i=0; i<app->num_species; ++i) {
    fluids[i] = moment_species_stage_state(&app->species[i], nstrang);

    if (app->species[i].app_accel_evolve) {
      gkyl_fv_proj_advance(app->species[i].app_accel_proj, tcurr, &app->local, app->species[i].app_accel);
    }
    app_accels[i] = app->species[i].app_accel;

    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      // Non-ideal variables are defined on an extended range with one additional "cell" in each direction.
      // This additional cell accounts for the fact that non-ideal variables are stored at cell vertices.
      stat = gkyl_ten_moment_grad_closure_advance(src->grad_closure_slvr[i],
        &src->non_ideal_local, &app->local,
        fluids[i], em_closures,
        src->non_ideal_cflrate[i], dt, src->non_ideal_vars[i], src->pr_rhs[i]);

      if (!stat.success)
        return (struct gkyl_update_status) {
          .success = false,
          .dt_suggested = stat.dt_suggested
        };

      dt_suggested = fmin(dt_suggested, stat.dt_suggested);
    }

    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_nn_closure) {
      // Non-ideal variables are defined on an extended range with one additional "cell" in each direction.
      // This additional cell accounts for the fact that non-ideal variables are stored at cell vertices.
      gkyl_ten_moment_nn_closure_advance(src->nn_closure_slvr[i], &src->non_ideal_local, &app->local, fluids[i],
        em_closures, src->non_ideal_vars[i], src->pr_rhs[i]);
    }
  }

  if (app->has_braginskii) {
    gkyl_moment_braginskii_advance(src->brag_slvr,
      src->non_ideal_local, app->local,
      fluids, em_closures,
      src->non_ideal_cflrate, src->non_ideal_vars, src->pr_rhs);
  }

  if (app->field.ext_em_evolve) {
    gkyl_fv_proj_advance(app->field.ext_em_proj, tcurr, &app->local, app->field.ext_em);
  }

  if (app->field.app_current_evolve) {
    if (app->field.use_explicit_em_coupling) {
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr, &app->local, app->field.app_current);
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr + dt*2, &app->local, app->field.app_current1);
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr + 2*dt/2.0, &app->local, app->field.app_current2);
    }
    else {
      gkyl_fv_proj_advance(app->field.app_current_proj, tcurr, &app->local, app->field.app_current);
    }
  }

  // Get the RHS pointer for accumulation during source update
  for (int i=0; i<app->num_species; ++i) {
    pr_rhs_const[i] = src->pr_rhs[i];
  }

  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].proj_nT_source
        && !(app->species[i].nT_source_set_only_once
             && app->species[i].nT_source_is_set))
    {
      gkyl_fv_proj_advance(app->species[i].proj_nT_source, tcurr,
                           &app->local, app->species[i].nT_source);
    }
    nT_sources[i] = app->species[i].nT_source;
    app->species[i].nT_source_is_set = true;
  }
  if (app->field.use_gr_em_coupling) {
    // Explicit SSP-RK3 relativistic multi-fluid + Maxwell coupling on modular
    // GR-Euler fluids (SI units). Conditionally stable: on under-resolution of
    // the plasma/cyclotron frequency it reports failure and a stable dt for the
    // app's UPDATE_REDO to retry with. The mutated sub-step is discarded from
    // the backup copy on redo.
    struct gkyl_moment_em_coupling_status gr_em_status =
      gkyl_moment_em_coupling_gr_em_explicit_advance(src->slvr, tcurr, dt, &app->local,
        fluids, em_state, app->field.ext_em);

    // Reduce locally computed CFL estimate from explicit sources across ranks. 
    double gr_red_local[2]  = { gr_em_status.success ? 1.0 : 0.0, gr_em_status.dt_suggested };
    double gr_red_global[2];
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 2, gr_red_local, gr_red_global);
    bool gr_em_success = gr_red_global[0] > 0.5;
    double gr_em_dt_suggested = gr_red_global[1];

    if (!gr_em_success) {
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = gr_em_dt_suggested
      };
    }
    dt_suggested = fmin(dt_suggested, gr_em_dt_suggested);
  }
  else if (app->field.use_explicit_em_coupling) {
    gkyl_moment_em_coupling_explicit_advance(src->slvr, tcurr, dt, &app->local,
      fluids, app_accels, pr_rhs_const,
      em_state, app->field.app_current, app->field.app_current1,
      app->field.app_current2, app->field.ext_em,
      nT_sources, app->field.app_current_proj,nstrang);
  }
  else {
    gkyl_moment_em_coupling_implicit_advance(src->slvr, tcurr, dt, &app->local,
      fluids, app_accels, pr_rhs_const,
      em_state, app_current, ext_em,
      nT_sources);
  }

  // Spacetime coupling: per-RK-stage GR source terms — gravity on the modular
  // GR fluids (reads the products array), and the vacuum-Einstein source on the
  // Einstein state for the dynamic backend (no-op for the static backend). The
  // Einstein state for this Strang half-step is f[0] (nstrang 0) or f[ndim]
  // (nstrang 1), matching the fluid Strang buffers; it is NULL for the static
  // backend, which leaves the Einstein source a no-op.
  if (src->spacetime_slvr) {
    gkyl_moment_spacetime_coupling_explicit_advance(src->spacetime_slvr,
      tcurr, dt, &app->local, fluids, app->spacetime.f[sidx[nstrang]]);
    // The source just modified the interior of the Einstein state; refresh its
    // ghost cells so the subsequent hyperbolic update (and the next Strang
    // half-step) read consistent boundary data. When the Einstein state is
    // evolved as a species (self-contained formulation) the per-species
    // apply_bc below covers this; the spacetime object has no species entry,
    // so apply its BC here. No-op for the static-analytic backend (no
    // Einstein state to evolve).
    moment_spacetime_apply_bc(app, tcurr, &app->spacetime,
      app->spacetime.f[sidx[nstrang]]);
  }

  for (int i=0; i<app->num_species; ++i) {
    moment_species_apply_bc(app, tcurr, &app->species[i], fluids[i]);
  }
  moment_field_apply_bc(app, tcurr, &app->field, em_state); // no-op without a field
  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

// free sources
void
moment_coupling_release(const struct gkyl_moment_app *app, const struct moment_coupling *src)
{
  gkyl_moment_em_coupling_release(src->slvr);
  if (src->nofield_em)
    gkyl_array_release(src->nofield_em);
  if (src->nofield_ext_em)
    gkyl_array_release(src->nofield_ext_em);
  if (src->nofield_app_current)
    gkyl_array_release(src->nofield_app_current);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_release(src->pr_rhs[i]);
    gkyl_array_release(src->non_ideal_cflrate[i]);
    gkyl_array_release(src->non_ideal_vars[i]);
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      gkyl_ten_moment_grad_closure_release(src->grad_closure_slvr[i]);
    }
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_nn_closure) {
      gkyl_ten_moment_nn_closure_release(src->nn_closure_slvr[i]);
    }
  }
  if (app->has_braginskii)
    gkyl_moment_braginskii_release(src->brag_slvr);
  if (src->spacetime_slvr)
    gkyl_moment_spacetime_coupling_release(src->spacetime_slvr);
  for (int i = 0; i < app->num_species; i++) {
    if (src->gr_euler_prim_status_wave_prop[i])
      gkyl_free(src->gr_euler_prim_status_wave_prop[i]);
    if (src->gr_euler_repair_status_wave_prop[i])
      gkyl_free(src->gr_euler_repair_status_wave_prop[i]);
    if (src->gr_euler_prim_status_source[i])
      gkyl_free(src->gr_euler_prim_status_source[i]);
    if (src->gr_euler_repair_status_source[i])
      gkyl_free(src->gr_euler_repair_status_source[i]);
  }
}

