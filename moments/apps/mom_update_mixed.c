#include <stdlib.h>

#include <gkyl_moment_priv.h>
#include <gkyl_sources_explicit_priv.h>

// Per-species mixed-scheme time step: 
// - the GR-Euler fluid is advanced by the MP (method-of-lines SSP-RK3) 
// - the vacuum-Einstein spacetime by the wave-propagation
// ,coupled each step.
//

// One SSP-RK3 forward-Euler stage of the spacetime source (h): gauge source (via the shared plus, if coupled, the frozen-fluid ADM matter source.
static void
spacetime_source_euler_stage(bool is_conformal, double excision_threshold,
  enum gkyl_spacetime_slicing slicing, enum gkyl_spacetime_evolution evolution,
  const double *qfluid, double gas_gamma, double tcurr, double h, int meqn,
  double *in, double *out)
{
  // Gauge source: out = in + h * S_gauge(in)
  if (is_conformal) {
    explicit_vacuum_einstein_conformal_source_update_euler(0, excision_threshold, slicing, evolution,
      tcurr, h, in, out);
  }
  else {
    explicit_vacuum_einstein_source_update_euler(0, excision_threshold, slicing, evolution,
      tcurr, h, in, out);
  }

  // Fluid -> spacetime matter source (only in live cells): out += h * M(in, fluid), M touching K_ij and V_k only
  if (qfluid && out[9] >= excision_threshold) {
    double mrhs[77] = { 0.0 };
    add_einstein_matter_source(in, qfluid, gas_gamma, excision_threshold, is_conformal, mrhs);
    for (int m = 0; m < meqn; ++m) {
      out[m] += h * mrhs[m];
    }
  }

  if (out[9] < excision_threshold) {
    for (int m = 0; m < meqn; ++m) {
      out[m] = 0.0;
    }
  }
}

// Integrate the spacetime source (gauge + optional frozen-fluid matter) over a half-step h, in place on the spacetime state array `f`, using SSP-RK3 with the per-cell forward-Euler stage above.
static void
spacetime_source_half_step(gkyl_moment_app *app, struct moment_species *einstein_sp,
  struct gkyl_array *f, const struct gkyl_array *fluid_f, double gas_gamma, double tcurr, double h)
{
  bool is_conformal = einstein_sp->has_vacuum_einstein_conformal;
  int meqn = einstein_sp->num_equations;

  double excision_threshold = is_conformal ? einstein_sp->vacuum_einstein_conformal_excision_threshold
    : einstein_sp->vacuum_einstein_excision_threshold;
  enum gkyl_spacetime_slicing slicing = is_conformal ? einstein_sp->vacuum_einstein_conformal_spacetime_slicing
    : einstein_sp->vacuum_einstein_spacetime_slicing;
  enum gkyl_spacetime_evolution evolution = is_conformal ? einstein_sp->vacuum_einstein_conformal_spacetime_evolution
    : einstein_sp->vacuum_einstein_spacetime_evolution;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double *q = gkyl_array_fetch(f, loc);
    const double *qfluid = fluid_f ? gkyl_array_cfetch(fluid_f, loc) : NULL;

    double f_old[77], f_stage[77], f_tmp[77];
    for (int m = 0; m < meqn; ++m) {
      f_old[m] = q[m];
    }

    // stage 1
    spacetime_source_euler_stage(is_conformal, excision_threshold, slicing, evolution, qfluid,
      gas_gamma, tcurr, h, meqn, f_old, f_stage);

    // stage 2: 3/4 old + 1/4 euler(stage1)
    spacetime_source_euler_stage(is_conformal, excision_threshold, slicing, evolution, qfluid,
      gas_gamma, tcurr + h, h, meqn, f_stage, f_tmp);
    for (int m = 0; m < meqn; ++m) {
      f_stage[m] = (0.75 * f_old[m]) + (0.25 * f_tmp[m]);
    }

    // stage 3: 1/3 old + 2/3 euler(stage2)
    spacetime_source_euler_stage(is_conformal, excision_threshold, slicing, evolution, qfluid,
      gas_gamma, tcurr + (0.5 * h), h, meqn, f_stage, f_tmp);
    for (int m = 0; m < meqn; ++m) {
      q[m] = ((1.0 / 3.0) * f_old[m]) + ((2.0 / 3.0) * f_tmp[m]);
    }
  }
}

// Copy the vacuum-Einstein metric into the fluid's metric slots
static void
mixed_sync_metric(gkyl_moment_app *app, int fluid_idx, int einstein_idx)
{
  bool is_conformal = app->species[einstein_idx].has_vacuum_einstein_conformal;
  double excision_threshold = is_conformal ? app->species[einstein_idx].vacuum_einstein_conformal_excision_threshold
    : app->species[einstein_idx].vacuum_einstein_excision_threshold;

  const struct gkyl_array *einstein_f = app->species[einstein_idx].fcurr;
  struct gkyl_array *fluid_f = app->species[fluid_idx].fcurr;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qe = gkyl_array_cfetch(einstein_f, loc);
    double *qf = gkyl_array_fetch(fluid_f, loc);
    sync_fluid_metric_from_einstein(qe, qf, excision_threshold, is_conformal);
  }
}

// Hyperbolicity guard: are all non-excised cell-average spacetime metrics valid inputs for the fluid solve?
// A cell is invalid if its spatial metric (slots 0-8) is (a) grossly asymmetric (numerical drift is fine - the SPD test below uses the symmetrized metric, but a large antisymmetric part signals breakdown); 
// (b) not SPD (positive leading principal minors of the symmetrized metric, scale-relative floors); 
// or, for the conformal system, (c) has a non-positive BSSN conformal factor chi (slot 64). 
static bool
mixed_spacetime_metric_valid(gkyl_moment_app *app, const struct moment_species *esp)
{
  bool is_conformal = esp->has_vacuum_einstein_conformal;
  double excision_threshold;
  if (is_conformal) {
    excision_threshold = esp->vacuum_einstein_conformal_excision_threshold;
  }
  else {
    excision_threshold = esp->vacuum_einstein_excision_threshold;
  }

  const struct gkyl_array *ef = esp->fcurr;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  int local_bad = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *e = gkyl_array_cfetch(ef, loc);
    if (e[9] < excision_threshold) {
      continue; // excised cell: handled by excision, not the SPD guard
    }
    double scale = fmax(fmax(fabs(e[0]), fabs(e[4])), fmax(fabs(e[8]), 1.0e-300));
    // (a) reject only gross asymmetry
    double asym = fmax(fmax(fabs(e[1] - e[3]), fabs(e[2] - e[6])), fabs(e[5] - e[7]));
    if (asym > 1.0e-2 * scale) {
      local_bad = 1;
      break;
    }
    // (b) SPD via leading principal minors (Sylvester) on the symmetrized metric (Sylvester needs symmetry), scale-relative floors (minors scale as scale^k)
    double g00 = e[0], g11 = e[4], g22 = e[8];
    double g01 = 0.5 * (e[1] + e[3]), g02 = 0.5 * (e[2] + e[6]), g12 = 0.5 * (e[5] + e[7]);
    double ep = 1.0e-12;
    double m1 = g00;
    double m2 = (g00 * g11) - (g01 * g01);
    double m3 = (g00 * ((g11 * g22) - (g12 * g12))) - (g01 * ((g01 * g22) - (g12 * g02))) +
      (g02 * ((g01 * g12) - (g11 * g02)));
    if (!(m1 > ep * scale && m2 > ep * scale * scale && m3 > ep * scale * scale * scale)) {
      local_bad = 1;
      break;
    }
    // (c) conformal: the BSSN conformal factor chi must be positive
    if (is_conformal && !(e[64] > 0.0)) {
      local_bad = 1;
      break;
    }
  }

  int global_bad = local_bad;
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_MAX, 1, &local_bad, &global_bad);
  return (global_bad == 0);
}

// Take a single mixed-scheme time-step.
struct gkyl_update_status
moment_update_mixed(gkyl_moment_app *app, double dt0)
{
  struct gkyl_update_status st = { .success = true, .dt_actual = dt0, .dt_suggested = DBL_MAX };
  double tcurr = app->tcurr;

  // Identify the coupled GR-Euler (MP) / vacuum-Einstein (wave_prop) pair.
  int fluid_idx = -1, einstein_idx = -1;
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].has_gr_euler && app->species[i].scheme_type == GKYL_MOMENT_MP) {
      fluid_idx = i;
    }
    if ((app->species[i].has_vacuum_einstein || app->species[i].has_vacuum_einstein_conformal)
      && app->species[i].scheme_type == GKYL_MOMENT_WAVE_PROP) {
      einstein_idx = i;
    }
  }
  bool coupled = (fluid_idx >= 0 && einstein_idx >= 0);

  // Cap dt by the spacetime (wave_prop) stable step up front
  double st_max_dt = (einstein_idx >= 0)
    ? moment_species_max_dt(app, &app->species[einstein_idx]) : DBL_MAX;
  double dt_use = fmin(dt0, st_max_dt);

  // 1.  Sync the live spacetime metric into the fluid state
  if (coupled) {
    mixed_sync_metric(app, fluid_idx, einstein_idx);
  }

  // 2. Advance the fluid (MP SSP-RK3)
  double dt = dt_use;
  if (fluid_idx >= 0) {
    struct gkyl_update_status fluid_st = moment_update_ssp_rk3(app, dt_use);
    if (!fluid_st.success) {
      return fluid_st;
    }
    dt = fluid_st.dt_actual;
    st.dt_suggested = fmin(st.dt_suggested, fluid_st.dt_suggested);
  }

  // 3. Advance the spacetime
  if (einstein_idx >= 0) {
    struct moment_species *esp = &app->species[einstein_idx];
    int ndim = esp->ndim;

    // Frozen (post-fluid-step) fluid state and adiabatic index for the matter source.
    const struct gkyl_array *fluid_f = coupled ? app->species[fluid_idx].fcurr : NULL;
    double gas_gamma = coupled ? app->species[fluid_idx].gr_euler_gas_gamma : 0.0;

    // 3a. Pre-transport half source on f[0] (= fcurr); refresh ghosts for the transport reconstruction.
    spacetime_source_half_step(app, esp, esp->f[0], fluid_f, gas_gamma, tcurr, 0.5 * dt);
    moment_species_apply_bc(app, tcurr, esp, esp->f[0]);

    // Hyperbolicity guard before the transport solve (so a non-SPD pre-source metric didn't enter wave_prop)
    if (!mixed_spacetime_metric_valid(app, esp)) {
      st.success = false;
      st.dt_actual = 0.0; // rejected step: do not advance app->tcurr (moment.c adds dt_actual unconditionally)
      st.dt_suggested = 0.5 * dt;
      return st;
    }

    // 3b. Transport: f[0] -> f[ndim] (pure wave_prop, no sources inside).
    struct gkyl_update_status ts = moment_species_update(app, esp, tcurr, dt);
    if (!ts.success) {
      return ts;
    }
    st.dt_suggested = fmin(st.dt_suggested, ts.dt_suggested);

    // Bring the transported state back into f[0] for the post source.
    gkyl_array_copy_range(esp->f[0], esp->f[ndim], &app->local_ext);

    // 3c. Post-transport half source on f[0]; reapply BC so ghosts are valid for the next step's sync.
    spacetime_source_half_step(app, esp, esp->f[0], fluid_f, gas_gamma, tcurr + dt, 0.5 * dt);
    moment_species_apply_bc(app, tcurr + dt, esp, esp->f[0]);

    // 3d. Hyperbolicity guard after the full spacetime step
    if (!mixed_spacetime_metric_valid(app, esp)) {
      st.success = false;
      st.dt_actual = 0.0;
      st.dt_suggested = 0.5 * dt;
      return st;
    }
  }

  st.dt_actual = dt;
  if (st.dt_suggested == DBL_MAX) {
    st.dt_suggested = dt;
  }
  return st;
}
