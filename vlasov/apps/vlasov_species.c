#include <gkyl_vlasov_priv.h>

#include <string.h>

// Vlasov species dispatcher. The concrete species type (kinetic vs fluid) is
// chosen at construction: vlasov_kinetic_species_new()/vlasov_fluid_species_new()
// allocate the aspect sub-object(s) and set the container's *_func pointers,
// and every vlasov_species_* wrapper below simply forwards to the appropriate
// pointer. This keeps species-type dispatch out of the time steppers and
// mirrors vlasov_field.c (vm_field_new/vp_field_new wiring). A future PKPM
// species wires implementations that drive both aspects.
//
// The RK-state arrays passed in (fin/fluidin/fout/fluidout) are indexed over
// the overall species count: fin[i] is the distribution of species i (NULL if
// it has no kinetic aspect) and fluidin[i] is its fluid moments (NULL if it
// has no fluid aspect).

// --- Shared no-op implementations -------------------------------------------
// Wired where a species type lacks the operation (e.g. a kinetic species has
// no fluid limiter); a PKPM species later replaces these with real methods.

static void
species_no_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr) { }

static void
species_no_calc_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin) { }

static void
species_no_limiter(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *fluid) { }

// No explicit field coupling: wired for the null field, and for fluid species
// whose EM coupling is the implicit op-split (see vm_fluid_em_coupling.c). An
// explicit-source fluid mode later wires a real method here instead.
static void
species_no_field_coupling(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target) { }

// --- Kinetic (dist-aspect) implementations ----------------------------------

static void
kinetic_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr)
{
  // evolve flag is data (set by vm_species_init), not type: checked at run time.
  if (sp->dist->collisionless.app_accel_evolve)
    vm_species_collisionless_app_accel(app, &sp->dist->collisionless, tcurr);
}

static void
kinetic_calc_self_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  vm_species_lbo_moms(app, sp->dist, &sp->dist->lbo, fin);
  vm_species_bgk_moms(app, sp->dist, &sp->dist->bgk, fin);
}

static void
kinetic_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  vm_species_lbo_cross_moms(app, sp->dist, &sp->dist->lbo, fin);
}

static double
kinetic_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, const struct gkyl_array *emin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  return vm_species_rhs(app, sp->dist, fin, emin, fout);
}

// Explicit coupling to Vlasov-Maxwell: accumulate this species' current onto
// the EM RHS ('target' is emout), J_s = q_s * m1i, scaled by -1/epsilon0.
static void
kinetic_accumulate_current(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target)
{
  struct vm_species *s = sp->dist;
  double qbyeps = sp->charge/app->field->info.epsilon0;

  vm_species_moment_calc(&s->m1i, s->local, app->local, fin);
  gkyl_array_accumulate_range(target, -qbyeps, s->m1i.marr, &app->local);

  if (app->field->use_ghost_current) {
    double avals_ghost_current[1], avals_ghost_current_global[1];
    // First set the scalar ghost current array to the cell average
    // current/(epsilon0*nx) where nx is the number of x cells.
    gkyl_array_set_range(app->field->ghost_current, qbyeps/app->grid.cells[0], s->m1i.marr, &app->local);
    // Integrate the current over the whole domain to find the globally averaged ghost current.
    if (app->use_gpu) {
      gkyl_array_reduce_range(app->field->red_ghost_current, app->field->ghost_current, GKYL_SUM, &app->local);
      gkyl_cu_memcpy(avals_ghost_current, app->field->red_ghost_current, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      gkyl_array_reduce_range(avals_ghost_current, app->field->ghost_current, GKYL_SUM, &app->local);
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, avals_ghost_current, avals_ghost_current_global);
    // Set the scalar ghost current array to the global average current and accumulate to the electric field.
    gkyl_array_clear(app->field->ghost_current, avals_ghost_current_global[0]);
    gkyl_array_accumulate_range(target, 1.0, app->field->ghost_current, &app->local);
  }
}

// Explicit coupling to Vlasov-Poisson: accumulate this species' charge density
// onto the Poisson source ('target' is the field's rho_c), rho_s = q_s * m0.
static void
kinetic_accumulate_charge_dens(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target)
{
  struct vm_species *s = sp->dist;

  vm_species_moment_calc(&s->m0, s->local, app->local, fin);
  gkyl_array_accumulate_range(target, sp->charge, s->m0.marr, &app->local);
}

static void
kinetic_calc_source_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  if (sp->dist->source_id)
    vm_species_source_adapt_moms(app, sp->dist, &sp->dist->src, fin);
}

static void
kinetic_source_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fluidout[])
{
  if (sp->dist->source_id) {
    if (sp->dist->src.evolve_source)
      vm_species_source_calc(app, sp->dist, &sp->dist->src, tcurr);
    vm_species_source_adapt(app, sp->dist, &sp->dist->src);
    vm_species_source_rhs(app, sp->dist, &sp->dist->src, fin, fout);
  }
}

static void
kinetic_step_f(struct vlasov_species *sp, double dt,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  vm_species_step_f(sp->dist, fout, dt, fin);
}

static void
kinetic_combine(gkyl_vlasov_app *app, struct vlasov_species *sp, double c1, double c2)
{
  struct vm_species *d = sp->dist;
  vm_species_combine(d, d->f1, c1, d->f, c2, d->fnew, &d->local_ext);
}

static void
kinetic_copy_range(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  struct vm_species *d = sp->dist;
  vm_species_copy_range(d, d->f, d->f1, &d->local_ext);
}

static void
kinetic_apply_bc(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *f, struct gkyl_array *fluid, double tcurr)
{
  vm_species_apply_bc(app, sp->dist, f, tcurr);
}

// --- Fluid (fluid-aspect) implementations -----------------------------------

static void
fluid_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  // Primitive moments; run in the cross-moms phase so they follow all self moments.
  vm_fluid_species_prim_vars(app, sp->fluid, fluidin);
}

static double
fluid_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, const struct gkyl_array *emin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  return vm_fluid_species_rhs(app, sp->fluid, fluidin, emin, fluidout);
}

static void
fluid_source_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fluidout[])
{
  if (sp->fluid->source_id)
    vm_fluid_species_source_rhs(app, sp->fluid, &sp->fluid->src, fluidin, fluidout);
}

static void
fluid_step_f(struct vlasov_species *sp, double dt,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  vm_fluid_species_step_f(sp->fluid, fluidout, dt, fluidin);
}

static void
fluid_combine(gkyl_vlasov_app *app, struct vlasov_species *sp, double c1, double c2)
{
  struct vm_fluid_species *f = sp->fluid;
  vm_fluid_species_combine(f, f->fluid1, c1, f->fluid, c2, f->fluidnew, &app->local_ext);
}

static void
fluid_copy_range(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  struct vm_fluid_species *f = sp->fluid;
  vm_fluid_species_copy_range(f, f->fluid, f->fluid1, &app->local_ext);
}

static void
fluid_apply_bc(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *f, struct gkyl_array *fluid, double tcurr)
{
  vm_fluid_species_apply_bc(app, sp->fluid, fluid);
}

static void
fluid_limiter(gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fluid)
{
  vm_fluid_species_limiter(app, sp->fluid, fluid);
}

// --- Constructors ------------------------------------------------------------

void
vlasov_kinetic_species_new(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_species *info, struct vlasov_species *sp)
{
  *sp = (struct vlasov_species) { };
  sp->type = GKYL_SPECIES_VLASOV;

  strcpy(sp->name, info->name);
  sp->charge = info->charge;
  sp->mass = info->mass;

  sp->dist = gkyl_malloc(sizeof(struct vm_species));
  *sp->dist = (struct vm_species) { };
  sp->dist->info = *info;

  sp->calc_app_accel_func = kinetic_calc_app_accel;
  sp->calc_self_moms_func = kinetic_calc_self_moms;
  sp->calc_cross_moms_func = kinetic_calc_cross_moms;
  sp->rhs_func = kinetic_rhs;
  // Explicit field coupling is wired per field type (the field object is
  // constructed before the species): Maxwell fields (E_B, GR_D_B) take the
  // current density, Poisson fields the charge density, the null field nothing.
  switch (app->field->field_id) {
    case GKYL_FIELD_NULL:
      sp->accumulate_field_coupling_func = species_no_field_coupling;
      break;
    case GKYL_FIELD_PHI:
    case GKYL_FIELD_PHI_EXT_POTENTIALS:
    case GKYL_FIELD_PHI_EXT_FIELDS:
      sp->accumulate_field_coupling_func = kinetic_accumulate_charge_dens;
      break;
    default:
      sp->accumulate_field_coupling_func = kinetic_accumulate_current;
      break;
  }
  sp->calc_source_moms_func = kinetic_calc_source_moms;
  sp->source_rhs_func = kinetic_source_rhs;
  sp->step_f_func = kinetic_step_f;
  sp->combine_func = kinetic_combine;
  sp->copy_range_func = kinetic_copy_range;
  sp->apply_bc_func = kinetic_apply_bc;
  sp->limiter_func = species_no_limiter;
}

void
vlasov_fluid_species_new(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_fluid_species *info, struct vlasov_species *sp)
{
  *sp = (struct vlasov_species) { };
  sp->type = GKYL_SPECIES_FLUID;

  strcpy(sp->name, info->name);
  sp->charge = info->charge;
  sp->mass = info->mass;

  sp->fluid = gkyl_malloc(sizeof(struct vm_fluid_species));
  *sp->fluid = (struct vm_fluid_species) { };
  sp->fluid->info = *info;

  sp->calc_app_accel_func = species_no_calc_app_accel;
  sp->calc_self_moms_func = species_no_calc_moms;
  sp->calc_cross_moms_func = fluid_calc_cross_moms;
  sp->rhs_func = fluid_rhs;
  sp->accumulate_field_coupling_func = species_no_field_coupling;
  sp->calc_source_moms_func = species_no_calc_moms;
  sp->source_rhs_func = fluid_source_rhs;
  sp->step_f_func = fluid_step_f;
  sp->combine_func = fluid_combine;
  sp->copy_range_func = fluid_copy_range;
  sp->apply_bc_func = fluid_apply_bc;
  sp->limiter_func = fluid_limiter;
}

// --- Type-agnostic wrappers (pure forwards) ----------------------------------

// Compute time-dependent applied acceleration.
void
vlasov_species_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr)
{
  sp->calc_app_accel_func(app, sp, tcurr);
}

// Self-collision moments and boundary corrections.
void
vlasov_species_calc_self_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  sp->calc_self_moms_func(app, sp, fin);
}

// Cross-collision moments (kinetic) and fluid primitive variables (fluid).
// Run after self moments for all species, so cross moments see them.
void
vlasov_species_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  sp->calc_cross_moms_func(app, sp, fin, fluidin);
}

// Compute the species RHS, returning the maximum stable time-step.
double
vlasov_species_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, const struct gkyl_array *emin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  return sp->rhs_func(app, sp, fin, fluidin, emin, fout, fluidout);
}

// Accumulate this species' explicit source contribution onto the field's
// target array (Maxwell: current onto emout; Poisson: charge density onto
// rho_c; no-op for the null field and for implicitly-coupled fluid species).
void
vlasov_species_accumulate_field_coupling(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, struct gkyl_array *target)
{
  sp->accumulate_field_coupling_func(app, sp, fin, fluidin, target);
}

// Pre-compute moments needed to rescale adaptive sources.
// Done in its own pass so the species RHS (and any boundary fluxes) are complete.
void
vlasov_species_calc_source_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  sp->calc_source_moms_func(app, sp, fin);
}

// Accumulate the source onto the RHS. Takes the full RK-state arrays since the
// source may read/write across species (e.g. boundary-flux sources).
void
vlasov_species_source_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fluidout[])
{
  sp->source_rhs_func(app, sp, tcurr, fin, fluidin, fout, fluidout);
}

// Forward-Euler accumulate: out = dt*out + in.
void
vlasov_species_step_f(struct vlasov_species *sp, double dt,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  sp->step_f_func(sp, dt, fin, fluidin, fout, fluidout);
}

// Combine RK stages into the first-stage buffer: f1 = c1*f + c2*fnew
// (the RK3 stepper's combine pattern).
void
vlasov_species_combine(gkyl_vlasov_app *app, struct vlasov_species *sp, double c1, double c2)
{
  sp->combine_func(app, sp, c1, c2);
}

// Copy the first-stage buffer back into the solution: f = f1.
void
vlasov_species_copy_range(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  sp->copy_range_func(app, sp);
}

// Apply boundary conditions.
void
vlasov_species_apply_bc(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *f, struct gkyl_array *fluid, double tcurr)
{
  sp->apply_bc_func(app, sp, f, fluid, tcurr);
}

// Apply the post-update slope limiter.
void
vlasov_species_limiter(gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fluid)
{
  sp->limiter_func(app, sp, fluid);
}
