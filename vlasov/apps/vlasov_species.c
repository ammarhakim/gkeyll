#include <gkyl_vlasov_priv.h>

#include <string.h>

// Vlasov species dispatcher over the unified container (struct vlasov_species).
//
// Two kinds of operations live here, deliberately organized differently:
//
// 1. Per-aspect operations (rhs, step/combine/copy, BCs, limiter, sources, the
//    implicit collision phases): structurally identical for every species
//    type -- "apply to whichever aspects are present" -- so the wrappers
//    NULL-test sp->dist / sp->fluid directly. A future PKPM species (both
//    aspects allocated) runs both halves sequentially with NO new dispatch
//    code here; per-aspect flavor differences (e.g. PKPM's dg_vlasov_pkpm vs
//    the default Vlasov equation) live on the aspect-level vtables wired by
//    the aspect inits.
//
// 2. Type-composed operations, wired once at construction through the
//    container vtable (mirroring vm_field_new/vp_field_new):
//    - accumulate_field_coupling_func: the species' explicit field source,
//      selected by (species type x field type) -- kinetic x Maxwell gives
//      current, kinetic x Poisson gives charge density, implicitly-coupled
//      fluids give a no-op.
//    - calc_self_moms_func / calc_cross_moms_func: the staging phases that
//      fill the pre-RHS auxiliary arrays. These genuinely differ by type
//      (kinetic: collision moments; fluid: primitive variables; PKPM later:
//      its own moments/vars machinery reading both aspects) rather than
//      composing per-aspect.
//
// The RK-state arrays passed in (fin/fluidin/fout/fluidout) are indexed over
// the overall species count: fin[i] is the distribution of species i (NULL if
// it has no kinetic aspect) and fluidin[i] is its fluid moments (NULL if it
// has no fluid aspect).

// --- Type-composed implementations (wired by the constructors) ---------------

// No-op moments phase: wired where a species type has no work in that phase.
static void
species_no_calc_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin) { }

// No explicit field coupling: wired for the null field, and for fluid species
// whose EM coupling is the implicit op-split (see vm_fluid_em_coupling.c). An
// explicit-source fluid mode later wires a real method here instead.
static void
species_no_field_coupling(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target) { }

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

static void
fluid_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  // Primitive moments; run in the cross-moms phase so they follow all self moments.
  vm_fluid_species_prim_vars(app, sp->fluid, fluidin);
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

// --- Constructors ------------------------------------------------------------

// Construct a species from a unified input: validate the declared type against
// the blocks, hoist the top-level identity into the selected block, and
// dispatch to the typed constructor.
void
vlasov_species_new(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_species *inp, struct vlasov_species *sp)
{
  switch (inp->type) {
    case GKYL_SPECIES_VLASOV: {
      assert(inp->kinetic.cells[0] > 0); // kinetic species require a velocity grid
      assert(!inp->fluid.equation); // ... and must not carry a fluid equation object

      struct gkyl_vlasov_kinetic_species ki = inp->kinetic;
      strcpy(ki.name, inp->name);
      ki.charge = inp->charge;
      ki.mass = inp->mass;
      vlasov_kinetic_species_new(app, &ki, sp);
      break;
    }
    case GKYL_SPECIES_FLUID: {
      assert(inp->fluid.equation); // fluid species require an equation object
      assert(inp->kinetic.cells[0] == 0); // ... and must not declare a velocity grid

      struct gkyl_vlasov_fluid_species fi = inp->fluid;
      strcpy(fi.name, inp->name);
      fi.charge = inp->charge;
      fi.mass = inp->mass;
      vlasov_fluid_species_new(app, &fi, sp);
      break;
    }
    default:
      assert(false); // GKYL_SPECIES_PKPM is not wired in the Vlasov app yet
      break;
  }
}

void
vlasov_kinetic_species_new(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_kinetic_species *info, struct vlasov_species *sp)
{
  *sp = (struct vlasov_species) { };
  sp->type = GKYL_SPECIES_VLASOV;

  strcpy(sp->name, info->name);
  sp->charge = info->charge;
  sp->mass = info->mass;

  sp->dist = gkyl_malloc(sizeof(struct vm_species));
  *sp->dist = (struct vm_species) { };
  sp->dist->info = *info;

  sp->calc_self_moms_func = kinetic_calc_self_moms;
  sp->calc_cross_moms_func = kinetic_calc_cross_moms;
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

  sp->calc_self_moms_func = species_no_calc_moms;
  sp->calc_cross_moms_func = fluid_calc_cross_moms;
  sp->accumulate_field_coupling_func = species_no_field_coupling;
}

// --- Type-composed wrappers (forward through the container vtable) -----------

// Self-collision moments and boundary corrections (staging phase 1).
void
vlasov_species_calc_self_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  sp->calc_self_moms_func(app, sp, fin);
}

// Cross-collision moments / fluid primitive variables (staging phase 2).
// Run after self moments for all species, so cross moments see them.
void
vlasov_species_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  sp->calc_cross_moms_func(app, sp, fin, fluidin);
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

// --- Per-aspect wrappers (apply to whichever aspects are present) -------------

// Compute time-dependent applied acceleration (kinetic aspect only).
void
vlasov_species_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr)
{
  if (sp->dist && sp->dist->collisionless.app_accel_evolve)
    vm_species_collisionless_app_accel(app, &sp->dist->collisionless, tcurr);
}

// Compute the species RHS, returning the maximum stable time-step across
// aspects. Each aspect advances through its own vtable (flavor dispatch); the
// coupling between a PKPM species' aspects is pre-staged in the moms/vars
// phases, so the aspect RHSs are independent here.
double
vlasov_species_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, const struct gkyl_array *emin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  double dt = DBL_MAX;
  if (sp->dist)  dt = fmin(dt, vm_species_rhs(app, sp->dist, fin, emin, fout));
  if (sp->fluid) dt = fmin(dt, vm_fluid_species_rhs(app, sp->fluid, fluidin, emin, fluidout));
  return dt;
}

// Pre-compute moments needed to rescale adaptive sources (kinetic aspect only).
// Done in its own pass so the species RHS (and any boundary fluxes) are complete.
void
vlasov_species_calc_source_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  if (sp->dist && sp->dist->source_id)
    vm_species_source_adapt_moms(app, sp->dist, &sp->dist->src, fin);
}

// Accumulate the source onto the RHS for each present aspect. Takes the full
// RK-state arrays since the source may read/write across species (e.g.
// boundary-flux sources).
void
vlasov_species_source_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fluidout[])
{
  if (sp->dist && sp->dist->source_id) {
    if (sp->dist->src.evolve_source)
      vm_species_source_calc(app, sp->dist, &sp->dist->src, tcurr);
    vm_species_source_adapt(app, sp->dist, &sp->dist->src);
    vm_species_source_rhs(app, sp->dist, &sp->dist->src, fin, fout);
  }
  if (sp->fluid && sp->fluid->source_id)
    vm_fluid_species_source_rhs(app, sp->fluid, &sp->fluid->src, fluidin, fluidout);
}

// Implicit (op-split) collision update, phase 1: moments of the input
// distribution (kinetic aspect only; the aspect-level BGK object no-ops for
// species without implicit BGK collisions).
void
vlasov_species_calc_implicit_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  if (sp->dist)
    vm_species_bgk_moms_implicit(app, sp->dist, &sp->dist->bgk, fin);
}

// Implicit collision update, phase 2: the implicit RHS. Run after phase 1 for
// all species (the RHS also computes cross moments).
void
vlasov_species_rhs_implicit(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *fout, double dt)
{
  if (sp->dist)
    vm_species_rhs_implicit(app, sp->dist, fin, fout, dt);
}

// Implicit collision update, phase 3: BCs and copy-back into the solution.
void
vlasov_species_finish_implicit_update(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *fout, double tcurr)
{
  if (sp->dist) {
    vm_species_apply_bc(app, sp->dist, fout, tcurr);
    gkyl_array_copy_range(sp->dist->f, fout, &sp->dist->local_ext);
  }
}

// Forward-Euler accumulate: out = dt*out + in, for each present aspect.
void
vlasov_species_step_f(struct vlasov_species *sp, double dt,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  if (sp->dist)  vm_species_step_f(sp->dist, fout, dt, fin);
  if (sp->fluid) vm_fluid_species_step_f(sp->fluid, fluidout, dt, fluidin);
}

// Combine RK stages into the first-stage buffer: f1 = c1*f + c2*fnew, for each
// present aspect (the RK3 stepper's combine pattern).
void
vlasov_species_combine(gkyl_vlasov_app *app, struct vlasov_species *sp, double c1, double c2)
{
  if (sp->dist) {
    struct vm_species *d = sp->dist;
    vm_species_combine(d, d->f1, c1, d->f, c2, d->fnew, &d->local_ext);
  }
  if (sp->fluid) {
    struct vm_fluid_species *f = sp->fluid;
    vm_fluid_species_combine(f, f->fluid1, c1, f->fluid, c2, f->fluidnew, &app->local_ext);
  }
}

// Copy the first-stage buffer back into the solution: f = f1, for each aspect.
void
vlasov_species_copy_range(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->dist) {
    struct vm_species *d = sp->dist;
    vm_species_copy_range(d, d->f, d->f1, &d->local_ext);
  }
  if (sp->fluid) {
    struct vm_fluid_species *f = sp->fluid;
    vm_fluid_species_copy_range(f, f->fluid, f->fluid1, &app->local_ext);
  }
}

// Apply boundary conditions to each present aspect.
void
vlasov_species_apply_bc(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *f, struct gkyl_array *fluid, double tcurr)
{
  if (sp->dist)  vm_species_apply_bc(app, sp->dist, f, tcurr);
  if (sp->fluid) vm_fluid_species_apply_bc(app, sp->fluid, fluid);
}

// Apply the post-update slope limiter (fluid aspect only).
void
vlasov_species_limiter(gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fluid)
{
  if (sp->fluid) vm_fluid_species_limiter(app, sp->fluid, fluid);
}
