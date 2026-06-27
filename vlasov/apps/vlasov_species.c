#include <gkyl_vlasov_priv.h>

// Vlasov species dispatcher: type-agnostic wrappers over the unified species
// container (struct vlasov_species). Each wrapper applies the operation to
// whichever aspects the species owns -- the kinetic distribution ('dist') and/or
// the fluid moments ('fluid') -- so the time steppers iterate the single
// app->species[] array (indexed over the overall species count) without
// branching on species type. A species missing an aspect simply skips it; a
// future PKPM species owns both and runs both halves. Modeled on vlasov_field.
//
// The RK-state arrays passed in (fin/fluidin/fout/fluidout) are indexed over the
// overall species count too: fin[i] is the distribution of species i (NULL if it
// has no kinetic aspect) and fluidin[i] is its fluid moments (NULL if it has no
// fluid aspect).

// Compute time-dependent applied acceleration (kinetic aspect only).
void
vlasov_species_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr)
{
  if (sp->dist && sp->dist->collisionless.app_accel_evolve)
    vm_species_collisionless_app_accel(app, &sp->dist->collisionless, tcurr);
}

// Self-collision moments and boundary corrections (kinetic aspect only).
void
vlasov_species_calc_self_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  if (sp->dist) {
    vm_species_lbo_moms(app, sp->dist, &sp->dist->lbo, fin);
    vm_species_bgk_moms(app, sp->dist, &sp->dist->bgk, fin);
  }
}

// Cross-collision moments (kinetic) and fluid primitive variables (fluid).
// Run after self moments for all species, so cross moments see them.
void
vlasov_species_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  if (sp->dist)
    vm_species_lbo_cross_moms(app, sp->dist, &sp->dist->lbo, fin);
  if (sp->fluid)
    vm_fluid_species_prim_vars(app, sp->fluid, fluidin);
}

// Compute the species RHS, returning the maximum stable time-step across aspects.
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

// Accumulate the source onto the RHS for both aspects. Takes the full RK-state
// arrays since the source may read/write across species (e.g. boundary-flux
// sources).
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
