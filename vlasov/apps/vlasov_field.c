#include <gkyl_vlasov_priv.h>

// Vlasov field dispatcher and type-agnostic wrappers. The concrete field type
// (Vlasov-Maxwell vs Vlasov-Poisson) is chosen here at construction; the
// constructors vm_field_new()/vp_field_new() set the field's *_func pointers,
// and every wrapper below simply forwards to the appropriate pointer. This
// keeps field-type dispatch out of vlasov.c / the time steppers, mirroring how
// vm_species and GK's gk_field are organized.

struct vm_field*
vlasov_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  return vm->is_electrostatic ? vp_field_new(vm, app) : vm_field_new(vm, app);
}

double
vlasov_field_update(gkyl_vlasov_app *app, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *emin, struct gkyl_array *emout)
{
  return app->field->update_func(app, tcurr, fin, emin, emout);
}

void
vlasov_field_combine(gkyl_vlasov_app *app, struct gkyl_array *out,
  double c1, const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2)
{
  app->field->combine_func(app, app->field, out, c1, arr1, c2, arr2);
}

void
vlasov_field_copy_range(gkyl_vlasov_app *app, struct gkyl_array *out, const struct gkyl_array *inp)
{
  app->field->copy_range_func(app, app->field, out, inp);
}

void
vlasov_field_apply_ic(gkyl_vlasov_app *app, const struct gkyl_array *fin[], double t0)
{
  app->field->apply_ic_func(app, app->field, fin, t0);
}

void
vlasov_field_apply_bc(gkyl_vlasov_app *app, struct gkyl_array *em)
{
  app->field->apply_bc_func(app, app->field, em);
}

void
vlasov_field_limiter(gkyl_vlasov_app *app, struct gkyl_array *em)
{
  app->field->limiter_func(app, app->field, em);
}

void
vlasov_field_complete_update(gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout)
{
  app->field->complete_update_func(app, dt, fin, fluidin, emin, emout);
}

void
vlasov_field_calc_ext_em(gkyl_vlasov_app *app, double tm)
{
  app->field->calc_ext_em_func(app, app->field, tm);
}

void
vlasov_field_calc_app_current(gkyl_vlasov_app *app, double tm)
{
  app->field->calc_app_current_func(app, app->field, tm);
}

void
vlasov_field_calc_ext_pot(gkyl_vlasov_app *app, double tm)
{
  app->field->calc_ext_pot_func(app, app->field, tm);
}

void
vlasov_field_calc_energy(gkyl_vlasov_app *app, double tm)
{
  app->field->calc_energy_func(app, tm, app->field);
}

void
vlasov_field_write(gkyl_vlasov_app *app, double tm, int frame)
{
  app->field->write_func(app, tm, frame);
}

void
vlasov_field_write_energy(gkyl_vlasov_app *app)
{
  app->field->write_energy_func(app);
}

void
vlasov_field_release(gkyl_vlasov_app *app)
{
  app->field->release_func(app, app->field);
}
