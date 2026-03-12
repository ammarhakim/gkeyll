#include <gkyl_gyrokinetic_priv.h>

// Take time-step using an implicit method for collisions.
// Use the actual timestep used to update explicit advection. 
void
gyrokinetic_update_implicit_coll(gkyl_gyrokinetic_app* app, double dt0)
{
  // Add contribution from implicit collision terms.
  double tcurr = app->tcurr, dt = dt0;

  struct gkyl_gyrokinetic_fdot_args *fdot_args = &app->fdot_args;
  const struct gkyl_array **fin = fdot_args->fin;
  const struct gkyl_array **fin_neut = fdot_args->fin_neut;
  const struct gkyl_array **fieldin = fdot_args->fieldin;
  struct gkyl_array **fout = fdot_args->fout;
  struct gkyl_array **fout_neut = fdot_args->fout_neut;
  struct gkyl_array ***bflux_in = fdot_args->bflux_in;
  struct gkyl_array ***bflux_in_neut = fdot_args->bflux_in_neut;
  struct gkyl_array ***bflux_out = fdot_args->bflux_out;
  struct gkyl_array ***bflux_out_neut = fdot_args->bflux_out_neut;
  struct gkyl_array **fieldout = fdot_args->fieldout;

  struct gkyl_update_status st = { .success = true };

  const int ns_charged = app->num_species;  
  const int ns_neut = app->num_neut_species;  
  const int nfields = app->field->num_fields;  

  // Fetch input and output arrays.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    fin[i] = gks->f;
    fout[i] = gks->f1;
    bflux_in[i] = gks->bflux.f;
    bflux_out[i] = gks->bflux.f;
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fin_neut[i] = gkns->f;
    fout_neut[i] = gkns->f1;
    bflux_in_neut[i] = gkns->bflux.f;
    bflux_out_neut[i] = gkns->bflux.f;
  }
  for (int i=0; i<nfields; ++i) {
    struct gk_field *gkf = app->field;
    fieldin[i] = gkf->f[i];
    fieldout[i] = gkf->f[i];
  }

  // Compute df/dt from implicit terms.
  gyrokinetic_rhs_implicit(app, tcurr, dt, fdot_args, &st);

  // Apply boundary conditions and copy solution.
  gyrokinetic_calc_field(app, tcurr, fout, bflux_out, fieldout);
  gyrokinetic_apply_bc(app, tcurr, fieldout, fout, fout_neut);

  // Copy the solution back into the main array.
  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    fin[i] = gks->f1;
    fout[i] = gks->f;
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fin_neut[i] = gkns->f1;
    fout_neut[i] = gkns->f;
  }

  for (int i=0; i<ns_charged; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_copy_range(gks, fout[i], fin[i], &gks->local_ext);
  }
  for (int i=0; i<ns_neut; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_copy_range(gkns, fout_neut[i], fin_neut[i], &gkns->local_ext);
  };
}
