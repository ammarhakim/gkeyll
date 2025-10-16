#include <gkyl_gyrokinetic_priv.h>

// Take time-step using an implicit method for collisions.
// Use the actual timestep used to update explicit advection. 
void
gyrokinetic_update_implicit_coll(gkyl_gyrokinetic_app* app, double dt0)
{
  // Update gyrokinetic species BGK collisions
  int ns = app->num_species;  
  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];

  // Fetch input and output arrays and compute necessary moments for collisions
  for (int i=0; i<ns; ++i) {
    struct gk_species *gks = &app->species[i];
    fin[i] = gks->f;
    fout[i] = gks->f1;
    gk_species_bgk_moms_implicit(app, gks, &gks->bgk, fin[i]);
  }

  // Compute necessary moments for cross-species collisions.
  // Needs to be done after self-collisions moments, so separate loop over species.
  for (int i=0; i<ns; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_bgk_cross_moms_implicit(app, gks, &gks->bgk, fin[i]);        
  }

  // implicit BGK contributions for gyrokinetic species
  for (int i=0; i<ns; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_rhs_implicit(app, gks, fin[i], fout[i], dt0);
  }

  // Update neutral species BGK collisions
  // If neutral species is static, collision_id = GKYL_NO_COLLISIONS by default
  int neuts = app->num_neut_species;  
  const struct gkyl_array *fin_neut[neuts];
  struct gkyl_array *fout_neut[neuts];

  // Fetch input and output neutral arrays and compute necessary neutral moments for collisions
  for (int i=0; i<neuts; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    fin_neut[i] = gkns->f;
    fout_neut[i] = gkns->f1;
    gk_neut_species_bgk_moms_implicit(app, gkns, &gkns->bgk, fin_neut[i]);
  }

  // implicit BGK contributions for neutral species
  for (int i=0; i<neuts; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_rhs_implicit(app, gkns, fin_neut[i], fout_neut[i], dt0);
  }

  // Apply boundary conditions and copy solution
  gyrokinetic_calc_field_and_apply_bc(app, app->tcurr, fout, fout_neut);
  for (int i=0; i<ns; ++i) {
    gkyl_array_copy_range(app->species[i].f, app->species[i].f1, &app->species[i].local_ext);
  };
  for (int i=0; i<neuts; ++i) {
    gkyl_array_copy_range(app->neut_species[i].f, app->neut_species[i].f1, &app->neut_species[i].local_ext);
  };
}
