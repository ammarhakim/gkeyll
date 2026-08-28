#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>

static void
wall_potential_advance(gkyl_gyrokinetic_app *app,
  struct gk_species_wall_potential *wall, double tm)
{
  if (wall->projector) {
    gkyl_eval_on_nodes_advance(wall->projector, tm, &app->local_ext, wall->phi_host);
    if (app->use_gpu)
      gkyl_array_copy(wall->phi, wall->phi_host);
  }
}

static void
wall_potential_init(gkyl_gyrokinetic_app *app,
  const struct gkyl_gyrokinetic_bc *bc, struct gk_species_wall_potential *wall)
{
  *wall = (struct gk_species_wall_potential) { };
  if (bc->type != GKYL_BC_GK_SPECIES_SHEATH)
    return;

  wall->phi = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(wall->phi, 0.0);

  wall->phi_host = wall->phi;
  if (bc->aux_profile) {
    if (app->use_gpu)
      wall->phi_host = mkarr(false, wall->phi->ncomp, wall->phi->size);

    wall->projector = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
      1, bc->aux_profile, bc->aux_ctx);
    wall_potential_advance(app, wall, 0.0);
  }
}

static void
wall_potential_release(const gkyl_gyrokinetic_app *app,
  const struct gk_species_wall_potential *wall)
{
  if (!wall->phi)
    return;

  gkyl_array_release(wall->phi);
  if (wall->projector) {
    gkyl_eval_on_nodes_release(wall->projector);
    if (app->use_gpu)
      gkyl_array_release(wall->phi_host);
  }
}

void
gk_species_wall_potential_init(gkyl_gyrokinetic_app *app, struct gk_species *species)
{
  int par_dir = app->cdim-1;
  wall_potential_init(app, &species->lower_bc[par_dir], &species->phi_wall_lo);
  wall_potential_init(app, &species->upper_bc[par_dir], &species->phi_wall_up);
}

void
gk_species_wall_potential_advance(gkyl_gyrokinetic_app *app,
  struct gk_species *species, double tm)
{
  wall_potential_advance(app, &species->phi_wall_lo, tm);
  wall_potential_advance(app, &species->phi_wall_up, tm);
}

void
gk_species_wall_potential_release(const gkyl_gyrokinetic_app *app,
  const struct gk_species *species)
{
  wall_potential_release(app, &species->phi_wall_lo);
  wall_potential_release(app, &species->phi_wall_up);
}
