#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <assert.h>
#include <math.h>

void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  assert(app->cdim > 1);
  f->invert_flr = gk_field_invert_flr;

  // Reference (squared) gyroradius in the operator A = 1 - rho^2*nabla_perp^2
  // used to retrieve phi from the modified potential Phi_0 (step 4 of the
  // algorithm in flr_effects.tex).
  double rhoSq_ref = 0.0;
  if (f->info.flr.type == GKYL_GK_FLR_PADE_CONST) {
    // Single reference species (e.g. the main ion).
    assert(f->info.flr.gyroradius > 0.0);
    rhoSq_ref = pow(f->info.flr.gyroradius, 2.0);
  }
  else {
    // GKYL_GK_FLR_PADE_CONST_SUM, or field-level FLR options not set (default):
    // polarization-weighted average of the species gyroradii,
    //   rho^2 = sum_s eps_s0*rho_s0^2 / sum_s eps_s0,  eps_s0 = n_s0*m_s/B^2,
    // which suppresses the electron contribution by m_e/m_i and reduces to
    // rho_i0^2 for a single ion species. Species without FLR enabled
    // contribute rho_s0^2 = 0.
    double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
    double eps_sum = 0.0;
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      double gyroradius_bmag = s->info.flr.bmag ? s->info.flr.bmag : app->bmag_ref;
      double rhoSq_s = s->info.flr.Tperp * s->info.mass / pow(s->info.charge * gyroradius_bmag, 2.0);
      double eps_s0 = s->info.polarization_density * s->info.mass / pow(polarization_bmag, 2.0);
      rhoSq_ref += eps_s0 * rhoSq_s;
      eps_sum += eps_s0;
    }
    rhoSq_ref /= eps_sum;
  }
  // Weight in the perpendicular Laplacian of A, staged for the forthcoming
  // Laplacian-apply infrastructure.
  f->flr_rhoSq = mkarr(app->use_gpu, (2 * (app->cdim - 1) - 1) * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_set_offset(f->flr_rhoSq, rhoSq_ref, app->gk_geom->geo_int.gxxj, 0 * app->basis.num_basis);
  if (app->cdim > 2) {
    gkyl_array_set_offset(f->flr_rhoSq, rhoSq_ref, app->gk_geom->geo_int.gxyj, 1 * app->basis.num_basis);
    gkyl_array_set_offset(f->flr_rhoSq, rhoSq_ref, app->gk_geom->geo_int.gyyj, 2 * app->basis.num_basis);
  }
}

void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  // Retrieve phi from the modified potential Phi_0 by applying
  //   A = 1 - rho_i0^2*nabla_perp^2
}

void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
}

void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->flr_rhoSq);
}
