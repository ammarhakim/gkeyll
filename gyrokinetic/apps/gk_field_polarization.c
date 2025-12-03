#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_position_map.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}

void
gk_field_polarization_potential_new(struct gk_field *f, struct gkyl_gyrokinetic_app *app)
{
  // Project the initial potential onto a p+1 tensor basis and compute the polarization
  // density to use use by species in calculating the initial ion density.
  f->init_phi_pol = true;
  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order + 1);

  f->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
  struct gkyl_array *phi_pol_ho = app->use_gpu ? mkarr(false, f->phi_pol->ncomp, f->phi_pol->size)
                                               : gkyl_array_acquire(f->phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &phi_pol_basis,
      .num_ret_vals = 1,
      .eval = f->info.polarization_potential,
      .ctx = f->info.polarization_potential_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
  });

  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &app->local, phi_pol_ho);
  gkyl_array_copy(f->phi_pol, phi_pol_ho);

  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_release(phi_pol_ho);
}

void
gk_field_polarization_potential_release(struct gk_field *f)
{
  gkyl_array_release(f->phi_pol);
}
