#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_return_flux_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_return_flux *reflux, struct gkyl_array *rhs)
{
}

static void
gk_species_return_flux_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_return_flux *reflux, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  for (int b=0; b<reflux->num_boundaries; ++b) {
    int dir = reflux->boundaries_dir[b];
    enum gkyl_edge_loc edge = reflux->boundaries_edge[b];
    const struct gkyl_array *bflux = gk_species_bflux_get_flux(&species->bflux, dir, edge);

    gkyl_bc_return_flux_gyrokinetic_advance(reflux->reflux_op[b], bflux, rhs);
  }

  app->stat.species_reflux_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_species_return_flux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_return_flux *reflux)
{
  reflux->rhs_func = gk_species_return_flux_rhs_disabled;

  // Check if any of the boundaries use RETURN_FLUX.
  reflux->num_boundaries = 0;

  int num_bound = 0;
  for (int d=0; d<app->cdim; ++d) {
    if ( gks->bc_is_np[d] ) {
      for (int e=0; e<2; ++e) {
        if ( ((e == 0 && gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_RETURN_FLUX_OPPOSITE) ||
              (e == 1 && gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_RETURN_FLUX_OPPOSITE)) ) {
          reflux->boundaries_dir[num_bound] = d;
          reflux->boundaries_edge[num_bound] = e==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
          reflux->boundaries_displace_dir[num_bound] = e==0? gks->lower_bc[d].index[0] : gks->upper_bc[d].index[0];

          // MF 2025/12/05: Restricted to displacing x flux at constant z (i.e. along y).
          assert(reflux->boundaries_dir[num_bound] == 0);
          assert(reflux->boundaries_displace_dir[num_bound] == 1);

          num_bound++;
        }
      }
    }
  }
  reflux->num_boundaries = num_bound;

  if (reflux->num_boundaries > 0) {
    for (int b=0; b<reflux->num_boundaries; ++b) {
      // Allocate return flux updater.
      int dir = reflux->boundaries_dir[b];
      enum gkyl_edge_loc edge = reflux->boundaries_edge[b];
      reflux->reflux_op[b] = gkyl_bc_return_flux_gyrokinetic_new(dir, edge,
        reflux->boundaries_displace_dir[b], &gks->grid, &gks->bflux.boundaries_phase_ghost_nosub[b],
        edge == GKYL_LOWER_EDGE? &gks->lower_skin[dir] : &gks->upper_skin[dir], app->use_gpu);
    }

    reflux->rhs_func = gk_species_return_flux_rhs_enabled;
  }
}

void
gk_species_return_flux_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_return_flux *reflux, struct gkyl_array *rhs)
{
  reflux->rhs_func(app, species, reflux, rhs);
}

void
gk_species_return_flux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_return_flux *reflux)
{
  if (reflux->num_boundaries > 0) {
    for (int b=0; b<reflux->num_boundaries; ++b) {
      gkyl_bc_return_flux_gyrokinetic_release(reflux->reflux_op[b]);
    }
  }
}
