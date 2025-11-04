#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_fv_proj.h>
#include <gkyl_wv_resistive_layer.h>

struct gkyl_wv_resistive_layer*
gkyl_wv_resistive_layer_new(void *resistivity_func, void *ctx)
{
  gkyl_wv_resistive_layer *res = gkyl_malloc(sizeof(gkyl_wv_resistive_layer));

  res->ctx = ctx;
  res->resistivity_func = resistivity_func;

  return res;
}

void
gkyl_wv_resistive_layer_new_profile(struct gkyl_wv_resistive_layer *res,
  struct gkyl_rect_grid *grid, struct gkyl_range *rng, struct gkyl_array *sigma)
{
  gkyl_fv_proj *proj = gkyl_fv_proj_new(grid, 1, 1, res->resistivity_func,
    res->ctx);
  gkyl_fv_proj_advance(proj, 0.0, rng, sigma);
  gkyl_fv_proj_release(proj);
}

void
gkyl_wv_resistive_layer_release(gkyl_wv_resistive_layer* res)
{
  gkyl_free(res);
}
