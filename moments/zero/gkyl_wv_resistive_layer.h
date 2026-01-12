#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_wv_resistive_layer gkyl_wv_resistive_layer;

struct gkyl_wv_resistive_layer {
  void *ctx;
  void *resistivity_func;
};

/**
 * Create new updater to apply set of boundary conditions.
 *
 * @param mask_func
 * @param type
 * @param ctx Context to pass to bcfunc.
 * @return New updater pointer.
 */
gkyl_wv_resistive_layer* gkyl_wv_resistive_layer_new(void *mask_func, void *ctx);

void gkyl_wv_resistive_layer_new_profile(struct gkyl_wv_resistive_layer *res,
  struct gkyl_rect_grid *grid, struct gkyl_range *rng, struct gkyl_array *sigma);

/**
 * Delete structure.
 *
 * @param embed_geo Structure to delete.
 */
void gkyl_wv_resistive_layer_release(gkyl_wv_resistive_layer *res);
