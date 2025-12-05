#pragma once

#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_bc_return_flux_gyrokinetic gkyl_bc_return_flux_gyrokinetic;

/**
 * Create a new updater to take the boundary flux from one skin cell, and put it
 * in another skin cell, displacing it along the `displacement_dir` direction.
 * 
 * MF 2025/12/05: At the moment it just shifts it to the other side of the domain.
 *
 * @param dir Direction perpendicular to the boundary.
 * @param edge Edge of the boundary (lower or upper).
 * @param displacement_dir Direction along which to move the boundary flux .
 * @param grid Grid object.
 * @param bflux_ghost_r Ghost range to index the boundary flux.
 * @param skin_r Skin range (to index df/dt).
 * @param use_gpu Whether to run on a GPU.
 */
struct gkyl_bc_return_flux_gyrokinetic* 
gkyl_bc_return_flux_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, int displacement_dir,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *bflux_ghost_r,
  const struct gkyl_range *skin_r, bool use_gpu);

/**
 * Displace the boundary flux and place it in the corresponding skin cell.
 *
 * @param up Return flux updater object.
 * @param bflux Boundary flux.
 * @param rhs Time rate of change (df/dt) to add displaced bflux to.
 */
void
gkyl_bc_return_flux_gyrokinetic_advance(const struct gkyl_bc_return_flux_gyrokinetic *up,
  const struct gkyl_array *bflux, struct gkyl_array *rhs);

/**
 * Free resources reserved for the return flux updater.
 *
 * @param up Return flux updater object.
 */
void
gkyl_bc_return_flux_gyrokinetic_release(struct gkyl_bc_return_flux_gyrokinetic *up);
