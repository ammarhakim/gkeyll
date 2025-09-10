#pragma once

#include <gkyl_boundary_flux.h>

struct gkyl_boundary_flux {
  int dir; // Direction perpendicular to the sheath boundary.
  enum gkyl_edge_loc edge; // Lower or upper boundary.
  struct gkyl_rect_grid grid; // Phase-space grid object.
  struct gkyl_range skin_r, ghost_r; // Skin and ghost ranges.
  const struct gkyl_dg_eqn *equation; // Equation object.
  bool use_gpu; // Whether to run on GPU.
  double skip_cell_threshold; // Threshold for skipping cells in the skin range.

  uint32_t flags;
  struct gkyl_boundary_flux *on_dev; // pointer to itself or device data
};


#ifdef GKYL_HAVE_CUDA
/**
 * Create new boundary_flux updater on the GPU.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param grid Grid object (on device)
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param equation Equation object
 * @param skip_cell_threshold Threshold for skipping cells.
 * @return New updater pointer.
 */
gkyl_boundary_flux*
gkyl_boundary_flux_cu_dev_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_dg_eqn *equation, double skip_cell_threshold);

/**
 * Compute the boundary flux on the GPU.
 *
 * @param up Boundary flux updater.
 * @param fIn Input distribution function.
 * @param fluxOut Output flux.
 */
void gkyl_boundary_flux_advance_cu(gkyl_boundary_flux *up,
  const struct gkyl_array *fIn, struct gkyl_array *fluxOut);
#endif
