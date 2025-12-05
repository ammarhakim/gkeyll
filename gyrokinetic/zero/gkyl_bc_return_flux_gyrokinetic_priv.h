#pragma once

#include <gkyl_bc_return_flux_gyrokinetic.h>

// Private header for bc_return_flux_gyrokinetic, not for direct use in user code.

// Primary struct in this updater.
struct gkyl_bc_return_flux_gyrokinetic {
  int dir; // Direction perpendicular to boundary.
  enum gkyl_edge_loc edge; // Edge of the boundary (lower or upper).
  int disp_dir; // Direction in which to displace flux.
  int disp_dir_num_cellsD2; // Number of cells in displacement direction / 2.
  const struct gkyl_range *bflux_ghost_r; // Ghost range to index bflux.
  const struct gkyl_range *skin_r; // Skin range.
  bool use_gpu; // Whether to apply the BC with the GPU.
};

#ifdef GKYL_HAVE_CUDA

/**
 * Apply the return flux BC on the GPU.
 *
 * @param up Return flux updater object.
 * @param bflux Boundary flux.
 * @param rhs Time rate of change (df/dt) to add displaced bflux to.
 */
void
gkyl_bc_return_flux_gyrokinetic_advance_cu(const struct gkyl_bc_return_flux_gyrokinetic *up,
  const struct gkyl_array *bflu, struct gkyl_array *rhs);

#endif
