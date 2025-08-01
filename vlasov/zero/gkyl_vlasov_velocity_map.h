#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_vlasov_velocity_map_inp {
  evalf_t eval_vmap; // Velocity mapping in each velocity-space dimension used by dg_basis_ops to construct C^1 mapping.
  void *ctx; // Context for function evaluation. Can be NULL.
};

/**
 * Construct a velocity map using a cubic C^1 representation based on input
 * mapping function. Returns vdim 1D mappings and vdim 1D quadratic 
 * representations of the inverse Jacobian. 
 *
 * @param vgrid Velocity-space grid object
 * @param vrange Velocity-space range
 * @param inp_vmap[GKYL_MAX_CDIM] Velocity mapping input (function and context) in each velocity-space dimension
 * @param vmap C^1 cubic representation of mapping in each velocity-space dimension.
 * @param jacob_vel Jacobian (derivative of vmap) in each velocity space dimension at 1V Gauss-Legendre quadrature points.
 * @param vmap_pgkyl C^1 cubic representation of mapping used for I/O (defined in the full 1V, 2V, or 3V).
 * @param jacob_vel_pgkyl Total Jacobian (derivative of vmap) used for I/O (defined in the full 1V, 2V, or 3V).
 * @param vmap_avg_pgkyl Cell average of the C^1 cubic representation of mapping used for I/O (defined in the full 1V, 2V, or 3V).
 * @param jacob_vel_avg_pgkyl Cell average of the Total Jacobian (derivative of vmap) used for I/O (defined in the full 1V, 2V, or 3V).
 * @param jacob_vel_gauss Total velocity space Jacobian evaluated at Gauss-Legendre quadrature points.
 */
void gkyl_vlasov_velocity_map_new(const struct gkyl_rect_grid *vgrid, 
  const struct gkyl_range *vrange, struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM], 
  struct gkyl_array *vmap, struct gkyl_array *jacob_vel, 
  struct gkyl_array *vmap_pgkyl, struct gkyl_array *jacob_vel_pgkyl, 
  struct gkyl_array *vmap_avg_pgkyl, struct gkyl_array *jacob_vel_avg_pgkyl, 
  struct gkyl_array *jacob_vel_gauss);