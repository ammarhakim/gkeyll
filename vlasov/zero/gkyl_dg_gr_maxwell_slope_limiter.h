#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>

// Object type
typedef struct gkyl_dg_gr_maxwell_slope_limiter gkyl_dg_gr_maxwell_slope_limiter;

// Inputs packaged in a struct.
struct gkyl_dg_gr_maxwell_slope_limiter_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  double limiter_fac; // Slope-to-cell-avg-jump scaling factor. Pass 0.0 to take the
                      // default 1/sqrt(3) (= ratio of orthonormal slope basis to
                      // adjacent-cell average jump). Decrease for more diffusion,
                      // increase for less.
  const int *theta_pole_lo; // (Lower side) 1 if direction has theta-pole BC. Default NULL = zeros.
  const int *theta_pole_up; // (Upper side) 1 if direction has theta-pole BC. Default NULL = zeros.
  bool use_gpu;
};

/**
 * Create a new tetrad-frame characteristic slope limiter updater for the DG
 * GR Maxwell scheme. Applies a 3-argument minmod limiter to slope modes
 * (q[1], q[2], q[3]) of each cell, projected into characteristic variables
 * via paper eq. 57 right eigenvectors in the locally flat tetrad. Cell
 * averages (q[0]) are preserved exactly.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_dg_gr_maxwell_slope_limiter*
gkyl_dg_gr_maxwell_slope_limiter_inew(const struct gkyl_dg_gr_maxwell_slope_limiter_inp *inp);

/**
 * Apply the slope limiter in-place over the configuration-space range. Runs
 * one r-direction pass, then one theta-direction pass. q[0] is preserved.
 *
 * @param up Slope limiter updater.
 * @param conf_range Configuration-space range (interior cells).
 * @param lapse Surface nodal lapse arrays.
 * @param shift Surface nodal shift arrays.
 * @param h_ij Surface nodal spatial 3-metric arrays.
 * @param det_h Surface nodal det(h) arrays.
 * @param field_con In/out contravariant (JD, JB) field array.
 */
void
gkyl_dg_gr_maxwell_slope_limiter_advance(struct gkyl_dg_gr_maxwell_slope_limiter *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *lapse,
  const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij,
  const struct gkyl_surf_and_vol_node_arrays *det_h,
  struct gkyl_array *field_con);

/**
 * Release updater.
 */
void
gkyl_dg_gr_maxwell_slope_limiter_release(struct gkyl_dg_gr_maxwell_slope_limiter *up);
