#pragma once

#include <gkyl_gr_spacetime.h>

// Forward declaration.
struct gkyl_tov;

/**
 * Create a TOV spacetime object.  For r >= R_star (from the star center),
 * all queries are delegated to an internal Schwarzschild (spin=0) black hole
 * spacetime with M = M_star.  For r < R_star, lapse, shift, spatial metric,
 * extrinsic curvature, and determinant come from the TOV solver in CKS coords.
 * The excision region function always returns false.
 *
 * @param use_gpu  Whether to use GPU.
 * @param tov      Pointer to an already-solved TOV object (borrowed, not owned).
 * @param pos_x    Star center x-coordinate.
 * @param pos_y    Star center y-coordinate.
 * @param pos_z    Star center z-coordinate.
 * @return New spacetime object.
 */
struct gkyl_gr_spacetime*
gkyl_gr_tov_spacetime_new(bool use_gpu, struct gkyl_tov *tov,
  double pos_x, double pos_y, double pos_z);
