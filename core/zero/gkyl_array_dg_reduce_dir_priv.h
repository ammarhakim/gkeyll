#pragma once

#include <float.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_dg_reduce_dir.h>
#include <gkyl_nodal_ops.h>

/**
 * Convert logical (reference) coordinates to computational (physical) coordinates.
 * xout[d] = xc[d] + 0.5*dx[d]*eta[d]
 */
static inline void
dg_reduce_dir_log_to_comp(int ndim, const double *eta,
  const double *GKYL_RESTRICT dx, const double *GKYL_RESTRICT xc,
  double *GKYL_RESTRICT xout)
{
  for (int d = 0; d < ndim; ++d)
    xout[d] = 0.5*dx[d]*eta[d] + xc[d];
}

/** Internal struct for dg_reduce_dir updater. */
struct gkyl_array_dg_reduce_dir {
  // Input parameters (copies).
  struct gkyl_rect_grid grid;       // Input grid (copy)
  struct gkyl_basis basis;          // Input basis (copy)
  struct gkyl_range range;          // Input local range (copy)
  struct gkyl_range range_ext;      // Input extended range (copy)
  int reduce_dir;                   // Direction to reduce
  enum gkyl_reduce_op op;           // Reduction operation (MAX or MIN)
  bool use_gpu;

  // Output grid/basis/range (owned).
  struct gkyl_rect_grid out_grid;   // Output grid (N-1 dim, or 1D 1-cell for 1D->0D)
  struct gkyl_basis out_basis;      // Output basis (N-1 dim, or p=0 1D for 1D->0D)
  struct gkyl_range out_range;      // Output range
  struct gkyl_range out_range_ext;  // Output extended range
  struct gkyl_range out_nrange;     // Nodal range for output

  // Output arrays (owned).
  struct gkyl_array *out_vals;          // Extremal values (modal DG)
  struct gkyl_array *out_coords;        // Extremal coordinates (modal DG)
  struct gkyl_array *out_vals_nodal;    // Nodal extremal values
  struct gkyl_array *out_coords_nodal;  // Nodal extremal coordinates
  struct gkyl_array *out_eval_nodal;    // Nodal array for eval_at_extremum

  // Internal working arrays.
  struct gkyl_array *nodes;         // Node locations in logical coords

  // Nodal-to-modal converter.
  struct gkyl_nodal_ops *n2m;
};
