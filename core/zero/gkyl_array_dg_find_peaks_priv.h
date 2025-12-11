#pragma once

#include <float.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_dg_find_peaks.h>
#include <gkyl_nodal_ops.h>

// Maximum number of peaks we can handle.
#define GKYL_DG_FIND_PEAKS_MAX 16

/**
 * Convert logical (reference) coordinates to computational (physical) coordinates.
 * xout[d] = xc[d] + 0.5*dx[d]*eta[d]
 */
static inline void
dg_find_peaks_log_to_comp(int ndim, const double *eta,
  const double *GKYL_RESTRICT dx, const double *GKYL_RESTRICT xc,
  double *GKYL_RESTRICT xout)
{
  for (int d = 0; d < ndim; ++d)
    xout[d] = 0.5*dx[d]*eta[d] + xc[d];
}

/** Internal struct for dg_find_peaks updater. */
struct gkyl_array_dg_find_peaks {
  // Input parameters (copies).
  struct gkyl_rect_grid grid;       // Input grid (copy)
  struct gkyl_basis basis;          // Input basis (copy)
  struct gkyl_range range;          // Input local range (copy)
  struct gkyl_range range_ext;      // Input extended range (copy)
  int search_dir;                   // Direction to search for peaks
  bool use_gpu;

  // Output grid/basis/range (owned).
  struct gkyl_rect_grid out_grid;   // Output grid (N-1 dim, or 1D 1-cell for 1D->0D)
  struct gkyl_basis out_basis;      // Output basis (N-1 dim, or p=0 1D for 1D->0D)
  struct gkyl_range out_range;      // Output range
  struct gkyl_range out_range_ext;  // Output extended range
  struct gkyl_range out_nrange;     // Nodal range for output

  // Peak information.
  int num_peaks;                              // Number of peaks detected
  enum gkyl_peak_type peak_types[GKYL_DG_FIND_PEAKS_MAX]; // Type of each peak

  // Output arrays (owned) - one per peak.
  struct gkyl_array *out_vals[GKYL_DG_FIND_PEAKS_MAX];        // Peak values (modal DG)
  struct gkyl_array *out_coords[GKYL_DG_FIND_PEAKS_MAX];      // Peak coordinates (modal DG)
  struct gkyl_array *out_vals_nodal[GKYL_DG_FIND_PEAKS_MAX];  // Nodal peak values
  struct gkyl_array *out_coords_nodal[GKYL_DG_FIND_PEAKS_MAX]; // Nodal peak coordinates

  // Internal working arrays.
  struct gkyl_array *nodes;         // Node locations in logical coords

  // Nodal-to-modal converter.
  struct gkyl_nodal_ops *n2m;
};

