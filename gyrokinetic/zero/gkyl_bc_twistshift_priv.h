#pragma once

// Private header for the bc_twistshift orchestrator, not for direct use in
// user code.

#include <gkyl_bc_twistshift.h>
#include <gkyl_twistshift_dg.h>
#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>

// Primary struct in this updater.
struct gkyl_bc_twistshift {
  bool use_gpu; // Whether to apply the BC on the GPU.

  struct gkyl_twistshift_dg *ts; // Pure twist-shift updater.
  int bc_dir; // Direction along which we treat the BC.
  struct gkyl_range *periodic_out_r; // Range of the periodic donor (output) field.
  struct gkyl_range *periodic_in_r; // Range of the periodic donor (input) field.
  
  struct gkyl_dg_lowpass_filter *filter; // Optional post-shift filter along shear_dir.
  int filter_half_width; // Filter stencil half-width M in cells (0 = no filter).
  double filter_cutoff_wavelength; // Filter cutoff wavelength.
  struct gkyl_range ghost_r; // Ghost plane the twist-shift fills.
  struct gkyl_array *filt_buff; // Buffer for the filter.
  void (*filter_func)(struct gkyl_dg_lowpass_filter *filt_up, struct gkyl_array *GKYL_RESTRICT finout, struct gkyl_array *GKYL_RESTRICT fbuff);

  struct gkyl_rect_grid ts_grid; // Grid refined along shear_dir.
  struct gkyl_range ts_update_r; // Update range on the ts_grid.
  int upsample_factor; // Supersampling factor along shear_dir.
  struct gkyl_dg_interpolate *refine; // Coarse ghost plane -> fine ghost plane.
  struct gkyl_dg_interpolate *coarsen; // Fine ghost plane -> coarse ghost plane.
  struct gkyl_array *ffine; // Field on the fine grid.
  struct gkyl_range coarse_ghost_r; // Ghost plane on the coarse grid (donor and target).
  void (*refine_func)(struct gkyl_dg_interpolate *refine, struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar);
  void (*coarsen_func)(struct gkyl_dg_interpolate *coarsen, struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar);
};
