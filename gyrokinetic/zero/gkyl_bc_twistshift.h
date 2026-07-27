#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_evalf_def.h>
#include <assert.h>

// Object type.
typedef struct gkyl_bc_twistshift gkyl_bc_twistshift;

struct gkyl_bc_twistshift_inp {
  int bc_dir; // Direction in which to apply this BC.
  struct gkyl_range *periodic_out_r; // Range of the periodic donor (output) field.
  struct gkyl_range *periodic_in_r; // Range of the periodic donor (input) field.
  int shift_dir; // Direction of the shift.
  int shear_dir; // Direction in which the shift varies (shear).
  enum gkyl_edge_loc edge; // Edge to apply this BC at (lower/upper).
  int cdim; // Configuration space dimensions.
  struct gkyl_range bcdir_ext_update_r; // Local range where to apply BC, extended in bc_dir.
  const int *num_ghost; // Number of ghost cells in each direction.
  struct gkyl_basis basis; // Basis of the field shifted.
  struct gkyl_rect_grid grid; // Grid the field shifted is defined on.
  evalf_t shift_func; // Function defining the shift.
  void *shift_func_ctx; // Context for shift_func.
  struct gkyl_array *shift_dg; // Discretized shift.
  bool use_gpu; // Whether to apply the BC using the GPU.
  // Optional inputs:
  int shift_poly_order; // Basis order for the DG representation of the shift.
  // Optional anti-aliasing: a low-pass filter along shear_dir (enabled when
  // filter_half_width > 0), applied on a grid supersampled by upsample_factor
  // along shear_dir (used only when the filter is on). Both act on the shifted
  // ghost plane only.
  int filter_half_width; // Filter stencil half-width M in cells (0 = no filter).
  double filter_cutoff_wavelength; // Filter cutoff wavelength (physical units).
  int upsample_factor; // Supersampling factor along shear_dir (0/1 = none).
};

/**
 * Create a new updater to apply the twist-shift BC, optionally with
 * supersampling and low-pass filtering to de-alias the shifted field. It
 * combines a gkyl_twistshift_dg updater with, when requested, a
 * gkyl_dg_lowpass_filter and gkyl_dg_interpolate operators.
 *
 * @param inp bc_twistshift_inp struct containing the inputs to the updater.
 * @return New updater pointer.
 */
struct gkyl_bc_twistshift* gkyl_bc_twistshift_new(const struct gkyl_bc_twistshift_inp *inp);

/**
 * Apply the twist-shift periodic BC. Can be used in-place.
 *
 * @param up Twist-shift BC updater object.
 * @param fdo Donor field.
 * @param ftar Target field.
 */
void gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Free memory associated with the bc_twistshift updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up);
