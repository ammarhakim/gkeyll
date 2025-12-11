#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/**
 * Reduce a DG array along one direction, computing max/min and argmax/argmin.
 * 
 * For a 2D input array f(psi, z), reducing along z (dir=1) gives:
 *   out_val(psi) = max_z f(psi, z)    or    min_z f(psi, z)
 *   out_coord(psi) = argmax_z f(psi, z)  or  argmin_z f(psi, z)
 * 
 * For a 1D input array f(z), reducing along z (dir=0) gives scalars:
 *   out_val = max f   or   min f
 *   out_coord = argmax f  or  argmin f
 * 
 * The reduction is computed by sampling the field at nodal points along the
 * reduction direction and finding the global extremum.
 */
typedef struct gkyl_array_dg_reduce_dir gkyl_array_dg_reduce_dir;

/** Type of reduction operation. */
enum gkyl_reduce_op {
  GKYL_REDUCE_OP_MAX,   // Find maximum
  GKYL_REDUCE_OP_MIN,   // Find minimum
};

/** Input parameters for dg_reduce_dir updater. */
struct gkyl_array_dg_reduce_dir_inp {
  const struct gkyl_basis *basis;       // Input basis (N-dimensional)
  const struct gkyl_rect_grid *grid;    // Input grid
  const struct gkyl_range *range;       // Input range (local)
  const struct gkyl_range *range_ext;   // Input extended range
  int reduce_dir;                       // Direction to reduce (0-indexed)
  enum gkyl_reduce_op op;               // Reduction operation (MAX or MIN)
  bool use_gpu;                         // Whether to run on GPU
};

/**
 * Create a new directional reduction updater.
 * 
 * @param inp Input parameters
 * @return New updater pointer
 */
struct gkyl_array_dg_reduce_dir* gkyl_array_dg_reduce_dir_new(
  const struct gkyl_array_dg_reduce_dir_inp *inp);

/**
 * Compute the reduction. For each point along the preserved dimensions,
 * find the maximum or minimum along the reduction direction.
 * 
 * @param up Updater object
 * @param in Input array (N-dimensional DG field)
 */
void gkyl_array_dg_reduce_dir_advance(struct gkyl_array_dg_reduce_dir *up,
  const struct gkyl_array *in);

/**
 * Evaluate a field at the extremum coordinates.
 * For example, if we found z_max(psi) = argmax_z f(psi,z), this evaluates
 * g(psi, z_max(psi)) for any field g.
 * 
 * @param up Updater object
 * @param in_array Input field to evaluate (N-dimensional)
 * @param out_val Output values at extremum coordinates (N-1 dimensional)
 */
void gkyl_array_dg_reduce_dir_eval_at_extremum(struct gkyl_array_dg_reduce_dir *up,
  const struct gkyl_array *in_array, struct gkyl_array *out_val);

/**
 * Get the output basis ((N-1)-dimensional, or p=0 1D for 1D->0D).
 * 
 * @param up Updater object
 * @return Pointer to output basis
 */
const struct gkyl_basis* gkyl_array_dg_reduce_dir_get_basis(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output grid.
 * 
 * @param up Updater object
 * @return Pointer to output grid
 */
const struct gkyl_rect_grid* gkyl_array_dg_reduce_dir_get_grid(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output range.
 * 
 * @param up Updater object
 * @return Pointer to output range
 */
const struct gkyl_range* gkyl_array_dg_reduce_dir_get_range(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output extended range.
 * 
 * @param up Updater object
 * @return Pointer to output extended range
 */
const struct gkyl_range* gkyl_array_dg_reduce_dir_get_range_ext(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output nodal range.
 * 
 * @param up Updater object
 * @return Pointer to output nodal range
 */
const struct gkyl_range* gkyl_array_dg_reduce_dir_get_nodal_range(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output array containing extremal values (modal DG expansion).
 * 
 * @param up Updater object
 * @return Pointer to output values array
 */
const struct gkyl_array* gkyl_array_dg_reduce_dir_get_vals(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output array containing extremal values (nodal representation).
 * 
 * @param up Updater object
 * @return Pointer to output nodal values array
 */
const struct gkyl_array* gkyl_array_dg_reduce_dir_get_vals_nodal(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output array containing coordinates of extrema (modal DG expansion).
 * 
 * @param up Updater object
 * @return Pointer to output coordinates array
 */
const struct gkyl_array* gkyl_array_dg_reduce_dir_get_coords(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Get the output array containing coordinates of extrema (nodal representation).
 * 
 * @param up Updater object
 * @return Pointer to output nodal coordinates array
 */
const struct gkyl_array* gkyl_array_dg_reduce_dir_get_coords_nodal(
  const struct gkyl_array_dg_reduce_dir *up);

/**
 * Release the updater and all internal arrays.
 * 
 * @param up Updater to delete
 */
void gkyl_array_dg_reduce_dir_release(struct gkyl_array_dg_reduce_dir *up);
