#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/**
 * Skip cell object: manages a boolean mask array that indicates which
 * phase-space cells should be skipped during simulation updates based
 * on distribution function values falling below a threshold.
 */

enum skip_cell_types {
  GKYL_GK_SKIP_CELL_NONE = 0, // No skip cell operation.
  GKYL_GK_SKIP_CELL_F_THRESHOLD, // Skip cells based on distribution function threshold.
};

/**
 * Input file parameters for skip cell object.
 */
struct gkyl_skip_cell_inp {
  enum skip_cell_types type; // Type of skip cell operation.
  double threshold; // Skips cells where |JBf| < threshold. JBf is the distribution function output from simulations.
};

/**
 * Skip cell object definition.
 */
struct gkyl_skip_cell {
  enum skip_cell_types type; // Type of skip cell operation.
  struct gkyl_array *mask; // Boolean mask array (GKYL_BOOL type).
  double f_threshold; // Threshold for marking cells as skippable.
  struct gkyl_range phase_rng; // Phase-space range.
  bool use_gpu; // Flag indicating GPU usage.
  
  uint32_t flags;
  struct gkyl_skip_cell *on_dev; // Pointer to device object.
  
  struct gkyl_ref_count ref_count; // Reference counter.
};

/**
 * Create a new skip cell object.
 *
 * @param skip_cell_inp Input parameters for skip cell.
 * @param phase_rng Phase-space range.
 * @param use_gpu Boolean to determine if running on GPU.
 * @return New skip cell object.
 */
struct gkyl_skip_cell*
gkyl_skip_cell_new(struct gkyl_skip_cell_inp skip_cell_inp,
  struct gkyl_range phase_rng, bool use_gpu);

/**
 * Update the boolean mask array based on distribution function values.
 * Marks cells as skippable if their distribution function falls below
 * the threshold.
 *
 * @param skip_cell Skip cell updater.
 * @param distf Distribution function array.
 */
void
gkyl_skip_cell_advance(struct gkyl_skip_cell *skip_cell, const struct gkyl_array *distf);

/**
 * Update the boolean mask array with inverted logic based on distribution 
 * function values. Marks cells as 1.0 where f >= threshold, 0.0 where f < threshold.
 * This avoids floating-point arithmetic issues with inversion.
 *
 * @param skip_cell Skip cell updater.
 * @param distf Distribution function array.
 */
void
gkyl_skip_cell_advance_inverse(struct gkyl_skip_cell *skip_cell, const struct gkyl_array *distf);

/**
 * Check if skip_cell object is on CUDA device.
 *
 * @param skip_cell Skip cell object.
 * @return True if on device, false otherwise.
 */
bool 
gkyl_skip_cell_is_cu_dev(const struct gkyl_skip_cell* skip_cell);

/**
 * Acquire a reference to the skip cell object.
 * 
 * @param skip_cell Skip cell object to acquire.
 */

struct gkyl_skip_cell*
gkyl_skip_cell_acquire(struct gkyl_skip_cell *skip_cell);

/**
 * Release memory associated with skip cell object.
 *
 * @param skip_cell Skip cell object to release.
 */
void
gkyl_skip_cell_release(struct gkyl_skip_cell *skip_cell);
