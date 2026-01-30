#pragma once

#include <gkyl_dg_array_mask.h>
#include <gkyl_gyrokinetic.h>

/**
 * Translate gyrokinetic skip cell type to DG array mask type.
 *
 * This function centralizes the translation from the gyrokinetic-specific
 * enum (used in simulation input) to the generic DG array mask type
 * (used internally by the mask updater).
 *
 * @param gk_type Gyrokinetic skip cell type from simulation input.
 * @return Corresponding DG array mask type.
 */
static inline enum gkyl_dg_array_mask_types
gkyl_gk_skip_cell_to_mask_type(enum gkyl_gyrokinetic_skip_cell_type gk_type)
{
  switch (gk_type) {
    case GKYL_GK_SKIP_CELL_BELOW:
      return GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_THRESHOLD;
    case GKYL_GK_SKIP_CELL_ABOVE:
      return GKYL_DG_ARRAY_MASK_C0_LESS_THAN_THRESHOLD;
    case GKYL_GK_SKIP_CELL_BELOW_FRAC_CONF:
      return GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD_SPATIAL;
    case GKYL_GK_SKIP_CELL_ABOVE_FRAC_CONF:
      return GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD_SPATIAL;
    case GKYL_GK_SKIP_CELL_BELOW_FRAC:
      return GKYL_DG_ARRAY_MASK_C0_GREATER_THAN_FRAC_THRESHOLD;
    case GKYL_GK_SKIP_CELL_ABOVE_FRAC:
      return GKYL_DG_ARRAY_MASK_C0_LESS_THAN_FRAC_THRESHOLD;
    case GKYL_GK_SKIP_CELL_NONE:
    default:
      return GKYL_DG_ARRAY_MASK_NONE;
  }
}
