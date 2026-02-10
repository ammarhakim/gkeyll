#pragma once

#include <gkyl_dg_array_mask.h>

enum gkyl_gyrokinetic_skip_cell_type {
  GKYL_GK_SKIP_CELL_NONE = 0, // No skipping of cells.
  GKYL_GK_SKIP_CELL_BELOW, // Skip cells where |J_tot f| < threshold. J_tot f is the distribution function output from simulations.
  GKYL_GK_SKIP_CELL_ABOVE, // Skip cells where |J_tot f| > threshold.
  GKYL_GK_SKIP_CELL_BELOW_FRAC, // Skip cells where |J_tot f| < threshold * max|J_tot f|. J_tot f is the distribution function output from simulations. Uses a global maximum of J_tot f.
  GKYL_GK_SKIP_CELL_ABOVE_FRAC, // Skip cells where |J_tot f| > threshold * max|J_tot f|.
  GKYL_GK_SKIP_CELL_BELOW_FRAC_CONF, // Skip cells where |J_tot f| < threshold * max|J_tot f(x)|. J_tot f is the distribution function output from simulations. This uses a spatially dependent mask which find the local maximum of J_tot f at each configuration space cell and applies the mask as a fraction of that.
  GKYL_GK_SKIP_CELL_ABOVE_FRAC_CONF // Skip cells where |J_tot f| > threshold * max|J_tot f(x)|.
};

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
      return GKYL_DG_ARRAY_MASK_C0_GREATER;
    case GKYL_GK_SKIP_CELL_ABOVE:
      return GKYL_DG_ARRAY_MASK_C0_LESS;
    case GKYL_GK_SKIP_CELL_BELOW_FRAC:
      return GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC;
    case GKYL_GK_SKIP_CELL_ABOVE_FRAC:
      return GKYL_DG_ARRAY_MASK_C0_LESS_FRAC;
    case GKYL_GK_SKIP_CELL_BELOW_FRAC_CONF:
      return GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF;
    case GKYL_GK_SKIP_CELL_ABOVE_FRAC_CONF:
      return GKYL_DG_ARRAY_MASK_C0_LESS_FRAC_CONF;
    case GKYL_GK_SKIP_CELL_NONE:
    default:
      return GKYL_DG_ARRAY_MASK_NONE;
  }
}
