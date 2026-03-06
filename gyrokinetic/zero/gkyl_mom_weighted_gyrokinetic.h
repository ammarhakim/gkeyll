#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

enum gkyl_mom_weight_type {
  GKYL_F_MOMENT_WEIGHT_NONE = 0, // No weight.
  GKYL_F_MOMENT_WEIGHT_CONF, // Weight is a function of conf-space.
  GKYL_F_MOMENT_WEIGHT_VEL, // Weight is a function of v-space.
  GKYL_F_MOMENT_WEIGHT_PHASE, // Weight is a function of phase-space.
};

// Object type.
typedef struct gkyl_mom_weighted_gyrokinetic gkyl_mom_weighted_gyrokinetic;

/**
 * Create new updated to compute weighted gyrokinetic moments.
 *
 * @param cbasis Configuration-space basis-functions.
 * @param pbasis Phase-space basis-functions.
 * @param mass Mass of species.
 * @param charge Charge of species.
 * @param vel_map Velocity space mapping object.
 * @param gk_geom Geometry object.
 * @param mom_type Name of moment to compute.
 * @param wgt_type Name of moment to compute.
 * @param is_integrated Whether to compute the volume integrated moment.
 * @param use_gpu bool to determine if on GPU.
 * @return A new weighted moment updater.
 */
struct gkyl_mom_weighted_gyrokinetic* gkyl_mom_weighted_gyrokinetic_new(double mass, double charge,
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, struct gkyl_rect_grid *phase_grid,
  const struct gkyl_velocity_map *vel_map, const struct gk_geometry *gk_geom,
  enum gkyl_distribution_moments mom_type, enum gkyl_mom_weight_type wgt_type,
  bool is_integrated, bool use_gpu);

/**
 * Compute a weighted moment of the gyrokinetic distribution function.
 * The phase_rng and conf_rng MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param up Weighted moment updater.
 * @param phase_rng Phase-space range.
 * @param conf_rng Config-space range.
 * @param wgt_rng Weight range.
 * @param phi Electrostatic potential (for Hamiltonian moment).
 * @param wgt Weight.
 * @param fin Input distribution function array.
 * @param mout Output moment array.
 */
void gkyl_mom_weighted_gyrokinetic_advance(struct gkyl_mom_weighted_gyrokinetic *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng, const struct gkyl_range *wgt_rng,
  struct gkyl_array *phi, struct gkyl_array *wgt, const struct gkyl_array *GKYL_RESTRICT fin,
  struct gkyl_array *GKYL_RESTRICT mout);

/**
 * Free memory associated with the weighted moment calculation.
 *
 * @param up Weighted moment updater.
 */
void gkyl_mom_weighted_gyrokinetic_release(struct gkyl_mom_weighted_gyrokinetic *up);
