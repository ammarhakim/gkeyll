#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

// Object type.
typedef struct gkyl_mom_gyrokinetic gkyl_mom_gyrokinetic;

/**
 * Create new updater to compute gyrokinetic moments.
 *
 * @param cbasis Configuration-space basis-functions.
 * @param pbasis Phase-space basis-functions.
 * @param mass Mass of species.
 * @param charge Charge of species.
 * @param vel_map Velocity space mapping object.
 * @param gk_geom Geometry object.
 * @param mom_type Name of moment to compute.
 * @param is_integrated Whether to compute the volume integrated moment.
 * @param use_gpu bool to determine if on GPU.
 * @return A new moment updater.
 */
struct gkyl_mom_gyrokinetic* gkyl_mom_gyrokinetic_new(double mass, double charge,
  const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_velocity_map *vel_map, const struct gk_geometry *gk_geom,
  enum gkyl_distribution_moments mom_type, bool is_integrated, bool use_gpu);

/**
 * Fetch the number of moments
 *
 * @param moment moment updater object
 * 
 * @return number of moments
 */
int gkyl_mom_gyrokinetic_num_mom(const struct gkyl_mom_gyrokinetic* up);

/**
 * Compute a moment of the gyrokinetic distribution function.
 * The phase_rng and conf_rng MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param up Moment updater.
 * @param phase_rng Phase-space range.
 * @param conf_rng Config-space range.
 * @param phi Electrostatic potential (for Hamiltonian moment).
 * @param fin Input distribution function array.
 * @param mout Output moment array.
 */
void gkyl_mom_gyrokinetic_advance(struct gkyl_mom_gyrokinetic *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT phi, const struct gkyl_array *GKYL_RESTRICT fin,
  struct gkyl_array *GKYL_RESTRICT mout);

/**
 * Free memory associated with the moment calculation.
 *
 * @param up Moment updater.
 */
void gkyl_mom_gyrokinetic_release(struct gkyl_mom_gyrokinetic *up);
