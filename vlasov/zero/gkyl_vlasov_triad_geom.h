#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>


enum gkyl_triad_preset_geom_type {
  GKYL_TRIAD_NONE = 0,
  GKYL_TRIAD_ANNULUS,
  GKYL_TRIAD_CYLINDRICAL_RZ,
  GKYL_TRIAD_GR_KERR_SCHILD_RPHI,
  GKYL_TRIAD_GR_KERR_SCHILD_3V,
};

struct gkyl_vlasov_triad_geom_inp {
  evalf_t eval_cov_tangent_basis; // The covariant tangent basis to be evaluated within each configuration space cell.
  evalf_t eval_triad_basis; // The triad basis to be evaluated within each configuration space cell.
  evalf_t eval_triad_basis_gradient; // The triad basis gradient to be evaluated within each configuration space cell.
  void *eval_cov_tangent_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_gradient_ctx; // Context for function evaluation. Can be NULL.

  evalf_t eval_vierbein; // The vierbein to be evaluated within each configuration space cell.
  evalf_t eval_vierbein_gradient; // The vierbein gradient to be evaluated within each configuration space cell.
  void *eval_vierbein_ctx; // Context for function evaluation. Can be NULL.
  void *eval_vierbein_gradient_ctx; // Context for function evaluation. Can be NULL.

  bool use_preset_geom; // bool for determining if we have specified a preset geometry.
  bool use_vierbein; // bool for determining which geometry convention we are constructing PT from
  enum gkyl_triad_preset_geom_type triad_preset_geom_type;  // geom type for preset geometries for triads
};

/**
 * Construct the triad geometry objects
 *
 * @param cgrid Configuration-space grid object
 * @param crange Configuration-space range
 * @param cbasis Configuration-space basis
 * @param pgrid Phase-space grid object
 * @param prange Phase-space range
 * @param pbasis Phase-space basis
 * @param inp_triad_geom Basis mapping input (function and context) cov_tangent_basis, traid_basis
 * @param conf_poisson_tensor The configuration component of the Poisson tensor
 */
void gkyl_vlasov_triad_geom_new(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis, 
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor);

/**
 * Preset function for the triad hamil 
 *
 * @param preset_geom_type Preset geometry type
 */
evalf_t gkyl_vlasov_triad_preset_hamil(const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type);