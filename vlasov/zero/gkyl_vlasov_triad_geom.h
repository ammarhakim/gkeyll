#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_vlasov_triad_geom_inp {
  evalf_t eval_cov_tangent_basis; // The covariant tangent basis to be evaluated within each configuration space cell.
  evalf_t eval_triad_basis; // The triad basis to be evaluated within each configuration space cell.
  evalf_t eval_triad_basis_gradient; // The triad basis gradient to be evaluated within each configuration space cell.
  void *eval_cov_tangent_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_ctx; // Context for function evaluation. Can be NULL.
  void *eval_triad_basis_gradient_ctx; // Context for function evaluation. Can be NULL.
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
 * @param inp_basis_vectors Basis mapping input (function and context) cov_tangent_basis, traid_basis
 * @param conf_poisson_tensor The configuration component of the Poisson tensor
 */
void gkyl_vlasov_triad_geom_new(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis, 
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors, struct gkyl_array *conf_poisson_tensor);