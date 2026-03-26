#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_fem_parproj_bcs.h>

// Object type
typedef struct gkyl_fem_parproj_couplex gkyl_fem_parproj_couplex;

/**
 * Create new updater to project a DG field onto the FEM (nodal) basis
 * in order to make the field continuous. That is, we solve
 *    wgtL*phi_{fem} \doteq wgtR*rho_{dg}
 * where wgt is the weight field (optional), phi_{fem} is the (continuous field)
 * we wish to compute, rho_{dg} is the (discontinuous) input source field,
 * and \doteq implies weak equality with respect to the FEM basis.
 *
 * The weak equality above is non-local, and couples cells along x
 * and z, i.e. along the radial direction and along the field line in
 * gyrokinetic. It doesn't enforce BCs along x, and z-boundaries are handled
 * according to 'bctype'.
 *
 * Free using gkyl_fem_parproj_couplex_release method.
 *
 * @param solve_range Range in which to perform the projection operation.
 * @param basis Basis functions of the DG field.
 * @param bctype Type of boundary condition (see gkyl_fem_parproj_bcs).
 * @param weight_left Weight on left-side of the operator (time-independent).
 * @param weight_right Weight on right-side of the operator (time-independent).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_parproj_couplex* gkyl_fem_parproj_couplex_new(const struct gkyl_range *solve_range,
  const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype,
  const struct gkyl_array *weight_left, const struct gkyl_array *weight_right,
  bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs (only use ghost cells).
 */
void gkyl_fem_parproj_couplex_set_rhs(struct gkyl_fem_parproj_couplex* up,
  const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_couplex_solve(struct gkyl_fem_parproj_couplex* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_parproj_couplex_release(struct gkyl_fem_parproj_couplex *up);
