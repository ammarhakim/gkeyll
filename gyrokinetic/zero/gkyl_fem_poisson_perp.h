#pragma once

#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_fem_poisson_perp gkyl_fem_poisson_perp;

/**
 * Create new updater to solve the Helmholtz problem
 *   - nabla . (epsilon * nabla phi) - kSq * phi = rho
 * using a FEM to ensure phi is continuous. This solver is also
 * used as a Poisson solver by passing a zero kSq. The input is the
 * DG field rho, which is translated to FEM. The output is the
 * DG field phi, after we've translated the FEM solution to DG.
 * Free using gkyl_fem_poisson_perp_release method.
 *
 * @param solve_range Range in which to perform the projection operation.
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param bcs Boundary conditions.
 * @param epsilon Spatially varying permittivity tensor.
 * @param kSq Squared wave number (factor multiplying phi in Helmholtz eq).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_poisson_perp* gkyl_fem_poisson_perp_new(
  const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid,
  const struct gkyl_basis basis, struct gkyl_poisson_bc *bcs, struct gkyl_array *epsilon,
  struct gkyl_array *kSq, bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM poisson updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void gkyl_fem_poisson_perp_set_rhs(gkyl_fem_poisson_perp* up, struct gkyl_array *rhsin);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_poisson_perp_solve(gkyl_fem_poisson_perp* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_poisson_perp_release(struct gkyl_fem_poisson_perp *up);
