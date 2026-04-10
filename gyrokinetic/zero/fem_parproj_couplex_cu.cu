/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fem_parproj_couplex.h>
#include <gkyl_fem_parproj_couplex_priv.h>
}

// CUDA kernel to set device pointers to l2g, RHS src and solution
// kernels. Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
fem_parproj_couplex_set_cu_ker_ptrs(const struct gkyl_basis basis, bool has_weight_lhs, bool has_weight_rhs,
  bool isperiodic, bool isdirichlet, struct gkyl_fem_parproj_couplex_kernels *kers)
{
  // Local-to-global kernel.
  int bckey[1] = {0};
  bckey[0] = isperiodic? 0 : 1;

  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<4; k++)
        kers->l2g[k] = CK(ser_loc2glob_list, basis.ndim, bckey[0], basis.poly_order, k);
      break;
    default:
      assert(false);
      break;
  }

  // RHS source kernel.
  bckey[0] = isdirichlet? 1 : 0;

  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<9; k++)
        kers->srcker[k] = has_weight_rhs? CK(ser_srcstencil_list_weighted, basis.ndim, bckey[0], basis.poly_order, k)
                                        : CK(ser_srcstencil_list_noweight, basis.ndim, bckey[0], basis.poly_order, k);
      break;
    default:
      assert(false);
      break;
  }

  // Solution kernel.
  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->solker = ser_solstencil_list[basis.ndim-1].kernels[basis.poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

void
fem_parproj_couplex_choose_kernels_cu(const struct gkyl_basis* basis, bool has_weight_lhs, bool has_weight_rhs,
  bool isperiodic, bool isdirichlet, struct gkyl_fem_parproj_couplex_kernels *kers)
{
  fem_parproj_couplex_set_cu_ker_ptrs<<<1,1>>>(*basis, has_weight_lhs, has_weight_rhs, isperiodic, isdirichlet, kers);
}

__global__ void
gkyl_fem_parproj_couplex_set_rhs_kernel(int pardir, double *rhs_global, const struct gkyl_array *rhsin,
  const struct gkyl_array *weight, const struct gkyl_array *phibc,
  struct gkyl_range range, struct gkyl_range perp_range, struct gkyl_range fem_range,
  struct gkyl_fem_parproj_couplex_kernels *kers, long numnodes_global)
{
  int idx1[GKYL_MAX_CDIM], idx0[GKYL_MAX_CDIM], fem_idx[GKYL_MAX_CDIM];;
  long globalidx[32];
  int ndim = range.ndim;
  int num_cells[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; d++)
    num_cells[d] = range.upper[d] - range.lower[d] + 1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx1);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linidx = gkyl_range_idx(&range, idx1);

    fem_idx[0] = idx1[0];
    fem_idx[GKYL_MIN2(pardir,1)] = idx1[pardir];

    for (size_t d=0; d<ndim; d++) idx0[d] = idx1[d]-range.lower[d];

    int keri = idx_to_inup_ker(fem_range.ndim, fem_range.upper, fem_idx);
    kers->l2g[keri](num_cells, idx0, globalidx);

    const double *wgt_p = weight? (const double *) gkyl_array_cfetch(weight, linidx) : NULL;
    const double *phibc_p = phibc? (const double *) gkyl_array_cfetch(phibc, linidx) : NULL;
    const double *rhsin_p = (const double*) gkyl_array_cfetch(rhsin, linidx);

    int idx_perp[] = {perp_range.lower[0]};
    for (int d=1; d<range.ndim-1; d++) idx_perp[d] = idx1[d];
    long linidx_perp = gkyl_range_idx(&perp_range, idx_perp);

    long perpProbOff = linidx_perp*numnodes_global;

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(fem_range.ndim, fem_range.lower, fem_range.upper, fem_idx);
    kers->srcker[keri](wgt_p, rhsin_p, phibc_p, perpProbOff, globalidx, rhs_global);
  }

}

void
gkyl_fem_parproj_couplex_set_rhs_cu(gkyl_fem_parproj_couplex *up, const struct gkyl_array *rhsin,
  const struct gkyl_array *phibc)
{
  gkyl_culinsolver_clear_rhs(up->prob_cu, 0);
  double *rhs_cu = gkyl_culinsolver_get_rhs_ptr(up->prob_cu, 0);

  const struct gkyl_array *phibc_cu = phibc? phibc->on_dev : NULL;
  const struct gkyl_array *wgt_cu = up->has_weight_rhs? up->weight_rhs->on_dev : NULL;

  gkyl_fem_parproj_couplex_set_rhs_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(up->pardir, rhs_cu,
    rhsin->on_dev, wgt_cu, phibc_cu, *up->solve_range, up->perp_range, up->fem_range,
    up->kernels_cu, up->numnodes_global);
}

__global__ void
gkyl_fem_parproj_couplex_get_sol_kernel(int pardir, struct gkyl_array *phiout, const double *x_global,
  struct gkyl_range range, struct gkyl_range perp_range, struct gkyl_range fem_range,
  struct gkyl_fem_parproj_couplex_kernels *kers, long numnodes_global)
{
  int idx1[GKYL_MAX_CDIM], idx0[GKYL_MAX_CDIM], fem_idx[GKYL_MAX_CDIM];;
  long globalidx[32];
  int ndim = range.ndim;
  int num_cells[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; d++)
    num_cells[d] = range.upper[d] - range.lower[d] + 1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx1);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linidx = gkyl_range_idx(&range, idx1);

    fem_idx[0] = idx1[0];
    fem_idx[GKYL_MIN2(pardir,1)] = idx1[pardir];

    for (size_t d=0; d<ndim; d++) idx0[d] = idx1[d]-range.lower[d];

    int keri = idx_to_inup_ker(fem_range.ndim, fem_range.upper, fem_idx);
    kers->l2g[keri](num_cells, idx0, globalidx);

    double *phiout_p = (double*) gkyl_array_cfetch(phiout, linidx);

    int idx_perp[] = {perp_range.lower[0]};
    for (int d=1; d<range.ndim-1; d++) idx_perp[d] = idx1[d];
    long linidx_perp = gkyl_range_idx(&perp_range, idx_perp);

    long perpProbOff = linidx_perp*numnodes_global;

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    kers->solker(x_global, perpProbOff, globalidx, phiout_p);
  }
}

void
gkyl_fem_parproj_couplex_solve_cu(gkyl_fem_parproj_couplex *up, struct gkyl_array *phiout)
{
  gkyl_culinsolver_solve(up->prob_cu);
  double *x_cu = gkyl_culinsolver_get_sol_ptr(up->prob_cu, 0);

  gkyl_fem_parproj_couplex_get_sol_kernel<<<phiout->nblocks, phiout->nthreads>>>(up->pardir, phiout->on_dev, x_cu,
    *up->solve_range, up->perp_range, up->fem_range, up->kernels_cu, up->numnodes_global);
}
