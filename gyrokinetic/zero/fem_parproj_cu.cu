/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_priv.h>
}

// CUDA kernel to set device pointers to l2g, RHS src and solution
// kernels. Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
fem_parproj_set_cu_ker_ptrs(struct gkyl_fem_parproj_kernels* kers, enum gkyl_basis_type b_type,
  int dim, int poly_order, bool isperiodic, bool isweighted, bool isdirichlet)
{
  // Set l2g kernels.
  int bckey_periodic = isperiodic? 0 : 1;
  const local2global_kern_list *local2global_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      local2global_kernels = ser_loc2glob_list;
      break;
    default:
      assert(false);
      break;
  }
  for (int k=0; k<2; k++)
    kers->l2g[k] = CK(local2global_kernels, dim, bckey_periodic, poly_order, k);

  // Set RHS stencil kernels.
  int bckey_dirichlet = isdirichlet? 1 : 0;
  const srcstencil_kern_list *srcstencil_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      srcstencil_kernels = isweighted? ser_srcstencil_list_weighted : ser_srcstencil_list_noweight;
      break;
    default:
      assert(false);
  }
  for (int k=0; k<3; k++)
    kers->srcker[k] = CK(srcstencil_kernels, dim, bckey_dirichlet, poly_order, k);

  // Set the get solution stencil kernel.
  const solstencil_kern_list *solstencil_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      solstencil_kernels = ser_solstencil_list;
      break;
    default:
      assert(false);
  }
  kers->solker = solstencil_kernels[dim-1].kernels[poly_order-1];

}

void
fem_parproj_choose_kernels_cu(const struct gkyl_basis *basis, bool has_weight_rhs,
  bool isperiodic, bool isdirichlet, struct gkyl_fem_parproj_kernels *kers)
{
  fem_parproj_set_cu_ker_ptrs<<<1,1>>>(kers, basis->b_type, basis->ndim,
    basis->poly_order, isperiodic, has_weight_rhs, isdirichlet);
}

__global__ void
gkyl_fem_parproj_set_rhs_kernel(double *rhs_global, const struct gkyl_array *rhsin,
  const struct gkyl_array *weight, const struct gkyl_array *phibc,
  struct gkyl_range range, struct gkyl_range perp_range2d, struct gkyl_range par_range1d,
  struct gkyl_fem_parproj_kernels *kers, long numnodes_global)
{
  int idx[GKYL_MAX_CDIM];
  long globalidx[32];
  int parnum_cells = range.upper[range.ndim-1]-range.lower[range.ndim-1]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linidx = gkyl_range_idx(&range, idx);

    const double *wgt_p = weight? (const double *) gkyl_array_cfetch(weight, linidx) : NULL;
    const double *phibc_p = phibc? (const double *) gkyl_array_cfetch(phibc, linidx) : NULL;
    const double *rhsin_p = (const double*) gkyl_array_cfetch(rhsin, linidx);

    int idx1d[] = {idx[range.ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    int keri = idx1d[0] == parnum_cells? 1 : 0;
    kers->l2g[keri](parnum_cells, paridx, globalidx);

    int idx2d[] = {perp_range2d.lower[0], perp_range2d.lower[0]};
    for (int d=0; d<range.ndim-1; d++) idx2d[d] = idx[d];
    long perpidx2d = gkyl_range_idx(&perp_range2d, idx2d);
    long perpProbOff = perpidx2d*numnodes_global;

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(parnum_cells, idx1d[0]);
    kers->srcker[keri](wgt_p, rhsin_p, phibc_p, perpProbOff, globalidx, rhs_global);
  }
}

void
gkyl_fem_parproj_set_rhs_cu(gkyl_fem_parproj *up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{
  gkyl_culinsolver_clear_rhs(up->prob_cu, 0);
  double *rhs_cu = gkyl_culinsolver_get_rhs_ptr(up->prob_cu, 0);

  const struct gkyl_array *phibc_cu = phibc? phibc->on_dev : NULL;
  const struct gkyl_array *wgt_cu = up->has_weight_rhs? up->weight_rhs->on_dev : NULL;

  gkyl_fem_parproj_set_rhs_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(rhs_cu, rhsin->on_dev, wgt_cu, phibc_cu,
    *up->solve_range, up->perp_range2d, up->par_range1d, up->kernels_cu, up->numnodes_global);
}

__global__ void
gkyl_fem_parproj_get_sol_kernel(struct gkyl_array *phiout, const double *x_global, struct gkyl_range range,
  struct gkyl_range perp_range2d, struct gkyl_range par_range1d, struct gkyl_fem_parproj_kernels *kers,
  long numnodes_global)
{
  int idx[GKYL_MAX_DIM];
  long globalidx[32];
  int parnum_cells = range.upper[range.ndim-1]-range.lower[range.ndim-1]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linidx = gkyl_range_idx(&range, idx);
    double *phiout_p = (double*) gkyl_array_cfetch(phiout, linidx);

    int idx1d[] = {idx[range.ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    int keri = idx1d[0] == parnum_cells? 1 : 0;
    kers->l2g[keri](parnum_cells, paridx, globalidx);

    int idx2d[] = {perp_range2d.lower[0], perp_range2d.lower[0]};
    for (int d=0; d<range.ndim-1; d++) idx2d[d] = idx[d];
    long perpidx2d = gkyl_range_idx(&perp_range2d, idx2d);
    long perpProbOff = perpidx2d*numnodes_global;

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    kers->solker(x_global, perpProbOff, globalidx, phiout_p);
  }
}

void
gkyl_fem_parproj_solve_cu(gkyl_fem_parproj *up, struct gkyl_array *phiout)
{
  gkyl_culinsolver_solve(up->prob_cu);
  double *x_cu = gkyl_culinsolver_get_sol_ptr(up->prob_cu, 0);

  gkyl_fem_parproj_get_sol_kernel<<<phiout->nblocks, phiout->nthreads>>>(phiout->on_dev, x_cu,
    *up->solve_range, up->perp_range2d, up->par_range1d, up->kernels_cu, up->numnodes_global);
}