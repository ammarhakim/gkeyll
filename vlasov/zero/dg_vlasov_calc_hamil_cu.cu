/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_dg_vlasov_calc_hamil_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_vlasov_calc_hamil_cu_kernel(struct gkyl_rect_grid vel_grid, 
  struct gkyl_basis vel_basis, struct gkyl_range vel_range, 
  enum gkyl_model_id model_id, struct gkyl_array *hamil)
{
  int vdim = vel_basis.ndim;
  int poly_order = vel_basis.poly_order;
  calc_hamil_t calc_hamil;
  switch (vel_basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      calc_hamil = choose_ser_hamil_kern(model_id, vdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      calc_hamil = choose_tensor_hamil_kern(model_id, vdim, poly_order);
      break;

    default:
      assert(false);
      break;    
  }

  int idx[GKYL_MAX_DIM];
  // Cell center array
  double xc[GKYL_MAX_DIM];  

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < vel_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&vel_range, linc1, idx);
    gkyl_rect_grid_cell_center(&vel_grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&vel_range, idx);

    double *hamil_d = (double*) gkyl_array_fetch(hamil, loc);
    calc_hamil(xc, vel_grid.dx, hamil_d);
  }
}

// Host-side wrapper for initialization of Hamiltonian. 
void
gkyl_dg_vlasov_calc_hamil_cu(const struct gkyl_rect_grid *vel_grid, 
  const struct gkyl_basis *vel_basis, const struct gkyl_range *vel_range, 
  enum gkyl_model_id model_id, struct gkyl_array *hamil)
{
  int nblocks = vel_range->nblocks;
  int nthreads = vel_range->nthreads;
  gkyl_dg_vlasov_calc_hamil_cu_kernel<<<nblocks, nthreads>>>(*vel_grid, *vel_basis, *vel_range, 
    model_id, hamil->on_dev);
}