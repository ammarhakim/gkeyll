/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_divide_Jv.h>
#include <gkyl_dg_vlasov_divide_Jv_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_vlasov_divide_Jv_cu_kernel(struct gkyl_basis conf_basis, struct gkyl_basis phase_basis, 
  struct gkyl_range vel_range, struct gkyl_range phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *Jf, 
  struct gkyl_array *f_no_J)
{
  int cdim = conf_basis.ndim; 
  int pdim = phase_basis.ndim;
  int vdim = pdim - cdim; 
  int poly_order = phase_basis.poly_order;
  divide_Jv_t divide_Jv;
  switch (conf_basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      divide_Jv = choose_ser_divide_Jv_kern(cdim, vdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      divide_Jv = choose_tensor_divide_Jv_kern(cdim, vdim, poly_order);
      break;

    default:
      assert(false);
      break;    
  }

  int idx[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);

    for (int i=0; i<vdim; ++i) {
      idx_vel[i] = idx[cdim+i];
    }
    long vidx = gkyl_range_idx(vel_range, idx_vel);
    long pidx = gkyl_range_idx(phase_range, idx);

    const double *Jf_d = (const double*) gkyl_array_fetch(Jf, pidx);
    double *f_no_J_d = (double*) gkyl_array_fetch(f_no_J, pidx);
    divide_Jv(jacob_vel_gauss ? (const double*) gkyl_array_cfetch(jacob_vel_gauss, vidx) : 0,
      Jf_d, f_no_J_d);
  }
}

// Host-side wrapper for dividing out velocity-space Jacobian from Jf. 
void
gkyl_dg_vlasov_divide_Jv_cu(const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *Jf, 
  struct gkyl_array *f_no_J)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_vlasov_divide_Jv_cu_kernel<<<nblocks, nthreads>>>(*conf_basis, *phase_basis, 
    *vel_range, *phase_range, jacob_vel_gauss ? jacob_vel_gauss->on_dev : 0, 
    Jf->on_dev, f_no_J->on_dev);
}