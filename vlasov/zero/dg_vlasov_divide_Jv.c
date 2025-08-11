#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_divide_Jv.h>
#include <gkyl_dg_vlasov_divide_Jv_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_vlasov_divide_Jv(const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *Jf, 
  struct gkyl_array *f_no_J, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    gkyl_dg_vlasov_divide_Jv_cu(conf_basis, phase_basis, vel_range, phase_range, 
      jacob_vel_gauss, Jf, f_no_J);
    return;
  } 
#endif 

  int cdim = conf_basis->ndim; 
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim; 
  int poly_order = phase_basis->poly_order;
  divide_Jv_t divide_Jv;
  switch (conf_basis->b_type) {
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

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  int idx_vel[GKYL_MAX_DIM];

  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<vdim; ++i) {
      idx_vel[i] = iter.idx[cdim+i];
    }
    long vidx = gkyl_range_idx(vel_range, idx_vel);
    long pidx = gkyl_range_idx(phase_range, iter.idx);

    const double *Jf_d = gkyl_array_cfetch(Jf, pidx);
    double *f_no_J_d = gkyl_array_fetch(f_no_J, pidx);
    divide_Jv(jacob_vel_gauss ? gkyl_array_cfetch(jacob_vel_gauss, vidx) : 0,
      Jf_d, f_no_J_d);
  }  
}

void
gkyl_dg_vlasov_rescale_Jv(const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *f_no_J, 
  struct gkyl_array *Jf, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    gkyl_dg_vlasov_rescale_Jv_cu(conf_basis, phase_basis, vel_range, phase_range, 
      jacob_vel_gauss, f_no_J, Jf);
    return;
  } 
#endif 

  int cdim = conf_basis->ndim; 
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim; 
  int poly_order = phase_basis->poly_order;
  rescale_Jv_t rescale_Jv;
  switch (conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      rescale_Jv = choose_ser_rescale_Jv_kern(cdim, vdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      rescale_Jv = choose_tensor_rescale_Jv_kern(cdim, vdim, poly_order);
      break;

    default:
      assert(false);
      break;    
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  int idx_vel[GKYL_MAX_DIM];

  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<vdim; ++i) {
      idx_vel[i] = iter.idx[cdim+i];
    }
    long vidx = gkyl_range_idx(vel_range, idx_vel);
    long pidx = gkyl_range_idx(phase_range, iter.idx);

    const double *f_no_J_d = gkyl_array_cfetch(f_no_J, pidx);   
    double *Jf_d = gkyl_array_fetch(Jf, pidx);
    rescale_Jv(jacob_vel_gauss ? gkyl_array_cfetch(jacob_vel_gauss, vidx) : 0,
      f_no_J_d, Jf_d);
  }  
}