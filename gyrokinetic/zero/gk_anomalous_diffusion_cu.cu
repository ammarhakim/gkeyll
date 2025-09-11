/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gk_anomalous_diffusion.h>    
#include <gkyl_gk_anomalous_diffusion_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_gk_anomalous_diffusion_set_auxfields_cu_kernel(const struct gkyl_dg_eqn* eqn,
  const struct gkyl_array* nu, const struct gkyl_array* jacobgeo_inv)
{
  struct gk_anomalous_diffusion* diffusion = container_of(eqn, struct gk_anomalous_diffusion, eqn);
  diffusion->auxfields.nu = nu;
  diffusion->auxfields.jacobgeo_inv = jacobgeo_inv;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_gk_anomalous_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn* eqn,
  struct gkyl_gk_anomalous_diffusion_auxfields auxin)
{
  gkyl_gk_anomalous_diffusion_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.nu->on_dev, auxin.jacobgeo_inv->on_dev);
}

__global__ void static
gk_anomalous_diffusion_set_cu_dev_ptrs(struct gk_anomalous_diffusion *diffusion, enum gkyl_basis_type b_type,
  int cdim, int vdim, int poly_order, bool use_bound_local)
{
  diffusion->auxfields.nu = 0; 
  diffusion->auxfields.jacobgeo_inv = 0; 

  const gkyl_gk_anomalous_diffusion_vol_kern_list *vol_kernels;
  const gkyl_gk_anomalous_diffusion_surf_kern_list *surfx_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_surfx_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_diagx_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels            = ser_vol_kernels;
      surfx_kernels          = ser_gyrokinetic_surfx_kernels;
      boundary_surfx_kernels = use_bound_local? ser_gyrokinetic_boundary_surfx_boundlocal_kernels
                                              : ser_gyrokinetic_boundary_surfx_zeroflux_kernels;
      boundary_diagx_kernels = use_bound_local? ser_gyrokinetic_boundary_diagx_boundlocal_kernels
                                              : ser_gyrokinetic_boundary_diagx_boundrecovery_kernels;
      break;

    default:
      assert(false);
      break;
  }

  diffusion->eqn.num_equations = 1;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;
  diffusion->eqn.boundary_diag_term = boundary_diag;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, cdim+vdim, poly_order);
  diffusion->surf = CKSURF(surfx_kernels, cdim+vdim, poly_order);
  diffusion->boundary_surf = CKSURF(boundary_surfx_kernels, cdim+vdim, poly_order);
  diffusion->boundary_diag = CKSURF(boundary_diagx_kernels, cdim+vdim, poly_order);
}

struct gkyl_dg_eqn*
gkyl_gk_anomalous_diffusion_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  const struct gkyl_range *diff_range, bool use_bound_local, double skip_cell_threshold)
{
  struct gk_anomalous_diffusion* diffusion = (struct gk_anomalous_diffusion*) gkyl_malloc(sizeof(struct gk_anomalous_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  if (skip_cell_threshold > 0.0)
    diffusion->skip_cell_thresh = skip_cell_threshold * pow(sqrt(2.0), cdim + vdim);
  else
    diffusion->skip_cell_thresh = -1.0;

  diffusion->diff_range = diff_range;

  diffusion->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion->eqn.flags);
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_gk_anomalous_diffusion_free);

  // Copy the host struct to device struct.
  struct gk_anomalous_diffusion* diffusion_cu = (struct gk_anomalous_diffusion*) gkyl_cu_malloc(sizeof(struct gk_anomalous_diffusion));
  gkyl_cu_memcpy(diffusion_cu, diffusion, sizeof(struct gk_anomalous_diffusion), GKYL_CU_MEMCPY_H2D);
  gk_anomalous_diffusion_set_cu_dev_ptrs<<<1,1>>>(diffusion_cu, cbasis->b_type, cdim, vdim, poly_order, use_bound_local);

  // Set parent on_dev pointer.
  diffusion->eqn.on_dev = &diffusion_cu->eqn;

  return &diffusion->eqn;
}
