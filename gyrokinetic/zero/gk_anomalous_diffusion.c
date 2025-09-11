#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gk_anomalous_diffusion.h>
#include <gkyl_gk_anomalous_diffusion_priv.h>
#include <gkyl_util.h>

void
gkyl_gk_anomalous_diffusion_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct gk_anomalous_diffusion *diffusion = container_of(base->on_dev, struct gk_anomalous_diffusion, eqn);
    gkyl_cu_free(diffusion);
  }
  
  struct gk_anomalous_diffusion *diffusion = container_of(base, struct gk_anomalous_diffusion, eqn);
  gkyl_free(diffusion);
}

void
gkyl_gk_anomalous_diffusion_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_gk_anomalous_diffusion_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.nu) && gkyl_array_is_cu_dev(auxin.jacobgeo_inv)) {
    gkyl_gk_anomalous_diffusion_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct gk_anomalous_diffusion *diffusion = container_of(eqn, struct gk_anomalous_diffusion, eqn);
  diffusion->auxfields.nu = auxin.nu;
  diffusion->auxfields.jacobgeo_inv = auxin.jacobgeo_inv;
}

struct gkyl_dg_eqn*
gkyl_gk_anomalous_diffusion_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  const struct gkyl_range *conf_range, double skip_cell_threshold, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_gk_anomalous_diffusion_cu_dev_new(basis, cbasis, conf_range, skip_cell_threshold);
#endif
  
  struct gk_anomalous_diffusion *diffusion = gkyl_malloc(sizeof(struct gk_anomalous_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  if (skip_cell_threshold > 0.0)
    diffusion->skip_cell_thresh = skip_cell_threshold * pow(sqrt(2.0), cdim + vdim);
  else
    diffusion->skip_cell_thresh = -1.0;

  const gkyl_gk_anomalous_diffusion_vol_kern_list *vol_kernels;
  const gkyl_gk_anomalous_diffusion_surf_kern_list *surfx_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_surfx_kernels;
  const gkyl_gk_anomalous_diffusion_boundary_surf_kern_list *boundary_diagx_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels            = ser_vol_kernels;
      surfx_kernels          = ser_gyrokinetic_surfx_kernels;
      boundary_surfx_kernels = ser_gyrokinetic_boundary_surfx_kernels;
      boundary_diagx_kernels = ser_gyrokinetic_boundary_diagx_kernels;
      break;
  
    default:
      assert(false);
      break;    
  } 

  diffusion->eqn.num_equations = 1;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;
  diffusion->eqn.boundary_diag_term = boundary_diag;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, cdim+vdim,  poly_order);
  diffusion->surf = CKSURF(surfx_kernels, cdim+vdim,  poly_order);
  diffusion->boundary_surf = CKSURF(boundary_surfx_kernels, cdim+vdim,  poly_order);
  diffusion->boundary_diag = CKSURF(boundary_diagx_kernels, cdim+vdim,  poly_order);

  // Ensure non-NULL pointers.
  assert(diffusion->surf);
  assert(diffusion->boundary_surf);
  assert(diffusion->boundary_diag);

  diffusion->auxfields.nu = 0;
  diffusion->auxfields.jacobgeo_inv = 0;
  diffusion->conf_range = conf_range;

  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_gk_anomalous_diffusion_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  
  return &diffusion->eqn;
}
