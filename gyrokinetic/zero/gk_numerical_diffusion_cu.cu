/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gk_numerical_diffusion.h>    
#include <gkyl_gk_numerical_diffusion_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_gk_numerical_diffusion_set_auxfields_cu_kernel(const struct gkyl_dg_eqn* eqn,
  const struct gkyl_array* nu, const struct gkyl_array* jacobgeo_inv)
{
  struct gk_numerical_diffusion* diffusion = container_of(eqn, struct gk_numerical_diffusion, eqn);
  diffusion->auxfields.nu = nu;
  diffusion->auxfields.jacobgeo_inv = jacobgeo_inv;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_gk_numerical_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn* eqn, struct gkyl_gk_numerical_diffusion_auxfields auxin)
{
  gkyl_gk_numerical_diffusion_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.nu->on_dev, auxin.jacobgeo_inv->on_dev);
}

__global__ void static
gk_numerical_diffusion_set_cu_dev_ptrs(struct gk_numerical_diffusion *diffusion, enum gkyl_basis_type b_type,
  int cdim, int vdim, int poly_order, const enum gkyl_gyrokinetic_bc_type *bc_lower, const enum gkyl_gyrokinetic_bc_type *bc_upper,
  const bool *has_diff_in_dir, int diff_order, int dirs_linidx)
{
  diffusion->auxfields.nu = 0; 
  diffusion->auxfields.jacobgeo_inv = 0; 

  const gkyl_gk_numerical_diffusion_vol_kern_list *vol_kernels;
  const gkyl_gk_numerical_diffusion_surf_kern_list *surf_kernels;
  const gkyl_gk_numerical_diffusion_boundary_surf_kern_list *boundary_surf_lower_kernels[GKYL_MAX_CDIM],
                                                            *boundary_surf_upper_kernels[GKYL_MAX_CDIM];
  const gkyl_gk_numerical_diffusion_boundary_surf_kern_list *boundary_diag_lower_kernels[GKYL_MAX_CDIM],
                                                            *boundary_diag_upper_kernels[GKYL_MAX_CDIM];

  // Choice of boundary_surf and boundary_diag kernels:
  //   boundary_surf: zero_flux or local
  //   boundary_diag: local or recovery
  // MF 2025/09/10: as of now these options are meant for (here
  // N/A means not applicable):
  //             bound_surf  bound_diag  hyper_dg-zero_flux
  // SKIP:       N/A         recovery    no
  // ABSORB:     N/A         N/A         no
  // PERIODIC:   N/A         N/A         no
  // FIXED_FUNC: N/A         N/A         no
  // ZERO_FLUX:  zero_flux   N/A         yes
  // ELSE:       local       local       yes

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels_varnu;

      surf_kernels  = ser_gyrokinetic_surf_kernels_varnu;

      for (int d=0; d<cdim; d++) {

        if ((bc_lower[d] == GKYL_BC_GK_SKIP) ||
            (bc_lower[d] == GKYL_BC_GK_SPECIES_PERIODIC)) {
          // Boundary surf and diag kernels not used.
          boundary_surf_lower_kernels[d] = &ser_gyrokinetic_boundary_surf_lower_zeroflux_kernels_varnu[d];
          boundary_diag_lower_kernels[d] = &ser_gyrokinetic_boundary_diag_lower_boundrecovery_kernels_varnu[d];
        }
        else if (bc_lower[d] == GKYL_BC_GK_SPECIES_ZERO_FLUX) {
          boundary_surf_lower_kernels[d] = &ser_gyrokinetic_boundary_surf_lower_zeroflux_kernels_varnu[d];
          // Boundary diag kernel not used.
          boundary_diag_lower_kernels[d] = &ser_gyrokinetic_boundary_diag_lower_boundrecovery_kernels_varnu[d];
        }
        else if ((bc_lower[d] == GKYL_BC_GK_SPECIES_ABSORB) ||
                 (bc_lower[d] == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
          // Boundary surf kernel not used.
          boundary_surf_lower_kernels[d] = &ser_gyrokinetic_boundary_surf_lower_zeroflux_kernels_varnu[d];
          boundary_diag_lower_kernels[d] = &ser_gyrokinetic_boundary_diag_lower_boundrecovery_kernels_varnu[d];
        }
        else {
          boundary_surf_lower_kernels[d] = &ser_gyrokinetic_boundary_surf_lower_boundlocal_kernels_varnu[d];
          boundary_diag_lower_kernels[d] = &ser_gyrokinetic_boundary_diag_lower_boundlocal_kernels_varnu[d];
        }

        if ((bc_upper[d] == GKYL_BC_GK_SKIP) ||
            (bc_upper[d] == GKYL_BC_GK_SPECIES_PERIODIC)) {
          // Boundary surf and diag kernels not used.
          boundary_surf_upper_kernels[d] = &ser_gyrokinetic_boundary_surf_upper_zeroflux_kernels_varnu[d];
          boundary_diag_upper_kernels[d] = &ser_gyrokinetic_boundary_diag_upper_boundrecovery_kernels_varnu[d];
        }
        else if (bc_upper[d] == GKYL_BC_GK_SPECIES_ZERO_FLUX) {
          boundary_surf_upper_kernels[d] = &ser_gyrokinetic_boundary_surf_upper_zeroflux_kernels_varnu[d];
          // Boundary diag kernel not used.
          boundary_diag_upper_kernels[d] = &ser_gyrokinetic_boundary_diag_upper_boundrecovery_kernels_varnu[d];
        }
        else if ((bc_upper[d] == GKYL_BC_GK_SPECIES_ABSORB) ||
                 (bc_upper[d] == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
          // Boundary surf kernel not used.
          boundary_surf_upper_kernels[d] = &ser_gyrokinetic_boundary_surf_upper_zeroflux_kernels_varnu[d];
          boundary_diag_upper_kernels[d] = &ser_gyrokinetic_boundary_diag_upper_boundrecovery_kernels_varnu[d];
        }
        else {
          boundary_surf_upper_kernels[d] = &ser_gyrokinetic_boundary_surf_upper_boundlocal_kernels_varnu[d];
          boundary_diag_upper_kernels[d] = &ser_gyrokinetic_boundary_diag_upper_boundlocal_kernels_varnu[d];
        }
      }

      break;
  
    default:
      assert(false);
      break;    
  } 

  diffusion->eqn.num_equations = 1;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;
  diffusion->eqn.boundary_diag_term = boundary_diag;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, cdim, vdim, diff_order, poly_order, dirs_linidx);

  for (int d=0; d<cdim; d++) {
    diffusion->surf[d] = CKSURF(surf_kernels, cdim, vdim, diff_order, poly_order, d);

    diffusion->boundary_surf_lower[d] = CKBSURF(boundary_surf_lower_kernels[d], cdim, vdim, diff_order, poly_order);
    diffusion->boundary_diag_lower[d] = CKBSURF(boundary_diag_lower_kernels[d], cdim, vdim, diff_order, poly_order);

    diffusion->boundary_surf_upper[d] = CKBSURF(boundary_surf_upper_kernels[d], cdim, vdim, diff_order, poly_order);
    diffusion->boundary_diag_upper[d] = CKBSURF(boundary_diag_upper_kernels[d], cdim, vdim, diff_order, poly_order);
  }

  // Ensure non-NULL pointers.
  assert(diffusion->eqn.vol_term);
  for (int i=0; i<cdim; ++i) {
    assert(diffusion->surf[i]);
    assert(diffusion->boundary_surf_lower[i]);
    assert(diffusion->boundary_diag_lower[i]);
    assert(diffusion->boundary_surf_upper[i]);
    assert(diffusion->boundary_diag_upper[i]);
  }

}

struct gkyl_dg_eqn*
gkyl_gk_numerical_diffusion_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  const struct gkyl_range *conf_range, const enum gkyl_gyrokinetic_bc_type *bc_lower, const enum gkyl_gyrokinetic_bc_type *bc_upper,
  const bool *has_diff_in_dir, int diff_order)
{
  struct gk_numerical_diffusion* diffusion = (struct gk_numerical_diffusion*) gkyl_malloc(sizeof(struct gk_numerical_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  diffusion->conf_range = *conf_range;

  for (int d=0; d<cdim; d++) diffusion->has_diff_in_dir[d] = 0;
  for (int d=0; d<cdim; d++) diffusion->has_diff_in_dir[d] = has_diff_in_dir[d];

  int dirs_linidx = diffdirs_linidx(has_diff_in_dir, cdim);

  diffusion->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion->eqn.flags);
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_gk_numerical_diffusion_free);

  size_t bc_sz = cdim*sizeof(enum gkyl_gyrokinetic_bc_type);
  enum gkyl_gyrokinetic_bc_type *bc_lower_cu = (enum gkyl_gyrokinetic_bc_type*) gkyl_cu_malloc(bc_sz);
  enum gkyl_gyrokinetic_bc_type *bc_upper_cu = (enum gkyl_gyrokinetic_bc_type*) gkyl_cu_malloc(bc_sz);
  gkyl_cu_memcpy(bc_lower_cu, bc_lower, bc_sz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(bc_upper_cu, bc_upper, bc_sz, GKYL_CU_MEMCPY_H2D);

  size_t dir_sz = GKYL_MAX_CDIM*sizeof(bool);
  bool *has_diff_in_dir_cu = (bool *) gkyl_cu_malloc(dir_sz);
  gkyl_cu_memcpy(has_diff_in_dir_cu, has_diff_in_dir, dir_sz, GKYL_CU_MEMCPY_H2D);

  // copy the host struct to device struct
  struct gk_numerical_diffusion* diffusion_cu = (struct gk_numerical_diffusion*) gkyl_cu_malloc(sizeof(struct gk_numerical_diffusion));
  gkyl_cu_memcpy(diffusion_cu, diffusion, sizeof(struct gk_numerical_diffusion), GKYL_CU_MEMCPY_H2D);
  gk_numerical_diffusion_set_cu_dev_ptrs<<<1,1>>>(diffusion_cu, cbasis->b_type, cdim, vdim, poly_order,
    bc_lower_cu, bc_upper_cu, has_diff_in_dir_cu, diff_order, dirs_linidx);

  // set parent on_dev pointer
  diffusion->eqn.on_dev = &diffusion_cu->eqn;

  gkyl_cu_free(has_diff_in_dir_cu);
  gkyl_cu_free(bc_lower_cu);
  gkyl_cu_free(bc_upper_cu);

  return &diffusion->eqn;
}
