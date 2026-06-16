#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gk_numerical_diffusion.h>
#include <gkyl_gk_numerical_diffusion_priv.h>
#include <gkyl_util.h>

void
gkyl_gk_numerical_diffusion_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct gk_numerical_diffusion *diffusion = container_of(base->on_dev, struct gk_numerical_diffusion, eqn);
    gkyl_cu_free(diffusion);
  }
  
  struct gk_numerical_diffusion *diffusion = container_of(base, struct gk_numerical_diffusion, eqn);
  gkyl_free(diffusion);
}

void
gkyl_gk_numerical_diffusion_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_gk_numerical_diffusion_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.nu) && gkyl_array_is_cu_dev(auxin.jacobgeo_inv)) {
    gkyl_gk_numerical_diffusion_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct gk_numerical_diffusion *diffusion = container_of(eqn, struct gk_numerical_diffusion, eqn);
  diffusion->auxfields.nu = auxin.nu;
  diffusion->auxfields.jacobgeo_inv = auxin.jacobgeo_inv;
}

struct gkyl_dg_eqn*
gkyl_gk_numerical_diffusion_new(const struct gkyl_basis *basis,
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  const enum gkyl_gyrokinetic_bc_type *bc_lower, const enum gkyl_gyrokinetic_bc_type *bc_upper,
  const bool *has_diff_in_dir, int diff_order, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_gk_numerical_diffusion_cu_dev_new(basis, cbasis,
      conf_range, bc_lower, bc_upper, has_diff_in_dir, diff_order);
#endif
  
  struct gk_numerical_diffusion *diffusion = gkyl_malloc(sizeof(struct gk_numerical_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  for (int d=0; d<cdim; d++) diffusion->has_diff_in_dir[d] = 0;
  for (int d=0; d<cdim; d++) diffusion->has_diff_in_dir[d] = has_diff_in_dir[d];

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

  switch (cbasis->b_type) {
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

  int dirs_linidx = diffdirs_linidx(has_diff_in_dir, cdim);

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

  diffusion->auxfields.nu = 0;
  diffusion->auxfields.jacobgeo_inv = 0;
  diffusion->conf_range = *conf_range;

  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_gk_numerical_diffusion_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  
  return &diffusion->eqn;
}
