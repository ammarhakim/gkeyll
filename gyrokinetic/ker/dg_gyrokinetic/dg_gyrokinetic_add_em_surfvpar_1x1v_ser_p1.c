#include <gkyl_dg_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double dg_gyrokinetic_add_em_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
    const double *flux_surf_l, const double *flux_surf_r, 
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // flux_surf_l: Surface expansion of phase space flux on the left.
  // flux_surf_r: Surface expansion of phase space flux on the right.
  // out: output increment in center cell.

  return 0.0; 

} 
