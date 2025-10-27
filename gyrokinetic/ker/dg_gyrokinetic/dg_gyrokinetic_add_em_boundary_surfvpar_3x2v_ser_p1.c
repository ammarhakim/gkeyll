#include <gkyl_dg_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double dg_gyrokinetic_add_em_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_edge, const double *vmap_prime_skin,
    const double *flux_surf_edge, const double *flux_surf_skin, 
    const int edge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // flux_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // flux_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // out: output increment in center cell.

  return 0.0; 

} 
