#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
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

  double rdvpar2 = 2.0/dxv[1];

  double GhatL[2]= {0.0}; 
  double GhatR[2]= {0.0}; 

  const double *fnodal_l = &flux_surf_skin[3]; 
  const double *fnodal_r = &flux_surf_edge[3]; 
  GhatL[0] = 0.7071067811865475*fnodal_l[1]+0.7071067811865475*fnodal_l[0]; 
  GhatL[1] = 0.7071067811865475*fnodal_l[1]-0.7071067811865475*fnodal_l[0]; 
  GhatR[0] = 0.7071067811865475*fnodal_r[1]+0.7071067811865475*fnodal_r[0]; 
  GhatR[1] = 0.7071067811865475*fnodal_r[1]-0.7071067811865475*fnodal_r[0]; 

  if (edge == -1) { 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[3] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[4] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[5] += -1.58113883008419*GhatR[1]*rdvpar2; 

  } else { 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[3] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[4] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[5] += 1.58113883008419*GhatL[1]*rdvpar2; 

  } 

  return 0.0; 

} 
