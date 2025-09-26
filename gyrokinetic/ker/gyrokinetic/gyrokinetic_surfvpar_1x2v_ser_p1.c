#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
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

  double rdvpar2 = 2.0/dxv[1];

  double GhatL[4]= {0.0}; 
  double GhatR[4]= {0.0}; 

  const double *fnodal_l = &flux_surf_l[6]; 
  const double *fnodal_r = &flux_surf_r[6]; 
  GhatL[0] = 0.5*fnodal_l[3]+0.5*fnodal_l[2]+0.5*fnodal_l[1]+0.5*fnodal_l[0]; 
  GhatL[1] = 0.5*fnodal_l[3]+0.5*fnodal_l[2]-0.5*fnodal_l[1]-0.5*fnodal_l[0]; 
  GhatL[2] = 0.5*fnodal_l[3]-0.5*fnodal_l[2]+0.5*fnodal_l[1]-0.5*fnodal_l[0]; 
  GhatL[3] = 0.5*fnodal_l[3]-0.5*fnodal_l[2]-0.5*fnodal_l[1]+0.5*fnodal_l[0]; 
  GhatR[0] = 0.5*fnodal_r[3]+0.5*fnodal_r[2]+0.5*fnodal_r[1]+0.5*fnodal_r[0]; 
  GhatR[1] = 0.5*fnodal_r[3]+0.5*fnodal_r[2]-0.5*fnodal_r[1]-0.5*fnodal_r[0]; 
  GhatR[2] = 0.5*fnodal_r[3]-0.5*fnodal_r[2]+0.5*fnodal_r[1]-0.5*fnodal_r[0]; 
  GhatR[3] = 0.5*fnodal_r[3]-0.5*fnodal_r[2]-0.5*fnodal_r[1]+0.5*fnodal_r[0]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdvpar2; 
  out[4] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[5] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdvpar2; 
  out[6] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdvpar2; 
  out[7] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdvpar2; 
  out[8] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdvpar2; 
  out[9] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdvpar2; 
  out[10] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdvpar2; 
  out[11] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdvpar2; 

  return 0.0; 

} 
