#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
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

  double rdvpar2 = 2.0/dxv[2];

  double GhatL[8]= {0.0}; 
  double GhatR[8]= {0.0}; 

  const double *fnodal_l = &flux_surf_skin[24]; 
  const double *fnodal_r = &flux_surf_edge[24]; 
  GhatL[0] = 0.3535533905932737*fnodal_l[7]+0.3535533905932737*fnodal_l[6]+0.3535533905932737*fnodal_l[5]+0.3535533905932737*fnodal_l[4]+0.3535533905932737*fnodal_l[3]+0.3535533905932737*fnodal_l[2]+0.3535533905932737*fnodal_l[1]+0.3535533905932737*fnodal_l[0]; 
  GhatL[1] = 0.3535533905932737*fnodal_l[7]+0.3535533905932737*fnodal_l[6]+0.3535533905932737*fnodal_l[5]+0.3535533905932737*fnodal_l[4]-0.3535533905932737*fnodal_l[3]-0.3535533905932737*fnodal_l[2]-0.3535533905932737*fnodal_l[1]-0.3535533905932737*fnodal_l[0]; 
  GhatL[2] = 0.3535533905932737*fnodal_l[7]+0.3535533905932737*fnodal_l[6]-0.3535533905932737*fnodal_l[5]-0.3535533905932737*fnodal_l[4]+0.3535533905932737*fnodal_l[3]+0.3535533905932737*fnodal_l[2]-0.3535533905932737*fnodal_l[1]-0.3535533905932737*fnodal_l[0]; 
  GhatL[3] = 0.3535533905932737*fnodal_l[7]-0.3535533905932737*fnodal_l[6]+0.3535533905932737*fnodal_l[5]-0.3535533905932737*fnodal_l[4]+0.3535533905932737*fnodal_l[3]-0.3535533905932737*fnodal_l[2]+0.3535533905932737*fnodal_l[1]-0.3535533905932737*fnodal_l[0]; 
  GhatL[4] = 0.3535533905932737*fnodal_l[7]+0.3535533905932737*fnodal_l[6]-0.3535533905932737*fnodal_l[5]-0.3535533905932737*fnodal_l[4]-0.3535533905932737*fnodal_l[3]-0.3535533905932737*fnodal_l[2]+0.3535533905932737*fnodal_l[1]+0.3535533905932737*fnodal_l[0]; 
  GhatL[5] = 0.3535533905932737*fnodal_l[7]-0.3535533905932737*fnodal_l[6]+0.3535533905932737*fnodal_l[5]-0.3535533905932737*fnodal_l[4]-0.3535533905932737*fnodal_l[3]+0.3535533905932737*fnodal_l[2]-0.3535533905932737*fnodal_l[1]+0.3535533905932737*fnodal_l[0]; 
  GhatL[6] = 0.3535533905932737*fnodal_l[7]-0.3535533905932737*fnodal_l[6]-0.3535533905932737*fnodal_l[5]+0.3535533905932737*fnodal_l[4]+0.3535533905932737*fnodal_l[3]-0.3535533905932737*fnodal_l[2]-0.3535533905932737*fnodal_l[1]+0.3535533905932737*fnodal_l[0]; 
  GhatL[7] = 0.3535533905932737*fnodal_l[7]-0.3535533905932737*fnodal_l[6]-0.3535533905932737*fnodal_l[5]+0.3535533905932737*fnodal_l[4]-0.3535533905932737*fnodal_l[3]+0.3535533905932737*fnodal_l[2]+0.3535533905932737*fnodal_l[1]-0.3535533905932737*fnodal_l[0]; 
  GhatR[0] = 0.3535533905932737*fnodal_r[7]+0.3535533905932737*fnodal_r[6]+0.3535533905932737*fnodal_r[5]+0.3535533905932737*fnodal_r[4]+0.3535533905932737*fnodal_r[3]+0.3535533905932737*fnodal_r[2]+0.3535533905932737*fnodal_r[1]+0.3535533905932737*fnodal_r[0]; 
  GhatR[1] = 0.3535533905932737*fnodal_r[7]+0.3535533905932737*fnodal_r[6]+0.3535533905932737*fnodal_r[5]+0.3535533905932737*fnodal_r[4]-0.3535533905932737*fnodal_r[3]-0.3535533905932737*fnodal_r[2]-0.3535533905932737*fnodal_r[1]-0.3535533905932737*fnodal_r[0]; 
  GhatR[2] = 0.3535533905932737*fnodal_r[7]+0.3535533905932737*fnodal_r[6]-0.3535533905932737*fnodal_r[5]-0.3535533905932737*fnodal_r[4]+0.3535533905932737*fnodal_r[3]+0.3535533905932737*fnodal_r[2]-0.3535533905932737*fnodal_r[1]-0.3535533905932737*fnodal_r[0]; 
  GhatR[3] = 0.3535533905932737*fnodal_r[7]-0.3535533905932737*fnodal_r[6]+0.3535533905932737*fnodal_r[5]-0.3535533905932737*fnodal_r[4]+0.3535533905932737*fnodal_r[3]-0.3535533905932737*fnodal_r[2]+0.3535533905932737*fnodal_r[1]-0.3535533905932737*fnodal_r[0]; 
  GhatR[4] = 0.3535533905932737*fnodal_r[7]+0.3535533905932737*fnodal_r[6]-0.3535533905932737*fnodal_r[5]-0.3535533905932737*fnodal_r[4]-0.3535533905932737*fnodal_r[3]-0.3535533905932737*fnodal_r[2]+0.3535533905932737*fnodal_r[1]+0.3535533905932737*fnodal_r[0]; 
  GhatR[5] = 0.3535533905932737*fnodal_r[7]-0.3535533905932737*fnodal_r[6]+0.3535533905932737*fnodal_r[5]-0.3535533905932737*fnodal_r[4]-0.3535533905932737*fnodal_r[3]+0.3535533905932737*fnodal_r[2]-0.3535533905932737*fnodal_r[1]+0.3535533905932737*fnodal_r[0]; 
  GhatR[6] = 0.3535533905932737*fnodal_r[7]-0.3535533905932737*fnodal_r[6]-0.3535533905932737*fnodal_r[5]+0.3535533905932737*fnodal_r[4]+0.3535533905932737*fnodal_r[3]-0.3535533905932737*fnodal_r[2]-0.3535533905932737*fnodal_r[1]+0.3535533905932737*fnodal_r[0]; 
  GhatR[7] = 0.3535533905932737*fnodal_r[7]-0.3535533905932737*fnodal_r[6]-0.3535533905932737*fnodal_r[5]+0.3535533905932737*fnodal_r[4]-0.3535533905932737*fnodal_r[3]+0.3535533905932737*fnodal_r[2]+0.3535533905932737*fnodal_r[1]-0.3535533905932737*fnodal_r[0]; 

  if (edge == -1) { 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[3] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdvpar2; 
  out[6] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[7] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[8] += -0.7071067811865475*GhatR[5]*rdvpar2; 
  out[9] += -0.7071067811865475*GhatR[6]*rdvpar2; 
  out[10] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[11] += -1.224744871391589*GhatR[4]*rdvpar2; 
  out[12] += -0.7071067811865475*GhatR[7]*rdvpar2; 
  out[13] += -1.224744871391589*GhatR[5]*rdvpar2; 
  out[14] += -1.224744871391589*GhatR[6]*rdvpar2; 
  out[15] += -1.224744871391589*GhatR[7]*rdvpar2; 
  out[16] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[17] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[18] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[19] += -1.58113883008419*GhatR[3]*rdvpar2; 
  out[20] += -1.58113883008419*GhatR[4]*rdvpar2; 
  out[21] += -1.58113883008419*GhatR[5]*rdvpar2; 
  out[22] += -1.58113883008419*GhatR[6]*rdvpar2; 
  out[23] += -1.58113883008419*GhatR[7]*rdvpar2; 

  } else { 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[3] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdvpar2; 
  out[6] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[7] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[8] += 0.7071067811865475*GhatL[5]*rdvpar2; 
  out[9] += 0.7071067811865475*GhatL[6]*rdvpar2; 
  out[10] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[11] += -1.224744871391589*GhatL[4]*rdvpar2; 
  out[12] += 0.7071067811865475*GhatL[7]*rdvpar2; 
  out[13] += -1.224744871391589*GhatL[5]*rdvpar2; 
  out[14] += -1.224744871391589*GhatL[6]*rdvpar2; 
  out[15] += -1.224744871391589*GhatL[7]*rdvpar2; 
  out[16] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[17] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[18] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[19] += 1.58113883008419*GhatL[3]*rdvpar2; 
  out[20] += 1.58113883008419*GhatL[4]*rdvpar2; 
  out[21] += 1.58113883008419*GhatL[5]*rdvpar2; 
  out[22] += 1.58113883008419*GhatL[6]*rdvpar2; 
  out[23] += 1.58113883008419*GhatL[7]*rdvpar2; 

  } 

  return 0.0; 

} 
