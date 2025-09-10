#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_edge, const double *vmap_prime_skin,
    const double *alpha_surf_edge, const double *alpha_surf_skin, 
    const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
    const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // alpha_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // alpha_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // sgn_alpha_surf_edge: sign(alpha_surf_edge) at quadrature points.
  // sgn_alpha_surf_skin: sign(alpha_surf_skin) at quadrature points.
  // const_sgn_alpha_edge: Boolean array true if sign(alpha_surf_edge) is only one sign, either +1 or -1.
  // const_sgn_alpha_skin: Boolean array true if sign(alpha_surf_skin) is only one sign, either +1 or -1.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double rdvpar2 = 2.0/dxv[1];

  const double *alphaL = &alpha_surf_skin[6];
  const double *alphaR = &alpha_surf_edge[6];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[6];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[6];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[1];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[1];

  if (edge == -1) { 

  double fUpR[4] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = (1.58113883008419*fskin[8]+1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  fUpR[1] = (1.58113883008419*fskin[9]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  fUpR[2] = (1.58113883008419*fskin[10]+1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  fUpR[3] = (1.58113883008419*fskin[11]+1.224744871391589*fskin[7]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 
    } else { 
  fUpR[0] = (1.58113883008419*fedge[8]-1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  fUpR[1] = (1.58113883008419*fedge[9]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  fUpR[2] = (1.58113883008419*fedge[10]-1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  fUpR[3] = (1.58113883008419*fedge[11]-1.224744871391589*fedge[7]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 
    } 
  } else { 
  double f_cr[4] = {0.};
  double f_rl[4] = {0.};
  double sgn_alphaUpR[4] = {0.};
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = (1.58113883008419*fskin[8]+1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  f_cr[1] = (1.58113883008419*fskin[9]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  f_cr[2] = (1.58113883008419*fskin[10]+1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  f_cr[3] = (1.58113883008419*fskin[11]+1.224744871391589*fskin[7]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 

  f_rl[0] = (1.58113883008419*fedge[8]-1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  f_rl[1] = (1.58113883008419*fedge[9]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  f_rl[2] = (1.58113883008419*fedge[10]-1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  f_rl[3] = (1.58113883008419*fedge[11]-1.224744871391589*fedge[7]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 

  fUpR[0] = (0.25*f_cr[3]-0.25*f_rl[3])*sgn_alphaUpR[3]+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[2]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[1]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[1]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[2]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = (0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[3]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[3]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[3]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.25*f_cr[2]-0.25*f_rl[2]); 

  } 
  double GhatR[4] = {0.};
  GhatR[0] = 0.5*(alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.5*(alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.5*(alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.5*(alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[4] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[6] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[7] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[8] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[9] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[10] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[11] += -1.58113883008419*GhatR[3]*rdvpar2; 

  } else { 

  double fUpL[4] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = (1.58113883008419*fedge[8]+1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  fUpL[1] = (1.58113883008419*fedge[9]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  fUpL[2] = (1.58113883008419*fedge[10]+1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  fUpL[3] = (1.58113883008419*fedge[11]+1.224744871391589*fedge[7]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 
    } else { 
  fUpL[0] = (1.58113883008419*fskin[8]-1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  fUpL[1] = (1.58113883008419*fskin[9]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  fUpL[2] = (1.58113883008419*fskin[10]-1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  fUpL[3] = (1.58113883008419*fskin[11]-1.224744871391589*fskin[7]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 
    } 
  } else { 
  double f_lr[4] = {0.};
  double f_cl[4] = {0.};
  double sgn_alphaUpL[4] = {0.};
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = (1.58113883008419*fedge[8]+1.224744871391589*fedge[2]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  f_lr[1] = (1.58113883008419*fedge[9]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  f_lr[2] = (1.58113883008419*fedge[10]+1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  f_lr[3] = (1.58113883008419*fedge[11]+1.224744871391589*fedge[7]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 

  f_cl[0] = (1.58113883008419*fskin[8]-1.224744871391589*fskin[2]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  f_cl[1] = (1.58113883008419*fskin[9]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  f_cl[2] = (1.58113883008419*fskin[10]-1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  f_cl[3] = (1.58113883008419*fskin[11]-1.224744871391589*fskin[7]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 

  fUpL[0] = (0.25*f_lr[3]-0.25*f_cl[3])*sgn_alphaUpL[3]+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[2]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[1]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[1]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[2]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = (0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[3]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[3]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.25*f_lr[2]-0.25*f_cl[2]); 

  } 
  double GhatL[4] = {0.};
  GhatL[0] = 0.5*(alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.5*(alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.5*(alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.5*(alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[4] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[6] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[7] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[8] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[9] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[10] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[11] += 1.58113883008419*GhatL[3]*rdvpar2; 

  } 

  double vmap_prime_min = fmin(fabs(vmap_prime_edge[0]),fabs(vmap_prime_skin[0]));
  double cflFreq = fmax(fabs(alphaL[0]/vmap_prime_min), fabs(alphaR[0]/vmap_prime_min)); 
  return 1.25*rdvpar2*cflFreq; 

} 
