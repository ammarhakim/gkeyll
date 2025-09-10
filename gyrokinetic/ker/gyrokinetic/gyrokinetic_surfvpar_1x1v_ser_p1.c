#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
    const double *alpha_surf_l, const double *alpha_surf_r, 
    const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r,
    const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // alpha_surf_l: Surface expansion of phase space flux on the left.
  // alpha_surf_r: Surface expansion of phase space flux on the right.
  // sgn_alpha_surf_l: sign(alpha_surf_l) at quadrature points.
  // sgn_alpha_surf_r: sign(alpha_surf_r) at quadrature points.
  // const_sgn_alpha_l: Boolean array true if sign(alpha_surf_l) is only one sign, either +1 or -1.
  // const_sgn_alpha_r: Boolean array true if sign(alpha_surf_r) is only one sign, either +1 or -1.
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  double rdvpar2 = 2.0/dxv[1];

  const double *alphaL = &alpha_surf_l[3];
  const double *alphaR = &alpha_surf_r[3];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_l[3];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_r[3];
  const int *const_sgn_alphaL = &const_sgn_alpha_l[1];
  const int *const_sgn_alphaR = &const_sgn_alpha_r[1];

  double fUpL[2] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = (1.58113883008419*fl[4]+1.224744871391589*fl[2]+0.7071067811865475*fl[0])/vmap_prime_l[0]; 
  fUpL[1] = (1.58113883008419*fl[5]+1.224744871391589*fl[3]+0.7071067811865475*fl[1])/vmap_prime_l[0]; 
    } else { 
  fUpL[0] = (1.58113883008419*fc[4]-1.224744871391589*fc[2]+0.7071067811865475*fc[0])/vmap_prime_c[0]; 
  fUpL[1] = (1.58113883008419*fc[5]-1.224744871391589*fc[3]+0.7071067811865475*fc[1])/vmap_prime_c[0]; 
    } 
  } else { 
  double f_lr[2] = {0.};
  double f_cl[2] = {0.};
  double sgn_alphaUpL[2] = {0.};
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = (1.58113883008419*fl[4]+1.224744871391589*fl[2]+0.7071067811865475*fl[0])/vmap_prime_l[0]; 
  f_lr[1] = (1.58113883008419*fl[5]+1.224744871391589*fl[3]+0.7071067811865475*fl[1])/vmap_prime_l[0]; 

  f_cl[0] = (1.58113883008419*fc[4]-1.224744871391589*fc[2]+0.7071067811865475*fc[0])/vmap_prime_c[0]; 
  f_cl[1] = (1.58113883008419*fc[5]-1.224744871391589*fc[3]+0.7071067811865475*fc[1])/vmap_prime_c[0]; 

  fUpL[0] = (0.3535533905932737*f_lr[1]-0.3535533905932737*f_cl[1])*sgn_alphaUpL[1]+(0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[1]+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[1]; 

  } 
  double fUpR[2] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = (1.58113883008419*fc[4]+1.224744871391589*fc[2]+0.7071067811865475*fc[0])/vmap_prime_c[0]; 
  fUpR[1] = (1.58113883008419*fc[5]+1.224744871391589*fc[3]+0.7071067811865475*fc[1])/vmap_prime_c[0]; 
    } else { 
  fUpR[0] = (1.58113883008419*fr[4]-1.224744871391589*fr[2]+0.7071067811865475*fr[0])/vmap_prime_r[0]; 
  fUpR[1] = (1.58113883008419*fr[5]-1.224744871391589*fr[3]+0.7071067811865475*fr[1])/vmap_prime_r[0]; 
    } 
  } else { 
  double f_cr[2] = {0.};
  double f_rl[2] = {0.};
  double sgn_alphaUpR[2] = {0.};
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = (1.58113883008419*fc[4]+1.224744871391589*fc[2]+0.7071067811865475*fc[0])/vmap_prime_c[0]; 
  f_cr[1] = (1.58113883008419*fc[5]+1.224744871391589*fc[3]+0.7071067811865475*fc[1])/vmap_prime_c[0]; 

  f_rl[0] = (1.58113883008419*fr[4]-1.224744871391589*fr[2]+0.7071067811865475*fr[0])/vmap_prime_r[0]; 
  f_rl[1] = (1.58113883008419*fr[5]-1.224744871391589*fr[3]+0.7071067811865475*fr[1])/vmap_prime_r[0]; 

  fUpR[0] = (0.3535533905932737*f_cr[1]-0.3535533905932737*f_rl[1])*sgn_alphaUpR[1]+(0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[1]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[1]; 

  } 
  double GhatL[2] = {0.};
  double GhatR[2] = {0.};
  GhatL[0] = 0.7071067811865475*(alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.7071067811865475*(alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 

  GhatR[0] = 0.7071067811865475*(alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.7071067811865475*(alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[3] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[4] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdvpar2; 
  out[5] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdvpar2; 

  double vmap_prime_min = fmin(fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_c[0])),fabs(vmap_prime_r[0]));
  double cflFreq = fmax(fabs(alphaL[0]/vmap_prime_min), fabs(alphaR[0]/vmap_prime_min)); 
  return 1.767766952966369*rdvpar2*cflFreq; 

} 
