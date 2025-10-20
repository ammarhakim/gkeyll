#include <gkyl_vlasov_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_dist_div_ppar_x_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar) 
{ 
  // w[NDIM]:         Cell-center coordinates.
  // dxv[NDIM]:       Cell spacing.
  // bvar_surf_l/c/r: Input surface magnetic field unit vector and tensor in left/center/right cells in each direction.
  // fl/fc/fr:        Input distribution functions [F_0, T_perp/m G = T_perp/m (F_0 - F_1)] in left/center/right cells.
  // bvar_c:          Input volume expansion of magnetic field unit vector and tensor in center cell.
  // pkpm_max_b:      Input surface expansion of max |b| for Lax penalization of streaming: lambda_i = |b_i|.
  // pkpm_div_ppar:   Increment to volume expansion of div(p_par b).
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double volFact = dxv[1]/2.0; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  double F_0_lr[3] = {0.0}; 
  double F_0_cl[3] = {0.0}; 
  double F_0_cr[3] = {0.0}; 
  double F_0_rl[3] = {0.0}; 
  double Ghat_F_0_vpar_l[2] = {0.0}; 
  double Ghat_F_0_vpar_r[2] = {0.0}; 

  F_0_lr[0] = 1.224744871391589*F_0l[1]+0.7071067811865475*F_0l[0]; 
  F_0_lr[1] = 1.224744871391589*F_0l[3]+0.7071067811865475*F_0l[2]; 
  F_0_lr[2] = 1.224744871391589*F_0l[5]+0.7071067811865475*F_0l[4]; 

  F_0_cl[0] = 0.7071067811865475*F_0c[0]-1.224744871391589*F_0c[1]; 
  F_0_cl[1] = 0.7071067811865475*F_0c[2]-1.224744871391589*F_0c[3]; 
  F_0_cl[2] = 0.7071067811865475*F_0c[4]-1.224744871391589*F_0c[5]; 

  F_0_cr[0] = 1.224744871391589*F_0c[1]+0.7071067811865475*F_0c[0]; 
  F_0_cr[1] = 1.224744871391589*F_0c[3]+0.7071067811865475*F_0c[2]; 
  F_0_cr[2] = 1.224744871391589*F_0c[5]+0.7071067811865475*F_0c[4]; 

  F_0_rl[0] = 0.7071067811865475*F_0r[0]-1.224744871391589*F_0r[1]; 
  F_0_rl[1] = 0.7071067811865475*F_0r[2]-1.224744871391589*F_0r[3]; 
  F_0_rl[2] = 0.7071067811865475*F_0r[4]-1.224744871391589*F_0r[5]; 

  const double *b_surf_lr = &bvar_surf_l[1]; 
  const double *b_surf_cl = &bvar_surf_c[0]; 
  const double *b_surf_cr = &bvar_surf_c[1]; 
  const double *b_surf_rl = &bvar_surf_r[0]; 

  const double *pkpm_max_b_l = &pkpm_max_b[0]; 
  const double *pkpm_max_b_r = &pkpm_max_b[1]; 

  const double *b_c = &bvar_c[0]; 
  double alpha_c[6] = {0.0}; 
  alpha_c[0] = 1.414213562373095*b_c[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*b_c[1]*wvpar; 
  alpha_c[2] = 0.408248290463863*b_c[0]*dvpar; 
  alpha_c[3] = 0.408248290463863*b_c[1]*dvpar; 

  double bl_r = b_surf_lr[0]; 
  double bc_l = b_surf_cl[0]; 
  double bc_r = b_surf_cr[0]; 
  double br_l = b_surf_rl[0]; 
  double avg_b_l = 0.5*(bl_r + bc_l); 
  double avg_b_r = 0.5*(bc_r + br_l); 

  double max_b_l = pkpm_max_b_l[0]; 
  double max_b_r = pkpm_max_b_r[0]; 

  if (wvpar>0) { 

  Ghat_F_0_vpar_l[0] = ((0.5*F_0_lr[0]-0.5*F_0_cl[0])*max_b_l+0.5*(F_0_lr[0]+F_0_cl[0])*avg_b_l)*wvpar+dvpar*((0.1443375672974065*F_0_lr[1]-0.1443375672974065*F_0_cl[1])*max_b_l+0.1443375672974065*(F_0_lr[1]+F_0_cl[1])*avg_b_l); 
  Ghat_F_0_vpar_l[1] = ((0.5*F_0_lr[1]-0.5*F_0_cl[1])*max_b_l+0.5*(F_0_lr[1]+F_0_cl[1])*avg_b_l)*wvpar+dvpar*((0.1290994448735806*F_0_lr[2]-0.1290994448735806*F_0_cl[2]+0.1443375672974065*F_0_lr[0]-0.1443375672974065*F_0_cl[0])*max_b_l+(0.1290994448735806*(F_0_lr[2]+F_0_cl[2])+0.1443375672974065*(F_0_lr[0]+F_0_cl[0]))*avg_b_l); 

  Ghat_F_0_vpar_r[0] = (0.5*(F_0_cr[0]*max_b_r+(F_0_rl[0]+F_0_cr[0])*avg_b_r)-0.5*F_0_rl[0]*max_b_r)*wvpar+dvpar*(0.1443375672974065*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.1443375672974065*F_0_rl[1]*max_b_r); 
  Ghat_F_0_vpar_r[1] = (0.5*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.5*F_0_rl[1]*max_b_r)*wvpar+dvpar*(((-0.1290994448735806*F_0_rl[2])+0.1290994448735806*F_0_cr[2]-0.1443375672974065*F_0_rl[0]+0.1443375672974065*F_0_cr[0])*max_b_r+(0.1290994448735806*(F_0_rl[2]+F_0_cr[2])+0.1443375672974065*(F_0_rl[0]+F_0_cr[0]))*avg_b_r); 

  } else { 

  Ghat_F_0_vpar_l[0] = (0.5*(F_0_cl[0]*max_b_l+(F_0_lr[0]+F_0_cl[0])*avg_b_l)-0.5*F_0_lr[0]*max_b_l)*wvpar+dvpar*(0.1443375672974065*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.1443375672974065*F_0_lr[1]*max_b_l); 
  Ghat_F_0_vpar_l[1] = (0.5*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.5*F_0_lr[1]*max_b_l)*wvpar+dvpar*(((-0.1290994448735806*F_0_lr[2])+0.1290994448735806*F_0_cl[2]-0.1443375672974065*F_0_lr[0]+0.1443375672974065*F_0_cl[0])*max_b_l+(0.1290994448735806*(F_0_lr[2]+F_0_cl[2])+0.1443375672974065*(F_0_lr[0]+F_0_cl[0]))*avg_b_l); 

  Ghat_F_0_vpar_r[0] = ((0.5*F_0_rl[0]-0.5*F_0_cr[0])*max_b_r+0.5*(F_0_rl[0]+F_0_cr[0])*avg_b_r)*wvpar+dvpar*((0.1443375672974065*F_0_rl[1]-0.1443375672974065*F_0_cr[1])*max_b_r+0.1443375672974065*(F_0_rl[1]+F_0_cr[1])*avg_b_r); 
  Ghat_F_0_vpar_r[1] = ((0.5*F_0_rl[1]-0.5*F_0_cr[1])*max_b_r+0.5*(F_0_rl[1]+F_0_cr[1])*avg_b_r)*wvpar+dvpar*((0.1290994448735806*F_0_rl[2]-0.1290994448735806*F_0_cr[2]+0.1443375672974065*F_0_rl[0]-0.1443375672974065*F_0_cr[0])*max_b_r+(0.1290994448735806*(F_0_rl[2]+F_0_cr[2])+0.1443375672974065*(F_0_rl[0]+F_0_cr[0]))*avg_b_r); 

  } 
  pkpm_div_ppar[0] += dx1*volFact*((Ghat_F_0_vpar_r[0]-1.0*Ghat_F_0_vpar_l[0])*wvpar+(0.2886751345948129*Ghat_F_0_vpar_r[1]-0.2886751345948129*Ghat_F_0_vpar_l[1])*dvpar); 
  pkpm_div_ppar[1] += dx1*volFact*((1.732050807568877*(Ghat_F_0_vpar_r[0]+Ghat_F_0_vpar_l[0])-1.224744871391589*(F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2]+F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0]))*wvpar+((-0.3162277660168379*alpha_c[3]*F_0c[5])-0.3162277660168379*alpha_c[2]*F_0c[4]-0.3535533905932737*(F_0c[1]*alpha_c[3]+alpha_c[1]*F_0c[3]+F_0c[0]*alpha_c[2]+alpha_c[0]*F_0c[2])+0.5*(Ghat_F_0_vpar_r[1]+Ghat_F_0_vpar_l[1]))*dvpar); 

} 
