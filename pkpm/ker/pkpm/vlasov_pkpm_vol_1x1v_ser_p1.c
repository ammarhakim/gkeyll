#include <gkyl_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *pkpm_prim, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *g_dist_source, const double *f, 
  double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:         Cell-center coordinates. 
  // dxv[NDIM]:       Cell spacing. 
  // bvar:            Input magnetic field unit vector and tensor (nine components; first three components, b_i, other six components, b_i b_j). 
  // pkpm_prim:       Input primitive variables [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp]. 
  // div_b:           Input volume expansion of div(b). 
  // pkpm_accel_vars: Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // g_dist_source:   Input [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
  //                  (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]. 
  //                  First input is mirror force source, second input is vperp characteristics source. 
  // f:               Input distribution functions [F_0, T_perp/m G = T_perp/m (F_0 - F_1)].
  // out:             Incremented output distribution functions. 
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &pkpm_prim[0]; 
  const double *uy = &pkpm_prim[2]; 
  const double *uz = &pkpm_prim[4]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  const double *F_0_source = &f[6]; 
  const double *G_1_source = &g_dist_source[0]; 
  const double *G_1_vperp = &g_dist_source[6]; 
  const double *p_perp_div_b = &pkpm_accel_vars[0]; 
  const double *bb_grad_u = &pkpm_accel_vars[2]; 
  const double *p_force = &pkpm_accel_vars[4]; 

  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[6] = {0.0}; 
  double alpha_vdim[6] = {0.0}; 
  double alpha_div_b[6] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*dx0*(bx[0]*wvpar+ux[0]); 
  alpha_cdim[1] = 1.414213562373095*dx0*(bx[1]*wvpar+ux[1]); 
  alpha_cdim[2] = 0.408248290463863*bx[0]*dvpar*dx0; 
  alpha_cdim[3] = 0.408248290463863*bx[1]*dvpar*dx0; 

  alpha_vdim[0] = 1.414213562373095*p_force[0]*dv1par-1.414213562373095*bb_grad_u[0]*dv1par*wvpar; 
  alpha_vdim[1] = 1.414213562373095*p_force[1]*dv1par-1.414213562373095*bb_grad_u[1]*dv1par*wvpar; 
  alpha_vdim[2] = -0.408248290463863*bb_grad_u[0]*dv1par*dvpar; 
  alpha_vdim[3] = -0.408248290463863*bb_grad_u[1]*dv1par*dvpar; 

  alpha_div_b[0] = 1.414213562373095*div_b[0]*dv1par; 
  alpha_div_b[1] = 1.414213562373095*div_b[1]*dv1par; 

  cflFreq_mid += 5.0*fabs(0.3535533905932737*p_perp_div_b[0])*dv1par; 

  out_F_0[1] += 0.8660254037844386*F_0[3]*alpha_cdim[3]+0.8660254037844386*F_0[2]*alpha_cdim[2]+0.8660254037844386*F_0[1]*alpha_cdim[1]+0.8660254037844386*F_0[0]*alpha_cdim[0]; 
  out_F_0[2] += 0.8660254037844386*F_0[3]*alpha_vdim[3]+0.8660254037844386*F_0[2]*alpha_vdim[2]+0.8660254037844386*F_0[1]*alpha_vdim[1]+0.8660254037844386*F_0_source[1]*alpha_div_b[1]+0.8660254037844386*F_0[0]*alpha_vdim[0]+0.8660254037844386*F_0_source[0]*alpha_div_b[0]; 
  out_F_0[3] += 0.7745966692414834*alpha_cdim[3]*F_0[5]+0.7745966692414833*alpha_cdim[2]*F_0[4]+0.8660254037844386*F_0[2]*alpha_vdim[3]+0.8660254037844386*F_0[1]*alpha_cdim[3]+0.8660254037844386*alpha_vdim[2]*F_0[3]+0.8660254037844386*alpha_cdim[1]*F_0[3]+0.8660254037844386*F_0[0]*alpha_cdim[2]+0.8660254037844386*alpha_cdim[0]*F_0[2]+0.8660254037844386*F_0[0]*alpha_vdim[1]+0.8660254037844386*F_0_source[0]*alpha_div_b[1]+0.8660254037844386*alpha_div_b[0]*F_0_source[1]+0.8660254037844386*alpha_vdim[0]*F_0[1]; 
  out_F_0[4] += 1.732050807568877*alpha_vdim[3]*F_0[5]+1.732050807568877*alpha_vdim[2]*F_0[4]+1.936491673103709*F_0[1]*alpha_vdim[3]+1.936491673103709*alpha_div_b[1]*F_0_source[3]+1.936491673103709*alpha_vdim[1]*F_0[3]+1.936491673103709*F_0[0]*alpha_vdim[2]+1.936491673103709*alpha_div_b[0]*F_0_source[2]+1.936491673103709*alpha_vdim[0]*F_0[2]; 
  out_F_0[5] += 1.732050807568877*alpha_vdim[2]*F_0[5]+0.8660254037844386*alpha_cdim[1]*F_0[5]+1.732050807568877*alpha_vdim[3]*F_0[4]+0.8660254037844387*alpha_cdim[0]*F_0[4]+1.936491673103709*F_0[0]*alpha_vdim[3]+0.7745966692414834*F_0[3]*alpha_cdim[3]+1.936491673103709*alpha_div_b[0]*F_0_source[3]+1.936491673103709*alpha_vdim[0]*F_0[3]+1.936491673103709*F_0[1]*alpha_vdim[2]+0.7745966692414834*F_0[2]*alpha_cdim[2]+1.936491673103709*alpha_div_b[1]*F_0_source[2]+1.936491673103709*alpha_vdim[1]*F_0[2]; 
  out_G_1[0] += G_1_vperp[0]; 
  out_G_1[1] += 0.8660254037844386*G_1[3]*alpha_cdim[3]+0.8660254037844386*G_1[2]*alpha_cdim[2]+0.8660254037844386*G_1[1]*alpha_cdim[1]+G_1_vperp[1]+0.8660254037844386*G_1[0]*alpha_cdim[0]; 
  out_G_1[2] += 0.8660254037844386*G_1[3]*alpha_vdim[3]+0.8660254037844386*G_1[2]*alpha_vdim[2]+G_1_vperp[2]+0.8660254037844386*G_1[1]*alpha_vdim[1]+0.8660254037844386*G_1_source[1]*alpha_div_b[1]+0.8660254037844386*G_1[0]*alpha_vdim[0]+0.8660254037844386*G_1_source[0]*alpha_div_b[0]; 
  out_G_1[3] += 0.7745966692414834*alpha_cdim[3]*G_1[5]+0.7745966692414833*alpha_cdim[2]*G_1[4]+0.8660254037844386*G_1[2]*alpha_vdim[3]+0.8660254037844386*G_1[1]*alpha_cdim[3]+G_1_vperp[3]+0.8660254037844386*alpha_vdim[2]*G_1[3]+0.8660254037844386*alpha_cdim[1]*G_1[3]+0.8660254037844386*G_1[0]*alpha_cdim[2]+0.8660254037844386*alpha_cdim[0]*G_1[2]+0.8660254037844386*G_1[0]*alpha_vdim[1]+0.8660254037844386*G_1_source[0]*alpha_div_b[1]+0.8660254037844386*alpha_div_b[0]*G_1_source[1]+0.8660254037844386*alpha_vdim[0]*G_1[1]; 
  out_G_1[4] += 1.732050807568877*alpha_vdim[3]*G_1[5]+G_1_vperp[4]+1.732050807568877*alpha_vdim[2]*G_1[4]+1.936491673103709*G_1[1]*alpha_vdim[3]+1.936491673103709*alpha_div_b[1]*G_1_source[3]+1.936491673103709*alpha_vdim[1]*G_1[3]+1.936491673103709*G_1[0]*alpha_vdim[2]+1.936491673103709*alpha_div_b[0]*G_1_source[2]+1.936491673103709*alpha_vdim[0]*G_1[2]; 
  out_G_1[5] += G_1_vperp[5]+1.732050807568877*alpha_vdim[2]*G_1[5]+0.8660254037844386*alpha_cdim[1]*G_1[5]+1.732050807568877*alpha_vdim[3]*G_1[4]+0.8660254037844387*alpha_cdim[0]*G_1[4]+1.936491673103709*G_1[0]*alpha_vdim[3]+0.7745966692414834*G_1[3]*alpha_cdim[3]+1.936491673103709*alpha_div_b[0]*G_1_source[3]+1.936491673103709*alpha_vdim[0]*G_1[3]+1.936491673103709*G_1[1]*alpha_vdim[2]+0.7745966692414834*G_1[2]*alpha_cdim[2]+1.936491673103709*alpha_div_b[1]*G_1_source[2]+1.936491673103709*alpha_vdim[1]*G_1[2]; 

  return cflFreq_mid; 
} 
