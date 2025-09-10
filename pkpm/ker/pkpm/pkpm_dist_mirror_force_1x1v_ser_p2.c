#include <gkyl_vlasov_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1) 
{ 
  // w[NDIM]:          Cell-center coordinates. 
  // dxv[NDIM]:        Cell spacing. 
  // pkpm_prim:        Input primitive variables [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp]. 
  // nu_prim_moms_sum: Input sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // div_b:            Input volume expansion of div(b). 
  // pkpm_accel_vars:  Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // f:                Input distribution functions [F_0, T_perp/m G = T_perp/m (F_0 - F_1)].
  // F_k_p_1:          Input k+1 distribution function. F_2 expansion is the first NP coefficients. 
  // g_dist_source:    Output [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)),  
  //                   (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ].
  //                   First output is mirror force source, second output is vperp characteristics source.
  // F_k_m_1:          Output k-1 distribution function. F_1 expansion is the first NP coefficients. 

  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *T_perp_over_m = &pkpm_prim[12]; 
  const double *T_perp_over_m_inv = &pkpm_prim[15]; 

  const double *p_perp_source = &pkpm_accel_vars[9]; 

  double alpha_G_1_source[8] = {0.0}; 

  alpha_G_1_source[0] = 1.414213562373095*p_perp_source[0]-1.414213562373095*div_b[0]*wvpar; 
  alpha_G_1_source[1] = 1.414213562373095*p_perp_source[1]-1.414213562373095*div_b[1]*wvpar; 
  alpha_G_1_source[2] = -0.408248290463863*div_b[0]*dvpar; 
  alpha_G_1_source[3] = -0.408248290463863*div_b[1]*dvpar; 
  alpha_G_1_source[4] = 1.414213562373095*p_perp_source[2]-1.414213562373095*div_b[2]*wvpar; 
  alpha_G_1_source[6] = -0.408248290463863*div_b[2]*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[8]; 
  const double *F_2 = &F_k_p_1[0]; 

  const double *nu_vtsq_sum = &nu_prim_moms_sum[3];

  double *out_G_1_mirror = &g_dist_source[0]; 
  double *out_G_1_vperp = &g_dist_source[8]; 
  double *out_F_1 = &F_k_m_1[0]; 

  double tmp_F_2_m_F_0[8] = {0.0}; 
  double tmp_T_perp_F_0[8] = {0.0}; 
  double tmp_T_perp_g_dist[8] = {0.0}; 
  double tmp_F_0_m_F_1[8] = {0.0}; 

  tmp_F_2_m_F_0[0] = F_2[0]-1.0*F_0[0]; 
  tmp_F_2_m_F_0[1] = F_2[1]-1.0*F_0[1]; 
  tmp_F_2_m_F_0[2] = F_2[2]-1.0*F_0[2]; 
  tmp_F_2_m_F_0[3] = F_2[3]-1.0*F_0[3]; 
  tmp_F_2_m_F_0[4] = F_2[4]-1.0*F_0[4]; 
  tmp_F_2_m_F_0[5] = F_2[5]-1.0*F_0[5]; 
  tmp_F_2_m_F_0[6] = F_2[6]-1.0*F_0[6]; 
  tmp_F_2_m_F_0[7] = F_2[7]-1.0*F_0[7]; 

  tmp_T_perp_F_0[0] = 0.7071067811865475*T_perp_over_m[2]*tmp_F_2_m_F_0[4]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[1]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[0]; 
  tmp_T_perp_F_0[1] = 0.6324555320336759*T_perp_over_m[1]*tmp_F_2_m_F_0[4]+0.6324555320336759*tmp_F_2_m_F_0[1]*T_perp_over_m[2]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[1]+0.7071067811865475*tmp_F_2_m_F_0[0]*T_perp_over_m[1]; 
  tmp_T_perp_F_0[2] = 0.7071067811865475*T_perp_over_m[2]*tmp_F_2_m_F_0[6]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[2]; 
  tmp_T_perp_F_0[3] = 0.632455532033676*T_perp_over_m[1]*tmp_F_2_m_F_0[6]+0.6324555320336759*T_perp_over_m[2]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[2]; 
  tmp_T_perp_F_0[4] = 0.4517539514526256*T_perp_over_m[2]*tmp_F_2_m_F_0[4]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[4]+0.7071067811865475*tmp_F_2_m_F_0[0]*T_perp_over_m[2]+0.6324555320336759*T_perp_over_m[1]*tmp_F_2_m_F_0[1]; 
  tmp_T_perp_F_0[5] = 0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[7]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[5]; 
  tmp_T_perp_F_0[6] = 0.4517539514526256*T_perp_over_m[2]*tmp_F_2_m_F_0[6]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[6]+0.632455532033676*T_perp_over_m[1]*tmp_F_2_m_F_0[3]+0.7071067811865475*T_perp_over_m[2]*tmp_F_2_m_F_0[2]; 
  tmp_T_perp_F_0[7] = 0.6324555320336759*T_perp_over_m[2]*tmp_F_2_m_F_0[7]+0.7071067811865475*T_perp_over_m[0]*tmp_F_2_m_F_0[7]+0.7071067811865475*T_perp_over_m[1]*tmp_F_2_m_F_0[5]; 

  tmp_T_perp_g_dist[0] = tmp_T_perp_F_0[0]+2.0*G_1[0]; 
  tmp_T_perp_g_dist[1] = tmp_T_perp_F_0[1]+2.0*G_1[1]; 
  tmp_T_perp_g_dist[2] = tmp_T_perp_F_0[2]+2.0*G_1[2]; 
  tmp_T_perp_g_dist[3] = tmp_T_perp_F_0[3]+2.0*G_1[3]; 
  tmp_T_perp_g_dist[4] = tmp_T_perp_F_0[4]+2.0*G_1[4]; 
  tmp_T_perp_g_dist[5] = tmp_T_perp_F_0[5]+2.0*G_1[5]; 
  tmp_T_perp_g_dist[6] = tmp_T_perp_F_0[6]+2.0*G_1[6]; 
  tmp_T_perp_g_dist[7] = tmp_T_perp_F_0[7]+2.0*G_1[7]; 

  out_G_1_mirror[0] = 1.414213562373095*T_perp_over_m[2]*tmp_T_perp_g_dist[4]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[1]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[0]; 
  out_G_1_mirror[1] = 1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[4]+1.264911064067352*tmp_T_perp_g_dist[1]*T_perp_over_m[2]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[1]+1.414213562373095*tmp_T_perp_g_dist[0]*T_perp_over_m[1]; 
  out_G_1_mirror[2] = 1.414213562373095*T_perp_over_m[2]*tmp_T_perp_g_dist[6]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[2]; 
  out_G_1_mirror[3] = 1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[6]+1.264911064067352*T_perp_over_m[2]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[2]; 
  out_G_1_mirror[4] = 0.9035079029052515*T_perp_over_m[2]*tmp_T_perp_g_dist[4]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[4]+1.414213562373095*tmp_T_perp_g_dist[0]*T_perp_over_m[2]+1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[1]; 
  out_G_1_mirror[5] = 1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[7]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[5]; 
  out_G_1_mirror[6] = 0.9035079029052515*T_perp_over_m[2]*tmp_T_perp_g_dist[6]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[6]+1.264911064067352*T_perp_over_m[1]*tmp_T_perp_g_dist[3]+1.414213562373095*T_perp_over_m[2]*tmp_T_perp_g_dist[2]; 
  out_G_1_mirror[7] = 1.264911064067352*T_perp_over_m[2]*tmp_T_perp_g_dist[7]+1.414213562373095*T_perp_over_m[0]*tmp_T_perp_g_dist[7]+1.414213562373095*T_perp_over_m[1]*tmp_T_perp_g_dist[5]; 

  out_G_1_vperp[0] = 0.5*G_1[6]*alpha_G_1_source[6]+0.5*G_1[4]*alpha_G_1_source[4]+1.414213562373095*nu_vtsq_sum[2]*F_0[4]+0.5*G_1[3]*alpha_G_1_source[3]+0.5*G_1[2]*alpha_G_1_source[2]+1.414213562373095*F_0[1]*nu_vtsq_sum[1]+0.5*G_1[1]*alpha_G_1_source[1]+1.414213562373095*F_0[0]*nu_vtsq_sum[0]+0.5*G_1[0]*alpha_G_1_source[0]; 
  out_G_1_vperp[1] = 0.447213595499958*G_1[3]*alpha_G_1_source[6]+0.447213595499958*alpha_G_1_source[3]*G_1[6]+0.4472135954999579*G_1[1]*alpha_G_1_source[4]+0.4472135954999579*alpha_G_1_source[1]*G_1[4]+1.264911064067352*nu_vtsq_sum[1]*F_0[4]+0.5*G_1[2]*alpha_G_1_source[3]+0.5*alpha_G_1_source[2]*G_1[3]+1.264911064067352*F_0[1]*nu_vtsq_sum[2]+1.414213562373095*F_0[0]*nu_vtsq_sum[1]+0.5*G_1[0]*alpha_G_1_source[1]+0.5*alpha_G_1_source[0]*G_1[1]+1.414213562373095*nu_vtsq_sum[0]*F_0[1]; 
  out_G_1_vperp[2] = 0.447213595499958*alpha_G_1_source[3]*G_1[7]+0.5000000000000001*G_1[4]*alpha_G_1_source[6]+0.5000000000000001*alpha_G_1_source[4]*G_1[6]+1.414213562373095*nu_vtsq_sum[2]*F_0[6]+0.4472135954999579*alpha_G_1_source[2]*G_1[5]+0.5*G_1[1]*alpha_G_1_source[3]+0.5*alpha_G_1_source[1]*G_1[3]+1.414213562373095*nu_vtsq_sum[1]*F_0[3]+0.5*G_1[0]*alpha_G_1_source[2]+0.5*alpha_G_1_source[0]*G_1[2]+1.414213562373095*nu_vtsq_sum[0]*F_0[2]; 
  out_G_1_vperp[3] = 0.4*alpha_G_1_source[6]*G_1[7]+0.447213595499958*alpha_G_1_source[2]*G_1[7]+0.447213595499958*G_1[1]*alpha_G_1_source[6]+0.447213595499958*alpha_G_1_source[1]*G_1[6]+1.264911064067352*nu_vtsq_sum[1]*F_0[6]+0.4472135954999579*alpha_G_1_source[3]*G_1[5]+0.4472135954999579*G_1[3]*alpha_G_1_source[4]+0.4472135954999579*alpha_G_1_source[3]*G_1[4]+0.5*G_1[0]*alpha_G_1_source[3]+0.5*alpha_G_1_source[0]*G_1[3]+1.264911064067352*nu_vtsq_sum[2]*F_0[3]+1.414213562373095*nu_vtsq_sum[0]*F_0[3]+0.5*G_1[1]*alpha_G_1_source[2]+0.5*alpha_G_1_source[1]*G_1[2]+1.414213562373095*nu_vtsq_sum[1]*F_0[2]; 
  out_G_1_vperp[4] = 0.31943828249997*G_1[6]*alpha_G_1_source[6]+0.5000000000000001*G_1[2]*alpha_G_1_source[6]+0.5000000000000001*alpha_G_1_source[2]*G_1[6]+0.31943828249997*G_1[4]*alpha_G_1_source[4]+0.5*G_1[0]*alpha_G_1_source[4]+0.5*alpha_G_1_source[0]*G_1[4]+0.9035079029052515*nu_vtsq_sum[2]*F_0[4]+1.414213562373095*nu_vtsq_sum[0]*F_0[4]+0.4472135954999579*G_1[3]*alpha_G_1_source[3]+1.414213562373095*F_0[0]*nu_vtsq_sum[2]+1.264911064067352*F_0[1]*nu_vtsq_sum[1]+0.4472135954999579*G_1[1]*alpha_G_1_source[1]; 
  out_G_1_vperp[5] = 0.5000000000000001*alpha_G_1_source[1]*G_1[7]+1.414213562373095*nu_vtsq_sum[1]*F_0[7]+0.4472135954999579*G_1[6]*alpha_G_1_source[6]+0.5*alpha_G_1_source[0]*G_1[5]+1.414213562373095*nu_vtsq_sum[0]*F_0[5]+0.4472135954999579*G_1[3]*alpha_G_1_source[3]+0.4472135954999579*G_1[2]*alpha_G_1_source[2]; 
  out_G_1_vperp[6] = 0.4*alpha_G_1_source[3]*G_1[7]+0.4472135954999579*G_1[5]*alpha_G_1_source[6]+0.31943828249997*G_1[4]*alpha_G_1_source[6]+0.5*G_1[0]*alpha_G_1_source[6]+0.31943828249997*alpha_G_1_source[4]*G_1[6]+0.5*alpha_G_1_source[0]*G_1[6]+0.9035079029052515*nu_vtsq_sum[2]*F_0[6]+1.414213562373095*nu_vtsq_sum[0]*F_0[6]+0.5000000000000001*G_1[2]*alpha_G_1_source[4]+0.5000000000000001*alpha_G_1_source[2]*G_1[4]+0.447213595499958*G_1[1]*alpha_G_1_source[3]+0.447213595499958*alpha_G_1_source[1]*G_1[3]+1.264911064067352*nu_vtsq_sum[1]*F_0[3]+1.414213562373095*F_0[2]*nu_vtsq_sum[2]; 
  out_G_1_vperp[7] = 0.4472135954999579*alpha_G_1_source[4]*G_1[7]+0.5*alpha_G_1_source[0]*G_1[7]+1.264911064067352*nu_vtsq_sum[2]*F_0[7]+1.414213562373095*nu_vtsq_sum[0]*F_0[7]+0.4*G_1[3]*alpha_G_1_source[6]+0.4*alpha_G_1_source[3]*G_1[6]+0.5000000000000001*alpha_G_1_source[1]*G_1[5]+1.414213562373095*nu_vtsq_sum[1]*F_0[5]+0.447213595499958*G_1[2]*alpha_G_1_source[3]+0.447213595499958*alpha_G_1_source[2]*G_1[3]; 

  tmp_F_0_m_F_1[0] = 0.7071067811865475*T_perp_over_m_inv[2]*G_1[4]+0.7071067811865475*G_1[1]*T_perp_over_m_inv[1]+0.7071067811865475*G_1[0]*T_perp_over_m_inv[0]; 
  tmp_F_0_m_F_1[1] = 0.6324555320336759*T_perp_over_m_inv[1]*G_1[4]+0.6324555320336759*G_1[1]*T_perp_over_m_inv[2]+0.7071067811865475*G_1[0]*T_perp_over_m_inv[1]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[1]; 
  tmp_F_0_m_F_1[2] = 0.7071067811865475*T_perp_over_m_inv[2]*G_1[6]+0.7071067811865475*T_perp_over_m_inv[1]*G_1[3]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[2]; 
  tmp_F_0_m_F_1[3] = 0.632455532033676*T_perp_over_m_inv[1]*G_1[6]+0.6324555320336759*T_perp_over_m_inv[2]*G_1[3]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[3]+0.7071067811865475*T_perp_over_m_inv[1]*G_1[2]; 
  tmp_F_0_m_F_1[4] = 0.4517539514526256*T_perp_over_m_inv[2]*G_1[4]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[4]+0.7071067811865475*G_1[0]*T_perp_over_m_inv[2]+0.6324555320336759*G_1[1]*T_perp_over_m_inv[1]; 
  tmp_F_0_m_F_1[5] = 0.7071067811865475*T_perp_over_m_inv[1]*G_1[7]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[5]; 
  tmp_F_0_m_F_1[6] = 0.4517539514526256*T_perp_over_m_inv[2]*G_1[6]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[6]+0.632455532033676*T_perp_over_m_inv[1]*G_1[3]+0.7071067811865475*G_1[2]*T_perp_over_m_inv[2]; 
  tmp_F_0_m_F_1[7] = 0.6324555320336759*T_perp_over_m_inv[2]*G_1[7]+0.7071067811865475*T_perp_over_m_inv[0]*G_1[7]+0.7071067811865475*T_perp_over_m_inv[1]*G_1[5]; 

  out_F_1[0] = F_0[0]-1.0*tmp_F_0_m_F_1[0]; 
  out_F_1[1] = F_0[1]-1.0*tmp_F_0_m_F_1[1]; 
  out_F_1[2] = F_0[2]-1.0*tmp_F_0_m_F_1[2]; 
  out_F_1[3] = F_0[3]-1.0*tmp_F_0_m_F_1[3]; 
  out_F_1[4] = F_0[4]-1.0*tmp_F_0_m_F_1[4]; 
  out_F_1[5] = F_0[5]-1.0*tmp_F_0_m_F_1[5]; 
  out_F_1[6] = F_0[6]-1.0*tmp_F_0_m_F_1[6]; 
  out_F_1[7] = F_0[7]-1.0*tmp_F_0_m_F_1[7]; 
} 
