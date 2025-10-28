#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_io_2x_tensor_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io) 
{ 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model. 
  // euler_pkpm:       Input [rho ux, rho uy, rho uz], Fluid state vector. 
  // p_ij:             Input p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij. 
  // prim:             Input [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp]. 
  // pkpm_accel:       Input volume expansion of pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // fluid_io:         Output fluid conserved variables. 
  // pkpm_vars_io:     Output pkpm variables (primitive and acceleration). 

  const double *rho = &vlasov_pkpm_moms[0]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[9]; 
  const double *rhouz = &euler_pkpm[18]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[9]; 
  const double *Pxz = &p_ij[18]; 
  const double *Pyy = &p_ij[27]; 
  const double *Pyz = &p_ij[36]; 
  const double *Pzz = &p_ij[45]; 

  const double *ux = &prim[0]; 
  const double *uy = &prim[9]; 
  const double *uz = &prim[18]; 
  const double *pkpm_div_ppar = &prim[27]; 
  const double *T_perp_over_m = &prim[36]; 
  const double *T_perp_over_m_inv = &prim[45]; 

  const double *p_perp_div_b = &pkpm_accel[0]; 
  const double *bb_grad_u = &pkpm_accel[9]; 

  double *fluid_io_rho = &fluid_io[0]; 
  double *fluid_io_rhoux = &fluid_io[9]; 
  double *fluid_io_rhouy = &fluid_io[18]; 
  double *fluid_io_rhouz = &fluid_io[27]; 
  double *fluid_io_Sxx = &fluid_io[36]; 
  double *fluid_io_Sxy = &fluid_io[45]; 
  double *fluid_io_Sxz = &fluid_io[54]; 
  double *fluid_io_Syy = &fluid_io[63]; 
  double *fluid_io_Syz = &fluid_io[72]; 
  double *fluid_io_Szz = &fluid_io[81]; 

  double *pkpm_vars_io_ux = &pkpm_vars_io[0]; 
  double *pkpm_vars_io_uy = &pkpm_vars_io[9]; 
  double *pkpm_vars_io_uz = &pkpm_vars_io[18]; 
  double *pkpm_vars_io_T_perp_over_m = &pkpm_vars_io[27]; 
  double *pkpm_vars_io_T_perp_over_m_inv = &pkpm_vars_io[36]; 
  double *pkpm_vars_io_pkpm_div_ppar = &pkpm_vars_io[45]; 
  double *pkpm_vars_io_p_perp_div_b = &pkpm_vars_io[54]; 
  double *pkpm_vars_io_bb_grad_u = &pkpm_vars_io[63]; 

  double rhouxux[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhoux, ux, rhouxux); 
 
  double rhouxuy[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhoux, uy, rhouxuy); 
 
  double rhouxuz[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhoux, uz, rhouxuz); 
 
  double rhouyuy[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhouy, uy, rhouyuy); 
 
  double rhouyuz[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhouy, uz, rhouyuz); 
 
  double rhouzuz[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhouz, uz, rhouzuz); 
 
  fluid_io_rho[0] = rho[0]; 
  fluid_io_rhoux[0] = rhoux[0]; 
  fluid_io_rhouy[0] = rhouy[0]; 
  fluid_io_rhouz[0] = rhouz[0]; 
  fluid_io_Sxx[0] = Pxx[0] + rhouxux[0]; 
  fluid_io_Sxy[0] = Pxy[0] + rhouxuy[0]; 
  fluid_io_Sxz[0] = Pxz[0] + rhouxuz[0]; 
  fluid_io_Syy[0] = Pyy[0] + rhouyuy[0]; 
  fluid_io_Syz[0] = Pyz[0] + rhouyuz[0]; 
  fluid_io_Szz[0] = Pzz[0] + rhouzuz[0]; 
 
  pkpm_vars_io_ux[0] = ux[0]; 
  pkpm_vars_io_uy[0] = uy[0]; 
  pkpm_vars_io_uz[0] = uz[0]; 
  pkpm_vars_io_T_perp_over_m[0] = T_perp_over_m[0]; 
  pkpm_vars_io_T_perp_over_m_inv[0] = T_perp_over_m_inv[0]; 
  pkpm_vars_io_pkpm_div_ppar[0] = pkpm_div_ppar[0]; 
  pkpm_vars_io_p_perp_div_b[0] = p_perp_div_b[0]; 
  pkpm_vars_io_bb_grad_u[0] = bb_grad_u[0]; 
 
  fluid_io_rho[1] = rho[1]; 
  fluid_io_rhoux[1] = rhoux[1]; 
  fluid_io_rhouy[1] = rhouy[1]; 
  fluid_io_rhouz[1] = rhouz[1]; 
  fluid_io_Sxx[1] = Pxx[1] + rhouxux[1]; 
  fluid_io_Sxy[1] = Pxy[1] + rhouxuy[1]; 
  fluid_io_Sxz[1] = Pxz[1] + rhouxuz[1]; 
  fluid_io_Syy[1] = Pyy[1] + rhouyuy[1]; 
  fluid_io_Syz[1] = Pyz[1] + rhouyuz[1]; 
  fluid_io_Szz[1] = Pzz[1] + rhouzuz[1]; 
 
  pkpm_vars_io_ux[1] = ux[1]; 
  pkpm_vars_io_uy[1] = uy[1]; 
  pkpm_vars_io_uz[1] = uz[1]; 
  pkpm_vars_io_T_perp_over_m[1] = T_perp_over_m[1]; 
  pkpm_vars_io_T_perp_over_m_inv[1] = T_perp_over_m_inv[1]; 
  pkpm_vars_io_pkpm_div_ppar[1] = pkpm_div_ppar[1]; 
  pkpm_vars_io_p_perp_div_b[1] = p_perp_div_b[1]; 
  pkpm_vars_io_bb_grad_u[1] = bb_grad_u[1]; 
 
  fluid_io_rho[2] = rho[2]; 
  fluid_io_rhoux[2] = rhoux[2]; 
  fluid_io_rhouy[2] = rhouy[2]; 
  fluid_io_rhouz[2] = rhouz[2]; 
  fluid_io_Sxx[2] = Pxx[2] + rhouxux[2]; 
  fluid_io_Sxy[2] = Pxy[2] + rhouxuy[2]; 
  fluid_io_Sxz[2] = Pxz[2] + rhouxuz[2]; 
  fluid_io_Syy[2] = Pyy[2] + rhouyuy[2]; 
  fluid_io_Syz[2] = Pyz[2] + rhouyuz[2]; 
  fluid_io_Szz[2] = Pzz[2] + rhouzuz[2]; 
 
  pkpm_vars_io_ux[2] = ux[2]; 
  pkpm_vars_io_uy[2] = uy[2]; 
  pkpm_vars_io_uz[2] = uz[2]; 
  pkpm_vars_io_T_perp_over_m[2] = T_perp_over_m[2]; 
  pkpm_vars_io_T_perp_over_m_inv[2] = T_perp_over_m_inv[2]; 
  pkpm_vars_io_pkpm_div_ppar[2] = pkpm_div_ppar[2]; 
  pkpm_vars_io_p_perp_div_b[2] = p_perp_div_b[2]; 
  pkpm_vars_io_bb_grad_u[2] = bb_grad_u[2]; 
 
  fluid_io_rho[3] = rho[3]; 
  fluid_io_rhoux[3] = rhoux[3]; 
  fluid_io_rhouy[3] = rhouy[3]; 
  fluid_io_rhouz[3] = rhouz[3]; 
  fluid_io_Sxx[3] = Pxx[3] + rhouxux[3]; 
  fluid_io_Sxy[3] = Pxy[3] + rhouxuy[3]; 
  fluid_io_Sxz[3] = Pxz[3] + rhouxuz[3]; 
  fluid_io_Syy[3] = Pyy[3] + rhouyuy[3]; 
  fluid_io_Syz[3] = Pyz[3] + rhouyuz[3]; 
  fluid_io_Szz[3] = Pzz[3] + rhouzuz[3]; 
 
  pkpm_vars_io_ux[3] = ux[3]; 
  pkpm_vars_io_uy[3] = uy[3]; 
  pkpm_vars_io_uz[3] = uz[3]; 
  pkpm_vars_io_T_perp_over_m[3] = T_perp_over_m[3]; 
  pkpm_vars_io_T_perp_over_m_inv[3] = T_perp_over_m_inv[3]; 
  pkpm_vars_io_pkpm_div_ppar[3] = pkpm_div_ppar[3]; 
  pkpm_vars_io_p_perp_div_b[3] = p_perp_div_b[3]; 
  pkpm_vars_io_bb_grad_u[3] = bb_grad_u[3]; 
 
  fluid_io_rho[4] = rho[4]; 
  fluid_io_rhoux[4] = rhoux[4]; 
  fluid_io_rhouy[4] = rhouy[4]; 
  fluid_io_rhouz[4] = rhouz[4]; 
  fluid_io_Sxx[4] = Pxx[4] + rhouxux[4]; 
  fluid_io_Sxy[4] = Pxy[4] + rhouxuy[4]; 
  fluid_io_Sxz[4] = Pxz[4] + rhouxuz[4]; 
  fluid_io_Syy[4] = Pyy[4] + rhouyuy[4]; 
  fluid_io_Syz[4] = Pyz[4] + rhouyuz[4]; 
  fluid_io_Szz[4] = Pzz[4] + rhouzuz[4]; 
 
  pkpm_vars_io_ux[4] = ux[4]; 
  pkpm_vars_io_uy[4] = uy[4]; 
  pkpm_vars_io_uz[4] = uz[4]; 
  pkpm_vars_io_T_perp_over_m[4] = T_perp_over_m[4]; 
  pkpm_vars_io_T_perp_over_m_inv[4] = T_perp_over_m_inv[4]; 
  pkpm_vars_io_pkpm_div_ppar[4] = pkpm_div_ppar[4]; 
  pkpm_vars_io_p_perp_div_b[4] = p_perp_div_b[4]; 
  pkpm_vars_io_bb_grad_u[4] = bb_grad_u[4]; 
 
  fluid_io_rho[5] = rho[5]; 
  fluid_io_rhoux[5] = rhoux[5]; 
  fluid_io_rhouy[5] = rhouy[5]; 
  fluid_io_rhouz[5] = rhouz[5]; 
  fluid_io_Sxx[5] = Pxx[5] + rhouxux[5]; 
  fluid_io_Sxy[5] = Pxy[5] + rhouxuy[5]; 
  fluid_io_Sxz[5] = Pxz[5] + rhouxuz[5]; 
  fluid_io_Syy[5] = Pyy[5] + rhouyuy[5]; 
  fluid_io_Syz[5] = Pyz[5] + rhouyuz[5]; 
  fluid_io_Szz[5] = Pzz[5] + rhouzuz[5]; 
 
  pkpm_vars_io_ux[5] = ux[5]; 
  pkpm_vars_io_uy[5] = uy[5]; 
  pkpm_vars_io_uz[5] = uz[5]; 
  pkpm_vars_io_T_perp_over_m[5] = T_perp_over_m[5]; 
  pkpm_vars_io_T_perp_over_m_inv[5] = T_perp_over_m_inv[5]; 
  pkpm_vars_io_pkpm_div_ppar[5] = pkpm_div_ppar[5]; 
  pkpm_vars_io_p_perp_div_b[5] = p_perp_div_b[5]; 
  pkpm_vars_io_bb_grad_u[5] = bb_grad_u[5]; 
 
  fluid_io_rho[6] = rho[6]; 
  fluid_io_rhoux[6] = rhoux[6]; 
  fluid_io_rhouy[6] = rhouy[6]; 
  fluid_io_rhouz[6] = rhouz[6]; 
  fluid_io_Sxx[6] = Pxx[6] + rhouxux[6]; 
  fluid_io_Sxy[6] = Pxy[6] + rhouxuy[6]; 
  fluid_io_Sxz[6] = Pxz[6] + rhouxuz[6]; 
  fluid_io_Syy[6] = Pyy[6] + rhouyuy[6]; 
  fluid_io_Syz[6] = Pyz[6] + rhouyuz[6]; 
  fluid_io_Szz[6] = Pzz[6] + rhouzuz[6]; 
 
  pkpm_vars_io_ux[6] = ux[6]; 
  pkpm_vars_io_uy[6] = uy[6]; 
  pkpm_vars_io_uz[6] = uz[6]; 
  pkpm_vars_io_T_perp_over_m[6] = T_perp_over_m[6]; 
  pkpm_vars_io_T_perp_over_m_inv[6] = T_perp_over_m_inv[6]; 
  pkpm_vars_io_pkpm_div_ppar[6] = pkpm_div_ppar[6]; 
  pkpm_vars_io_p_perp_div_b[6] = p_perp_div_b[6]; 
  pkpm_vars_io_bb_grad_u[6] = bb_grad_u[6]; 
 
  fluid_io_rho[7] = rho[7]; 
  fluid_io_rhoux[7] = rhoux[7]; 
  fluid_io_rhouy[7] = rhouy[7]; 
  fluid_io_rhouz[7] = rhouz[7]; 
  fluid_io_Sxx[7] = Pxx[7] + rhouxux[7]; 
  fluid_io_Sxy[7] = Pxy[7] + rhouxuy[7]; 
  fluid_io_Sxz[7] = Pxz[7] + rhouxuz[7]; 
  fluid_io_Syy[7] = Pyy[7] + rhouyuy[7]; 
  fluid_io_Syz[7] = Pyz[7] + rhouyuz[7]; 
  fluid_io_Szz[7] = Pzz[7] + rhouzuz[7]; 
 
  pkpm_vars_io_ux[7] = ux[7]; 
  pkpm_vars_io_uy[7] = uy[7]; 
  pkpm_vars_io_uz[7] = uz[7]; 
  pkpm_vars_io_T_perp_over_m[7] = T_perp_over_m[7]; 
  pkpm_vars_io_T_perp_over_m_inv[7] = T_perp_over_m_inv[7]; 
  pkpm_vars_io_pkpm_div_ppar[7] = pkpm_div_ppar[7]; 
  pkpm_vars_io_p_perp_div_b[7] = p_perp_div_b[7]; 
  pkpm_vars_io_bb_grad_u[7] = bb_grad_u[7]; 
 
  fluid_io_rho[8] = rho[8]; 
  fluid_io_rhoux[8] = rhoux[8]; 
  fluid_io_rhouy[8] = rhouy[8]; 
  fluid_io_rhouz[8] = rhouz[8]; 
  fluid_io_Sxx[8] = Pxx[8] + rhouxux[8]; 
  fluid_io_Sxy[8] = Pxy[8] + rhouxuy[8]; 
  fluid_io_Sxz[8] = Pxz[8] + rhouxuz[8]; 
  fluid_io_Syy[8] = Pyy[8] + rhouyuy[8]; 
  fluid_io_Syz[8] = Pyz[8] + rhouyuz[8]; 
  fluid_io_Szz[8] = Pzz[8] + rhouzuz[8]; 
 
  pkpm_vars_io_ux[8] = ux[8]; 
  pkpm_vars_io_uy[8] = uy[8]; 
  pkpm_vars_io_uz[8] = uz[8]; 
  pkpm_vars_io_T_perp_over_m[8] = T_perp_over_m[8]; 
  pkpm_vars_io_T_perp_over_m_inv[8] = T_perp_over_m_inv[8]; 
  pkpm_vars_io_pkpm_div_ppar[8] = pkpm_div_ppar[8]; 
  pkpm_vars_io_p_perp_div_b[8] = p_perp_div_b[8]; 
  pkpm_vars_io_bb_grad_u[8] = bb_grad_u[8]; 
 
} 
