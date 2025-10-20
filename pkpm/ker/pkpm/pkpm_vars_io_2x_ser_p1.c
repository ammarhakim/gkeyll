#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_io_2x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
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
  const double *rhouy = &euler_pkpm[4]; 
  const double *rhouz = &euler_pkpm[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[4]; 
  const double *Pxz = &p_ij[8]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pyz = &p_ij[16]; 
  const double *Pzz = &p_ij[20]; 

  const double *ux = &prim[0]; 
  const double *uy = &prim[4]; 
  const double *uz = &prim[8]; 
  const double *pkpm_div_ppar = &prim[12]; 
  const double *T_perp_over_m = &prim[16]; 
  const double *T_perp_over_m_inv = &prim[20]; 

  const double *p_perp_div_b = &pkpm_accel[0]; 
  const double *bb_grad_u = &pkpm_accel[4]; 

  double *fluid_io_rho = &fluid_io[0]; 
  double *fluid_io_rhoux = &fluid_io[4]; 
  double *fluid_io_rhouy = &fluid_io[8]; 
  double *fluid_io_rhouz = &fluid_io[12]; 
  double *fluid_io_Sxx = &fluid_io[16]; 
  double *fluid_io_Sxy = &fluid_io[20]; 
  double *fluid_io_Sxz = &fluid_io[24]; 
  double *fluid_io_Syy = &fluid_io[28]; 
  double *fluid_io_Syz = &fluid_io[32]; 
  double *fluid_io_Szz = &fluid_io[36]; 

  double *pkpm_vars_io_ux = &pkpm_vars_io[0]; 
  double *pkpm_vars_io_uy = &pkpm_vars_io[4]; 
  double *pkpm_vars_io_uz = &pkpm_vars_io[8]; 
  double *pkpm_vars_io_T_perp_over_m = &pkpm_vars_io[12]; 
  double *pkpm_vars_io_T_perp_over_m_inv = &pkpm_vars_io[16]; 
  double *pkpm_vars_io_pkpm_div_ppar = &pkpm_vars_io[20]; 
  double *pkpm_vars_io_p_perp_div_b = &pkpm_vars_io[24]; 
  double *pkpm_vars_io_bb_grad_u = &pkpm_vars_io[28]; 

  double rhouxux[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhoux, ux, rhouxux); 
 
  double rhouxuy[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhoux, uy, rhouxuy); 
 
  double rhouxuz[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhoux, uz, rhouxuz); 
 
  double rhouyuy[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhouy, uy, rhouyuy); 
 
  double rhouyuz[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhouy, uz, rhouyuz); 
 
  double rhouzuz[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhouz, uz, rhouzuz); 
 
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
 
} 
