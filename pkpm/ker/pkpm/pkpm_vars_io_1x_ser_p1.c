#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_io_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
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
  const double *rhouy = &euler_pkpm[2]; 
  const double *rhouz = &euler_pkpm[4]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[2]; 
  const double *Pxz = &p_ij[4]; 
  const double *Pyy = &p_ij[6]; 
  const double *Pyz = &p_ij[8]; 
  const double *Pzz = &p_ij[10]; 

  const double *ux = &prim[0]; 
  const double *uy = &prim[2]; 
  const double *uz = &prim[4]; 
  const double *pkpm_div_ppar = &prim[6]; 
  const double *T_perp_over_m = &prim[8]; 
  const double *T_perp_over_m_inv = &prim[10]; 

  const double *p_perp_div_b = &pkpm_accel[0]; 
  const double *bb_grad_u = &pkpm_accel[2]; 

  double *fluid_io_rho = &fluid_io[0]; 
  double *fluid_io_rhoux = &fluid_io[2]; 
  double *fluid_io_rhouy = &fluid_io[4]; 
  double *fluid_io_rhouz = &fluid_io[6]; 
  double *fluid_io_Sxx = &fluid_io[8]; 
  double *fluid_io_Sxy = &fluid_io[10]; 
  double *fluid_io_Sxz = &fluid_io[12]; 
  double *fluid_io_Syy = &fluid_io[14]; 
  double *fluid_io_Syz = &fluid_io[16]; 
  double *fluid_io_Szz = &fluid_io[18]; 

  double *pkpm_vars_io_ux = &pkpm_vars_io[0]; 
  double *pkpm_vars_io_uy = &pkpm_vars_io[2]; 
  double *pkpm_vars_io_uz = &pkpm_vars_io[4]; 
  double *pkpm_vars_io_T_perp_over_m = &pkpm_vars_io[6]; 
  double *pkpm_vars_io_T_perp_over_m_inv = &pkpm_vars_io[8]; 
  double *pkpm_vars_io_pkpm_div_ppar = &pkpm_vars_io[10]; 
  double *pkpm_vars_io_p_perp_div_b = &pkpm_vars_io[12]; 
  double *pkpm_vars_io_bb_grad_u = &pkpm_vars_io[14]; 

  double rhouxux[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhoux, ux, rhouxux); 
 
  double rhouxuy[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhoux, uy, rhouxuy); 
 
  double rhouxuz[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhoux, uz, rhouxuz); 
 
  double rhouyuy[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhouy, uy, rhouyuy); 
 
  double rhouyuz[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhouy, uz, rhouyuz); 
 
  double rhouzuz[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhouz, uz, rhouzuz); 
 
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
 
} 
