#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, 
    double* GKYL_RESTRICT p_ij) 
{ 
  // bvar:             Input volume expansion of magnetic field unit vector and tensor.
  //                   [bx, by, bz, bxbx, bxby, bxbz, byby, bybz, bzbz] 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // p_ij:             Output volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.

  // Parallel/Perp pressure are first/second component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *p_parallel = &vlasov_pkpm_moms[3]; 
  const double *p_perp = &vlasov_pkpm_moms[6]; 
  const double *bxbx = &bvar[9]; 
  const double *bxby = &bvar[12]; 
  const double *bxbz = &bvar[15]; 
  const double *byby = &bvar[18]; 
  const double *bybz = &bvar[21]; 
  const double *bzbz = &bvar[24]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[3]; 
  double *Pxz = &p_ij[6]; 
  double *Pyy = &p_ij[9]; 
  double *Pyz = &p_ij[12]; 
  double *Pzz = &p_ij[15]; 

  double DP[3] = {0.0}; 
  DP[0] = p_parallel[0] - p_perp[0]; 
  DP[1] = p_parallel[1] - p_perp[1]; 
  DP[2] = p_parallel[2] - p_perp[2]; 
  // DP b_i b_j. 
  double DP_bxbx[3] = {0.0}; 
  binop_mul_1d_ser_p2(DP, bxbx, DP_bxbx); 
 
  double DP_bxby[3] = {0.0}; 
  binop_mul_1d_ser_p2(DP, bxby, DP_bxby); 
 
  double DP_bxbz[3] = {0.0}; 
  binop_mul_1d_ser_p2(DP, bxbz, DP_bxbz); 
 
  double DP_byby[3] = {0.0}; 
  binop_mul_1d_ser_p2(DP, byby, DP_byby); 
 
  double DP_bybz[3] = {0.0}; 
  binop_mul_1d_ser_p2(DP, bybz, DP_bybz); 
 
  double DP_bzbz[3] = {0.0}; 
  binop_mul_1d_ser_p2(DP, bzbz, DP_bzbz); 
 
  Pxx[0] = DP_bxbx[0] + p_perp[0]; 
  Pxy[0] = DP_bxby[0]; 
  Pxz[0] = DP_bxbz[0]; 
  Pyy[0] = DP_byby[0] + p_perp[0]; 
  Pyz[0] = DP_bybz[0]; 
  Pzz[0] = DP_bzbz[0] + p_perp[0]; 
  Pxx[1] = DP_bxbx[1] + p_perp[1]; 
  Pxy[1] = DP_bxby[1]; 
  Pxz[1] = DP_bxbz[1]; 
  Pyy[1] = DP_byby[1] + p_perp[1]; 
  Pyz[1] = DP_bybz[1]; 
  Pzz[1] = DP_bzbz[1] + p_perp[1]; 
  Pxx[2] = DP_bxbx[2] + p_perp[2]; 
  Pxy[2] = DP_bxby[2]; 
  Pxz[2] = DP_bxbz[2]; 
  Pyy[2] = DP_byby[2] + p_perp[2]; 
  Pyz[2] = DP_bybz[2]; 
  Pzz[2] = DP_bzbz[2] + p_perp[2]; 
 
} 
