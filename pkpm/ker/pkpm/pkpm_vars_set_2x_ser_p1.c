#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double *p_ij, const double *pkpm_div_ppar) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       [rho ux, rho uy, rho uz], Fluid input state vector.
  // p_ij:             p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.
  // pkpm_div_ppar:    div(p_par b) computed from kinetic equation for consistency.

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_pkpm_div_ppar = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_T_perp_over_m = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_T_perp_over_m_inv = gkyl_nmat_get(rhs, count+5); 
  struct gkyl_mat rhs_Txx = gkyl_nmat_get(rhs, count+6); 
  struct gkyl_mat rhs_Tyy = gkyl_nmat_get(rhs, count+7); 
  struct gkyl_mat rhs_Tzz = gkyl_nmat_get(rhs, count+8); 
  // Clear rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&rhs_uz, 0.0); 
  gkyl_mat_clear(&rhs_pkpm_div_ppar, 0.0); 
  gkyl_mat_clear(&rhs_T_perp_over_m, 0.0); 
  gkyl_mat_clear(&rhs_T_perp_over_m_inv, 0.0); 
  gkyl_mat_clear(&rhs_Txx, 0.0); 
  gkyl_mat_clear(&rhs_Tyy, 0.0); 
  gkyl_mat_clear(&rhs_Tzz, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[4]; 
  const double *rhouz = &euler_pkpm[8]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_par = &vlasov_pkpm_moms[4]; 
  const double *p_perp = &vlasov_pkpm_moms[8]; 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pzz = &p_ij[20]; 

  int cell_avg = 0;
  // Check if rho, p_par, or p_perp < 0 at control points. 
  // *THIS IS ONLY A CHECK RIGHT NOW AND UNUSED* 
  if (1.5*rho[3]-0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (1.5*p_par[3]-0.8660254037844386*p_par[2]-0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if (1.5*p_perp[3]-0.8660254037844386*p_perp[2]-0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if ((-1.5*rho[3])-0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.5*p_par[3])-0.8660254037844386*p_par[2]+0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if ((-1.5*p_perp[3])-0.8660254037844386*p_perp[2]+0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if ((-1.5*rho[3])+0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.5*p_par[3])+0.8660254037844386*p_par[2]-0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if ((-1.5*p_perp[3])+0.8660254037844386*p_perp[2]-0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if (1.5*rho[3]+0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (1.5*p_par[3]+0.8660254037844386*p_par[2]+0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if (1.5*p_perp[3]+0.8660254037844386*p_perp[2]+0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  double rho_inv[4] = {0.0}; 
  double p_perp_inv[4] = {0.0}; 
  ser_2x_p1_inv(rho, rho_inv); 
  ser_2x_p1_inv(p_perp, p_perp_inv); 
  // Calculate expansions of primitive variables, which can be calculated free of aliasing errors. 
  double ux[4] = {0.0}; 
  double uy[4] = {0.0}; 
  double uz[4] = {0.0}; 
  double p_force[4] = {0.0}; 
  double T_perp_over_m[4] = {0.0}; 
  double T_perp_over_m_inv[4] = {0.0}; 
  double Txx[4] = {0.0}; 
  double Tyy[4] = {0.0}; 
  double Tzz[4] = {0.0}; 
 
  binop_mul_2d_ser_p1(rho_inv, rhoux, ux); 
  binop_mul_2d_ser_p1(rho_inv, rhouy, uy); 
  binop_mul_2d_ser_p1(rho_inv, rhouz, uz); 
  binop_mul_2d_ser_p1(rho_inv, pkpm_div_ppar, p_force); 
  binop_mul_2d_ser_p1(rho_inv, p_perp, T_perp_over_m); 
  binop_mul_2d_ser_p1(p_perp_inv, rho, T_perp_over_m_inv); 
  binop_mul_2d_ser_p1(rho_inv, Pxx, Txx); 
  binop_mul_2d_ser_p1(rho_inv, Pyy, Tyy); 
  binop_mul_2d_ser_p1(rho_inv, Pzz, Tzz); 
 
    gkyl_mat_set(&rhs_ux,0,0,ux[0]); 
    gkyl_mat_set(&rhs_uy,0,0,uy[0]); 
    gkyl_mat_set(&rhs_uz,0,0,uz[0]); 
    gkyl_mat_set(&rhs_pkpm_div_ppar,0,0,p_force[0]); 
    gkyl_mat_set(&rhs_T_perp_over_m,0,0,T_perp_over_m[0]); 
    gkyl_mat_set(&rhs_T_perp_over_m_inv,0,0,T_perp_over_m_inv[0]); 
    gkyl_mat_set(&rhs_Txx,0,0,Txx[0]); 
    gkyl_mat_set(&rhs_Tyy,0,0,Tyy[0]); 
    gkyl_mat_set(&rhs_Tzz,0,0,Tzz[0]); 
    gkyl_mat_set(&rhs_ux,1,0,ux[1]); 
    gkyl_mat_set(&rhs_uy,1,0,uy[1]); 
    gkyl_mat_set(&rhs_uz,1,0,uz[1]); 
    gkyl_mat_set(&rhs_pkpm_div_ppar,1,0,p_force[1]); 
    gkyl_mat_set(&rhs_T_perp_over_m,1,0,T_perp_over_m[1]); 
    gkyl_mat_set(&rhs_T_perp_over_m_inv,1,0,T_perp_over_m_inv[1]); 
    gkyl_mat_set(&rhs_Txx,1,0,Txx[1]); 
    gkyl_mat_set(&rhs_Tyy,1,0,Tyy[1]); 
    gkyl_mat_set(&rhs_Tzz,1,0,Tzz[1]); 
    gkyl_mat_set(&rhs_ux,2,0,ux[2]); 
    gkyl_mat_set(&rhs_uy,2,0,uy[2]); 
    gkyl_mat_set(&rhs_uz,2,0,uz[2]); 
    gkyl_mat_set(&rhs_pkpm_div_ppar,2,0,p_force[2]); 
    gkyl_mat_set(&rhs_T_perp_over_m,2,0,T_perp_over_m[2]); 
    gkyl_mat_set(&rhs_T_perp_over_m_inv,2,0,T_perp_over_m_inv[2]); 
    gkyl_mat_set(&rhs_Txx,2,0,Txx[2]); 
    gkyl_mat_set(&rhs_Tyy,2,0,Tyy[2]); 
    gkyl_mat_set(&rhs_Tzz,2,0,Tzz[2]); 
    gkyl_mat_set(&rhs_ux,3,0,ux[3]); 
    gkyl_mat_set(&rhs_uy,3,0,uy[3]); 
    gkyl_mat_set(&rhs_uz,3,0,uz[3]); 
    gkyl_mat_set(&rhs_pkpm_div_ppar,3,0,p_force[3]); 
    gkyl_mat_set(&rhs_T_perp_over_m,3,0,T_perp_over_m[3]); 
    gkyl_mat_set(&rhs_T_perp_over_m_inv,3,0,T_perp_over_m_inv[3]); 
    gkyl_mat_set(&rhs_Txx,3,0,Txx[3]); 
    gkyl_mat_set(&rhs_Tyy,3,0,Tyy[3]); 
    gkyl_mat_set(&rhs_Tzz,3,0,Tzz[3]); 
 
  return cell_avg;
} 
