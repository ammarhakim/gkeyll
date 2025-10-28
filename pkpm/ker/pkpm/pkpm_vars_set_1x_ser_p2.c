#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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

  struct gkyl_mat A_ux = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_uy = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_uz = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat A_pkpm_div_ppar = gkyl_nmat_get(A, count+3); 
  struct gkyl_mat A_T_perp_over_m = gkyl_nmat_get(A, count+4); 
  struct gkyl_mat A_T_perp_over_m_inv = gkyl_nmat_get(A, count+5); 
  struct gkyl_mat A_Txx = gkyl_nmat_get(A, count+6); 
  struct gkyl_mat A_Tyy = gkyl_nmat_get(A, count+7); 
  struct gkyl_mat A_Tzz = gkyl_nmat_get(A, count+8); 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_pkpm_div_ppar = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_T_perp_over_m = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_T_perp_over_m_inv = gkyl_nmat_get(rhs, count+5); 
  struct gkyl_mat rhs_Txx = gkyl_nmat_get(rhs, count+6); 
  struct gkyl_mat rhs_Tyy = gkyl_nmat_get(rhs, count+7); 
  struct gkyl_mat rhs_Tzz = gkyl_nmat_get(rhs, count+8); 
  // Clear matrix and rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&A_ux, 0.0); gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&A_uy, 0.0); gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&A_uz, 0.0); gkyl_mat_clear(&rhs_uz, 0.0); 
  gkyl_mat_clear(&A_pkpm_div_ppar, 0.0); gkyl_mat_clear(&rhs_pkpm_div_ppar, 0.0); 
  gkyl_mat_clear(&A_T_perp_over_m, 0.0); gkyl_mat_clear(&rhs_T_perp_over_m, 0.0); 
  gkyl_mat_clear(&A_T_perp_over_m_inv, 0.0); gkyl_mat_clear(&rhs_T_perp_over_m_inv, 0.0); 
  gkyl_mat_clear(&A_Txx, 0.0); gkyl_mat_clear(&rhs_Txx, 0.0); 
  gkyl_mat_clear(&A_Tyy, 0.0); gkyl_mat_clear(&rhs_Tyy, 0.0); 
  gkyl_mat_clear(&A_Tzz, 0.0); gkyl_mat_clear(&rhs_Tzz, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[3]; 
  const double *rhouz = &euler_pkpm[6]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_par = &vlasov_pkpm_moms[3]; 
  const double *p_perp = &vlasov_pkpm_moms[6]; 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[9]; 
  const double *Pzz = &p_ij[15]; 

  int cell_avg = 0;
  // Check if rho, p_par, or p_perp < 0 at control points. 
  // *THIS IS ONLY A CHECK RIGHT NOW AND UNUSED* 
  if (1.58113883008419*rho[2]-1.224744871391589*rho[1]+0.7071067811865475*rho[0] < 0.0) cell_avg = 1; 
  if (1.58113883008419*p_par[2]-1.224744871391589*p_par[1]+0.7071067811865475*p_par[0] < 0.0) cell_avg = 1; 
  if (1.58113883008419*p_perp[2]-1.224744871391589*p_perp[1]+0.7071067811865475*p_perp[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*rho[0]-0.7905694150420947*rho[2] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*p_par[0]-0.7905694150420947*p_par[2] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*p_perp[0]-0.7905694150420947*p_perp[2] < 0.0) cell_avg = 1; 
  if (1.58113883008419*rho[2]+1.224744871391589*rho[1]+0.7071067811865475*rho[0] < 0.0) cell_avg = 1; 
  if (1.58113883008419*p_par[2]+1.224744871391589*p_par[1]+0.7071067811865475*p_par[0] < 0.0) cell_avg = 1; 
  if (1.58113883008419*p_perp[2]+1.224744871391589*p_perp[1]+0.7071067811865475*p_perp[0] < 0.0) cell_avg = 1; 
 
  gkyl_mat_set(&rhs_ux,0,0,rhoux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,rhouy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,rhouz[0]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,0,0,pkpm_div_ppar[0]); 
  gkyl_mat_set(&rhs_T_perp_over_m,0,0,p_perp[0]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,0,0,rho[0]); 
  gkyl_mat_set(&rhs_Txx,0,0,Pxx[0]); 
  gkyl_mat_set(&rhs_Tyy,0,0,Pyy[0]); 
  gkyl_mat_set(&rhs_Tzz,0,0,Pzz[0]); 
  gkyl_mat_set(&rhs_ux,1,0,rhoux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,rhouy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,rhouz[1]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,1,0,pkpm_div_ppar[1]); 
  gkyl_mat_set(&rhs_T_perp_over_m,1,0,p_perp[1]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,1,0,rho[1]); 
  gkyl_mat_set(&rhs_Txx,1,0,Pxx[1]); 
  gkyl_mat_set(&rhs_Tyy,1,0,Pyy[1]); 
  gkyl_mat_set(&rhs_Tzz,1,0,Pzz[1]); 
  gkyl_mat_set(&rhs_ux,2,0,rhoux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,rhouy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,rhouz[2]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,2,0,pkpm_div_ppar[2]); 
  gkyl_mat_set(&rhs_T_perp_over_m,2,0,p_perp[2]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,2,0,rho[2]); 
  gkyl_mat_set(&rhs_Txx,2,0,Pxx[2]); 
  gkyl_mat_set(&rhs_Tyy,2,0,Pyy[2]); 
  gkyl_mat_set(&rhs_Tzz,2,0,Pzz[2]); 
 
  double temp_rho = 0.0; 
  double temp_p_perp = 0.0; 
  temp_rho = 0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,temp_rho); 
  gkyl_mat_set(&A_uy,0,0,temp_rho); 
  gkyl_mat_set(&A_uz,0,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,0,temp_rho); 
  gkyl_mat_set(&A_Txx,0,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,0,temp_rho); 
 
  temp_p_perp = 0.7071067811865475*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,0,temp_p_perp); 
 
  temp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,temp_rho); 
  gkyl_mat_set(&A_uy,0,1,temp_rho); 
  gkyl_mat_set(&A_uz,0,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,1,temp_rho); 
  gkyl_mat_set(&A_Txx,0,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,1,temp_rho); 
 
  temp_p_perp = 0.7071067811865475*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,1,temp_p_perp); 
 
  temp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,temp_rho); 
  gkyl_mat_set(&A_uy,0,2,temp_rho); 
  gkyl_mat_set(&A_uz,0,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,2,temp_rho); 
  gkyl_mat_set(&A_Txx,0,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,2,temp_rho); 
 
  temp_p_perp = 0.7071067811865475*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,2,temp_p_perp); 
 
  temp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,temp_rho); 
  gkyl_mat_set(&A_uy,1,0,temp_rho); 
  gkyl_mat_set(&A_uz,1,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,0,temp_rho); 
  gkyl_mat_set(&A_Txx,1,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,0,temp_rho); 
 
  temp_p_perp = 0.7071067811865475*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,0,temp_p_perp); 
 
  temp_rho = 0.6324555320336759*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,temp_rho); 
  gkyl_mat_set(&A_uy,1,1,temp_rho); 
  gkyl_mat_set(&A_uz,1,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,1,temp_rho); 
  gkyl_mat_set(&A_Txx,1,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,1,temp_rho); 
 
  temp_p_perp = 0.6324555320336759*p_perp[2]+0.7071067811865475*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,1,temp_p_perp); 
 
  temp_rho = 0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,1,2,temp_rho); 
  gkyl_mat_set(&A_uy,1,2,temp_rho); 
  gkyl_mat_set(&A_uz,1,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,2,temp_rho); 
  gkyl_mat_set(&A_Txx,1,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,2,temp_rho); 
 
  temp_p_perp = 0.6324555320336759*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,2,temp_p_perp); 
 
  temp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,temp_rho); 
  gkyl_mat_set(&A_uy,2,0,temp_rho); 
  gkyl_mat_set(&A_uz,2,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,0,temp_rho); 
  gkyl_mat_set(&A_Txx,2,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,0,temp_rho); 
 
  temp_p_perp = 0.7071067811865475*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,0,temp_p_perp); 
 
  temp_rho = 0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,2,1,temp_rho); 
  gkyl_mat_set(&A_uy,2,1,temp_rho); 
  gkyl_mat_set(&A_uz,2,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,1,temp_rho); 
  gkyl_mat_set(&A_Txx,2,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,1,temp_rho); 
 
  temp_p_perp = 0.6324555320336759*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,1,temp_p_perp); 
 
  temp_rho = 0.4517539514526256*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,temp_rho); 
  gkyl_mat_set(&A_uy,2,2,temp_rho); 
  gkyl_mat_set(&A_uz,2,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,2,temp_rho); 
  gkyl_mat_set(&A_Txx,2,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,2,temp_rho); 
 
  temp_p_perp = 0.4517539514526256*p_perp[2]+0.7071067811865475*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,2,temp_p_perp); 
 
  return cell_avg;
} 
