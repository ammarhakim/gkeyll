#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH int pkpm_vars_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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
  const double *rhouy = &euler_pkpm[9]; 
  const double *rhouz = &euler_pkpm[18]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_par = &vlasov_pkpm_moms[9]; 
  const double *p_perp = &vlasov_pkpm_moms[18]; 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[27]; 
  const double *Pzz = &p_ij[45]; 

  int cell_avg = 0;
  // Check if rho, p_par, or p_perp < 0 at control points. 
  // *THIS IS ONLY A CHECK RIGHT NOW AND UNUSED* 
  if (2.5*rho[8]-1.936491673103709*rho[7]-1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]+1.5*rho[3]-0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*p_par[8]-1.936491673103709*p_par[7]-1.936491673103709*p_par[6]+1.118033988749895*p_par[5]+1.118033988749895*p_par[4]+1.5*p_par[3]-0.8660254037844386*p_par[2]-0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if (2.5*p_perp[8]-1.936491673103709*p_perp[7]-1.936491673103709*p_perp[6]+1.118033988749895*p_perp[5]+1.118033988749895*p_perp[4]+1.5*p_perp[3]-0.8660254037844386*p_perp[2]-0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])+0.9682458365518543*rho[6]+1.118033988749895*rho[5]-0.5590169943749475*rho[4]-0.8660254037844386*rho[2]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_par[8])+0.9682458365518543*p_par[6]+1.118033988749895*p_par[5]-0.5590169943749475*p_par[4]-0.8660254037844386*p_par[2]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_perp[8])+0.9682458365518543*p_perp[6]+1.118033988749895*p_perp[5]-0.5590169943749475*p_perp[4]-0.8660254037844386*p_perp[2]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]+1.936491673103709*rho[7]-1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]-1.5*rho[3]-0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*p_par[8]+1.936491673103709*p_par[7]-1.936491673103709*p_par[6]+1.118033988749895*p_par[5]+1.118033988749895*p_par[4]-1.5*p_par[3]-0.8660254037844386*p_par[2]+0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if (2.5*p_perp[8]+1.936491673103709*p_perp[7]-1.936491673103709*p_perp[6]+1.118033988749895*p_perp[5]+1.118033988749895*p_perp[4]-1.5*p_perp[3]-0.8660254037844386*p_perp[2]+0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])+0.9682458365518543*rho[7]-0.5590169943749475*rho[5]+1.118033988749895*rho[4]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_par[8])+0.9682458365518543*p_par[7]-0.5590169943749475*p_par[5]+1.118033988749895*p_par[4]-0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_perp[8])+0.9682458365518543*p_perp[7]-0.5590169943749475*p_perp[5]+1.118033988749895*p_perp[4]-0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if (0.625*rho[8]-0.5590169943749475*rho[5]-0.5590169943749475*rho[4]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (0.625*p_par[8]-0.5590169943749475*p_par[5]-0.5590169943749475*p_par[4]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if (0.625*p_perp[8]-0.5590169943749475*p_perp[5]-0.5590169943749475*p_perp[4]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])-0.9682458365518543*rho[7]-0.5590169943749475*rho[5]+1.118033988749895*rho[4]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_par[8])-0.9682458365518543*p_par[7]-0.5590169943749475*p_par[5]+1.118033988749895*p_par[4]+0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_perp[8])-0.9682458365518543*p_perp[7]-0.5590169943749475*p_perp[5]+1.118033988749895*p_perp[4]+0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]-1.936491673103709*rho[7]+1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]-1.5*rho[3]+0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*p_par[8]-1.936491673103709*p_par[7]+1.936491673103709*p_par[6]+1.118033988749895*p_par[5]+1.118033988749895*p_par[4]-1.5*p_par[3]+0.8660254037844386*p_par[2]-0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if (2.5*p_perp[8]-1.936491673103709*p_perp[7]+1.936491673103709*p_perp[6]+1.118033988749895*p_perp[5]+1.118033988749895*p_perp[4]-1.5*p_perp[3]+0.8660254037844386*p_perp[2]-0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if ((-1.25*rho[8])-0.9682458365518543*rho[6]+1.118033988749895*rho[5]-0.5590169943749475*rho[4]+0.8660254037844386*rho[2]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_par[8])-0.9682458365518543*p_par[6]+1.118033988749895*p_par[5]-0.5590169943749475*p_par[4]+0.8660254037844386*p_par[2]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if ((-1.25*p_perp[8])-0.9682458365518543*p_perp[6]+1.118033988749895*p_perp[5]-0.5590169943749475*p_perp[4]+0.8660254037844386*p_perp[2]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]+1.936491673103709*rho[7]+1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]+1.5*rho[3]+0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*p_par[8]+1.936491673103709*p_par[7]+1.936491673103709*p_par[6]+1.118033988749895*p_par[5]+1.118033988749895*p_par[4]+1.5*p_par[3]+0.8660254037844386*p_par[2]+0.8660254037844386*p_par[1]+0.5*p_par[0] < 0.0) cell_avg = 1; 
  if (2.5*p_perp[8]+1.936491673103709*p_perp[7]+1.936491673103709*p_perp[6]+1.118033988749895*p_perp[5]+1.118033988749895*p_perp[4]+1.5*p_perp[3]+0.8660254037844386*p_perp[2]+0.8660254037844386*p_perp[1]+0.5*p_perp[0] < 0.0) cell_avg = 1; 
 
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
  gkyl_mat_set(&rhs_ux,3,0,rhoux[3]); 
  gkyl_mat_set(&rhs_uy,3,0,rhouy[3]); 
  gkyl_mat_set(&rhs_uz,3,0,rhouz[3]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,3,0,pkpm_div_ppar[3]); 
  gkyl_mat_set(&rhs_T_perp_over_m,3,0,p_perp[3]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,3,0,rho[3]); 
  gkyl_mat_set(&rhs_Txx,3,0,Pxx[3]); 
  gkyl_mat_set(&rhs_Tyy,3,0,Pyy[3]); 
  gkyl_mat_set(&rhs_Tzz,3,0,Pzz[3]); 
  gkyl_mat_set(&rhs_ux,4,0,rhoux[4]); 
  gkyl_mat_set(&rhs_uy,4,0,rhouy[4]); 
  gkyl_mat_set(&rhs_uz,4,0,rhouz[4]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,4,0,pkpm_div_ppar[4]); 
  gkyl_mat_set(&rhs_T_perp_over_m,4,0,p_perp[4]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,4,0,rho[4]); 
  gkyl_mat_set(&rhs_Txx,4,0,Pxx[4]); 
  gkyl_mat_set(&rhs_Tyy,4,0,Pyy[4]); 
  gkyl_mat_set(&rhs_Tzz,4,0,Pzz[4]); 
  gkyl_mat_set(&rhs_ux,5,0,rhoux[5]); 
  gkyl_mat_set(&rhs_uy,5,0,rhouy[5]); 
  gkyl_mat_set(&rhs_uz,5,0,rhouz[5]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,5,0,pkpm_div_ppar[5]); 
  gkyl_mat_set(&rhs_T_perp_over_m,5,0,p_perp[5]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,5,0,rho[5]); 
  gkyl_mat_set(&rhs_Txx,5,0,Pxx[5]); 
  gkyl_mat_set(&rhs_Tyy,5,0,Pyy[5]); 
  gkyl_mat_set(&rhs_Tzz,5,0,Pzz[5]); 
  gkyl_mat_set(&rhs_ux,6,0,rhoux[6]); 
  gkyl_mat_set(&rhs_uy,6,0,rhouy[6]); 
  gkyl_mat_set(&rhs_uz,6,0,rhouz[6]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,6,0,pkpm_div_ppar[6]); 
  gkyl_mat_set(&rhs_T_perp_over_m,6,0,p_perp[6]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,6,0,rho[6]); 
  gkyl_mat_set(&rhs_Txx,6,0,Pxx[6]); 
  gkyl_mat_set(&rhs_Tyy,6,0,Pyy[6]); 
  gkyl_mat_set(&rhs_Tzz,6,0,Pzz[6]); 
  gkyl_mat_set(&rhs_ux,7,0,rhoux[7]); 
  gkyl_mat_set(&rhs_uy,7,0,rhouy[7]); 
  gkyl_mat_set(&rhs_uz,7,0,rhouz[7]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,7,0,pkpm_div_ppar[7]); 
  gkyl_mat_set(&rhs_T_perp_over_m,7,0,p_perp[7]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,7,0,rho[7]); 
  gkyl_mat_set(&rhs_Txx,7,0,Pxx[7]); 
  gkyl_mat_set(&rhs_Tyy,7,0,Pyy[7]); 
  gkyl_mat_set(&rhs_Tzz,7,0,Pzz[7]); 
  gkyl_mat_set(&rhs_ux,8,0,rhoux[8]); 
  gkyl_mat_set(&rhs_uy,8,0,rhouy[8]); 
  gkyl_mat_set(&rhs_uz,8,0,rhouz[8]); 
  gkyl_mat_set(&rhs_pkpm_div_ppar,8,0,pkpm_div_ppar[8]); 
  gkyl_mat_set(&rhs_T_perp_over_m,8,0,p_perp[8]); 
  gkyl_mat_set(&rhs_T_perp_over_m_inv,8,0,rho[8]); 
  gkyl_mat_set(&rhs_Txx,8,0,Pxx[8]); 
  gkyl_mat_set(&rhs_Tyy,8,0,Pyy[8]); 
  gkyl_mat_set(&rhs_Tzz,8,0,Pzz[8]); 
 
  double temp_rho = 0.0; 
  double temp_p_perp = 0.0; 
  temp_rho = 0.5*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,temp_rho); 
  gkyl_mat_set(&A_uy,0,0,temp_rho); 
  gkyl_mat_set(&A_uz,0,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,0,temp_rho); 
  gkyl_mat_set(&A_Txx,0,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,0,temp_p_perp); 
 
  temp_rho = 0.5*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,temp_rho); 
  gkyl_mat_set(&A_uy,0,1,temp_rho); 
  gkyl_mat_set(&A_uz,0,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,1,temp_rho); 
  gkyl_mat_set(&A_Txx,0,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,1,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,1,temp_p_perp); 
 
  temp_rho = 0.5*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,temp_rho); 
  gkyl_mat_set(&A_uy,0,2,temp_rho); 
  gkyl_mat_set(&A_uz,0,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,2,temp_rho); 
  gkyl_mat_set(&A_Txx,0,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,2,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,2,temp_p_perp); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,0,3,temp_rho); 
  gkyl_mat_set(&A_uy,0,3,temp_rho); 
  gkyl_mat_set(&A_uz,0,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,3,temp_rho); 
  gkyl_mat_set(&A_Txx,0,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,3,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,3,temp_p_perp); 
 
  temp_rho = 0.5*rho[4]; 
  gkyl_mat_set(&A_ux,0,4,temp_rho); 
  gkyl_mat_set(&A_uy,0,4,temp_rho); 
  gkyl_mat_set(&A_uz,0,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,4,temp_rho); 
  gkyl_mat_set(&A_Txx,0,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,4,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[4]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,4,temp_p_perp); 
 
  temp_rho = 0.5*rho[5]; 
  gkyl_mat_set(&A_ux,0,5,temp_rho); 
  gkyl_mat_set(&A_uy,0,5,temp_rho); 
  gkyl_mat_set(&A_uz,0,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,5,temp_rho); 
  gkyl_mat_set(&A_Txx,0,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,5,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[5]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,5,temp_p_perp); 
 
  temp_rho = 0.5*rho[6]; 
  gkyl_mat_set(&A_ux,0,6,temp_rho); 
  gkyl_mat_set(&A_uy,0,6,temp_rho); 
  gkyl_mat_set(&A_uz,0,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,6,temp_rho); 
  gkyl_mat_set(&A_Txx,0,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,6,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,6,temp_p_perp); 
 
  temp_rho = 0.5*rho[7]; 
  gkyl_mat_set(&A_ux,0,7,temp_rho); 
  gkyl_mat_set(&A_uy,0,7,temp_rho); 
  gkyl_mat_set(&A_uz,0,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,7,temp_rho); 
  gkyl_mat_set(&A_Txx,0,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,7,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,7,temp_p_perp); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,0,8,temp_rho); 
  gkyl_mat_set(&A_uy,0,8,temp_rho); 
  gkyl_mat_set(&A_uz,0,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,0,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,0,8,temp_rho); 
  gkyl_mat_set(&A_Txx,0,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,0,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,0,8,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[8]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,0,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,temp_rho); 
  gkyl_mat_set(&A_uy,1,0,temp_rho); 
  gkyl_mat_set(&A_uz,1,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,0,temp_rho); 
  gkyl_mat_set(&A_Txx,1,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,0,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,temp_rho); 
  gkyl_mat_set(&A_uy,1,1,temp_rho); 
  gkyl_mat_set(&A_uz,1,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,1,temp_rho); 
  gkyl_mat_set(&A_Txx,1,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,1,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[4]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,1,temp_p_perp); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,1,2,temp_rho); 
  gkyl_mat_set(&A_uy,1,2,temp_rho); 
  gkyl_mat_set(&A_uz,1,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,2,temp_rho); 
  gkyl_mat_set(&A_Txx,1,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,2,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,2,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[6]+0.5*rho[2]; 
  gkyl_mat_set(&A_ux,1,3,temp_rho); 
  gkyl_mat_set(&A_uy,1,3,temp_rho); 
  gkyl_mat_set(&A_uz,1,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,3,temp_rho); 
  gkyl_mat_set(&A_Txx,1,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,3,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[6]+0.5*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,3,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[1]; 
  gkyl_mat_set(&A_ux,1,4,temp_rho); 
  gkyl_mat_set(&A_uy,1,4,temp_rho); 
  gkyl_mat_set(&A_uz,1,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,4,temp_rho); 
  gkyl_mat_set(&A_Txx,1,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,4,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,4,temp_p_perp); 
 
  temp_rho = 0.5000000000000001*rho[7]; 
  gkyl_mat_set(&A_ux,1,5,temp_rho); 
  gkyl_mat_set(&A_uy,1,5,temp_rho); 
  gkyl_mat_set(&A_uz,1,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,5,temp_rho); 
  gkyl_mat_set(&A_Txx,1,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,5,temp_rho); 
 
  temp_p_perp = 0.5000000000000001*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,5,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,1,6,temp_rho); 
  gkyl_mat_set(&A_uy,1,6,temp_rho); 
  gkyl_mat_set(&A_uz,1,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,6,temp_rho); 
  gkyl_mat_set(&A_Txx,1,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,6,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,6,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[5]; 
  gkyl_mat_set(&A_ux,1,7,temp_rho); 
  gkyl_mat_set(&A_uy,1,7,temp_rho); 
  gkyl_mat_set(&A_uz,1,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,7,temp_rho); 
  gkyl_mat_set(&A_Txx,1,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,7,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[8]+0.5000000000000001*p_perp[5]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,7,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[7]; 
  gkyl_mat_set(&A_ux,1,8,temp_rho); 
  gkyl_mat_set(&A_uy,1,8,temp_rho); 
  gkyl_mat_set(&A_uz,1,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,1,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,1,8,temp_rho); 
  gkyl_mat_set(&A_Txx,1,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,1,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,1,8,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,1,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,temp_rho); 
  gkyl_mat_set(&A_uy,2,0,temp_rho); 
  gkyl_mat_set(&A_uz,2,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,0,temp_rho); 
  gkyl_mat_set(&A_Txx,2,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,0,temp_p_perp); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,2,1,temp_rho); 
  gkyl_mat_set(&A_uy,2,1,temp_rho); 
  gkyl_mat_set(&A_uz,2,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,1,temp_rho); 
  gkyl_mat_set(&A_Txx,2,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,1,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,1,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[5]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,temp_rho); 
  gkyl_mat_set(&A_uy,2,2,temp_rho); 
  gkyl_mat_set(&A_uz,2,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,2,temp_rho); 
  gkyl_mat_set(&A_Txx,2,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,2,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[5]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,2,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[7]+0.5*rho[1]; 
  gkyl_mat_set(&A_ux,2,3,temp_rho); 
  gkyl_mat_set(&A_uy,2,3,temp_rho); 
  gkyl_mat_set(&A_uz,2,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,3,temp_rho); 
  gkyl_mat_set(&A_Txx,2,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,3,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[7]+0.5*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,3,temp_p_perp); 
 
  temp_rho = 0.5000000000000001*rho[6]; 
  gkyl_mat_set(&A_ux,2,4,temp_rho); 
  gkyl_mat_set(&A_uy,2,4,temp_rho); 
  gkyl_mat_set(&A_uz,2,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,4,temp_rho); 
  gkyl_mat_set(&A_Txx,2,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,4,temp_rho); 
 
  temp_p_perp = 0.5000000000000001*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,4,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[2]; 
  gkyl_mat_set(&A_ux,2,5,temp_rho); 
  gkyl_mat_set(&A_uy,2,5,temp_rho); 
  gkyl_mat_set(&A_uz,2,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,5,temp_rho); 
  gkyl_mat_set(&A_Txx,2,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,5,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,5,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[4]; 
  gkyl_mat_set(&A_ux,2,6,temp_rho); 
  gkyl_mat_set(&A_uy,2,6,temp_rho); 
  gkyl_mat_set(&A_uz,2,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,6,temp_rho); 
  gkyl_mat_set(&A_Txx,2,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,6,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[8]+0.5000000000000001*p_perp[4]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,6,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,2,7,temp_rho); 
  gkyl_mat_set(&A_uy,2,7,temp_rho); 
  gkyl_mat_set(&A_uz,2,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,7,temp_rho); 
  gkyl_mat_set(&A_Txx,2,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,7,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,7,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[6]; 
  gkyl_mat_set(&A_ux,2,8,temp_rho); 
  gkyl_mat_set(&A_uy,2,8,temp_rho); 
  gkyl_mat_set(&A_uz,2,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,2,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,2,8,temp_rho); 
  gkyl_mat_set(&A_Txx,2,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,2,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,2,8,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,2,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,3,0,temp_rho); 
  gkyl_mat_set(&A_uy,3,0,temp_rho); 
  gkyl_mat_set(&A_uz,3,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,0,temp_rho); 
  gkyl_mat_set(&A_Txx,3,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,0,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[6]+0.5*rho[2]; 
  gkyl_mat_set(&A_ux,3,1,temp_rho); 
  gkyl_mat_set(&A_uy,3,1,temp_rho); 
  gkyl_mat_set(&A_uz,3,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,1,temp_rho); 
  gkyl_mat_set(&A_Txx,3,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,1,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[6]+0.5*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,1,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[7]+0.5*rho[1]; 
  gkyl_mat_set(&A_ux,3,2,temp_rho); 
  gkyl_mat_set(&A_uy,3,2,temp_rho); 
  gkyl_mat_set(&A_uz,3,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,2,temp_rho); 
  gkyl_mat_set(&A_Txx,3,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,2,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[7]+0.5*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,2,temp_p_perp); 
 
  temp_rho = 0.4*rho[8]+0.4472135954999579*rho[5]+0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,3,3,temp_rho); 
  gkyl_mat_set(&A_uy,3,3,temp_rho); 
  gkyl_mat_set(&A_uz,3,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,3,temp_rho); 
  gkyl_mat_set(&A_Txx,3,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,3,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[8]+0.4472135954999579*p_perp[5]+0.4472135954999579*p_perp[4]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,3,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,3,4,temp_rho); 
  gkyl_mat_set(&A_uy,3,4,temp_rho); 
  gkyl_mat_set(&A_uz,3,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,4,temp_rho); 
  gkyl_mat_set(&A_Txx,3,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,4,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,4,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,3,5,temp_rho); 
  gkyl_mat_set(&A_uy,3,5,temp_rho); 
  gkyl_mat_set(&A_uz,3,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,5,temp_rho); 
  gkyl_mat_set(&A_Txx,3,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,5,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,5,temp_p_perp); 
 
  temp_rho = 0.4*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,3,6,temp_rho); 
  gkyl_mat_set(&A_uy,3,6,temp_rho); 
  gkyl_mat_set(&A_uz,3,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,6,temp_rho); 
  gkyl_mat_set(&A_Txx,3,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,6,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[7]+0.447213595499958*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,6,temp_p_perp); 
 
  temp_rho = 0.4*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,3,7,temp_rho); 
  gkyl_mat_set(&A_uy,3,7,temp_rho); 
  gkyl_mat_set(&A_uz,3,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,7,temp_rho); 
  gkyl_mat_set(&A_Txx,3,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,7,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[6]+0.447213595499958*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,7,temp_p_perp); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,3,8,temp_rho); 
  gkyl_mat_set(&A_uy,3,8,temp_rho); 
  gkyl_mat_set(&A_uz,3,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,3,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,3,8,temp_rho); 
  gkyl_mat_set(&A_Txx,3,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,3,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,3,8,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,3,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[4]; 
  gkyl_mat_set(&A_ux,4,0,temp_rho); 
  gkyl_mat_set(&A_uy,4,0,temp_rho); 
  gkyl_mat_set(&A_uz,4,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,0,temp_rho); 
  gkyl_mat_set(&A_Txx,4,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[4]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,0,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[1]; 
  gkyl_mat_set(&A_ux,4,1,temp_rho); 
  gkyl_mat_set(&A_uy,4,1,temp_rho); 
  gkyl_mat_set(&A_uz,4,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,1,temp_rho); 
  gkyl_mat_set(&A_Txx,4,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,1,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,1,temp_p_perp); 
 
  temp_rho = 0.5000000000000001*rho[6]; 
  gkyl_mat_set(&A_ux,4,2,temp_rho); 
  gkyl_mat_set(&A_uy,4,2,temp_rho); 
  gkyl_mat_set(&A_uz,4,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,2,temp_rho); 
  gkyl_mat_set(&A_Txx,4,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,2,temp_rho); 
 
  temp_p_perp = 0.5000000000000001*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,2,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,4,3,temp_rho); 
  gkyl_mat_set(&A_uy,4,3,temp_rho); 
  gkyl_mat_set(&A_uz,4,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,3,temp_rho); 
  gkyl_mat_set(&A_Txx,4,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,3,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,3,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,4,4,temp_rho); 
  gkyl_mat_set(&A_uy,4,4,temp_rho); 
  gkyl_mat_set(&A_uz,4,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,4,temp_rho); 
  gkyl_mat_set(&A_Txx,4,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,4,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[4]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,4,temp_p_perp); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,4,5,temp_rho); 
  gkyl_mat_set(&A_uy,4,5,temp_rho); 
  gkyl_mat_set(&A_uz,4,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,5,temp_rho); 
  gkyl_mat_set(&A_Txx,4,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,5,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[8]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,5,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[6]+0.5000000000000001*rho[2]; 
  gkyl_mat_set(&A_ux,4,6,temp_rho); 
  gkyl_mat_set(&A_uy,4,6,temp_rho); 
  gkyl_mat_set(&A_uz,4,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,6,temp_rho); 
  gkyl_mat_set(&A_Txx,4,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,6,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[6]+0.5000000000000001*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,6,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[7]; 
  gkyl_mat_set(&A_ux,4,7,temp_rho); 
  gkyl_mat_set(&A_uy,4,7,temp_rho); 
  gkyl_mat_set(&A_uz,4,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,7,temp_rho); 
  gkyl_mat_set(&A_Txx,4,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,7,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,7,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[5]; 
  gkyl_mat_set(&A_ux,4,8,temp_rho); 
  gkyl_mat_set(&A_uy,4,8,temp_rho); 
  gkyl_mat_set(&A_uz,4,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,4,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,4,8,temp_rho); 
  gkyl_mat_set(&A_Txx,4,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,4,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,4,8,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[8]+0.5*p_perp[5]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,4,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[5]; 
  gkyl_mat_set(&A_ux,5,0,temp_rho); 
  gkyl_mat_set(&A_uy,5,0,temp_rho); 
  gkyl_mat_set(&A_uz,5,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,0,temp_rho); 
  gkyl_mat_set(&A_Txx,5,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[5]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,0,temp_p_perp); 
 
  temp_rho = 0.5000000000000001*rho[7]; 
  gkyl_mat_set(&A_ux,5,1,temp_rho); 
  gkyl_mat_set(&A_uy,5,1,temp_rho); 
  gkyl_mat_set(&A_uz,5,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,1,temp_rho); 
  gkyl_mat_set(&A_Txx,5,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,1,temp_rho); 
 
  temp_p_perp = 0.5000000000000001*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,1,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[2]; 
  gkyl_mat_set(&A_ux,5,2,temp_rho); 
  gkyl_mat_set(&A_uy,5,2,temp_rho); 
  gkyl_mat_set(&A_uz,5,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,2,temp_rho); 
  gkyl_mat_set(&A_Txx,5,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,2,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,2,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[3]; 
  gkyl_mat_set(&A_ux,5,3,temp_rho); 
  gkyl_mat_set(&A_uy,5,3,temp_rho); 
  gkyl_mat_set(&A_uz,5,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,3,temp_rho); 
  gkyl_mat_set(&A_Txx,5,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,3,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,3,temp_p_perp); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,5,4,temp_rho); 
  gkyl_mat_set(&A_uy,5,4,temp_rho); 
  gkyl_mat_set(&A_uz,5,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,4,temp_rho); 
  gkyl_mat_set(&A_Txx,5,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,4,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[8]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,4,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[5]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,5,5,temp_rho); 
  gkyl_mat_set(&A_uy,5,5,temp_rho); 
  gkyl_mat_set(&A_uz,5,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,5,temp_rho); 
  gkyl_mat_set(&A_Txx,5,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,5,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[5]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,5,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[6]; 
  gkyl_mat_set(&A_ux,5,6,temp_rho); 
  gkyl_mat_set(&A_uy,5,6,temp_rho); 
  gkyl_mat_set(&A_uz,5,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,6,temp_rho); 
  gkyl_mat_set(&A_Txx,5,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,6,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,6,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[7]+0.5000000000000001*rho[1]; 
  gkyl_mat_set(&A_ux,5,7,temp_rho); 
  gkyl_mat_set(&A_uy,5,7,temp_rho); 
  gkyl_mat_set(&A_uz,5,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,7,temp_rho); 
  gkyl_mat_set(&A_Txx,5,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,7,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[7]+0.5000000000000001*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,7,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[4]; 
  gkyl_mat_set(&A_ux,5,8,temp_rho); 
  gkyl_mat_set(&A_uy,5,8,temp_rho); 
  gkyl_mat_set(&A_uz,5,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,5,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,5,8,temp_rho); 
  gkyl_mat_set(&A_Txx,5,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,5,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,5,8,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[8]+0.5*p_perp[4]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,5,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[6]; 
  gkyl_mat_set(&A_ux,6,0,temp_rho); 
  gkyl_mat_set(&A_uy,6,0,temp_rho); 
  gkyl_mat_set(&A_uz,6,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,0,temp_rho); 
  gkyl_mat_set(&A_Txx,6,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,0,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,6,1,temp_rho); 
  gkyl_mat_set(&A_uy,6,1,temp_rho); 
  gkyl_mat_set(&A_uz,6,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,1,temp_rho); 
  gkyl_mat_set(&A_Txx,6,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,1,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,1,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[4]; 
  gkyl_mat_set(&A_ux,6,2,temp_rho); 
  gkyl_mat_set(&A_uy,6,2,temp_rho); 
  gkyl_mat_set(&A_uz,6,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,2,temp_rho); 
  gkyl_mat_set(&A_Txx,6,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,2,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[8]+0.5000000000000001*p_perp[4]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,2,temp_p_perp); 
 
  temp_rho = 0.4*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,6,3,temp_rho); 
  gkyl_mat_set(&A_uy,6,3,temp_rho); 
  gkyl_mat_set(&A_uz,6,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,3,temp_rho); 
  gkyl_mat_set(&A_Txx,6,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,3,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[7]+0.447213595499958*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,3,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[6]+0.5000000000000001*rho[2]; 
  gkyl_mat_set(&A_ux,6,4,temp_rho); 
  gkyl_mat_set(&A_uy,6,4,temp_rho); 
  gkyl_mat_set(&A_uz,6,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,4,temp_rho); 
  gkyl_mat_set(&A_Txx,6,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,4,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[6]+0.5000000000000001*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,4,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[6]; 
  gkyl_mat_set(&A_ux,6,5,temp_rho); 
  gkyl_mat_set(&A_uy,6,5,temp_rho); 
  gkyl_mat_set(&A_uz,6,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,5,temp_rho); 
  gkyl_mat_set(&A_Txx,6,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,5,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,5,temp_p_perp); 
 
  temp_rho = 0.2857142857142857*rho[8]+0.4472135954999579*rho[5]+0.31943828249997*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,6,6,temp_rho); 
  gkyl_mat_set(&A_uy,6,6,temp_rho); 
  gkyl_mat_set(&A_uz,6,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,6,temp_rho); 
  gkyl_mat_set(&A_Txx,6,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,6,temp_rho); 
 
  temp_p_perp = 0.2857142857142857*p_perp[8]+0.4472135954999579*p_perp[5]+0.31943828249997*p_perp[4]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,6,temp_p_perp); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,6,7,temp_rho); 
  gkyl_mat_set(&A_uy,6,7,temp_rho); 
  gkyl_mat_set(&A_uz,6,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,7,temp_rho); 
  gkyl_mat_set(&A_Txx,6,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,7,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,7,temp_p_perp); 
 
  temp_rho = 0.2857142857142857*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,6,8,temp_rho); 
  gkyl_mat_set(&A_uy,6,8,temp_rho); 
  gkyl_mat_set(&A_uz,6,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,6,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,6,8,temp_rho); 
  gkyl_mat_set(&A_Txx,6,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,6,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,6,8,temp_rho); 
 
  temp_p_perp = 0.2857142857142857*p_perp[6]+0.447213595499958*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,6,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[7]; 
  gkyl_mat_set(&A_ux,7,0,temp_rho); 
  gkyl_mat_set(&A_uy,7,0,temp_rho); 
  gkyl_mat_set(&A_uz,7,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,0,temp_rho); 
  gkyl_mat_set(&A_Txx,7,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,0,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[8]+0.5000000000000001*rho[5]; 
  gkyl_mat_set(&A_ux,7,1,temp_rho); 
  gkyl_mat_set(&A_uy,7,1,temp_rho); 
  gkyl_mat_set(&A_uz,7,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,1,temp_rho); 
  gkyl_mat_set(&A_Txx,7,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,1,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[8]+0.5000000000000001*p_perp[5]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,1,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[3]; 
  gkyl_mat_set(&A_ux,7,2,temp_rho); 
  gkyl_mat_set(&A_uy,7,2,temp_rho); 
  gkyl_mat_set(&A_uz,7,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,2,temp_rho); 
  gkyl_mat_set(&A_Txx,7,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,2,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,2,temp_p_perp); 
 
  temp_rho = 0.4*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,7,3,temp_rho); 
  gkyl_mat_set(&A_uy,7,3,temp_rho); 
  gkyl_mat_set(&A_uz,7,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,3,temp_rho); 
  gkyl_mat_set(&A_Txx,7,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,3,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[6]+0.447213595499958*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,3,temp_p_perp); 
 
  temp_rho = 0.4472135954999579*rho[7]; 
  gkyl_mat_set(&A_ux,7,4,temp_rho); 
  gkyl_mat_set(&A_uy,7,4,temp_rho); 
  gkyl_mat_set(&A_uz,7,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,4,temp_rho); 
  gkyl_mat_set(&A_Txx,7,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,4,temp_rho); 
 
  temp_p_perp = 0.4472135954999579*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,4,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[7]+0.5000000000000001*rho[1]; 
  gkyl_mat_set(&A_ux,7,5,temp_rho); 
  gkyl_mat_set(&A_uy,7,5,temp_rho); 
  gkyl_mat_set(&A_uz,7,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,5,temp_rho); 
  gkyl_mat_set(&A_Txx,7,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,5,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[7]+0.5000000000000001*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,5,temp_p_perp); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,7,6,temp_rho); 
  gkyl_mat_set(&A_uy,7,6,temp_rho); 
  gkyl_mat_set(&A_uz,7,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,6,temp_rho); 
  gkyl_mat_set(&A_Txx,7,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,6,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,6,temp_p_perp); 
 
  temp_rho = 0.2857142857142857*rho[8]+0.31943828249997*rho[5]+0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,7,7,temp_rho); 
  gkyl_mat_set(&A_uy,7,7,temp_rho); 
  gkyl_mat_set(&A_uz,7,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,7,temp_rho); 
  gkyl_mat_set(&A_Txx,7,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,7,temp_rho); 
 
  temp_p_perp = 0.2857142857142857*p_perp[8]+0.31943828249997*p_perp[5]+0.4472135954999579*p_perp[4]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,7,temp_p_perp); 
 
  temp_rho = 0.2857142857142857*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,7,8,temp_rho); 
  gkyl_mat_set(&A_uy,7,8,temp_rho); 
  gkyl_mat_set(&A_uz,7,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,7,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,7,8,temp_rho); 
  gkyl_mat_set(&A_Txx,7,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,7,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,7,8,temp_rho); 
 
  temp_p_perp = 0.2857142857142857*p_perp[7]+0.447213595499958*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,7,8,temp_p_perp); 
 
  temp_rho = 0.5*rho[8]; 
  gkyl_mat_set(&A_ux,8,0,temp_rho); 
  gkyl_mat_set(&A_uy,8,0,temp_rho); 
  gkyl_mat_set(&A_uz,8,0,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,0,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,0,temp_rho); 
  gkyl_mat_set(&A_Txx,8,0,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,0,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,0,temp_rho); 
 
  temp_p_perp = 0.5*p_perp[8]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,0,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[7]; 
  gkyl_mat_set(&A_ux,8,1,temp_rho); 
  gkyl_mat_set(&A_uy,8,1,temp_rho); 
  gkyl_mat_set(&A_uz,8,1,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,1,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,1,temp_rho); 
  gkyl_mat_set(&A_Txx,8,1,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,1,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,1,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[7]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,1,temp_p_perp); 
 
  temp_rho = 0.447213595499958*rho[6]; 
  gkyl_mat_set(&A_ux,8,2,temp_rho); 
  gkyl_mat_set(&A_uy,8,2,temp_rho); 
  gkyl_mat_set(&A_uz,8,2,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,2,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,2,temp_rho); 
  gkyl_mat_set(&A_Txx,8,2,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,2,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,2,temp_rho); 
 
  temp_p_perp = 0.447213595499958*p_perp[6]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,2,temp_p_perp); 
 
  temp_rho = 0.4*rho[3]; 
  gkyl_mat_set(&A_ux,8,3,temp_rho); 
  gkyl_mat_set(&A_uy,8,3,temp_rho); 
  gkyl_mat_set(&A_uz,8,3,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,3,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,3,temp_rho); 
  gkyl_mat_set(&A_Txx,8,3,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,3,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,3,temp_rho); 
 
  temp_p_perp = 0.4*p_perp[3]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,3,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[5]; 
  gkyl_mat_set(&A_ux,8,4,temp_rho); 
  gkyl_mat_set(&A_uy,8,4,temp_rho); 
  gkyl_mat_set(&A_uz,8,4,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,4,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,4,temp_rho); 
  gkyl_mat_set(&A_Txx,8,4,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,4,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,4,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[8]+0.5*p_perp[5]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,4,temp_p_perp); 
 
  temp_rho = 0.31943828249997*rho[8]+0.5*rho[4]; 
  gkyl_mat_set(&A_ux,8,5,temp_rho); 
  gkyl_mat_set(&A_uy,8,5,temp_rho); 
  gkyl_mat_set(&A_uz,8,5,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,5,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,5,temp_rho); 
  gkyl_mat_set(&A_Txx,8,5,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,5,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,5,temp_rho); 
 
  temp_p_perp = 0.31943828249997*p_perp[8]+0.5*p_perp[4]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,5,temp_p_perp); 
 
  temp_rho = 0.2857142857142857*rho[6]+0.447213595499958*rho[2]; 
  gkyl_mat_set(&A_ux,8,6,temp_rho); 
  gkyl_mat_set(&A_uy,8,6,temp_rho); 
  gkyl_mat_set(&A_uz,8,6,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,6,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,6,temp_rho); 
  gkyl_mat_set(&A_Txx,8,6,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,6,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,6,temp_rho); 
 
  temp_p_perp = 0.2857142857142857*p_perp[6]+0.447213595499958*p_perp[2]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,6,temp_p_perp); 
 
  temp_rho = 0.2857142857142857*rho[7]+0.447213595499958*rho[1]; 
  gkyl_mat_set(&A_ux,8,7,temp_rho); 
  gkyl_mat_set(&A_uy,8,7,temp_rho); 
  gkyl_mat_set(&A_uz,8,7,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,7,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,7,temp_rho); 
  gkyl_mat_set(&A_Txx,8,7,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,7,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,7,temp_rho); 
 
  temp_p_perp = 0.2857142857142857*p_perp[7]+0.447213595499958*p_perp[1]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,7,temp_p_perp); 
 
  temp_rho = 0.2040816326530612*rho[8]+0.31943828249997*rho[5]+0.31943828249997*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,8,8,temp_rho); 
  gkyl_mat_set(&A_uy,8,8,temp_rho); 
  gkyl_mat_set(&A_uz,8,8,temp_rho); 
  gkyl_mat_set(&A_pkpm_div_ppar,8,8,temp_rho); 
  gkyl_mat_set(&A_T_perp_over_m,8,8,temp_rho); 
  gkyl_mat_set(&A_Txx,8,8,temp_rho); 
  gkyl_mat_set(&A_Tyy,8,8,temp_rho); 
  gkyl_mat_set(&A_Tzz,8,8,temp_rho); 
 
  temp_p_perp = 0.2040816326530612*p_perp[8]+0.31943828249997*p_perp[5]+0.31943828249997*p_perp[4]+0.5*p_perp[0]; 
  gkyl_mat_set(&A_T_perp_over_m_inv,8,8,temp_p_perp); 
 
  return cell_avg;
} 
