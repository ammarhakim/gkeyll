#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms, double gas_gamma, double mass) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: Moments [rho, rho ux, rho uy, rho uz, totalE].
  // gas_gamma: Adiabatic index. 
  // mass: Species mass. 

  struct gkyl_mat A_ux = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_uy = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_uz = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat A_temp = gkyl_nmat_get(A, count+3); 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_temp = gkyl_nmat_get(rhs, count+3); 
  // Clear matrix and rhs for each component of flow velocity being solved for 
  gkyl_mat_clear(&A_ux, 0.0); gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&A_uy, 0.0); gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&A_uz, 0.0); gkyl_mat_clear(&rhs_uz, 0.0); 
  gkyl_mat_clear(&A_temp, 0.0); gkyl_mat_clear(&rhs_temp, 0.0); 
  const double *rho   = &moms[0]; 
  const double *rhoux = &moms[3]; 
  const double *rhouy = &moms[6]; 
  const double *rhouz = &moms[9]; 
  const double *totE  = &moms[12]; 

  double rhouxSq[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouz, rhouz, rhouzSq); 
 
  double rho_temp[3]; 
  rho_temp[0] = (gas_gamma - 1.0)*(mass * totE[0] - 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0])); 
  rho_temp[1] = (gas_gamma - 1.0)*(mass * totE[1] - 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1])); 
  rho_temp[2] = (gas_gamma - 1.0)*(mass * totE[2] - 0.5*(rhouxSq[2] + rhouySq[2] + rhouzSq[2])); 

  gkyl_mat_set(&rhs_ux,0,0,rhoux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,rhouy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,rhouz[0]); 
  gkyl_mat_set(&rhs_temp,0,0,rho_temp[0]); 
  gkyl_mat_set(&rhs_ux,1,0,rhoux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,rhouy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,rhouz[1]); 
  gkyl_mat_set(&rhs_temp,1,0,rho_temp[1]); 
  gkyl_mat_set(&rhs_ux,2,0,rhoux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,rhouy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,rhouz[2]); 
  gkyl_mat_set(&rhs_temp,2,0,rho_temp[2]); 
 
  double tmp_rho = 0.0; 
  tmp_rho = 0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,tmp_rho); 
  gkyl_mat_set(&A_uy,0,0,tmp_rho); 
  gkyl_mat_set(&A_uz,0,0,tmp_rho); 
  gkyl_mat_set(&A_temp,0,0,tmp_rho); 
 
  tmp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,tmp_rho); 
  gkyl_mat_set(&A_uy,0,1,tmp_rho); 
  gkyl_mat_set(&A_uz,0,1,tmp_rho); 
  gkyl_mat_set(&A_temp,0,1,tmp_rho); 
 
  tmp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,tmp_rho); 
  gkyl_mat_set(&A_uy,0,2,tmp_rho); 
  gkyl_mat_set(&A_uz,0,2,tmp_rho); 
  gkyl_mat_set(&A_temp,0,2,tmp_rho); 
 
  tmp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,tmp_rho); 
  gkyl_mat_set(&A_uy,1,0,tmp_rho); 
  gkyl_mat_set(&A_uz,1,0,tmp_rho); 
  gkyl_mat_set(&A_temp,1,0,tmp_rho); 
 
  tmp_rho = 0.6324555320336759*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,tmp_rho); 
  gkyl_mat_set(&A_uy,1,1,tmp_rho); 
  gkyl_mat_set(&A_uz,1,1,tmp_rho); 
  gkyl_mat_set(&A_temp,1,1,tmp_rho); 
 
  tmp_rho = 0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,1,2,tmp_rho); 
  gkyl_mat_set(&A_uy,1,2,tmp_rho); 
  gkyl_mat_set(&A_uz,1,2,tmp_rho); 
  gkyl_mat_set(&A_temp,1,2,tmp_rho); 
 
  tmp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,tmp_rho); 
  gkyl_mat_set(&A_uy,2,0,tmp_rho); 
  gkyl_mat_set(&A_uz,2,0,tmp_rho); 
  gkyl_mat_set(&A_temp,2,0,tmp_rho); 
 
  tmp_rho = 0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,2,1,tmp_rho); 
  gkyl_mat_set(&A_uy,2,1,tmp_rho); 
  gkyl_mat_set(&A_uz,2,1,tmp_rho); 
  gkyl_mat_set(&A_temp,2,1,tmp_rho); 
 
  tmp_rho = 0.45175395145262565*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,tmp_rho); 
  gkyl_mat_set(&A_uy,2,2,tmp_rho); 
  gkyl_mat_set(&A_uz,2,2,tmp_rho); 
  gkyl_mat_set(&A_temp,2,2,tmp_rho); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, 
    double* GKYL_RESTRICT out) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // xsol: Input solution vector. 
  // out: Output volume expansion of flow velocity and temperature;
 
  struct gkyl_mat x_ux = gkyl_nmat_get(xsol, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(xsol, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(xsol, count+2); 
  struct gkyl_mat x_temp = gkyl_nmat_get(xsol, count+3); 
  double *ux = &out[0]; 
  double *uy = &out[3]; 
  double *uz = &out[6]; 
  double *temp = &out[9]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  temp[0] = gkyl_mat_get(&x_temp,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  temp[1] = gkyl_mat_get(&x_temp,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  temp[2] = gkyl_mat_get(&x_temp,2,0); 

} 
 
