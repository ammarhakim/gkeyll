#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms, double gas_gamma, double mass) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: Moments [rho, rho ux, rho uy, rho uz, totalE].
  // gas_gamma: Adiabatic index. 
  // mass: Species mass. 

  struct gkyl_mat A_temp = gkyl_nmat_get(A, count); 
  struct gkyl_mat rhs_temp = gkyl_nmat_get(rhs, count); 
  // Clear matrix and rhs. 
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

  gkyl_mat_set(&rhs_temp,0,0,rho_temp[0]); 
  gkyl_mat_set(&rhs_temp,1,0,rho_temp[1]); 
  gkyl_mat_set(&rhs_temp,2,0,rho_temp[2]); 
 
  gkyl_mat_set(&A_temp,0,0,0.0); 
 
  gkyl_mat_set(&A_temp,0,1,0.0); 
 
  gkyl_mat_set(&A_temp,0,2,0.0); 
 
  gkyl_mat_set(&A_temp,1,0,0.0); 
 
  gkyl_mat_set(&A_temp,1,1,0.0); 
 
  gkyl_mat_set(&A_temp,1,2,0.0); 
 
  gkyl_mat_set(&A_temp,2,0,0.0); 
 
  gkyl_mat_set(&A_temp,2,1,0.0); 
 
  gkyl_mat_set(&A_temp,2,2,0.0); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, 
    double* GKYL_RESTRICT out) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // xsol: Input solution vector. 
  // out: Output volume expansion of temperaure. 
 
  struct gkyl_mat x_temp = gkyl_nmat_get(xsol, count); 

  out[0] = gkyl_mat_get(&x_temp,0,0); 
  out[1] = gkyl_mat_get(&x_temp,1,0); 
  out[2] = gkyl_mat_get(&x_temp,2,0); 

} 
 
