#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms, double gas_gamma, double mass) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: moments (rho, rho ux, rho uy, rho uz, totalE).
  // gas_gamma: Adiabatic index. 
  // mass: Species mass. 

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_temp = gkyl_nmat_get(rhs, count); 
  // Clear rhs for each component of flow velocity being solved for 
  gkyl_mat_clear(&rhs_temp, 0.0); 
  const double *rho   = &moms[0]; 
  const double *rhoux = &moms[8]; 
  const double *rhouy = &moms[16]; 
  const double *rhouz = &moms[24]; 
  const double *totE  = &moms[32]; 

  double rhoSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rho, rho, rhoSq); 
 
  double rho_totE[8] = {0.0}; 
  binop_mul_3d_ser_p1(rho, totE, rho_totE); 
 
  double rhouxSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouz, rhouz, rhouzSq); 
 
  double rhoSq_temp[8]; 
  rhoSq_temp[0] = mass*(gas_gamma-1.0)*(rho_totE[0] - 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0])); 
  rhoSq_temp[1] = mass*(gas_gamma-1.0)*(rho_totE[1] - 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1])); 
  rhoSq_temp[2] = mass*(gas_gamma-1.0)*(rho_totE[2] - 0.5*(rhouxSq[2] + rhouySq[2] + rhouzSq[2])); 
  rhoSq_temp[3] = mass*(gas_gamma-1.0)*(rho_totE[3] - 0.5*(rhouxSq[3] + rhouySq[3] + rhouzSq[3])); 
  rhoSq_temp[4] = mass*(gas_gamma-1.0)*(rho_totE[4] - 0.5*(rhouxSq[4] + rhouySq[4] + rhouzSq[4])); 
  rhoSq_temp[5] = mass*(gas_gamma-1.0)*(rho_totE[5] - 0.5*(rhouxSq[5] + rhouySq[5] + rhouzSq[5])); 
  rhoSq_temp[6] = mass*(gas_gamma-1.0)*(rho_totE[6] - 0.5*(rhouxSq[6] + rhouySq[6] + rhouzSq[6])); 
  rhoSq_temp[7] = mass*(gas_gamma-1.0)*(rho_totE[7] - 0.5*(rhouxSq[7] + rhouySq[7] + rhouzSq[7])); 

  double rhoSq_inv[8] = {0.0}; 
  ser_3x_p1_inv(rhoSq, rhoSq_inv); 
  // Calculate expansions of temperature. 
  double temp[8] = {0.0}; 
 
  binop_mul_3d_ser_p1(rhoSq_inv, rhoSq_temp, temp); 
 
  gkyl_mat_set(&rhs_temp,0,0,temp[0]); 
  gkyl_mat_set(&rhs_temp,1,0,temp[1]); 
  gkyl_mat_set(&rhs_temp,2,0,temp[2]); 
  gkyl_mat_set(&rhs_temp,3,0,temp[3]); 
  gkyl_mat_set(&rhs_temp,4,0,temp[4]); 
  gkyl_mat_set(&rhs_temp,5,0,temp[5]); 
  gkyl_mat_set(&rhs_temp,6,0,temp[6]); 
  gkyl_mat_set(&rhs_temp,7,0,temp[7]); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, 
    double* GKYL_RESTRICT out) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // xsol: Input solution vector. 
  // out: Output volume expansion of temperaure. 
 
  struct gkyl_mat x_temp = gkyl_nmat_get(xsol, count); 

  out[0] = gkyl_mat_get(&x_temp,0,0); 
  out[1] = gkyl_mat_get(&x_temp,1,0); 
  out[2] = gkyl_mat_get(&x_temp,2,0); 
  out[3] = gkyl_mat_get(&x_temp,3,0); 
  out[4] = gkyl_mat_get(&x_temp,4,0); 
  out[5] = gkyl_mat_get(&x_temp,5,0); 
  out[6] = gkyl_mat_get(&x_temp,6,0); 
  out[7] = gkyl_mat_get(&x_temp,7,0); 

} 
 
