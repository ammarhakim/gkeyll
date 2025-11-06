#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_set_prob_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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
  const double *rhoux = &moms[4]; 
  const double *rhouy = &moms[8]; 
  const double *rhouz = &moms[12]; 
  const double *totE  = &moms[16]; 

  double rhouxSq[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[4] = {0.0}; 
  binop_mul_2d_ser_p1(rhouz, rhouz, rhouzSq); 
 
  double rho_temp[4]; 
  rho_temp[0] = (gas_gamma - 1.0)*(mass * totE[0] - 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0])); 
  rho_temp[1] = (gas_gamma - 1.0)*(mass * totE[1] - 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1])); 
  rho_temp[2] = (gas_gamma - 1.0)*(mass * totE[2] - 0.5*(rhouxSq[2] + rhouySq[2] + rhouzSq[2])); 
  rho_temp[3] = (gas_gamma - 1.0)*(mass * totE[3] - 0.5*(rhouxSq[3] + rhouySq[3] + rhouzSq[3])); 

  double rho_inv[4] = {0.0}; 
  ser_2x_p1_inv(rho, rho_inv); 
  // Calculate expansions of temperature. 
  double temp[4] = {0.0}; 
 
  binop_mul_2d_ser_p1(rho_inv, rho_temp, temp); 
 
  gkyl_mat_set(&rhs_temp,0,0,temp[0]); 
  gkyl_mat_set(&rhs_temp,1,0,temp[1]); 
  gkyl_mat_set(&rhs_temp,2,0,temp[2]); 
  gkyl_mat_set(&rhs_temp,3,0,temp[3]); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_temp_get_sol_2x_ser_p1(int count, struct gkyl_nmat *xsol, 
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

} 
 
