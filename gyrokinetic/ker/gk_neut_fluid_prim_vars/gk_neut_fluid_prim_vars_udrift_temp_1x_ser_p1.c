#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_set_prob_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms, double gas_gamma, double mass) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: moments (rho, rho ux, rho uy, rho uz, totalE).
  // gas_gamma: Adiabatic index. 
  // mass: Species mass. 

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_temp = gkyl_nmat_get(rhs, count+3); 
  // Clear rhs for each component of flow velocity being solved for 
  gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&rhs_uz, 0.0); 
  gkyl_mat_clear(&rhs_temp, 0.0); 
  const double *rho   = &moms[0]; 
  const double *rhoux = &moms[2]; 
  const double *rhouy = &moms[4]; 
  const double *rhouz = &moms[6]; 
  const double *totE  = &moms[8]; 

  double rhoSq[2] = {0.0}; 
  binop_mul_1d_ser_p1(rho, rho, rhoSq); 
 
  double rho_totE[2] = {0.0}; 
  binop_mul_1d_ser_p1(rho, totE, rho_totE); 
 
  double rhouxSq[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhouz, rhouz, rhouzSq); 
 
  double rho_inv[2] = {0.0}; 
  ser_1x_p1_inv(rho, rho_inv); 
  // Calculate expansions of flow velocity. 
  double ux[2] = {0.0}; 
  double uy[2] = {0.0}; 
  double uz[2] = {0.0}; 
 
  binop_mul_1d_ser_p1(rho_inv, rhoux, ux); 
  binop_mul_1d_ser_p1(rho_inv, rhouy, uy); 
  binop_mul_1d_ser_p1(rho_inv, rhouz, uz); 
 
  double rhoSq_temp[2]; 
  rhoSq_temp[0] = mass*(gas_gamma-1.0)*(rho_totE[0] - 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0])); 
  rhoSq_temp[1] = mass*(gas_gamma-1.0)*(rho_totE[1] - 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1])); 

  double rhoSq_inv[2] = {0.0}; 
  ser_1x_p1_inv(rhoSq, rhoSq_inv); 
  // Calculate expansions of temperature. 
  double temp[2] = {0.0}; 
 
  binop_mul_1d_ser_p1(rhoSq_inv, rhoSq_temp, temp); 
  gkyl_mat_set(&rhs_ux,0,0,ux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,uy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,uz[0]); 
  gkyl_mat_set(&rhs_temp,0,0,temp[0]); 
  gkyl_mat_set(&rhs_ux,1,0,ux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,uy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,uz[1]); 
  gkyl_mat_set(&rhs_temp,1,0,temp[1]); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_get_sol_1x_ser_p1(int count, struct gkyl_nmat *xsol, 
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
  double *uy = &out[2]; 
  double *uz = &out[4]; 
  double *temp = &out[6]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  temp[0] = gkyl_mat_get(&x_temp,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  temp[1] = gkyl_mat_get(&x_temp,1,0); 

} 
 
