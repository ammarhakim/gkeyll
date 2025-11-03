#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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
  const double *rhoux = &moms[8]; 
  const double *rhouy = &moms[16]; 
  const double *rhouz = &moms[24]; 
  const double *totE  = &moms[32]; 

  double rhouxSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouz, rhouz, rhouzSq); 
 
  double rho_inv[8] = {0.0}; 
  ser_3x_p1_inv(rho, rho_inv); 
  // Calculate expansions of flow velocity. 
  double ux[8] = {0.0}; 
  double uy[8] = {0.0}; 
  double uz[8] = {0.0}; 
 
  binop_mul_3d_ser_p1(rho_inv, rhoux, ux); 
  binop_mul_3d_ser_p1(rho_inv, rhouy, uy); 
  binop_mul_3d_ser_p1(rho_inv, rhouz, uz); 
 
  double rho_temp[8]; 
  rho_temp[0] = (gas_gamma - 1.0)*(mass * totE[0] - 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0])); 
  rho_temp[1] = (gas_gamma - 1.0)*(mass * totE[1] - 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1])); 
  rho_temp[2] = (gas_gamma - 1.0)*(mass * totE[2] - 0.5*(rhouxSq[2] + rhouySq[2] + rhouzSq[2])); 
  rho_temp[3] = (gas_gamma - 1.0)*(mass * totE[3] - 0.5*(rhouxSq[3] + rhouySq[3] + rhouzSq[3])); 
  rho_temp[4] = (gas_gamma - 1.0)*(mass * totE[4] - 0.5*(rhouxSq[4] + rhouySq[4] + rhouzSq[4])); 
  rho_temp[5] = (gas_gamma - 1.0)*(mass * totE[5] - 0.5*(rhouxSq[5] + rhouySq[5] + rhouzSq[5])); 
  rho_temp[6] = (gas_gamma - 1.0)*(mass * totE[6] - 0.5*(rhouxSq[6] + rhouySq[6] + rhouzSq[6])); 
  rho_temp[7] = (gas_gamma - 1.0)*(mass * totE[7] - 0.5*(rhouxSq[7] + rhouySq[7] + rhouzSq[7])); 

  // Calculate expansions of temperature. 
  double temp[8] = {0.0}; 
 
  binop_mul_3d_ser_p1(rho_inv, rho_temp, temp); 
  gkyl_mat_set(&rhs_ux,0,0,ux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,uy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,uz[0]); 
  gkyl_mat_set(&rhs_temp,0,0,temp[0]); 
  gkyl_mat_set(&rhs_ux,1,0,ux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,uy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,uz[1]); 
  gkyl_mat_set(&rhs_temp,1,0,temp[1]); 
  gkyl_mat_set(&rhs_ux,2,0,ux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,uy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,uz[2]); 
  gkyl_mat_set(&rhs_temp,2,0,temp[2]); 
  gkyl_mat_set(&rhs_ux,3,0,ux[3]); 
  gkyl_mat_set(&rhs_uy,3,0,uy[3]); 
  gkyl_mat_set(&rhs_uz,3,0,uz[3]); 
  gkyl_mat_set(&rhs_temp,3,0,temp[3]); 
  gkyl_mat_set(&rhs_ux,4,0,ux[4]); 
  gkyl_mat_set(&rhs_uy,4,0,uy[4]); 
  gkyl_mat_set(&rhs_uz,4,0,uz[4]); 
  gkyl_mat_set(&rhs_temp,4,0,temp[4]); 
  gkyl_mat_set(&rhs_ux,5,0,ux[5]); 
  gkyl_mat_set(&rhs_uy,5,0,uy[5]); 
  gkyl_mat_set(&rhs_uz,5,0,uz[5]); 
  gkyl_mat_set(&rhs_temp,5,0,temp[5]); 
  gkyl_mat_set(&rhs_ux,6,0,ux[6]); 
  gkyl_mat_set(&rhs_uy,6,0,uy[6]); 
  gkyl_mat_set(&rhs_uz,6,0,uz[6]); 
  gkyl_mat_set(&rhs_temp,6,0,temp[6]); 
  gkyl_mat_set(&rhs_ux,7,0,ux[7]); 
  gkyl_mat_set(&rhs_uy,7,0,uy[7]); 
  gkyl_mat_set(&rhs_uz,7,0,uz[7]); 
  gkyl_mat_set(&rhs_temp,7,0,temp[7]); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_temp_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, 
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
  double *uy = &out[8]; 
  double *uz = &out[16]; 
  double *temp = &out[24]; 

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
  ux[3] = gkyl_mat_get(&x_ux,3,0); 
  uy[3] = gkyl_mat_get(&x_uy,3,0); 
  uz[3] = gkyl_mat_get(&x_uz,3,0); 
  temp[3] = gkyl_mat_get(&x_temp,3,0); 
  ux[4] = gkyl_mat_get(&x_ux,4,0); 
  uy[4] = gkyl_mat_get(&x_uy,4,0); 
  uz[4] = gkyl_mat_get(&x_uz,4,0); 
  temp[4] = gkyl_mat_get(&x_temp,4,0); 
  ux[5] = gkyl_mat_get(&x_ux,5,0); 
  uy[5] = gkyl_mat_get(&x_uy,5,0); 
  uz[5] = gkyl_mat_get(&x_uz,5,0); 
  temp[5] = gkyl_mat_get(&x_temp,5,0); 
  ux[6] = gkyl_mat_get(&x_ux,6,0); 
  uy[6] = gkyl_mat_get(&x_uy,6,0); 
  uz[6] = gkyl_mat_get(&x_uz,6,0); 
  temp[6] = gkyl_mat_get(&x_temp,6,0); 
  ux[7] = gkyl_mat_get(&x_ux,7,0); 
  uy[7] = gkyl_mat_get(&x_uy,7,0); 
  uz[7] = gkyl_mat_get(&x_uz,7,0); 
  temp[7] = gkyl_mat_get(&x_temp,7,0); 

} 
 
