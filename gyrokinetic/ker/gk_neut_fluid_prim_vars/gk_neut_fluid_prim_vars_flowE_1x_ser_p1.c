#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_set_prob_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: moments (rho, rho ux, rho uy, rho uz, totalE).

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_flowE = gkyl_nmat_get(rhs, count); 
  // Clear rhs for each component of flow velocity being solved for 
  gkyl_mat_clear(&rhs_flowE, 0.0); 
  const double *rho   = &moms[0]; 
  const double *rhoux = &moms[2]; 
  const double *rhouy = &moms[4]; 
  const double *rhouz = &moms[6]; 

  double rhouxSq[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[2] = {0.0}; 
  binop_mul_1d_ser_p1(rhouz, rhouz, rhouzSq); 
 
  double rhouSqD2[2]; 
  rhouSqD2[0] = 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0]); 
  rhouSqD2[1] = 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1]); 

  double rho_inv[2] = {0.0}; 
  ser_1x_p1_inv(rho, rho_inv); 
  // Calculate expansions of flow energy. 
  double flowE[2] = {0.0}; 
  binop_mul_1d_ser_p1(rho_inv, rhouSqD2, flowE); 
 
  gkyl_mat_set(&rhs_flowE,0,0,flowE[0]); 
  gkyl_mat_set(&rhs_flowE,1,0,flowE[1]); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_get_sol_1x_ser_p1(int count, struct gkyl_nmat *xsol, 
    double* GKYL_RESTRICT out) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // xsol: Input solution vector. 
  // out: Output volume expansion of flow velocity and temperature;
 
  struct gkyl_mat x_flowE = gkyl_nmat_get(xsol, count); 
  double *flowE = &out[0]; 

  flowE[0] = gkyl_mat_get(&x_flowE,0,0); 
  flowE[1] = gkyl_mat_get(&x_flowE,1,0); 

} 
 
