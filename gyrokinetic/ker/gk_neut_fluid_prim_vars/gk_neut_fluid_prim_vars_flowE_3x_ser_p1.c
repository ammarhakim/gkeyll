#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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
  const double *rhoux = &moms[8]; 
  const double *rhouy = &moms[16]; 
  const double *rhouz = &moms[24]; 

  double rhouxSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(rhouz, rhouz, rhouzSq); 
 
  double rhouSqD2[8]; 
  rhouSqD2[0] = 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0]); 
  rhouSqD2[1] = 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1]); 
  rhouSqD2[2] = 0.5*(rhouxSq[2] + rhouySq[2] + rhouzSq[2]); 
  rhouSqD2[3] = 0.5*(rhouxSq[3] + rhouySq[3] + rhouzSq[3]); 
  rhouSqD2[4] = 0.5*(rhouxSq[4] + rhouySq[4] + rhouzSq[4]); 
  rhouSqD2[5] = 0.5*(rhouxSq[5] + rhouySq[5] + rhouzSq[5]); 
  rhouSqD2[6] = 0.5*(rhouxSq[6] + rhouySq[6] + rhouzSq[6]); 
  rhouSqD2[7] = 0.5*(rhouxSq[7] + rhouySq[7] + rhouzSq[7]); 

  double rho_inv[8] = {0.0}; 
  ser_3x_p1_inv(rho, rho_inv); 
  // Calculate expansions of flow energy. 
  double flowE[8] = {0.0}; 
  binop_mul_3d_ser_p1(rho_inv, rhouSqD2, flowE); 
 
  gkyl_mat_set(&rhs_flowE,0,0,flowE[0]); 
  gkyl_mat_set(&rhs_flowE,1,0,flowE[1]); 
  gkyl_mat_set(&rhs_flowE,2,0,flowE[2]); 
  gkyl_mat_set(&rhs_flowE,3,0,flowE[3]); 
  gkyl_mat_set(&rhs_flowE,4,0,flowE[4]); 
  gkyl_mat_set(&rhs_flowE,5,0,flowE[5]); 
  gkyl_mat_set(&rhs_flowE,6,0,flowE[6]); 
  gkyl_mat_set(&rhs_flowE,7,0,flowE[7]); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, 
    double* GKYL_RESTRICT out) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // xsol: Input solution vector. 
  // out: Output volume expansion of flow velocity and temperature;
 
  struct gkyl_mat x_flowE = gkyl_nmat_get(xsol, count); 
  double *flowE = &out[0]; 

  flowE[0] = gkyl_mat_get(&x_flowE,0,0); 
  flowE[1] = gkyl_mat_get(&x_flowE,1,0); 
  flowE[2] = gkyl_mat_get(&x_flowE,2,0); 
  flowE[3] = gkyl_mat_get(&x_flowE,3,0); 
  flowE[4] = gkyl_mat_get(&x_flowE,4,0); 
  flowE[5] = gkyl_mat_get(&x_flowE,5,0); 
  flowE[6] = gkyl_mat_get(&x_flowE,6,0); 
  flowE[7] = gkyl_mat_get(&x_flowE,7,0); 

} 
 
