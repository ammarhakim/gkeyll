#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: Moments [rho, rho ux, rho uy, rho uz, totalE].

  struct gkyl_mat A_flowE = gkyl_nmat_get(A, count); 
  struct gkyl_mat rhs_flowE = gkyl_nmat_get(rhs, count); 
  // Clear matrix and rhs.
  gkyl_mat_clear(&A_flowE, 0.0); gkyl_mat_clear(&rhs_flowE, 0.0); 
  const double *rho   = &moms[0]; 
  const double *rhoux = &moms[3]; 
  const double *rhouy = &moms[6]; 
  const double *rhouz = &moms[9]; 

  double rhouxSq[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhoux, rhoux, rhouxSq); 
 
  double rhouySq[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouy, rhouy, rhouySq); 
 
  double rhouzSq[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouz, rhouz, rhouzSq); 
 
  double rhouSqD2[3]; 
  rhouSqD2[0] = 0.5*(rhouxSq[0] + rhouySq[0] + rhouzSq[0]); 
  rhouSqD2[1] = 0.5*(rhouxSq[1] + rhouySq[1] + rhouzSq[1]); 
  rhouSqD2[2] = 0.5*(rhouxSq[2] + rhouySq[2] + rhouzSq[2]); 

  gkyl_mat_set(&rhs_flowE,0,0,rhouSqD2[0]); 
  gkyl_mat_set(&rhs_flowE,1,0,rhouSqD2[1]); 
  gkyl_mat_set(&rhs_flowE,2,0,rhouSqD2[2]); 
 
  gkyl_mat_set(&A_flowE,0,0,0.7071067811865475*rho[0]); 
 
  gkyl_mat_set(&A_flowE,0,1,0.7071067811865475*rho[1]); 
 
  gkyl_mat_set(&A_flowE,0,2,0.7071067811865475*rho[2]); 
 
  gkyl_mat_set(&A_flowE,1,0,0.7071067811865475*rho[1]); 
 
  gkyl_mat_set(&A_flowE,1,1,0.6324555320336759*rho[2]+0.7071067811865475*rho[0]); 
 
  gkyl_mat_set(&A_flowE,1,2,0.6324555320336759*rho[1]); 
 
  gkyl_mat_set(&A_flowE,2,0,0.7071067811865475*rho[2]); 
 
  gkyl_mat_set(&A_flowE,2,1,0.6324555320336759*rho[1]); 
 
  gkyl_mat_set(&A_flowE,2,2,0.45175395145262565*rho[2]+0.7071067811865475*rho[0]); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_flowE_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, 
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

} 
 
