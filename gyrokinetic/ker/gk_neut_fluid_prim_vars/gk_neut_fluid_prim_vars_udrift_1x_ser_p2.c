#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_set_prob_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: Moments [rho, rho ux, rho uy, rho uz, totalE].

  struct gkyl_mat A_ux = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_uy = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_uz = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  // Clear matrix and rhs for each component of flow velocity being solved for 
  gkyl_mat_clear(&A_ux, 0.0); gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&A_uy, 0.0); gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&A_uz, 0.0); gkyl_mat_clear(&rhs_uz, 0.0); 
  const double *rho   = &moms[0]; 
  const double *rhoux = &moms[3]; 
  const double *rhouy = &moms[6]; 
  const double *rhouz = &moms[9]; 

  gkyl_mat_set(&rhs_ux,0,0,rhoux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,rhouy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,rhouz[0]); 
  gkyl_mat_set(&rhs_ux,1,0,rhoux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,rhouy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,rhouz[1]); 
  gkyl_mat_set(&rhs_ux,2,0,rhoux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,rhouy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,rhouz[2]); 
 
  double temp_rho = 0.0; 
  temp_rho = 0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,temp_rho); 
  gkyl_mat_set(&A_uy,0,0,temp_rho); 
  gkyl_mat_set(&A_uz,0,0,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,temp_rho); 
  gkyl_mat_set(&A_uy,0,1,temp_rho); 
  gkyl_mat_set(&A_uz,0,1,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,temp_rho); 
  gkyl_mat_set(&A_uy,0,2,temp_rho); 
  gkyl_mat_set(&A_uz,0,2,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,temp_rho); 
  gkyl_mat_set(&A_uy,1,0,temp_rho); 
  gkyl_mat_set(&A_uz,1,0,temp_rho); 
 
  temp_rho = 0.6324555320336759*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,temp_rho); 
  gkyl_mat_set(&A_uy,1,1,temp_rho); 
  gkyl_mat_set(&A_uz,1,1,temp_rho); 
 
  temp_rho = 0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,1,2,temp_rho); 
  gkyl_mat_set(&A_uy,1,2,temp_rho); 
  gkyl_mat_set(&A_uz,1,2,temp_rho); 
 
  temp_rho = 0.7071067811865475*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,temp_rho); 
  gkyl_mat_set(&A_uy,2,0,temp_rho); 
  gkyl_mat_set(&A_uz,2,0,temp_rho); 
 
  temp_rho = 0.6324555320336759*rho[1]; 
  gkyl_mat_set(&A_ux,2,1,temp_rho); 
  gkyl_mat_set(&A_uy,2,1,temp_rho); 
  gkyl_mat_set(&A_uz,2,1,temp_rho); 
 
  temp_rho = 0.45175395145262565*rho[2]+0.7071067811865475*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,temp_rho); 
  gkyl_mat_set(&A_uy,2,2,temp_rho); 
  gkyl_mat_set(&A_uz,2,2,temp_rho); 
 
} 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_get_sol_1x_ser_p2(int count, struct gkyl_nmat *xsol, 
    double* GKYL_RESTRICT out) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // xsol: Input solution vector. 
  // out: Output volume expansion of flow velocity;
 
  struct gkyl_mat x_ux = gkyl_nmat_get(xsol, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(xsol, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(xsol, count+2); 
  double *ux = &out[0]; 
  double *uy = &out[3]; 
  double *uz = &out[6]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 

} 
 
