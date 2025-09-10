#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH int pkpm_vars_u_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       [rho ux, rho uy, rho uz], Fluid input state vector.

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  // Clear rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&rhs_uz, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[4]; 
  const double *rhouz = &euler_pkpm[8]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  int cell_avg = 0;
  // Check if rho < 0 at control points. 
  // *THIS IS ONLY A CHECK RIGHT NOW AND UNUSED* 
  if (1.5*rho[3]-0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.5*rho[3])-0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if ((-1.5*rho[3])+0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (1.5*rho[3]+0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  double rho_inv[4] = {0.0}; 
  ser_2x_p1_inv(rho, rho_inv); 
  // Calculate expansions of primitive variables, which can be calculated free of aliasing errors. 
  double ux[4] = {0.0}; 
  double uy[4] = {0.0}; 
  double uz[4] = {0.0}; 
 
  binop_mul_2d_ser_p1(rho_inv, rhoux, ux); 
  binop_mul_2d_ser_p1(rho_inv, rhouy, uy); 
  binop_mul_2d_ser_p1(rho_inv, rhouz, uz); 
 
  gkyl_mat_set(&rhs_ux,0,0,ux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,uy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,uz[0]); 
  gkyl_mat_set(&rhs_ux,1,0,ux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,uy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,uz[1]); 
  gkyl_mat_set(&rhs_ux,2,0,ux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,uy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,uz[2]); 
  gkyl_mat_set(&rhs_ux,3,0,ux[3]); 
  gkyl_mat_set(&rhs_uy,3,0,uy[3]); 
  gkyl_mat_set(&rhs_uz,3,0,uz[3]); 
 
  return cell_avg;
} 
