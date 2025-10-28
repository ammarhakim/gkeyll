#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH int em_set_ExB_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB) 
{ 
  // count:   integer to indicate which matrix being fetched. 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // num_ExB: Input (E x B)_i (numerator of E x B velocity) and B_i^2 (|B|^2 is E x B velocity denominator). 
 
  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ExB_x = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_ExB_y = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_ExB_z = gkyl_nmat_get(rhs, count+2); 
  // Clear rhs for each component of E x B 
  gkyl_mat_clear(&rhs_ExB_x, 0.0); 
  gkyl_mat_clear(&rhs_ExB_y, 0.0); 
  gkyl_mat_clear(&rhs_ExB_z, 0.0); 
  const double *num_ExB_x = &num_ExB[0]; 
  const double *num_ExB_y = &num_ExB[4]; 
  const double *num_ExB_z = &num_ExB[8]; 
  const double *B_x_sq  = &num_ExB[12]; 
  const double *B_y_sq  = &num_ExB[16]; 
  const double *B_z_sq  = &num_ExB[20]; 
 
  double magB2[4] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB2[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB2[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if (1.5*magB2[3]-0.8660254037844386*magB2[2]-0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.5*magB2[3])-0.8660254037844386*magB2[2]+0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.5*magB2[3])+0.8660254037844386*magB2[2]-0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  if (1.5*magB2[3]+0.8660254037844386*magB2[2]+0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  double magB2_inv[4] = {0.0}; 

  if (cell_avg) { 
  // If |B|^2 < 0 at control points, only use cell average to get 1/|B|^2. 
  magB2_inv[0] = 4.0/magB2[0]; 
  } else { 
  ser_2x_p1_inv(magB2, magB2_inv); 
  } 
  // Calculate expansions of E x B/|B|^2, which can be calculated free of aliasing errors. 
  double ExB_x[4] = {0.0}; 
  double ExB_y[4] = {0.0}; 
  double ExB_z[4] = {0.0}; 
 
  binop_mul_2d_ser_p1(magB2_inv, num_ExB_x, ExB_x); 
  binop_mul_2d_ser_p1(magB2_inv, num_ExB_y, ExB_y); 
  binop_mul_2d_ser_p1(magB2_inv, num_ExB_z, ExB_z); 
 
  if (cell_avg) { 
    gkyl_mat_set(&rhs_ExB_x,0,0,ExB_x[0]); 
    gkyl_mat_set(&rhs_ExB_y,0,0,ExB_y[0]); 
    gkyl_mat_set(&rhs_ExB_z,0,0,ExB_z[0]); 
    gkyl_mat_set(&rhs_ExB_x,1,0,0.0); 
    gkyl_mat_set(&rhs_ExB_y,1,0,0.0); 
    gkyl_mat_set(&rhs_ExB_z,1,0,0.0); 
    gkyl_mat_set(&rhs_ExB_x,2,0,0.0); 
    gkyl_mat_set(&rhs_ExB_y,2,0,0.0); 
    gkyl_mat_set(&rhs_ExB_z,2,0,0.0); 
    gkyl_mat_set(&rhs_ExB_x,3,0,0.0); 
    gkyl_mat_set(&rhs_ExB_y,3,0,0.0); 
    gkyl_mat_set(&rhs_ExB_z,3,0,0.0); 
  } else { 
    gkyl_mat_set(&rhs_ExB_x,0,0,ExB_x[0]); 
    gkyl_mat_set(&rhs_ExB_y,0,0,ExB_y[0]); 
    gkyl_mat_set(&rhs_ExB_z,0,0,ExB_z[0]); 
    gkyl_mat_set(&rhs_ExB_x,1,0,ExB_x[1]); 
    gkyl_mat_set(&rhs_ExB_y,1,0,ExB_y[1]); 
    gkyl_mat_set(&rhs_ExB_z,1,0,ExB_z[1]); 
    gkyl_mat_set(&rhs_ExB_x,2,0,ExB_x[2]); 
    gkyl_mat_set(&rhs_ExB_y,2,0,ExB_y[2]); 
    gkyl_mat_set(&rhs_ExB_z,2,0,ExB_z[2]); 
    gkyl_mat_set(&rhs_ExB_x,3,0,ExB_x[3]); 
    gkyl_mat_set(&rhs_ExB_y,3,0,ExB_y[3]); 
    gkyl_mat_set(&rhs_ExB_z,3,0,ExB_z[3]); 
  } 
 
  return cell_avg;
} 
