#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH int em_set_bvar_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // BB:    Input magnetic field tensor B_i B_j. 
 
  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_bxbx = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_bxby = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_bxbz = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_byby = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_bybz = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_bzbz = gkyl_nmat_get(rhs, count+5); 
  // Clear rhs for each component of b_i b_j 
  gkyl_mat_clear(&rhs_bxbx, 0.0); 
  gkyl_mat_clear(&rhs_bxby, 0.0); 
  gkyl_mat_clear(&rhs_bxbz, 0.0); 
  gkyl_mat_clear(&rhs_byby, 0.0); 
  gkyl_mat_clear(&rhs_bybz, 0.0); 
  gkyl_mat_clear(&rhs_bzbz, 0.0); 
  const double *B_x_sq  = &BB[0]; 
  const double *B_x_B_y = &BB[4]; 
  const double *B_x_B_z = &BB[8]; 
  const double *B_y_sq  = &BB[12]; 
  const double *B_y_B_z = &BB[16]; 
  const double *B_z_sq  = &BB[20]; 
 
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
  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  double bxbx[4] = {0.0}; 
  double bxby[4] = {0.0}; 
  double bxbz[4] = {0.0}; 
  double byby[4] = {0.0}; 
  double bybz[4] = {0.0}; 
  double bzbz[4] = {0.0}; 
 
  binop_mul_2d_ser_p1(magB2_inv, B_x_sq, bxbx); 
  binop_mul_2d_ser_p1(magB2_inv, B_y_sq, byby); 
  binop_mul_2d_ser_p1(magB2_inv, B_z_sq, bzbz); 
  binop_mul_2d_ser_p1(magB2_inv, B_x_B_y, bxby); 
  binop_mul_2d_ser_p1(magB2_inv, B_x_B_z, bxbz); 
  binop_mul_2d_ser_p1(magB2_inv, B_y_B_z, bybz); 
 
  if (cell_avg) { 
    gkyl_mat_set(&rhs_bxbx,0,0,bxbx[0]); 
    gkyl_mat_set(&rhs_bxby,0,0,bxby[0]); 
    gkyl_mat_set(&rhs_bxbz,0,0,bxbz[0]); 
    gkyl_mat_set(&rhs_byby,0,0,byby[0]); 
    gkyl_mat_set(&rhs_bybz,0,0,bybz[0]); 
    gkyl_mat_set(&rhs_bzbz,0,0,bzbz[0]); 
    gkyl_mat_set(&rhs_bxbx,1,0,0.0); 
    gkyl_mat_set(&rhs_bxby,1,0,0.0); 
    gkyl_mat_set(&rhs_bxbz,1,0,0.0); 
    gkyl_mat_set(&rhs_byby,1,0,0.0); 
    gkyl_mat_set(&rhs_bybz,1,0,0.0); 
    gkyl_mat_set(&rhs_bzbz,1,0,0.0); 
    gkyl_mat_set(&rhs_bxbx,2,0,0.0); 
    gkyl_mat_set(&rhs_bxby,2,0,0.0); 
    gkyl_mat_set(&rhs_bxbz,2,0,0.0); 
    gkyl_mat_set(&rhs_byby,2,0,0.0); 
    gkyl_mat_set(&rhs_bybz,2,0,0.0); 
    gkyl_mat_set(&rhs_bzbz,2,0,0.0); 
    gkyl_mat_set(&rhs_bxbx,3,0,0.0); 
    gkyl_mat_set(&rhs_bxby,3,0,0.0); 
    gkyl_mat_set(&rhs_bxbz,3,0,0.0); 
    gkyl_mat_set(&rhs_byby,3,0,0.0); 
    gkyl_mat_set(&rhs_bybz,3,0,0.0); 
    gkyl_mat_set(&rhs_bzbz,3,0,0.0); 
  } else { 
    gkyl_mat_set(&rhs_bxbx,0,0,bxbx[0]); 
    gkyl_mat_set(&rhs_bxby,0,0,bxby[0]); 
    gkyl_mat_set(&rhs_bxbz,0,0,bxbz[0]); 
    gkyl_mat_set(&rhs_byby,0,0,byby[0]); 
    gkyl_mat_set(&rhs_bybz,0,0,bybz[0]); 
    gkyl_mat_set(&rhs_bzbz,0,0,bzbz[0]); 
    gkyl_mat_set(&rhs_bxbx,1,0,bxbx[1]); 
    gkyl_mat_set(&rhs_bxby,1,0,bxby[1]); 
    gkyl_mat_set(&rhs_bxbz,1,0,bxbz[1]); 
    gkyl_mat_set(&rhs_byby,1,0,byby[1]); 
    gkyl_mat_set(&rhs_bybz,1,0,bybz[1]); 
    gkyl_mat_set(&rhs_bzbz,1,0,bzbz[1]); 
    gkyl_mat_set(&rhs_bxbx,2,0,bxbx[2]); 
    gkyl_mat_set(&rhs_bxby,2,0,bxby[2]); 
    gkyl_mat_set(&rhs_bxbz,2,0,bxbz[2]); 
    gkyl_mat_set(&rhs_byby,2,0,byby[2]); 
    gkyl_mat_set(&rhs_bybz,2,0,bybz[2]); 
    gkyl_mat_set(&rhs_bzbz,2,0,bzbz[2]); 
    gkyl_mat_set(&rhs_bxbx,3,0,bxbx[3]); 
    gkyl_mat_set(&rhs_bxby,3,0,bxby[3]); 
    gkyl_mat_set(&rhs_bxbz,3,0,bxbz[3]); 
    gkyl_mat_set(&rhs_byby,3,0,byby[3]); 
    gkyl_mat_set(&rhs_bybz,3,0,bybz[3]); 
    gkyl_mat_set(&rhs_bzbz,3,0,bzbz[3]); 
  } 
 
  return cell_avg;
} 
