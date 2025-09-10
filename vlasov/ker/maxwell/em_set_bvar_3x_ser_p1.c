#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH int em_set_bvar_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB) 
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
  const double *B_x_B_y = &BB[8]; 
  const double *B_x_B_z = &BB[16]; 
  const double *B_y_sq  = &BB[24]; 
  const double *B_y_B_z = &BB[32]; 
  const double *B_z_sq  = &BB[40]; 
 
  double magB2[8] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB2[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB2[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 
  magB2[4] = B_z_sq[4]+B_y_sq[4]+B_x_sq[4]; 
  magB2[5] = B_z_sq[5]+B_y_sq[5]+B_x_sq[5]; 
  magB2[6] = B_z_sq[6]+B_y_sq[6]+B_x_sq[6]; 
  magB2[7] = B_z_sq[7]+B_y_sq[7]+B_x_sq[7]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if ((-1.837117307087383*magB2[7])+1.060660171779821*magB2[6]+1.060660171779821*magB2[5]+1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]+1.060660171779821*magB2[6]-1.060660171779821*magB2[5]-1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]-1.060660171779821*magB2[6]+1.060660171779821*magB2[5]-1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.837117307087383*magB2[7])-1.060660171779821*magB2[6]-1.060660171779821*magB2[5]+1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]-1.060660171779821*magB2[6]-1.060660171779821*magB2[5]+1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.837117307087383*magB2[7])-1.060660171779821*magB2[6]+1.060660171779821*magB2[5]-1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.837117307087383*magB2[7])+1.060660171779821*magB2[6]-1.060660171779821*magB2[5]-1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]+1.060660171779821*magB2[6]+1.060660171779821*magB2[5]+1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  double magB2_inv[8] = {0.0}; 

  if (cell_avg) { 
  // If |B|^2 < 0 at control points, only use cell average to get 1/|B|^2. 
  magB2_inv[0] = 8.0/magB2[0]; 
  } else { 
  ser_3x_p1_inv(magB2, magB2_inv); 
  } 
  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  double bxbx[8] = {0.0}; 
  double bxby[8] = {0.0}; 
  double bxbz[8] = {0.0}; 
  double byby[8] = {0.0}; 
  double bybz[8] = {0.0}; 
  double bzbz[8] = {0.0}; 
 
  binop_mul_3d_ser_p1(magB2_inv, B_x_sq, bxbx); 
  binop_mul_3d_ser_p1(magB2_inv, B_y_sq, byby); 
  binop_mul_3d_ser_p1(magB2_inv, B_z_sq, bzbz); 
  binop_mul_3d_ser_p1(magB2_inv, B_x_B_y, bxby); 
  binop_mul_3d_ser_p1(magB2_inv, B_x_B_z, bxbz); 
  binop_mul_3d_ser_p1(magB2_inv, B_y_B_z, bybz); 
 
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
    gkyl_mat_set(&rhs_bxbx,4,0,0.0); 
    gkyl_mat_set(&rhs_bxby,4,0,0.0); 
    gkyl_mat_set(&rhs_bxbz,4,0,0.0); 
    gkyl_mat_set(&rhs_byby,4,0,0.0); 
    gkyl_mat_set(&rhs_bybz,4,0,0.0); 
    gkyl_mat_set(&rhs_bzbz,4,0,0.0); 
    gkyl_mat_set(&rhs_bxbx,5,0,0.0); 
    gkyl_mat_set(&rhs_bxby,5,0,0.0); 
    gkyl_mat_set(&rhs_bxbz,5,0,0.0); 
    gkyl_mat_set(&rhs_byby,5,0,0.0); 
    gkyl_mat_set(&rhs_bybz,5,0,0.0); 
    gkyl_mat_set(&rhs_bzbz,5,0,0.0); 
    gkyl_mat_set(&rhs_bxbx,6,0,0.0); 
    gkyl_mat_set(&rhs_bxby,6,0,0.0); 
    gkyl_mat_set(&rhs_bxbz,6,0,0.0); 
    gkyl_mat_set(&rhs_byby,6,0,0.0); 
    gkyl_mat_set(&rhs_bybz,6,0,0.0); 
    gkyl_mat_set(&rhs_bzbz,6,0,0.0); 
    gkyl_mat_set(&rhs_bxbx,7,0,0.0); 
    gkyl_mat_set(&rhs_bxby,7,0,0.0); 
    gkyl_mat_set(&rhs_bxbz,7,0,0.0); 
    gkyl_mat_set(&rhs_byby,7,0,0.0); 
    gkyl_mat_set(&rhs_bybz,7,0,0.0); 
    gkyl_mat_set(&rhs_bzbz,7,0,0.0); 
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
    gkyl_mat_set(&rhs_bxbx,4,0,bxbx[4]); 
    gkyl_mat_set(&rhs_bxby,4,0,bxby[4]); 
    gkyl_mat_set(&rhs_bxbz,4,0,bxbz[4]); 
    gkyl_mat_set(&rhs_byby,4,0,byby[4]); 
    gkyl_mat_set(&rhs_bybz,4,0,bybz[4]); 
    gkyl_mat_set(&rhs_bzbz,4,0,bzbz[4]); 
    gkyl_mat_set(&rhs_bxbx,5,0,bxbx[5]); 
    gkyl_mat_set(&rhs_bxby,5,0,bxby[5]); 
    gkyl_mat_set(&rhs_bxbz,5,0,bxbz[5]); 
    gkyl_mat_set(&rhs_byby,5,0,byby[5]); 
    gkyl_mat_set(&rhs_bybz,5,0,bybz[5]); 
    gkyl_mat_set(&rhs_bzbz,5,0,bzbz[5]); 
    gkyl_mat_set(&rhs_bxbx,6,0,bxbx[6]); 
    gkyl_mat_set(&rhs_bxby,6,0,bxby[6]); 
    gkyl_mat_set(&rhs_bxbz,6,0,bxbz[6]); 
    gkyl_mat_set(&rhs_byby,6,0,byby[6]); 
    gkyl_mat_set(&rhs_bybz,6,0,bybz[6]); 
    gkyl_mat_set(&rhs_bzbz,6,0,bzbz[6]); 
    gkyl_mat_set(&rhs_bxbx,7,0,bxbx[7]); 
    gkyl_mat_set(&rhs_bxby,7,0,bxby[7]); 
    gkyl_mat_set(&rhs_bxbz,7,0,bxbz[7]); 
    gkyl_mat_set(&rhs_byby,7,0,byby[7]); 
    gkyl_mat_set(&rhs_bybz,7,0,bybz[7]); 
    gkyl_mat_set(&rhs_bzbz,7,0,bzbz[7]); 
  } 
 
  return cell_avg;
} 
