#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH int em_set_ExB_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *num_ExB) 
{ 
  // count:   integer to indicate which matrix being fetched. 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // num_ExB: Input (E x B)_i (numerator of E x B velocity) and B_i^2 (|B|^2 is E x B velocity denominator). 
 
  struct gkyl_mat A_ExB_x = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_ExB_y = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_ExB_z = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat rhs_ExB_x = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_ExB_y = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_ExB_z = gkyl_nmat_get(rhs, count+2); 
  // Clear matrix and rhs for each component of E x B 
  gkyl_mat_clear(&A_ExB_x, 0.0); gkyl_mat_clear(&rhs_ExB_x, 0.0); 
  gkyl_mat_clear(&A_ExB_y, 0.0); gkyl_mat_clear(&rhs_ExB_y, 0.0); 
  gkyl_mat_clear(&A_ExB_z, 0.0); gkyl_mat_clear(&rhs_ExB_z, 0.0); 
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
  if ((-1.870828693386971*magB2[3])+1.58113883008419*magB2[2]-1.224744871391589*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if (0.7621894676761731*magB2[3]-0.5270462766947298*magB2[2]-0.408248290463863*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if ((-0.7621894676761731*magB2[3])-0.5270462766947298*magB2[2]+0.408248290463863*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if (1.870828693386971*magB2[3]+1.58113883008419*magB2[2]+1.224744871391589*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    magB2[1] = 0.0; 
    magB2[2] = 0.0; 
    magB2[3] = 0.0; 
    gkyl_mat_set(&rhs_ExB_x,0,0,num_ExB_x[0]); 
    gkyl_mat_set(&rhs_ExB_y,0,0,num_ExB_y[0]); 
    gkyl_mat_set(&rhs_ExB_z,0,0,num_ExB_z[0]); 
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
    gkyl_mat_set(&rhs_ExB_x,0,0,num_ExB_x[0]); 
    gkyl_mat_set(&rhs_ExB_y,0,0,num_ExB_y[0]); 
    gkyl_mat_set(&rhs_ExB_z,0,0,num_ExB_z[0]); 
    gkyl_mat_set(&rhs_ExB_x,1,0,num_ExB_x[1]); 
    gkyl_mat_set(&rhs_ExB_y,1,0,num_ExB_y[1]); 
    gkyl_mat_set(&rhs_ExB_z,1,0,num_ExB_z[1]); 
    gkyl_mat_set(&rhs_ExB_x,2,0,num_ExB_x[2]); 
    gkyl_mat_set(&rhs_ExB_y,2,0,num_ExB_y[2]); 
    gkyl_mat_set(&rhs_ExB_z,2,0,num_ExB_z[2]); 
    gkyl_mat_set(&rhs_ExB_x,3,0,num_ExB_x[3]); 
    gkyl_mat_set(&rhs_ExB_y,3,0,num_ExB_y[3]); 
    gkyl_mat_set(&rhs_ExB_z,3,0,num_ExB_z[3]); 
  } 
 
  double temp = 0.0; 
  temp = 0.7071067811865475*magB2[0]; 
  gkyl_mat_set(&A_ExB_x,0,0,temp); 
  gkyl_mat_set(&A_ExB_y,0,0,temp); 
  gkyl_mat_set(&A_ExB_z,0,0,temp); 
  temp = 0.7071067811865475*magB2[1]; 
  gkyl_mat_set(&A_ExB_x,0,1,temp); 
  gkyl_mat_set(&A_ExB_y,0,1,temp); 
  gkyl_mat_set(&A_ExB_z,0,1,temp); 
  temp = 0.7071067811865475*magB2[2]; 
  gkyl_mat_set(&A_ExB_x,0,2,temp); 
  gkyl_mat_set(&A_ExB_y,0,2,temp); 
  gkyl_mat_set(&A_ExB_z,0,2,temp); 
  temp = 0.7071067811865475*magB2[3]; 
  gkyl_mat_set(&A_ExB_x,0,3,temp); 
  gkyl_mat_set(&A_ExB_y,0,3,temp); 
  gkyl_mat_set(&A_ExB_z,0,3,temp); 
  temp = 0.7071067811865475*magB2[1]; 
  gkyl_mat_set(&A_ExB_x,1,0,temp); 
  gkyl_mat_set(&A_ExB_y,1,0,temp); 
  gkyl_mat_set(&A_ExB_z,1,0,temp); 
  temp = 0.6324555320336759*magB2[2]+0.7071067811865475*magB2[0]; 
  gkyl_mat_set(&A_ExB_x,1,1,temp); 
  gkyl_mat_set(&A_ExB_y,1,1,temp); 
  gkyl_mat_set(&A_ExB_z,1,1,temp); 
  temp = 0.6210590034081186*magB2[3]+0.6324555320336759*magB2[1]; 
  gkyl_mat_set(&A_ExB_x,1,2,temp); 
  gkyl_mat_set(&A_ExB_y,1,2,temp); 
  gkyl_mat_set(&A_ExB_z,1,2,temp); 
  temp = 0.6210590034081186*magB2[2]; 
  gkyl_mat_set(&A_ExB_x,1,3,temp); 
  gkyl_mat_set(&A_ExB_y,1,3,temp); 
  gkyl_mat_set(&A_ExB_z,1,3,temp); 
  temp = 0.7071067811865475*magB2[2]; 
  gkyl_mat_set(&A_ExB_x,2,0,temp); 
  gkyl_mat_set(&A_ExB_y,2,0,temp); 
  gkyl_mat_set(&A_ExB_z,2,0,temp); 
  temp = 0.6210590034081186*magB2[3]+0.6324555320336759*magB2[1]; 
  gkyl_mat_set(&A_ExB_x,2,1,temp); 
  gkyl_mat_set(&A_ExB_y,2,1,temp); 
  gkyl_mat_set(&A_ExB_z,2,1,temp); 
  temp = 0.4517539514526256*magB2[2]+0.7071067811865475*magB2[0]; 
  gkyl_mat_set(&A_ExB_x,2,2,temp); 
  gkyl_mat_set(&A_ExB_y,2,2,temp); 
  gkyl_mat_set(&A_ExB_z,2,2,temp); 
  temp = 0.421637021355784*magB2[3]+0.6210590034081186*magB2[1]; 
  gkyl_mat_set(&A_ExB_x,2,3,temp); 
  gkyl_mat_set(&A_ExB_y,2,3,temp); 
  gkyl_mat_set(&A_ExB_z,2,3,temp); 
  temp = 0.7071067811865475*magB2[3]; 
  gkyl_mat_set(&A_ExB_x,3,0,temp); 
  gkyl_mat_set(&A_ExB_y,3,0,temp); 
  gkyl_mat_set(&A_ExB_z,3,0,temp); 
  temp = 0.6210590034081186*magB2[2]; 
  gkyl_mat_set(&A_ExB_x,3,1,temp); 
  gkyl_mat_set(&A_ExB_y,3,1,temp); 
  gkyl_mat_set(&A_ExB_z,3,1,temp); 
  temp = 0.421637021355784*magB2[3]+0.6210590034081186*magB2[1]; 
  gkyl_mat_set(&A_ExB_x,3,2,temp); 
  gkyl_mat_set(&A_ExB_y,3,2,temp); 
  gkyl_mat_set(&A_ExB_z,3,2,temp); 
  temp = 0.421637021355784*magB2[2]+0.7071067811865475*magB2[0]; 
  gkyl_mat_set(&A_ExB_x,3,3,temp); 
  gkyl_mat_set(&A_ExB_y,3,3,temp); 
  gkyl_mat_set(&A_ExB_z,3,3,temp); 
  return cell_avg;
} 
