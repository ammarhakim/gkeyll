#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH int em_set_bvar_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // BB:    Input magnetic field tensor B_i B_j. 
 
  struct gkyl_mat A_bxbx = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_bxby = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_bxbz = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat A_byby = gkyl_nmat_get(A, count+3); 
  struct gkyl_mat A_bybz = gkyl_nmat_get(A, count+4); 
  struct gkyl_mat A_bzbz = gkyl_nmat_get(A, count+5); 
  struct gkyl_mat rhs_bxbx = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_bxby = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_bxbz = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_byby = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_bybz = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_bzbz = gkyl_nmat_get(rhs, count+5); 
  // Clear matrix and rhs for each component of b_i b_j 
  gkyl_mat_clear(&A_bxbx, 0.0); gkyl_mat_clear(&rhs_bxbx, 0.0); 
  gkyl_mat_clear(&A_bxby, 0.0); gkyl_mat_clear(&rhs_bxby, 0.0); 
  gkyl_mat_clear(&A_bxbz, 0.0); gkyl_mat_clear(&rhs_bxbz, 0.0); 
  gkyl_mat_clear(&A_byby, 0.0); gkyl_mat_clear(&rhs_byby, 0.0); 
  gkyl_mat_clear(&A_bybz, 0.0); gkyl_mat_clear(&rhs_bybz, 0.0); 
  gkyl_mat_clear(&A_bzbz, 0.0); gkyl_mat_clear(&rhs_bzbz, 0.0); 
  const double *B_x_sq  = &BB[0]; 
  const double *B_x_B_y = &BB[3]; 
  const double *B_x_B_z = &BB[6]; 
  const double *B_y_sq  = &BB[9]; 
  const double *B_y_B_z = &BB[12]; 
  const double *B_z_sq  = &BB[15]; 
 
  double magB2[3] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB2[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if (1.58113883008419*magB2[2]-1.224744871391589*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*magB2[0]-0.7905694150420947*magB2[2] < 0.0) cell_avg = 1; 
  if (1.58113883008419*magB2[2]+1.224744871391589*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    magB2[1] = 0.0; 
    magB2[2] = 0.0; 
    gkyl_mat_set(&rhs_bxbx,0,0,B_x_sq[0]); 
    gkyl_mat_set(&rhs_bxby,0,0,B_x_B_y[0]); 
    gkyl_mat_set(&rhs_bxbz,0,0,B_x_B_z[0]); 
    gkyl_mat_set(&rhs_byby,0,0,B_y_sq[0]); 
    gkyl_mat_set(&rhs_bybz,0,0,B_y_B_z[0]); 
    gkyl_mat_set(&rhs_bzbz,0,0,B_z_sq[0]); 
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
  } else { 
    gkyl_mat_set(&rhs_bxbx,0,0,B_x_sq[0]); 
    gkyl_mat_set(&rhs_bxby,0,0,B_x_B_y[0]); 
    gkyl_mat_set(&rhs_bxbz,0,0,B_x_B_z[0]); 
    gkyl_mat_set(&rhs_byby,0,0,B_y_sq[0]); 
    gkyl_mat_set(&rhs_bybz,0,0,B_y_B_z[0]); 
    gkyl_mat_set(&rhs_bzbz,0,0,B_z_sq[0]); 
    gkyl_mat_set(&rhs_bxbx,1,0,B_x_sq[1]); 
    gkyl_mat_set(&rhs_bxby,1,0,B_x_B_y[1]); 
    gkyl_mat_set(&rhs_bxbz,1,0,B_x_B_z[1]); 
    gkyl_mat_set(&rhs_byby,1,0,B_y_sq[1]); 
    gkyl_mat_set(&rhs_bybz,1,0,B_y_B_z[1]); 
    gkyl_mat_set(&rhs_bzbz,1,0,B_z_sq[1]); 
    gkyl_mat_set(&rhs_bxbx,2,0,B_x_sq[2]); 
    gkyl_mat_set(&rhs_bxby,2,0,B_x_B_y[2]); 
    gkyl_mat_set(&rhs_bxbz,2,0,B_x_B_z[2]); 
    gkyl_mat_set(&rhs_byby,2,0,B_y_sq[2]); 
    gkyl_mat_set(&rhs_bybz,2,0,B_y_B_z[2]); 
    gkyl_mat_set(&rhs_bzbz,2,0,B_z_sq[2]); 
  } 
 
  double temp = 0.0; 
  temp = 0.7071067811865475*magB2[0]; 
  gkyl_mat_set(&A_bxbx,0,0,temp); 
  gkyl_mat_set(&A_bxby,0,0,temp); 
  gkyl_mat_set(&A_bxbz,0,0,temp); 
  gkyl_mat_set(&A_byby,0,0,temp); 
  gkyl_mat_set(&A_bybz,0,0,temp); 
  gkyl_mat_set(&A_bzbz,0,0,temp); 
  temp = 0.7071067811865475*magB2[1]; 
  gkyl_mat_set(&A_bxbx,0,1,temp); 
  gkyl_mat_set(&A_bxby,0,1,temp); 
  gkyl_mat_set(&A_bxbz,0,1,temp); 
  gkyl_mat_set(&A_byby,0,1,temp); 
  gkyl_mat_set(&A_bybz,0,1,temp); 
  gkyl_mat_set(&A_bzbz,0,1,temp); 
  temp = 0.7071067811865475*magB2[2]; 
  gkyl_mat_set(&A_bxbx,0,2,temp); 
  gkyl_mat_set(&A_bxby,0,2,temp); 
  gkyl_mat_set(&A_bxbz,0,2,temp); 
  gkyl_mat_set(&A_byby,0,2,temp); 
  gkyl_mat_set(&A_bybz,0,2,temp); 
  gkyl_mat_set(&A_bzbz,0,2,temp); 
  temp = 0.7071067811865475*magB2[1]; 
  gkyl_mat_set(&A_bxbx,1,0,temp); 
  gkyl_mat_set(&A_bxby,1,0,temp); 
  gkyl_mat_set(&A_bxbz,1,0,temp); 
  gkyl_mat_set(&A_byby,1,0,temp); 
  gkyl_mat_set(&A_bybz,1,0,temp); 
  gkyl_mat_set(&A_bzbz,1,0,temp); 
  temp = 0.6324555320336759*magB2[2]+0.7071067811865475*magB2[0]; 
  gkyl_mat_set(&A_bxbx,1,1,temp); 
  gkyl_mat_set(&A_bxby,1,1,temp); 
  gkyl_mat_set(&A_bxbz,1,1,temp); 
  gkyl_mat_set(&A_byby,1,1,temp); 
  gkyl_mat_set(&A_bybz,1,1,temp); 
  gkyl_mat_set(&A_bzbz,1,1,temp); 
  temp = 0.6324555320336759*magB2[1]; 
  gkyl_mat_set(&A_bxbx,1,2,temp); 
  gkyl_mat_set(&A_bxby,1,2,temp); 
  gkyl_mat_set(&A_bxbz,1,2,temp); 
  gkyl_mat_set(&A_byby,1,2,temp); 
  gkyl_mat_set(&A_bybz,1,2,temp); 
  gkyl_mat_set(&A_bzbz,1,2,temp); 
  temp = 0.7071067811865475*magB2[2]; 
  gkyl_mat_set(&A_bxbx,2,0,temp); 
  gkyl_mat_set(&A_bxby,2,0,temp); 
  gkyl_mat_set(&A_bxbz,2,0,temp); 
  gkyl_mat_set(&A_byby,2,0,temp); 
  gkyl_mat_set(&A_bybz,2,0,temp); 
  gkyl_mat_set(&A_bzbz,2,0,temp); 
  temp = 0.6324555320336759*magB2[1]; 
  gkyl_mat_set(&A_bxbx,2,1,temp); 
  gkyl_mat_set(&A_bxby,2,1,temp); 
  gkyl_mat_set(&A_bxbz,2,1,temp); 
  gkyl_mat_set(&A_byby,2,1,temp); 
  gkyl_mat_set(&A_bybz,2,1,temp); 
  gkyl_mat_set(&A_bzbz,2,1,temp); 
  temp = 0.4517539514526256*magB2[2]+0.7071067811865475*magB2[0]; 
  gkyl_mat_set(&A_bxbx,2,2,temp); 
  gkyl_mat_set(&A_bxby,2,2,temp); 
  gkyl_mat_set(&A_bxbz,2,2,temp); 
  gkyl_mat_set(&A_byby,2,2,temp); 
  gkyl_mat_set(&A_bybz,2,2,temp); 
  gkyl_mat_set(&A_bzbz,2,2,temp); 
  return cell_avg;
} 
