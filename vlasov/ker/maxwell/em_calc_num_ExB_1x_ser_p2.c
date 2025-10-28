#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_calc_num_ExB_1x_ser_p2(const double *em, double* GKYL_RESTRICT out) 
{ 
  // em:  Input electromagnetic fields. 
  // out: Output (E x B)_i (E x B velocity numerator) and B_i^2 (for denominator). 
 
  double *num_ExB_x = &out[0]; 
  double *num_ExB_y = &out[3]; 
  double *num_ExB_z = &out[6]; 
  double *B_x_sq  = &out[9]; 
  double *B_y_sq  = &out[12]; 
  double *B_z_sq  = &out[15]; 
 
  const double *E_x = &em[0]; 
  const double *E_y = &em[3]; 
  const double *E_z = &em[6]; 
  const double *B_x = &em[9]; 
  const double *B_y = &em[12]; 
  const double *B_z = &em[15]; 
 
  // Calculate E_i B_j. 
  double E_x_B_y[3] = {0.0}; 
  binop_mul_1d_ser_p2(E_x, B_y, E_x_B_y); 
 
  double E_x_B_z[3] = {0.0}; 
  binop_mul_1d_ser_p2(E_x, B_z, E_x_B_z); 
 
  double E_y_B_x[3] = {0.0}; 
  binop_mul_1d_ser_p2(E_y, B_x, E_y_B_x); 
 
  double E_y_B_z[3] = {0.0}; 
  binop_mul_1d_ser_p2(E_y, B_z, E_y_B_z); 
 
  double E_z_B_x[3] = {0.0}; 
  binop_mul_1d_ser_p2(E_z, B_x, E_z_B_x); 
 
  double E_z_B_y[3] = {0.0}; 
  binop_mul_1d_ser_p2(E_z, B_y, E_z_B_y); 
 
  num_ExB_x[0] = E_y_B_z[0] - E_z_B_y[0]; 
  num_ExB_y[0] = E_z_B_x[0] - E_x_B_z[0]; 
  num_ExB_z[0] = E_x_B_y[0] - E_y_B_x[0]; 
  num_ExB_x[1] = E_y_B_z[1] - E_z_B_y[1]; 
  num_ExB_y[1] = E_z_B_x[1] - E_x_B_z[1]; 
  num_ExB_z[1] = E_x_B_y[1] - E_y_B_x[1]; 
  num_ExB_x[2] = E_y_B_z[2] - E_z_B_y[2]; 
  num_ExB_y[2] = E_z_B_x[2] - E_x_B_z[2]; 
  num_ExB_z[2] = E_x_B_y[2] - E_y_B_x[2]; 
 
  // Calculate B_i^2. 
  binop_mul_1d_ser_p2(B_x, B_x, B_x_sq); 
  binop_mul_1d_ser_p2(B_y, B_y, B_y_sq); 
  binop_mul_1d_ser_p2(B_z, B_z, B_z_sq); 
 
} 
 
