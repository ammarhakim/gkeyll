#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_calc_BB_1x_ser_p1(const double *em, double* GKYL_RESTRICT out) 
{ 
  // em:  Input electromagnetic fields. 
  // out: Output B_i B_j tensor. 
 
  double *B_x_sq  = &out[0]; 
  double *B_x_B_y = &out[2]; 
  double *B_x_B_z = &out[4]; 
  double *B_y_sq  = &out[6]; 
  double *B_y_B_z = &out[8]; 
  double *B_z_sq  = &out[10]; 
 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  // Calculate B_i B_j. 
  binop_mul_1d_ser_p1(B_x, B_x, B_x_sq); 
  binop_mul_1d_ser_p1(B_x, B_y, B_x_B_y); 
  binop_mul_1d_ser_p1(B_x, B_z, B_x_B_z); 
  binop_mul_1d_ser_p1(B_y, B_y, B_y_sq); 
  binop_mul_1d_ser_p1(B_y, B_z, B_y_B_z); 
  binop_mul_1d_ser_p1(B_z, B_z, B_z_sq); 
 
} 
 
