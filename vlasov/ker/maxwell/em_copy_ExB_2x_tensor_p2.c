#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_copy_ExB_2x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
  double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf) 
{ 
  // count:          Integer to indicate which matrix being fetched. 
  // x:              Input solution vector. 
  // em:             Input electromagnetic fields. 
  // cell_avg_magB2: Input flag for cell average if 1/|B|^2 only used cell averages. 
  // ExB:            E x B velocity = E x B/|B|^2. 
 
  struct gkyl_mat x_ExBx = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_ExBy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_ExBz = gkyl_nmat_get(x, count+2); 
  double *ExB_x = &ExB[0]; 
  double *ExB_y = &ExB[9]; 
  double *ExB_z = &ExB[18]; 
 
  ExB_x[0] = gkyl_mat_get(&x_ExBx,0,0); 
  ExB_y[0] = gkyl_mat_get(&x_ExBy,0,0); 
  ExB_z[0] = gkyl_mat_get(&x_ExBz,0,0); 
  ExB_x[1] = gkyl_mat_get(&x_ExBx,1,0); 
  ExB_y[1] = gkyl_mat_get(&x_ExBy,1,0); 
  ExB_z[1] = gkyl_mat_get(&x_ExBz,1,0); 
  ExB_x[2] = gkyl_mat_get(&x_ExBx,2,0); 
  ExB_y[2] = gkyl_mat_get(&x_ExBy,2,0); 
  ExB_z[2] = gkyl_mat_get(&x_ExBz,2,0); 
  ExB_x[3] = gkyl_mat_get(&x_ExBx,3,0); 
  ExB_y[3] = gkyl_mat_get(&x_ExBy,3,0); 
  ExB_z[3] = gkyl_mat_get(&x_ExBz,3,0); 
  ExB_x[4] = gkyl_mat_get(&x_ExBx,4,0); 
  ExB_y[4] = gkyl_mat_get(&x_ExBy,4,0); 
  ExB_z[4] = gkyl_mat_get(&x_ExBz,4,0); 
  ExB_x[5] = gkyl_mat_get(&x_ExBx,5,0); 
  ExB_y[5] = gkyl_mat_get(&x_ExBy,5,0); 
  ExB_z[5] = gkyl_mat_get(&x_ExBz,5,0); 
  ExB_x[6] = gkyl_mat_get(&x_ExBx,6,0); 
  ExB_y[6] = gkyl_mat_get(&x_ExBy,6,0); 
  ExB_z[6] = gkyl_mat_get(&x_ExBz,6,0); 
  ExB_x[7] = gkyl_mat_get(&x_ExBx,7,0); 
  ExB_y[7] = gkyl_mat_get(&x_ExBy,7,0); 
  ExB_z[7] = gkyl_mat_get(&x_ExBz,7,0); 
  ExB_x[8] = gkyl_mat_get(&x_ExBx,8,0); 
  ExB_y[8] = gkyl_mat_get(&x_ExBy,8,0); 
  ExB_z[8] = gkyl_mat_get(&x_ExBz,8,0); 
} 
 
