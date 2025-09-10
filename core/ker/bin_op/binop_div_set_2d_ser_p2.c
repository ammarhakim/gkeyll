#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
GKYL_CU_DH void binop_div_set_2d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g) 
{ 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // f:       numerator field (must be a scalar). 
  // g:       denominator field (must be a scalar). 
 
  // If a corner value is below zero, use cell average g.
  bool avgg = false;
  if ((-1.936491673103709*(g[7]+g[6]))+1.118033988749895*(g[5]+g[4])+1.5*g[3]-0.8660254037844386*(g[2]+g[1])+0.5*g[0] < 0.0) { 
    avgg = true;
  }
  if (1.936491673103709*g[7]-1.936491673103709*g[6]+1.118033988749895*(g[5]+g[4])-1.5*g[3]-0.8660254037844386*g[2]+0.8660254037844386*g[1]+0.5*g[0] < 0.0) { 
    avgg = true;
  }
  if ((-1.936491673103709*g[7])+1.936491673103709*g[6]+1.118033988749895*(g[5]+g[4])-1.5*g[3]+0.8660254037844386*g[2]-0.8660254037844386*g[1]+0.5*g[0] < 0.0) { 
    avgg = true;
  }
  if (1.936491673103709*(g[7]+g[6])+1.118033988749895*(g[5]+g[4])+1.5*g[3]+0.8660254037844386*(g[2]+g[1])+0.5*g[0] < 0.0) { 
    avgg = true;
  }
 
  double lhs[8]; 
  if (avgg) { 
    lhs[0] = g[0]; 
    lhs[1] = 0.0; 
    lhs[2] = 0.0; 
    lhs[3] = 0.0; 
    lhs[4] = 0.0; 
    lhs[5] = 0.0; 
    lhs[6] = 0.0; 
    lhs[7] = 0.0; 
    gkyl_mat_set(rhs,0,0,f[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
  } else { 
    lhs[0] = g[0]; 
    lhs[1] = g[1]; 
    lhs[2] = g[2]; 
    lhs[3] = g[3]; 
    lhs[4] = g[4]; 
    lhs[5] = g[5]; 
    lhs[6] = g[6]; 
    lhs[7] = g[7]; 
    gkyl_mat_set(rhs,0,0,f[0]); 
    gkyl_mat_set(rhs,1,0,f[1]); 
    gkyl_mat_set(rhs,2,0,f[2]); 
    gkyl_mat_set(rhs,3,0,f[3]); 
    gkyl_mat_set(rhs,4,0,f[4]); 
    gkyl_mat_set(rhs,5,0,f[5]); 
    gkyl_mat_set(rhs,6,0,f[6]); 
    gkyl_mat_set(rhs,7,0,f[7]); 
  } 
 
  // Fill LHS matrix. 
  gkyl_mat_set(A,0,0,0.5*lhs[0]); 
  gkyl_mat_set(A,0,1,0.5*lhs[1]); 
  gkyl_mat_set(A,0,2,0.5*lhs[2]); 
  gkyl_mat_set(A,0,3,0.5*lhs[3]); 
  gkyl_mat_set(A,0,4,0.5*lhs[4]); 
  gkyl_mat_set(A,0,5,0.5*lhs[5]); 
  gkyl_mat_set(A,0,6,0.5*lhs[6]); 
  gkyl_mat_set(A,0,7,0.5*lhs[7]); 
  gkyl_mat_set(A,1,0,0.5*lhs[1]); 
  gkyl_mat_set(A,1,1,0.4472135954999579*lhs[4]+0.5*lhs[0]); 
  gkyl_mat_set(A,1,2,0.5*lhs[3]); 
  gkyl_mat_set(A,1,3,0.447213595499958*lhs[6]+0.5*lhs[2]); 
  gkyl_mat_set(A,1,4,0.4472135954999579*lhs[1]); 
  gkyl_mat_set(A,1,5,0.5000000000000001*lhs[7]); 
  gkyl_mat_set(A,1,6,0.447213595499958*lhs[3]); 
  gkyl_mat_set(A,1,7,0.5000000000000001*lhs[5]); 
  gkyl_mat_set(A,2,0,0.5*lhs[2]); 
  gkyl_mat_set(A,2,1,0.5*lhs[3]); 
  gkyl_mat_set(A,2,2,0.4472135954999579*lhs[5]+0.5*lhs[0]); 
  gkyl_mat_set(A,2,3,0.447213595499958*lhs[7]+0.5*lhs[1]); 
  gkyl_mat_set(A,2,4,0.5000000000000001*lhs[6]); 
  gkyl_mat_set(A,2,5,0.4472135954999579*lhs[2]); 
  gkyl_mat_set(A,2,6,0.5000000000000001*lhs[4]); 
  gkyl_mat_set(A,2,7,0.447213595499958*lhs[3]); 
  gkyl_mat_set(A,3,0,0.5*lhs[3]); 
  gkyl_mat_set(A,3,1,0.447213595499958*lhs[6]+0.5*lhs[2]); 
  gkyl_mat_set(A,3,2,0.447213595499958*lhs[7]+0.5*lhs[1]); 
  gkyl_mat_set(A,3,3,0.4472135954999579*lhs[5]+0.4472135954999579*lhs[4]+0.5*lhs[0]); 
  gkyl_mat_set(A,3,4,0.4472135954999579*lhs[3]); 
  gkyl_mat_set(A,3,5,0.4472135954999579*lhs[3]); 
  gkyl_mat_set(A,3,6,0.4*lhs[7]+0.447213595499958*lhs[1]); 
  gkyl_mat_set(A,3,7,0.4*lhs[6]+0.447213595499958*lhs[2]); 
  gkyl_mat_set(A,4,0,0.5*lhs[4]); 
  gkyl_mat_set(A,4,1,0.4472135954999579*lhs[1]); 
  gkyl_mat_set(A,4,2,0.5000000000000001*lhs[6]); 
  gkyl_mat_set(A,4,3,0.4472135954999579*lhs[3]); 
  gkyl_mat_set(A,4,4,0.31943828249997*lhs[4]+0.5*lhs[0]); 
  gkyl_mat_set(A,4,6,0.31943828249997*lhs[6]+0.5000000000000001*lhs[2]); 
  gkyl_mat_set(A,4,7,0.4472135954999579*lhs[7]); 
  gkyl_mat_set(A,5,0,0.5*lhs[5]); 
  gkyl_mat_set(A,5,1,0.5000000000000001*lhs[7]); 
  gkyl_mat_set(A,5,2,0.4472135954999579*lhs[2]); 
  gkyl_mat_set(A,5,3,0.4472135954999579*lhs[3]); 
  gkyl_mat_set(A,5,5,0.31943828249997*lhs[5]+0.5*lhs[0]); 
  gkyl_mat_set(A,5,6,0.4472135954999579*lhs[6]); 
  gkyl_mat_set(A,5,7,0.31943828249997*lhs[7]+0.5000000000000001*lhs[1]); 
  gkyl_mat_set(A,6,0,0.5*lhs[6]); 
  gkyl_mat_set(A,6,1,0.447213595499958*lhs[3]); 
  gkyl_mat_set(A,6,2,0.5000000000000001*lhs[4]); 
  gkyl_mat_set(A,6,3,0.4*lhs[7]+0.447213595499958*lhs[1]); 
  gkyl_mat_set(A,6,4,0.31943828249997*lhs[6]+0.5000000000000001*lhs[2]); 
  gkyl_mat_set(A,6,5,0.4472135954999579*lhs[6]); 
  gkyl_mat_set(A,6,6,0.4472135954999579*lhs[5]+0.31943828249997*lhs[4]+0.5*lhs[0]); 
  gkyl_mat_set(A,6,7,0.4*lhs[3]); 
  gkyl_mat_set(A,7,0,0.5*lhs[7]); 
  gkyl_mat_set(A,7,1,0.5000000000000001*lhs[5]); 
  gkyl_mat_set(A,7,2,0.447213595499958*lhs[3]); 
  gkyl_mat_set(A,7,3,0.4*lhs[6]+0.447213595499958*lhs[2]); 
  gkyl_mat_set(A,7,4,0.4472135954999579*lhs[7]); 
  gkyl_mat_set(A,7,5,0.31943828249997*lhs[7]+0.5000000000000001*lhs[1]); 
  gkyl_mat_set(A,7,6,0.4*lhs[3]); 
  gkyl_mat_set(A,7,7,0.31943828249997*lhs[5]+0.4472135954999579*lhs[4]+0.5*lhs[0]); 

} 
