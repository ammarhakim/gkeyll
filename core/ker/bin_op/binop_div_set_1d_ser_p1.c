#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
GKYL_CU_DH void binop_div_set_1d_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g) 
{ 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // f:       numerator field (must be a scalar). 
  // g:       denominator field (must be a scalar). 
 
  // If a corner value is below zero, use cell average g.
  bool avgg = false;
  if (0.7071067811865475*g[0]-1.224744871391589*g[1] < 0.0) { 
    avgg = true;
  }
  if (1.224744871391589*g[1]+0.7071067811865475*g[0] < 0.0) { 
    avgg = true;
  }
 
  double lhs[2]; 
  if (avgg) { 
    lhs[0] = g[0]; 
    lhs[1] = 0.0; 
    gkyl_mat_set(rhs,0,0,f[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
  } else { 
    lhs[0] = g[0]; 
    lhs[1] = g[1]; 
    gkyl_mat_set(rhs,0,0,f[0]); 
    gkyl_mat_set(rhs,1,0,f[1]); 
  } 
 
  // Fill LHS matrix. 
  gkyl_mat_set(A,0,0,0.7071067811865475*lhs[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*lhs[1]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*lhs[1]); 
  gkyl_mat_set(A,1,1,0.7071067811865475*lhs[0]); 

} 
