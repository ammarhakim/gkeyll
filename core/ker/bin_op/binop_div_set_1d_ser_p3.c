#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
GKYL_CU_DH void binop_div_set_1d_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g) 
{ 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // f:       numerator field (must be a scalar). 
  // g:       denominator field (must be a scalar). 
 
  // If a corner value is below zero, use cell average g.
  bool avgg = false;
  if ((-1.870828693386971*g[3])+1.58113883008419*g[2]-1.224744871391589*g[1]+0.7071067811865475*g[0] < 0.0) { 
    avgg = true;
  }
  if (1.870828693386971*g[3]+1.58113883008419*g[2]+1.224744871391589*g[1]+0.7071067811865475*g[0] < 0.0) { 
    avgg = true;
  }
 
  double lhs[4]; 
  if (avgg) { 
    lhs[0] = g[0]; 
    lhs[1] = 0.0; 
    lhs[2] = 0.0; 
    lhs[3] = 0.0; 
    gkyl_mat_set(rhs,0,0,f[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
  } else { 
    lhs[0] = g[0]; 
    lhs[1] = g[1]; 
    lhs[2] = g[2]; 
    lhs[3] = g[3]; 
    gkyl_mat_set(rhs,0,0,f[0]); 
    gkyl_mat_set(rhs,1,0,f[1]); 
    gkyl_mat_set(rhs,2,0,f[2]); 
    gkyl_mat_set(rhs,3,0,f[3]); 
  } 
 
  // Fill LHS matrix. 
  gkyl_mat_set(A,0,0,0.7071067811865475*lhs[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*lhs[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*lhs[2]); 
  gkyl_mat_set(A,0,3,0.7071067811865475*lhs[3]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*lhs[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*lhs[2]+0.7071067811865475*lhs[0]); 
  gkyl_mat_set(A,1,2,0.6210590034081186*lhs[3]+0.6324555320336759*lhs[1]); 
  gkyl_mat_set(A,1,3,0.6210590034081186*lhs[2]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*lhs[2]); 
  gkyl_mat_set(A,2,1,0.6210590034081186*lhs[3]+0.6324555320336759*lhs[1]); 
  gkyl_mat_set(A,2,2,0.4517539514526256*lhs[2]+0.7071067811865475*lhs[0]); 
  gkyl_mat_set(A,2,3,0.421637021355784*lhs[3]+0.6210590034081186*lhs[1]); 
  gkyl_mat_set(A,3,0,0.7071067811865475*lhs[3]); 
  gkyl_mat_set(A,3,1,0.6210590034081186*lhs[2]); 
  gkyl_mat_set(A,3,2,0.421637021355784*lhs[3]+0.6210590034081186*lhs[1]); 
  gkyl_mat_set(A,3,3,0.421637021355784*lhs[2]+0.7071067811865475*lhs[0]); 

} 
