#include <gkyl_binop_cross_mul_ser.h> 
 
GKYL_CU_DH
void
binop_cross_mul_accumulate_1d_2d_ser_p1(double a, const double *f, const double *g, double *fg) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  fg[0] += a*(0.7071067811865475*f[1]*g[1]+0.7071067811865475*f[0]*g[0]); 
  fg[1] += a*(0.7071067811865475*f[0]*g[1]+0.7071067811865475*g[0]*f[1]); 
  fg[2] += a*(0.7071067811865475*f[1]*g[3]+0.7071067811865475*f[0]*g[2]); 
  fg[3] += a*(0.7071067811865475*f[0]*g[3]+0.7071067811865475*f[1]*g[2]); 
} 
