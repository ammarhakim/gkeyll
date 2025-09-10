#include <gkyl_binop_cross_mul_gkhyb.h> 
 
GKYL_CU_DH
void
binop_cross_mul_1x1v_gkhyb_p1(const double *f, const double *g, double *fg) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  double tmp[6] = {0.0}; 
  tmp[0] = 0.7071067811865475*f[1]*g[1]+0.7071067811865475*f[0]*g[0]; 
  tmp[1] = 0.7071067811865475*f[0]*g[1]+0.7071067811865475*g[0]*f[1]; 
  tmp[2] = 0.7071067811865475*f[1]*g[3]+0.7071067811865475*f[0]*g[2]; 
  tmp[3] = 0.7071067811865475*f[0]*g[3]+0.7071067811865475*f[1]*g[2]; 
  tmp[4] = 0.7071067811865475*f[1]*g[5]+0.7071067811865475*f[0]*g[4]; 
  tmp[5] = 0.7071067811865475*f[0]*g[5]+0.7071067811865475*f[1]*g[4]; 
 
  fg[0] = tmp[0]; 
  fg[1] = tmp[1]; 
  fg[2] = tmp[2]; 
  fg[3] = tmp[3]; 
  fg[4] = tmp[4]; 
  fg[5] = tmp[5]; 
} 
