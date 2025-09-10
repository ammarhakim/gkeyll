#include <gkyl_binop_cross_mul_ser.h> 
 
GKYL_CU_DH
void
binop_cross_mul_2d_3d_ser_p1(const double *f, const double *g, double *fg) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  double tmp[8] = {0.0}; 
  tmp[0] = 0.5*f[3]*g[4]+0.5*f[2]*g[2]+0.5*f[1]*g[1]+0.5*f[0]*g[0]; 
  tmp[1] = 0.5*f[2]*g[4]+0.5*g[2]*f[3]+0.5*f[0]*g[1]+0.5*g[0]*f[1]; 
  tmp[2] = 0.5*f[1]*g[4]+0.5*g[1]*f[3]+0.5*f[0]*g[2]+0.5*g[0]*f[2]; 
  tmp[3] = 0.5*f[3]*g[7]+0.5*f[2]*g[6]+0.5*f[1]*g[5]+0.5*f[0]*g[3]; 
  tmp[4] = 0.5*f[0]*g[4]+0.5*g[0]*f[3]+0.5*f[1]*g[2]+0.5*g[1]*f[2]; 
  tmp[5] = 0.5*f[2]*g[7]+0.5*f[3]*g[6]+0.5*f[0]*g[5]+0.5*f[1]*g[3]; 
  tmp[6] = 0.5*f[1]*g[7]+0.5*f[0]*g[6]+0.5*f[3]*g[5]+0.5*f[2]*g[3]; 
  tmp[7] = 0.5*f[0]*g[7]+0.5*f[1]*g[6]+0.5*f[2]*g[5]+0.5*f[3]*g[3]; 
 
  fg[0] = tmp[0]; 
  fg[1] = tmp[1]; 
  fg[2] = tmp[2]; 
  fg[3] = tmp[3]; 
  fg[4] = tmp[4]; 
  fg[5] = tmp[5]; 
  fg[6] = tmp[6]; 
  fg[7] = tmp[7]; 
} 
