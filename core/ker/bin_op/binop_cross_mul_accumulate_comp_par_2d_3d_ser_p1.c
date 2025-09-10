#include <gkyl_binop_cross_mul_ser.h> 
 
GKYL_CU_DH
void
binop_cross_mul_accumulate_comp_par_2d_3d_ser_p1(double a, const double *f, const double *g, double *fg, int linc2) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  switch (linc2) { 
    case 0: 
      fg[0] += a*(0.5*f[3]*g[4]+0.5*f[2]*g[2]+0.5*f[1]*g[1]+0.5*f[0]*g[0]); 
      break; 
    case 1: 
      fg[1] += a*(0.5*f[2]*g[4]+0.5*g[2]*f[3]+0.5*f[0]*g[1]+0.5*g[0]*f[1]); 
      break; 
    case 2: 
      fg[2] += a*(0.5*f[1]*g[4]+0.5*g[1]*f[3]+0.5*f[0]*g[2]+0.5*g[0]*f[2]); 
      break; 
    case 3: 
      fg[3] += a*(0.5*f[3]*g[7]+0.5*f[2]*g[6]+0.5*f[1]*g[5]+0.5*f[0]*g[3]); 
      break; 
    case 4: 
      fg[4] += a*(0.5*f[0]*g[4]+0.5*g[0]*f[3]+0.5*f[1]*g[2]+0.5*g[1]*f[2]); 
      break; 
    case 5: 
      fg[5] += a*(0.5*f[2]*g[7]+0.5*f[3]*g[6]+0.5*f[0]*g[5]+0.5*f[1]*g[3]); 
      break; 
    case 6: 
      fg[6] += a*(0.5*f[1]*g[7]+0.5*f[0]*g[6]+0.5*f[3]*g[5]+0.5*f[2]*g[3]); 
      break; 
    case 7: 
      fg[7] += a*(0.5*f[0]*g[7]+0.5*f[1]*g[6]+0.5*f[2]*g[5]+0.5*f[3]*g[3]); 
      break; 
  } 
} 
