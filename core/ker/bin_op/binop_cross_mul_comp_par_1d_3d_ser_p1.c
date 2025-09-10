#include <gkyl_binop_cross_mul_ser.h> 
 
GKYL_CU_DH
void
binop_cross_mul_comp_par_1d_3d_ser_p1(const double *f, const double *g, double *fg, int linc2) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  switch (linc2) { 
    case 0: 
      fg[0] += 0.7071067811865475*f[1]*g[1]+0.7071067811865475*f[0]*g[0]; 
    case 1: 
      fg[1] += 0.7071067811865475*f[0]*g[1]+0.7071067811865475*g[0]*f[1]; 
    case 2: 
      fg[2] += 0.7071067811865475*f[1]*g[4]+0.7071067811865475*f[0]*g[2]; 
    case 3: 
      fg[3] += 0.7071067811865475*f[1]*g[5]+0.7071067811865475*f[0]*g[3]; 
    case 4: 
      fg[4] += 0.7071067811865475*f[0]*g[4]+0.7071067811865475*f[1]*g[2]; 
    case 5: 
      fg[5] += 0.7071067811865475*f[0]*g[5]+0.7071067811865475*f[1]*g[3]; 
    case 6: 
      fg[6] += 0.7071067811865475*f[1]*g[7]+0.7071067811865475*f[0]*g[6]; 
    case 7: 
      fg[7] += 0.7071067811865475*f[0]*g[7]+0.7071067811865475*f[1]*g[6]; 
  } 
} 
