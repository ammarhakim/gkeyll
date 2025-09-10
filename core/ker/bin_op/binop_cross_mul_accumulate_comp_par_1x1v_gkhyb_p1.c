#include <gkyl_binop_cross_mul_gkhyb.h> 
 
GKYL_CU_DH
void
binop_cross_mul_accumulate_comp_par_1x1v_gkhyb_p1(double a, const double *f, const double *g, double *fg, int linc2) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  switch (linc2) { 
    case 0: 
      fg[0] += a*(0.7071067811865475*f[1]*g[1]+0.7071067811865475*f[0]*g[0]); 
      break; 
    case 1: 
      fg[1] += a*(0.7071067811865475*f[0]*g[1]+0.7071067811865475*g[0]*f[1]); 
      break; 
    case 2: 
      fg[2] += a*(0.7071067811865475*f[1]*g[3]+0.7071067811865475*f[0]*g[2]); 
      break; 
    case 3: 
      fg[3] += a*(0.7071067811865475*f[0]*g[3]+0.7071067811865475*f[1]*g[2]); 
      break; 
    case 4: 
      fg[4] += a*(0.7071067811865475*f[1]*g[5]+0.7071067811865475*f[0]*g[4]); 
      break; 
    case 5: 
      fg[5] += a*(0.7071067811865475*f[0]*g[5]+0.7071067811865475*f[1]*g[4]); 
      break; 
  } 
} 
