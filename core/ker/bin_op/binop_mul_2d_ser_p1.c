// Thu Jul 28 13:04:18 2022
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_2d_ser_p1(const double *f, const double *g, double *fg )
{
  double tmp[4] = {0.};
  tmp[0] =  5.0000000000000000e-01*f[3]*g[3]+5.0000000000000000e-01*g[0]*f[0]+5.0000000000000000e-01*f[2]*g[2]+5.0000000000000000e-01*g[1]*f[1];
  tmp[1] =  5.0000000000000000e-01*f[3]*g[2]+5.0000000000000000e-01*g[0]*f[1]+5.0000000000000000e-01*f[2]*g[3]+5.0000000000000000e-01*g[1]*f[0];
  tmp[2] =  5.0000000000000000e-01*g[2]*f[0]+5.0000000000000000e-01*g[1]*f[3]+5.0000000000000000e-01*f[2]*g[0]+5.0000000000000000e-01*g[3]*f[1];
  tmp[3] =  5.0000000000000000e-01*g[0]*f[3]+5.0000000000000000e-01*g[2]*f[1]+5.0000000000000000e-01*f[2]*g[1]+5.0000000000000000e-01*g[3]*f[0];
 
  fg[0] = tmp[0];
  fg[1] = tmp[1];
  fg[2] = tmp[2];
  fg[3] = tmp[3];
  // nsum = 12, nprod = 32
}

struct gkyl_kern_op_count op_count_binop_mul_2d_ser_p1(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 12, .num_prod = 32 };
}
