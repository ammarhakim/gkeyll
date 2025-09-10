// Thu Jul 28 13:04:18 2022
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_1d_ser_p1(const double *f, const double *g, double *fg )
{
  double tmp[2] = {0.};
  tmp[0] =  7.0710678118654757e-01*f[1]*g[1]+7.0710678118654757e-01*f[0]*g[0];
  tmp[1] =  7.0710678118654757e-01*f[1]*g[0]+7.0710678118654757e-01*f[0]*g[1];
 
  fg[0] = tmp[0];
  fg[1] = tmp[1];
  // nsum = 2, nprod = 8
}

struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p1(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 2, .num_prod = 8 };
}
