#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
GKYL_CU_DH void binop_div_copy_sol(const struct gkyl_mat *rhs, double* GKYL_RESTRICT fdivg)
{
  for(size_t i=0; i<rhs->nr; ++i) 
    fdivg[i] = gkyl_mat_get(rhs,i,0); 
}