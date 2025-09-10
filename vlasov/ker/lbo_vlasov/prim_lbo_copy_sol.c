#include <gkyl_mat.h> 
#include <gkyl_prim_lbo_kernels.h>
 
GKYL_CU_DH void prim_lbo_copy_sol(const struct gkyl_mat *rhs, const int nc, const int udim, double* GKYL_RESTRICT out)
{
  for(size_t i=0; i<(udim+1)*nc; i++) 
    out[i] = gkyl_mat_get(rhs, i, 0);
}
