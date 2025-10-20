#include <gkyl_inflate_surf_kernels.h>

GKYL_CU_DH void inflate_surfx_2x_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[2] = 1.4142135623730951*fld_deflated[1]; 
}
 
GKYL_CU_DH void inflate_surfy_2x_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[1] = 1.4142135623730951*fld_deflated[1]; 
}
 
