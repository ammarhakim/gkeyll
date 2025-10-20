#include <gkyl_deflate_surf_kernels.h>

GKYL_CU_DH void deflate_surfx_lower_1x_ser_p2(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.5811388300841895*fld[2]-1.224744871391589*fld[1]+0.7071067811865475*fld[0]; 
}
 
GKYL_CU_DH void deflate_surfx_upper_1x_ser_p2(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.5811388300841895*fld[2]+1.224744871391589*fld[1]+0.7071067811865475*fld[0]; 
}
 
