#include "gkyl_deflate_geo_surf_kernels.h"

GKYL_CU_DH void deflate_geo_surfx_2x_ser_p1(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]; 
  fld_deflated[1] = 0.7071067811865475*fld[2]; 
 
}
GKYL_CU_DH void deflate_geo_surfx_1x_ser_p1(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.5*fld[0]; 
 
}
GKYL_CU_DH void deflate_geo_surfy_2x_ser_p1(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]; 
  fld_deflated[1] = 0.7071067811865475*fld[1]; 
 
}
