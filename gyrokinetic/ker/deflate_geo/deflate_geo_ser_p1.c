#include "gkyl_deflate_geo_kernels.h"

GKYL_CU_DH void deflate_geo_2x_ser_p1_remx(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]; 
  fld_deflated[1] = 0.7071067811865475*fld[2]; 
  fld_deflated[2] = 0.7071067811865475*fld[3]; 
  fld_deflated[3] = 0.7071067811865475*fld[6]; 
 
}
GKYL_CU_DH void deflate_geo_1x_ser_p1_remxy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.5*fld[0]; 
  fld_deflated[1] = 0.5*fld[3]; 
 
}
GKYL_CU_DH void deflate_geo_2x_ser_p1_remy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]; 
  fld_deflated[1] = 0.7071067811865475*fld[1]; 
  fld_deflated[2] = 0.7071067811865475*fld[3]; 
  fld_deflated[3] = 0.7071067811865475*fld[5]; 
 
}
