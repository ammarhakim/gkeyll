#include "gkyl_deflate_geo_kernels.h"

GKYL_CU_DH void deflate_geo_2x_ser_p2_remx(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]-0.7905694150420947*fld[7]; 
  fld_deflated[1] = 0.7071067811865475*fld[2]-0.7905694150420948*fld[11]; 
  fld_deflated[2] = 0.7071067811865475*fld[3]-0.7905694150420948*fld[13]; 
  fld_deflated[3] = 0.7071067811865475*fld[6]-0.7905694150420947*fld[17]; 
  fld_deflated[4] = 0.7071067811865475*fld[8]; 
  fld_deflated[5] = 0.7071067811865475*fld[9]; 
  fld_deflated[6] = 0.7071067811865475*fld[14]; 
  fld_deflated[7] = 0.7071067811865475*fld[16]; 
 
}
GKYL_CU_DH void deflate_geo_1x_ser_p2_remxy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = (-0.5590169943749475*fld[8])-0.5590169943749475*fld[7]+0.5*fld[0]; 
  fld_deflated[1] = (-0.5590169943749476*fld[14])-0.5590169943749476*fld[13]+0.5*fld[3]; 
  fld_deflated[2] = 0.5*fld[9]; 
 
}
GKYL_CU_DH void deflate_geo_2x_ser_p2_remy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]-0.7905694150420947*fld[8]; 
  fld_deflated[1] = 0.7071067811865475*fld[1]-0.7905694150420948*fld[12]; 
  fld_deflated[2] = 0.7071067811865475*fld[3]-0.7905694150420948*fld[14]; 
  fld_deflated[3] = 0.7071067811865475*fld[5]-0.7905694150420947*fld[18]; 
  fld_deflated[4] = 0.7071067811865475*fld[7]; 
  fld_deflated[5] = 0.7071067811865475*fld[9]; 
  fld_deflated[6] = 0.7071067811865475*fld[13]; 
  fld_deflated[7] = 0.7071067811865475*fld[15]; 
 
}
