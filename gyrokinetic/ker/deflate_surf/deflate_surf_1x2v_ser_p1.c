#include <gkyl_deflate_surf_kernels.h>

GKYL_CU_DH void deflate_surfx_lower_1x2v_ser_p1(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]-1.224744871391589*fld[1]; 
  fld_deflated[1] = 0.7071067811865475*fld[2]-1.224744871391589*fld[4]; 
  fld_deflated[2] = 0.7071067811865475*fld[3]-1.224744871391589*fld[5]; 
  fld_deflated[3] = 0.7071067811865475*fld[6]-1.224744871391589*fld[7]; 
  fld_deflated[4] = 0.7071067811865475*fld[8]-1.224744871391589*fld[9]; 
  fld_deflated[5] = 0.7071067811865475*fld[10]-1.224744871391589*fld[11]; 
}
 
GKYL_CU_DH void deflate_surfx_upper_1x2v_ser_p1(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.224744871391589*fld[1]+0.7071067811865475*fld[0]; 
  fld_deflated[1] = 1.224744871391589*fld[4]+0.7071067811865475*fld[2]; 
  fld_deflated[2] = 1.224744871391589*fld[5]+0.7071067811865475*fld[3]; 
  fld_deflated[3] = 1.224744871391589*fld[7]+0.7071067811865475*fld[6]; 
  fld_deflated[4] = 1.224744871391589*fld[9]+0.7071067811865475*fld[8]; 
  fld_deflated[5] = 1.224744871391589*fld[11]+0.7071067811865475*fld[10]; 
}
 
