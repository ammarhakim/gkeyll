#include <gkyl_inflate_surf_kernels.h>

GKYL_CU_DH void inflate_surfx_3x2v_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[2] = 1.4142135623730951*fld_deflated[1]; 
  fld[3] = 1.4142135623730951*fld_deflated[2]; 
  fld[4] = 1.4142135623730951*fld_deflated[3]; 
  fld[5] = 1.4142135623730951*fld_deflated[4]; 
  fld[8] = 1.4142135623730951*fld_deflated[5]; 
  fld[10] = 1.4142135623730951*fld_deflated[6]; 
  fld[11] = 1.4142135623730951*fld_deflated[7]; 
  fld[13] = 1.4142135623730951*fld_deflated[8]; 
  fld[14] = 1.4142135623730951*fld_deflated[9]; 
  fld[15] = 1.4142135623730951*fld_deflated[10]; 
  fld[19] = 1.4142135623730951*fld_deflated[11]; 
  fld[22] = 1.4142135623730951*fld_deflated[12]; 
  fld[24] = 1.4142135623730951*fld_deflated[13]; 
  fld[25] = 1.4142135623730951*fld_deflated[14]; 
  fld[30] = 1.4142135623730951*fld_deflated[15]; 
  fld[32] = 1.4142135623730951*fld_deflated[16]; 
  fld[34] = 1.4142135623730951*fld_deflated[17]; 
  fld[35] = 1.4142135623730951*fld_deflated[18]; 
  fld[36] = 1.4142135623730951*fld_deflated[19]; 
  fld[39] = 1.4142135623730951*fld_deflated[20]; 
  fld[41] = 1.4142135623730951*fld_deflated[21]; 
  fld[42] = 1.4142135623730951*fld_deflated[22]; 
  fld[46] = 1.4142135623730951*fld_deflated[23]; 
}
 
GKYL_CU_DH void inflate_surfy_3x2v_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[1] = 1.4142135623730951*fld_deflated[1]; 
  fld[3] = 1.4142135623730951*fld_deflated[2]; 
  fld[4] = 1.4142135623730951*fld_deflated[3]; 
  fld[5] = 1.4142135623730951*fld_deflated[4]; 
  fld[7] = 1.4142135623730951*fld_deflated[5]; 
  fld[9] = 1.4142135623730951*fld_deflated[6]; 
  fld[11] = 1.4142135623730951*fld_deflated[7]; 
  fld[12] = 1.4142135623730951*fld_deflated[8]; 
  fld[14] = 1.4142135623730951*fld_deflated[9]; 
  fld[15] = 1.4142135623730951*fld_deflated[10]; 
  fld[18] = 1.4142135623730951*fld_deflated[11]; 
  fld[21] = 1.4142135623730951*fld_deflated[12]; 
  fld[23] = 1.4142135623730951*fld_deflated[13]; 
  fld[25] = 1.4142135623730951*fld_deflated[14]; 
  fld[29] = 1.4142135623730951*fld_deflated[15]; 
  fld[32] = 1.4142135623730951*fld_deflated[16]; 
  fld[33] = 1.4142135623730951*fld_deflated[17]; 
  fld[35] = 1.4142135623730951*fld_deflated[18]; 
  fld[36] = 1.4142135623730951*fld_deflated[19]; 
  fld[38] = 1.4142135623730951*fld_deflated[20]; 
  fld[40] = 1.4142135623730951*fld_deflated[21]; 
  fld[42] = 1.4142135623730951*fld_deflated[22]; 
  fld[45] = 1.4142135623730951*fld_deflated[23]; 
}
 
GKYL_CU_DH void inflate_surfz_3x2v_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[1] = 1.4142135623730951*fld_deflated[1]; 
  fld[2] = 1.4142135623730951*fld_deflated[2]; 
  fld[4] = 1.4142135623730951*fld_deflated[3]; 
  fld[5] = 1.4142135623730951*fld_deflated[4]; 
  fld[6] = 1.4142135623730951*fld_deflated[5]; 
  fld[9] = 1.4142135623730951*fld_deflated[6]; 
  fld[10] = 1.4142135623730951*fld_deflated[7]; 
  fld[12] = 1.4142135623730951*fld_deflated[8]; 
  fld[13] = 1.4142135623730951*fld_deflated[9]; 
  fld[15] = 1.4142135623730951*fld_deflated[10]; 
  fld[17] = 1.4142135623730951*fld_deflated[11]; 
  fld[20] = 1.4142135623730951*fld_deflated[12]; 
  fld[23] = 1.4142135623730951*fld_deflated[13]; 
  fld[24] = 1.4142135623730951*fld_deflated[14]; 
  fld[28] = 1.4142135623730951*fld_deflated[15]; 
  fld[32] = 1.4142135623730951*fld_deflated[16]; 
  fld[33] = 1.4142135623730951*fld_deflated[17]; 
  fld[34] = 1.4142135623730951*fld_deflated[18]; 
  fld[36] = 1.4142135623730951*fld_deflated[19]; 
  fld[37] = 1.4142135623730951*fld_deflated[20]; 
  fld[40] = 1.4142135623730951*fld_deflated[21]; 
  fld[41] = 1.4142135623730951*fld_deflated[22]; 
  fld[44] = 1.4142135623730951*fld_deflated[23]; 
}
 
