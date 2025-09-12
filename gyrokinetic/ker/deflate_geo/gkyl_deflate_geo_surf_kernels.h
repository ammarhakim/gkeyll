#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void deflate_geo_surfx_2x_ser_p1(const double *fld, double *fld_deflated);
GKYL_CU_DH void deflate_geo_surfx_1x_ser_p1(const double *fld, double *fld_deflated);
GKYL_CU_DH void deflate_geo_surfy_2x_ser_p1(const double *fld, double *fld_deflated);

EXTERN_C_END 
