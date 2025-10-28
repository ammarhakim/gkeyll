#include <gkyl_translate_dim_kernels.h> 

GKYL_CU_DH void translate_dim_1x_ser_p1_to_2x_p1(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 1.4142135623730951*fdo[0]; 
  ftar[1] = 1.4142135623730951*fdo[1]; 
  ftar[2] = 0.0; 
  ftar[3] = 0.0; 
}

