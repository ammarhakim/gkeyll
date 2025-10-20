#include <gkyl_translate_dim_kernels.h> 

GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_dirx_lo(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]-1.224744871391589*fdo[1]; 
  ftar[1] = 0.7071067811865475*fdo[2]-1.224744871391589*fdo[3]; 
}

GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_dirx_mid(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[2]; 
}

GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_dirx_up(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 1.224744871391589*fdo[1]+0.7071067811865475*fdo[0]; 
  ftar[1] = 1.224744871391589*fdo[3]+0.7071067811865475*fdo[2]; 
}

GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_diry_lo(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]-1.224744871391589*fdo[2]; 
  ftar[1] = 0.7071067811865475*fdo[1]-1.224744871391589*fdo[3]; 
}

GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_diry_mid(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[1]; 
}

GKYL_CU_DH void translate_dim_2x_ser_p1_to_1x_p1_diry_up(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 1.224744871391589*fdo[2]+0.7071067811865475*fdo[0]; 
  ftar[1] = 1.224744871391589*fdo[3]+0.7071067811865475*fdo[1]; 
}

GKYL_CU_DH void translate_dim_2x_ser_p1_to_3x_p1(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 1.4142135623730951*fdo[0]; 
  ftar[1] = 1.4142135623730951*fdo[1]; 
  ftar[2] = 1.4142135623730951*fdo[2]; 
  ftar[3] = 0.0; 
  ftar[4] = 1.4142135623730951*fdo[3]; 
  ftar[5] = 0.0; 
  ftar[6] = 0.0; 
  ftar[7] = 0.0; 
}

