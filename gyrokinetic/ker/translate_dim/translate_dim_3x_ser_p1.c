#include <gkyl_translate_dim_kernels.h> 

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirx_lo(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]-1.224744871391589*fdo[1]; 
  ftar[1] = 0.7071067811865475*fdo[2]-1.224744871391589*fdo[4]; 
  ftar[2] = 0.7071067811865475*fdo[3]-1.224744871391589*fdo[5]; 
  ftar[3] = 0.7071067811865475*fdo[6]-1.224744871391589*fdo[7]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirx_mid(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[2]; 
  ftar[2] = 0.7071067811865475*fdo[3]; 
  ftar[3] = 0.7071067811865475*fdo[6]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirx_up(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 1.224744871391589*fdo[1]+0.7071067811865475*fdo[0]; 
  ftar[1] = 1.224744871391589*fdo[4]+0.7071067811865475*fdo[2]; 
  ftar[2] = 1.224744871391589*fdo[5]+0.7071067811865475*fdo[3]; 
  ftar[3] = 1.224744871391589*fdo[7]+0.7071067811865475*fdo[6]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_diry_lo(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]-1.224744871391589*fdo[2]; 
  ftar[1] = 0.7071067811865475*fdo[1]-1.224744871391589*fdo[4]; 
  ftar[2] = 0.7071067811865475*fdo[3]-1.224744871391589*fdo[6]; 
  ftar[3] = 0.7071067811865475*fdo[5]-1.224744871391589*fdo[7]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_diry_mid(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[1]; 
  ftar[2] = 0.7071067811865475*fdo[3]; 
  ftar[3] = 0.7071067811865475*fdo[5]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_diry_up(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 1.224744871391589*fdo[2]+0.7071067811865475*fdo[0]; 
  ftar[1] = 1.224744871391589*fdo[4]+0.7071067811865475*fdo[1]; 
  ftar[2] = 1.224744871391589*fdo[6]+0.7071067811865475*fdo[3]; 
  ftar[3] = 1.224744871391589*fdo[7]+0.7071067811865475*fdo[5]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirz_lo(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]-1.224744871391589*fdo[3]; 
  ftar[1] = 0.7071067811865475*fdo[1]-1.224744871391589*fdo[5]; 
  ftar[2] = 0.7071067811865475*fdo[2]-1.224744871391589*fdo[6]; 
  ftar[3] = 0.7071067811865475*fdo[4]-1.224744871391589*fdo[7]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirz_mid(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 0.7071067811865475*fdo[0]; 
  ftar[1] = 0.7071067811865475*fdo[1]; 
  ftar[2] = 0.7071067811865475*fdo[2]; 
  ftar[3] = 0.7071067811865475*fdo[4]; 
}

GKYL_CU_DH void translate_dim_3x_ser_p1_to_2x_p1_dirz_up(const double *fdo, double *ftar) 
{ 
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  ftar[0] = 1.224744871391589*fdo[3]+0.7071067811865475*fdo[0]; 
  ftar[1] = 1.224744871391589*fdo[5]+0.7071067811865475*fdo[1]; 
  ftar[2] = 1.224744871391589*fdo[6]+0.7071067811865475*fdo[2]; 
  ftar[3] = 1.224744871391589*fdo[7]+0.7071067811865475*fdo[4]; 
}

