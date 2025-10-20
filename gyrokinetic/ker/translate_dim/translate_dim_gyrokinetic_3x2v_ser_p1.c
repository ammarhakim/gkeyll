#include <gkyl_translate_dim_kernels.h> 

GKYL_CU_DH void translate_dim_gyrokinetic_3x2v_ser_p1_from_1x2v_p1(const double *flow, double *fout) 
{ 
  // flow: lower dimensional field to get DG coefficients from.
  // fout: field whose DG coefficients to populate.

  fout[0] = 2.0*flow[0]; 
  fout[1] = 0.0; 
  fout[2] = 0.0; 
  fout[3] = 2.0*flow[1]; 
  fout[4] = 2.0*flow[2]; 
  fout[5] = 2.0*flow[3]; 
  fout[6] = 0.0; 
  fout[7] = 0.0; 
  fout[8] = 0.0; 
  fout[9] = 0.0; 
  fout[10] = 0.0; 
  fout[11] = 2.0*flow[4]; 
  fout[12] = 0.0; 
  fout[13] = 0.0; 
  fout[14] = 2.0*flow[5]; 
  fout[15] = 2.0*flow[6]; 
  fout[16] = 0.0; 
  fout[17] = 0.0; 
  fout[18] = 0.0; 
  fout[19] = 0.0; 
  fout[20] = 0.0; 
  fout[21] = 0.0; 
  fout[22] = 0.0; 
  fout[23] = 0.0; 
  fout[24] = 0.0; 
  fout[25] = 2.0*flow[7]; 
  fout[26] = 0.0; 
  fout[27] = 0.0; 
  fout[28] = 0.0; 
  fout[29] = 0.0; 
  fout[30] = 0.0; 
  fout[31] = 0.0; 
  fout[32] = 2.0*flow[8]; 
  fout[33] = 0.0; 
  fout[34] = 0.0; 
  fout[35] = 2.0*flow[9]; 
  fout[36] = 2.0*flow[10]; 
  fout[37] = 0.0; 
  fout[38] = 0.0; 
  fout[39] = 0.0; 
  fout[40] = 0.0; 
  fout[41] = 0.0; 
  fout[42] = 2.0*flow[11]; 
  fout[43] = 0.0; 
  fout[44] = 0.0; 
  fout[45] = 0.0; 
  fout[46] = 0.0; 
  fout[47] = 0.0; 
}

GKYL_CU_DH void translate_dim_gyrokinetic_3x2v_ser_p1_from_2x2v_p1(const double *flow, double *fout) 
{ 
  // flow: lower dimensional field to get DG coefficients from.
  // fout: field whose DG coefficients to populate.

  fout[0] = 1.4142135623730951*flow[0]; 
  fout[1] = 1.4142135623730951*flow[1]; 
  fout[2] = 0.0; 
  fout[3] = 1.4142135623730951*flow[2]; 
  fout[4] = 1.4142135623730951*flow[3]; 
  fout[5] = 1.4142135623730951*flow[4]; 
  fout[6] = 0.0; 
  fout[7] = 1.4142135623730951*flow[5]; 
  fout[8] = 0.0; 
  fout[9] = 1.4142135623730951*flow[6]; 
  fout[10] = 0.0; 
  fout[11] = 1.4142135623730951*flow[7]; 
  fout[12] = 1.4142135623730951*flow[8]; 
  fout[13] = 0.0; 
  fout[14] = 1.4142135623730951*flow[9]; 
  fout[15] = 1.4142135623730951*flow[10]; 
  fout[16] = 0.0; 
  fout[17] = 0.0; 
  fout[18] = 1.4142135623730951*flow[11]; 
  fout[19] = 0.0; 
  fout[20] = 0.0; 
  fout[21] = 1.4142135623730951*flow[12]; 
  fout[22] = 0.0; 
  fout[23] = 1.4142135623730951*flow[13]; 
  fout[24] = 0.0; 
  fout[25] = 1.4142135623730951*flow[14]; 
  fout[26] = 0.0; 
  fout[27] = 0.0; 
  fout[28] = 0.0; 
  fout[29] = 1.4142135623730951*flow[15]; 
  fout[30] = 0.0; 
  fout[31] = 0.0; 
  fout[32] = 1.4142135623730951*flow[16]; 
  fout[33] = 1.4142135623730951*flow[17]; 
  fout[34] = 0.0; 
  fout[35] = 1.4142135623730951*flow[18]; 
  fout[36] = 1.4142135623730951*flow[19]; 
  fout[37] = 0.0; 
  fout[38] = 1.4142135623730951*flow[20]; 
  fout[39] = 0.0; 
  fout[40] = 1.4142135623730951*flow[21]; 
  fout[41] = 0.0; 
  fout[42] = 1.4142135623730951*flow[22]; 
  fout[43] = 0.0; 
  fout[44] = 0.0; 
  fout[45] = 1.4142135623730951*flow[23]; 
  fout[46] = 0.0; 
  fout[47] = 0.0; 
}

