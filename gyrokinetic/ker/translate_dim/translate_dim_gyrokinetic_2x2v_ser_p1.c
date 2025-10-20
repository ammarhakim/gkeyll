#include <gkyl_translate_dim_kernels.h> 

GKYL_CU_DH void translate_dim_gyrokinetic_2x2v_ser_p1_from_1x2v_p1(const double *flow, double *fout) 
{ 
  // flow: lower dimensional field to get DG coefficients from.
  // fout: field whose DG coefficients to populate.

  fout[0] = 1.4142135623730951*flow[0]; 
  fout[1] = 0.0; 
  fout[2] = 1.4142135623730951*flow[1]; 
  fout[3] = 1.4142135623730951*flow[2]; 
  fout[4] = 1.4142135623730951*flow[3]; 
  fout[5] = 0.0; 
  fout[6] = 0.0; 
  fout[7] = 1.4142135623730951*flow[4]; 
  fout[8] = 0.0; 
  fout[9] = 1.4142135623730951*flow[5]; 
  fout[10] = 1.4142135623730951*flow[6]; 
  fout[11] = 0.0; 
  fout[12] = 0.0; 
  fout[13] = 0.0; 
  fout[14] = 1.4142135623730951*flow[7]; 
  fout[15] = 0.0; 
  fout[16] = 1.4142135623730951*flow[8]; 
  fout[17] = 0.0; 
  fout[18] = 1.4142135623730951*flow[9]; 
  fout[19] = 1.4142135623730951*flow[10]; 
  fout[20] = 0.0; 
  fout[21] = 0.0; 
  fout[22] = 1.4142135623730951*flow[11]; 
  fout[23] = 0.0; 
}

