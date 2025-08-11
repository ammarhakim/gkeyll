#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rescale_Jv_1x2v_ser_p2(const double *jacob_vel_gauss, const double *f_no_J, double* GKYL_RESTRICT Jf) 
{ 
  Jf[0] = f_no_J[0]; 
  Jf[1] = f_no_J[1]; 
  Jf[2] = f_no_J[2]; 
  Jf[3] = f_no_J[3]; 
  Jf[4] = f_no_J[4]; 
  Jf[5] = f_no_J[5]; 
  Jf[6] = f_no_J[6]; 
  Jf[7] = f_no_J[7]; 
  Jf[8] = f_no_J[8]; 
  Jf[9] = f_no_J[9]; 
  Jf[10] = f_no_J[10]; 
  Jf[11] = f_no_J[11]; 
  Jf[12] = f_no_J[12]; 
  Jf[13] = f_no_J[13]; 
  Jf[14] = f_no_J[14]; 
  Jf[15] = f_no_J[15]; 
  Jf[16] = f_no_J[16]; 
  Jf[17] = f_no_J[17]; 
  Jf[18] = f_no_J[18]; 
  Jf[19] = f_no_J[19]; 
} 
