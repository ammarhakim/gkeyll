#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rescale_Jv_2x1v_ser_p1(const double *jacob_vel_gauss, const double *f_no_J, double* GKYL_RESTRICT Jf) 
{ 
  Jf[0] = f_no_J[0]; 
  Jf[1] = f_no_J[1]; 
  Jf[2] = f_no_J[2]; 
  Jf[3] = f_no_J[3]; 
  Jf[4] = f_no_J[4]; 
  Jf[5] = f_no_J[5]; 
  Jf[6] = f_no_J[6]; 
  Jf[7] = f_no_J[7]; 
} 
