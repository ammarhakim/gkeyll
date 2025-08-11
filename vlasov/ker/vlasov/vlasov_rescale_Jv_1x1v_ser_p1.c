#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_ser_p1(const double *jacob_vel_gauss, const double *f_no_J, double* GKYL_RESTRICT Jf) 
{ 
  Jf[0] = f_no_J[0]; 
  Jf[1] = f_no_J[1]; 
  Jf[2] = f_no_J[2]; 
  Jf[3] = f_no_J[3]; 
} 
