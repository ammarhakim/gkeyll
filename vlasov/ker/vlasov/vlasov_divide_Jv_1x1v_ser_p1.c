#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_divide_Jv_1x1v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J) 
{ 
  f_no_J[0] = Jf[0]; 
  f_no_J[1] = Jf[1]; 
  f_no_J[2] = Jf[2]; 
  f_no_J[3] = Jf[3]; 
} 
