#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vx_1x1v_tensor_p2(const double *dxv, 
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*(-(6.458522336955395*hamil[8])+1.9364916731037085*hamil[7]+5.002749890427371*hamil[6]-2.888338995926613*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[1] += -dx10*(-(2.549857826886422*hamil[8])+1.9364916731037085*hamil[7]+1.9751113797455495*hamil[6]-1.140331086775586*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[2] += -dx10*(2.549857826886422*hamil[8]+1.9364916731037085*hamil[7]-1.9751113797455495*hamil[6]+1.140331086775586*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[3] += -dx10*(6.458522336955395*hamil[8]+1.9364916731037085*hamil[7]-5.002749890427371*hamil[6]+2.888338995926613*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 

} 
