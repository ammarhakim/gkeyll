#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vx_1x1v_ser_p1(const double *dxv, 
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*(0.8660254037844386*hamil[1]-1.5*hamil[3]); 
  alpha_quad[1] += -dx10*(0.8660254037844386*hamil[1]-1.5*hamil[3]); 

} 
