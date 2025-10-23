#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*(1.936491673103709*hamil[7]+4.499999999999999*hamil[6]-2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[1] += -dx10*(1.936491673103709*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[2] += -dx10*(1.936491673103709*hamil[7]-4.499999999999999*hamil[6]+2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 

} 
