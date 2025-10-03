#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*((-2.29128784747792*hamil[11])-9.306450766252633*hamil[10]+5.373081855095957*hamil[8]+1.936491673103709*hamil[7]+5.002749890427371*hamil[6]-2.888338995926613*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[1] += -dx10*((-2.29128784747792*hamil[11])+1.450606717756908*hamil[10]-0.8375081789852307*hamil[8]+1.936491673103709*hamil[7]+1.97511137974555*hamil[6]-1.140331086775586*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[2] += -dx10*((-2.29128784747792*hamil[11])+1.450606717756908*hamil[10]-0.8375081789852307*hamil[8]+1.936491673103709*hamil[7]-1.97511137974555*hamil[6]+1.140331086775586*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[3] += -dx10*((-2.29128784747792*hamil[11])-9.306450766252633*hamil[10]+5.373081855095957*hamil[8]+1.936491673103709*hamil[7]-5.002749890427371*hamil[6]+2.888338995926613*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 

} 
