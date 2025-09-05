#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*((-2.29128784747792*hamil[11])-10.67445557248417*hamil[10]+6.162899798226436*hamil[8]+1.936491673103709*hamil[7]+5.264429177983873*hamil[6]-3.039419603038709*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[1] += -dx10*((-2.29128784747792*hamil[11])-1.545746280731402*hamil[10]+0.8924370312791381*hamil[8]+1.936491673103709*hamil[7]+3.128224005724662*hamil[6]-1.806080971790583*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[2] += -dx10*((-2.29128784747792*hamil[11])+3.43693177121688*hamil[10]-1.984313483298443*hamil[8]+1.936491673103709*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[3] += -dx10*((-2.29128784747792*hamil[11])-1.545746280731402*hamil[10]+0.8924370312791381*hamil[8]+1.936491673103709*hamil[7]-3.128224005724662*hamil[6]+1.806080971790583*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 
  alpha_quad[4] += -dx10*((-2.29128784747792*hamil[11])-10.67445557248417*hamil[10]+6.162899798226436*hamil[8]+1.936491673103709*hamil[7]-5.264429177983873*hamil[6]+3.039419603038709*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1]); 

} 
