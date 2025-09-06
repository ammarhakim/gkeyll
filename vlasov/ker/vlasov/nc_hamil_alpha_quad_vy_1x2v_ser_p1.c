#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_alpha_quad_vy_1x2v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double wx1 = w[1]; 
  const double wx2 = w[2]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[10]; 
  alpha_quad[0] += -(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[0]-0.7071067811865475*poisson_tensor_conf_1[1])*(wx2-0.5*dv11))-1.0*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7071067811865475*poisson_tensor_conf_0[1])*(wx1-0.2886751345948129*dv10)); 
  alpha_quad[1] += -(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[0]-0.7071067811865475*poisson_tensor_conf_1[1])*(wx2-0.5*dv11))-1.0*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7071067811865475*poisson_tensor_conf_0[1])*(wx1+0.2886751345948129*dv10)); 
  alpha_quad[2] += -(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv11))-1.0*(0.7071067811865475*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.2886751345948129*dv10)); 
  alpha_quad[3] += -(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv11))-1.0*(0.7071067811865475*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1+0.2886751345948129*dv10)); 

} 
