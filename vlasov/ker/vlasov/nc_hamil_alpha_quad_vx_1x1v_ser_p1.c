#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_alpha_quad_vx_1x1v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double wx1 = w[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[2]; 
  alpha_quad[0] += -dx10*(1.224744871391589*hamil[1]*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7071067811865475*poisson_tensor_conf_0[1])*(wx1-0.5*dv10)); 
  alpha_quad[1] += -dx10*(1.224744871391589*hamil[1]*(0.7071067811865475*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.5*dv10)); 

} 
