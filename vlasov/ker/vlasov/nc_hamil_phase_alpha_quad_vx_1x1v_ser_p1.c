#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double wx1 = w[1]; 
  const double dv0 = dxv[1]; 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  alpha_quad[0] += -1.0*(0.7071067811865475*poisson_tensor_conf_x0[0]-0.7071067811865475*poisson_tensor_conf_x0[1])*(0.8660254037844386*hamil[1]-1.5*hamil[3])*dx10; 
  alpha_quad[1] += -1.0*(0.7071067811865475*poisson_tensor_conf_x0[1]+0.7071067811865475*poisson_tensor_conf_x0[0])*(0.8660254037844386*hamil[1]-1.5*hamil[3])*dx10; 

} 
