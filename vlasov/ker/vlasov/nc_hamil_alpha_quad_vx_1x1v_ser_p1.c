#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_alpha_quad_vx_1x1v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double wx1 = w[1]; 
  const double dv10 = 2.0/dxv[1]; 
  alpha_quad[0] += -dv10*(0.0); 
  alpha_quad[1] += -dv10*(0.0); 

} 
