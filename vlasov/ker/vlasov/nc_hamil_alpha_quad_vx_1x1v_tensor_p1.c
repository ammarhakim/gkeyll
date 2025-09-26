#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double wx1 = w[1]; 
  const double dv0 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  alpha_quad[0] += 0.0; 
  alpha_quad[1] += 0.0; 
  alpha_quad[2] += 0.0; 
  alpha_quad[3] += 0.0; 
  alpha_quad[4] += 0.0; 
  alpha_quad[5] += 0.0; 
  alpha_quad[6] += 0.0; 
  alpha_quad[7] += 0.0; 

} 
