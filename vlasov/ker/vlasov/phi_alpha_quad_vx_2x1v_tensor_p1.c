#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
  const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*(0.8660254037844386*phi[1]-0.8660254037844386*phi[3]); 
  alpha_quad[1] += -dx10*(0.8660254037844386*phi[3]+0.8660254037844386*phi[1]); 
  alpha_quad[2] += -dx10*(0.8660254037844386*phi[1]-0.8660254037844386*phi[3]); 
  alpha_quad[3] += -dx10*(0.8660254037844386*phi[3]+0.8660254037844386*phi[1]); 

} 
