#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_tensor_p1(const double *dxv, 
  const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  alpha_quad[0] += -dx11*(0.8660254037844386*phi[2]-0.8660254037844386*phi[3]); 
  alpha_quad[1] += -dx11*(0.8660254037844386*phi[2]-0.8660254037844386*phi[3]); 
  alpha_quad[2] += -dx11*(0.8660254037844386*phi[2]-0.8660254037844386*phi[3]); 
  alpha_quad[3] += -dx11*(0.8660254037844386*phi[2]-0.8660254037844386*phi[3]); 
  alpha_quad[4] += -dx11*(0.8660254037844386*phi[3]+0.8660254037844386*phi[2]); 
  alpha_quad[5] += -dx11*(0.8660254037844386*phi[3]+0.8660254037844386*phi[2]); 
  alpha_quad[6] += -dx11*(0.8660254037844386*phi[3]+0.8660254037844386*phi[2]); 
  alpha_quad[7] += -dx11*(0.8660254037844386*phi[3]+0.8660254037844386*phi[2]); 

} 
