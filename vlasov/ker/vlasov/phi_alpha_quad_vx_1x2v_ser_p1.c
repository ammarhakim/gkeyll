#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vx_1x2v_ser_p1(const double *dxv, 
  const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*(1.224744871391589*phi[1]); 
  alpha_quad[1] += -dx10*(1.224744871391589*phi[1]); 
  alpha_quad[2] += -dx10*(1.224744871391589*phi[1]); 
  alpha_quad[3] += -dx10*(1.224744871391589*phi[1]); 

} 
