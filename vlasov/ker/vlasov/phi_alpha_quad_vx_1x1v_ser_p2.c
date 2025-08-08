#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vx_1x1v_ser_p2(const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double force_quad = 0.0;
  force_quad = 1.224744871391589*phi[1]-4.084728180770504*phi[2];
  alpha_quad[0] -= dx10*force_quad;

  force_quad = 1.224744871391589*phi[1]-1.6126716885136845*phi[2];
  alpha_quad[1] -= dx10*force_quad;

  force_quad = 1.6126716885136845*phi[2]+1.224744871391589*phi[1];
  alpha_quad[2] -= dx10*force_quad;

  force_quad = 4.084728180770504*phi[2]+1.224744871391589*phi[1];
  alpha_quad[3] -= dx10*force_quad;

} 
