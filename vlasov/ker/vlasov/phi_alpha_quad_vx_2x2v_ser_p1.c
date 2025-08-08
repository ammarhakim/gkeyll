#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vx_2x2v_ser_p1(const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double force_quad = 0.0;
  force_quad = 0.8660254037844386*phi[1]-0.8660254037844386*phi[3];
  alpha_quad[0] -= dx10*force_quad;
  alpha_quad[1] -= dx10*force_quad;

  force_quad = 0.8660254037844386*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[2] -= dx10*force_quad;
  alpha_quad[3] -= dx10*force_quad;

  force_quad = 0.8660254037844386*phi[1]-0.8660254037844386*phi[3];
  alpha_quad[4] -= dx10*force_quad;
  alpha_quad[5] -= dx10*force_quad;

  force_quad = 0.8660254037844386*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[6] -= dx10*force_quad;
  alpha_quad[7] -= dx10*force_quad;

} 
