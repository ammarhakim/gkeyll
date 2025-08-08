#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void E_alpha_quad_vx_2x2v_ser_p1(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  double force_quad = 0.0;
  force_quad = 0.5*Ex[3]-0.5*Ex[2]-0.5*Ex[1]+0.5*Ex[0];
  alpha_quad[0] += force_quad;
  alpha_quad[1] += force_quad;

  force_quad = -(0.5*Ex[3])+0.5*Ex[2]-0.5*Ex[1]+0.5*Ex[0];
  alpha_quad[2] += force_quad;
  alpha_quad[3] += force_quad;

  force_quad = -(0.5*Ex[3])-0.5*Ex[2]+0.5*Ex[1]+0.5*Ex[0];
  alpha_quad[4] += force_quad;
  alpha_quad[5] += force_quad;

  force_quad = 0.5*Ex[3]+0.5*Ex[2]+0.5*Ex[1]+0.5*Ex[0];
  alpha_quad[6] += force_quad;
  alpha_quad[7] += force_quad;

} 
