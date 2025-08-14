#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  double force_quad = 0.0;
  force_quad = -(0.9373434910716291*Ex[3])+1.1569870650434417*Ex[2]-1.1098391188717989*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[0] += force_quad;

  force_quad = 0.7808504157580675*Ex[3]-0.10289451165398211*Ex[2]-0.6594875259537024*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[1] += force_quad;

  force_quad = 0.7071067811865475*Ex[0]-0.7905694150420947*Ex[2];
  alpha_quad[2] += force_quad;

  force_quad = -(0.7808504157580675*Ex[3])-0.10289451165398211*Ex[2]+0.6594875259537024*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[3] += force_quad;

  force_quad = 0.9373434910716291*Ex[3]+1.1569870650434417*Ex[2]+1.1098391188717989*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[4] += force_quad;

} 
