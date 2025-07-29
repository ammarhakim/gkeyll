#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void EB_alpha_quad_vx_2x1v_tensor_p1(const double *dxv, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  alpha_quad[0] += 0.5*Ex[3]-0.5*Ex[2]-0.5*Ex[1]+0.5*Ex[0]; 
  alpha_quad[1] += -(0.5*Ex[3])+0.5*Ex[2]-0.5*Ex[1]+0.5*Ex[0]; 
  alpha_quad[2] += -(0.5*Ex[3])-0.5*Ex[2]+0.5*Ex[1]+0.5*Ex[0]; 
  alpha_quad[3] += 0.5*Ex[3]+0.5*Ex[2]+0.5*Ex[1]+0.5*Ex[0]; 

} 
