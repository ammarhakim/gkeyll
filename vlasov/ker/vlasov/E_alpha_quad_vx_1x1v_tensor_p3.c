#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void E_alpha_quad_vx_1x1v_tensor_p3(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  double force_quad = 0.0;
  force_quad = (-0.5701294036773671*Ex[3])+0.9681844646844028*Ex[2]-1.054672281193885*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[0] += force_quad;

  force_quad = 0.7702725556588816*Ex[3]-0.5164305132317774*Ex[2]-0.416390039500913*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[1] += force_quad;

  force_quad = (-0.7702725556588816*Ex[3])-0.5164305132317774*Ex[2]+0.416390039500913*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[2] += force_quad;

  force_quad = 0.5701294036773671*Ex[3]+0.9681844646844028*Ex[2]+1.054672281193885*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[3] += force_quad;

} 
