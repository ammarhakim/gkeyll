#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void EB_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  alpha_quad[0] += 0.9681844646844028*Ex[2]-1.0546722811938851*Ex[1]+0.7071067811865475*Ex[0]; 
  alpha_quad[1] += -(0.5164305132317774*Ex[2])-0.41639003950091297*Ex[1]+0.7071067811865475*Ex[0]; 
  alpha_quad[2] += -(0.5164305132317774*Ex[2])+0.41639003950091297*Ex[1]+0.7071067811865475*Ex[0]; 
  alpha_quad[3] += 0.9681844646844028*Ex[2]+1.0546722811938851*Ex[1]+0.7071067811865475*Ex[0]; 

} 
