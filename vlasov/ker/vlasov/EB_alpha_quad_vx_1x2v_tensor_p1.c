#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void EB_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *Bz = &qmem[10]; 
  alpha_quad[0] += (0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1])*(0.8660254037844386*hamil[2]-1.5*hamil[3])*dv11-0.7071067811865475*Ex[1]+0.7071067811865475*Ex[0]; 
  alpha_quad[1] += (0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1])*(0.8660254037844386*hamil[2]-1.5*hamil[3])*dv11-0.7071067811865475*Ex[1]+0.7071067811865475*Ex[0]; 
  alpha_quad[2] += (0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0])*(0.8660254037844386*hamil[2]-1.5*hamil[3])*dv11+0.7071067811865475*Ex[1]+0.7071067811865475*Ex[0]; 
  alpha_quad[3] += (0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0])*(0.8660254037844386*hamil[2]-1.5*hamil[3])*dv11+0.7071067811865475*Ex[1]+0.7071067811865475*Ex[0]; 

} 
