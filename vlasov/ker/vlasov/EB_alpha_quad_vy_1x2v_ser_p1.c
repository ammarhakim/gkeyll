#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void EB_alpha_quad_vy_1x2v_ser_p1(const double *dxv, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ey = &qmem[2]; 

  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bz = &qmem[10]; 
  alpha_quad[0] += -(1.0*(0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1])*(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10)-0.7071067811865475*Ey[1]+0.7071067811865475*Ey[0]; 
  alpha_quad[1] += -(1.0*(0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1])*(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10)-0.7071067811865475*Ey[1]+0.7071067811865475*Ey[0]; 
  alpha_quad[2] += -(1.0*(0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0])*(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10)+0.7071067811865475*Ey[1]+0.7071067811865475*Ey[0]; 
  alpha_quad[3] += -(1.0*(0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0])*(0.8660254037844386*hamil[1]-1.5*hamil[3])*dv10)+0.7071067811865475*Ey[1]+0.7071067811865475*Ey[0]; 

} 
