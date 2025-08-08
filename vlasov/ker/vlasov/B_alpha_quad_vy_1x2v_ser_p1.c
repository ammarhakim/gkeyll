#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_alpha_quad_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *Bz = &qmem[10]; 
  double Bz_quad = 0.0;
  Bz_quad = 0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1];
  alpha_quad[0] -= 2.0/dxv[1]*(0.8660254037844386*hamil[1]-1.5*hamil[3])*Bz_quad;
  alpha_quad[1] -= 2.0/dxv[1]*(0.8660254037844386*hamil[1]-1.5*hamil[3])*Bz_quad;

  Bz_quad = 0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[2] -= 2.0/dxv[1]*(0.8660254037844386*hamil[1]-1.5*hamil[3])*Bz_quad;
  alpha_quad[3] -= 2.0/dxv[1]*(0.8660254037844386*hamil[1]-1.5*hamil[3])*Bz_quad;

} 
