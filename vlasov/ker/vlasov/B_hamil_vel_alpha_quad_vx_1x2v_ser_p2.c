#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_hamil_vel_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *Bz = &qmem[15]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[4]; 
  double Bz_quad = 0.0;
  Bz_quad = 0.6324555320336759*Bz[2]-0.9486832980505137*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[0] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(4.499999999999999*hamil[7]+1.9364916731037085*hamil[6]-2.5980762113533156*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[1] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(1.9364916731037085*hamil[6]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[2] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(-(4.499999999999999*hamil[7])+1.9364916731037085*hamil[6]+2.5980762113533156*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;

  Bz_quad = 0.7071067811865475*Bz[0]-0.7905694150420947*Bz[2];
  alpha_quad[3] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(4.499999999999999*hamil[7]+1.9364916731037085*hamil[6]-2.5980762113533156*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[4] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(1.9364916731037085*hamil[6]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[5] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(-(4.499999999999999*hamil[7])+1.9364916731037085*hamil[6]+2.5980762113533156*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;

  Bz_quad = 0.6324555320336759*Bz[2]+0.9486832980505137*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[6] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(4.499999999999999*hamil[7]+1.9364916731037085*hamil[6]-2.5980762113533156*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[7] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(1.9364916731037085*hamil[6]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[8] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(-(4.499999999999999*hamil[7])+1.9364916731037085*hamil[6]+2.5980762113533156*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;

} 
